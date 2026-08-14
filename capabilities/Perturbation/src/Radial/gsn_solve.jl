# :GSN backend — the unified GSN radial solve pipeline (real ω). Collapses the upstream
# Solutions.jl/ComplexFrequencies.jl fork (GSN_PLAN §1): the real-ω path solved here is the
# β=0, single-segment case of the contour-parametrized pipeline; the complex-ω contour
# (Im ω ≠ 0, QNM use) is explicitly QUARANTINED per ruling OQ-4 — see `_gsn_solve_X`.
#
# GSN master equation (arXiv:2306.16469 Eq. (13); upstream Solutions.jl l.48):
#     d²X/dr*² − sF(r)·dX/dr* − sU(r)·X = 0
# Riccati form (X = e^{iΦ}; upstream Solutions.jl l.34):
#     Φ″ = −i·sU + sF·Φ′ − i·(Φ′)²
#
# FIREWALL: integration is driven by the VENDORED NumericsCore integrators (never OrdinaryDiffEq;
# GSN_PLAN §3a, user-authorized redesign of upstream's AutoVern9(Rosenbrock23)):
#   Float64 (and generic T): adaptive dp5_step! mesh with cubic-Hermite dense output — the exact
#     pattern proven in-tree by radial_sn.jl (_sn_integrate), complex state ↔ real 4-vector.
#   BigFloat: gbs_integrate_to_targets! (high-order Gragg-Bulirsch-Stoer; dp5's order-5 step
#     count is infeasible at 2^(20−bits) tolerances), on-demand per evaluation point.
# Step tolerance defaults to the precision-inherited NumericsCore.rk_tol(T).

"""
    _GSNXSolution

Carrier for one GSN homogeneous solution (In or Up) at real ω: the adaptive numerical mesh of
(X, X′, X″) on [rsin, rsout] (the type-generic DP5 mesh is stored for EVERY precision — under
deviation D-GBS the BigFloat tier also solves on this mesh; the GBS on-demand `demand` cache
remains only as a defensive fallback for an empty mesh), the SN-frame asymptotic amplitudes
(unit transmission), and the mode/domain metadata needed for the semianalytic continuation
beyond the numerical domain. Evaluate via `_gsn_X_at`.
"""
struct _GSNXSolution{T<:Real}
    s::Int; l::Int; m::Int
    a::T; omega::T; lambda::Complex{T}
    bc::Symbol                    # :In | :Up
    rsin::T; rsout::T
    hor_order::Int; inf_order::Int
    form::Symbol                  # :riccati | :linear (what the mesh stores)
    mesh_rs::Vector{T}            # ascending r*
    mesh_u::Vector{Complex{T}}    # X (linear) or Φ (riccati)
    mesh_du::Vector{Complex{T}}
    mesh_d2u::Vector{Complex{T}}
    inc_SN::Complex{T}            # Binc (In) / Cinc (Up); SN frame, transmission ≡ 1
    ref_SN::Complex{T}            # Bref (In) / Cref (Up)
    demand::Dict{T,Tuple{Complex{T},Complex{T}}}  # BigFloat on-demand (X, X′) cache
end

# ------------------------------------------------------------------------------------------------
# RHS closures (complex 2-state carried as real 4-vector [Re u1, Im u1, Re u2, Im u2] — the
# radial_sn.jl pattern that lets the real-typed NumericsCore primitives drive a complex ODE)
# ------------------------------------------------------------------------------------------------

"Linear-form RHS: u=[X, X′], u′=[X′, sF·X′ + sU·X]. [GSN Solutions.jl GSN_linear_eqn]"
function _gsn_rhs_linear!(g, rs, y, s::Int, m::Int, a, omega, lambda)
    r = _gsn_r_from_rstar(a, rs)
    F = _gsn_sF(s, m, a, omega, lambda, r)
    U = _gsn_sU(s, m, a, omega, lambda, r)
    X = complex(y[1], y[2]); dX = complex(y[3], y[4])
    d2X = F*dX + U*X
    g[1] = y[3]; g[2] = y[4]; g[3] = real(d2X); g[4] = imag(d2X)
    return g
end

"Riccati-form RHS: u=[Φ, Φ′], Φ″ = −i·sU + sF·Φ′ − i·Φ′². [GSN Solutions.jl GSN_Riccati_eqn]"
function _gsn_rhs_riccati!(g, rs, y, s::Int, m::Int, a, omega, lambda)
    r = _gsn_r_from_rstar(a, rs)
    F = _gsn_sF(s, m, a, omega, lambda, r)
    U = _gsn_sU(s, m, a, omega, lambda, r)
    dP = complex(y[3], y[4])
    d2P = -im*U + F*dP - im*dP*dP
    g[1] = y[3]; g[2] = y[4]; g[3] = real(d2P); g[4] = imag(d2P)
    return g
end

# ------------------------------------------------------------------------------------------------
# Integrators
# ------------------------------------------------------------------------------------------------

"""
    _gsn_integrate_mesh(f!, u0, du0, rs0, rs1; tol) -> (rs, u, du, d2u) or nothing

Adaptive DP5 integration of the complex 2-state ODE from rs0 to rs1 (either direction), storing
the full mesh for cubic-Hermite dense output (radial_sn.jl `_sn_integrate` pattern, made
direction- and type-generic). Returns `nothing` on breakdown (non-finite state / step-budget
exhaustion) — the signal the `:auto` driver uses to fall back from Riccati to linear form,
replacing the upstream LoggingExtras warning-filter mechanism (GSN_PLAN §1 gsn_solve row).
"""
function _gsn_integrate_mesh(f!, u0::Complex{T}, du0::Complex{T}, rs0::T, rs1::T; tol::T) where {T<:Real}
    dirp = rs1 > rs0
    span = abs(rs1 - rs0)
    h = (dirp ? one(T) : -one(T)) * min(span, max(span/200, T(1)/100))
    hmin = span*eps(T)*16 + eps(T)
    y = T[real(u0), imag(u0), real(du0), imag(du0)]
    g = similar(y); f!(g, rs0, y); k1 = copy(g)
    rcur = rs0
    rsm = T[rcur]; um = Complex{T}[u0]; dum = Complex{T}[du0]
    d2m = Complex{T}[complex(g[3], g[4])]
    maxsteps = 5_000_000; nsteps = 0
    while (dirp ? rcur < rs1 : rcur > rs1)
        nsteps += 1
        nsteps > maxsteps && return nothing
        if (dirp && rcur + h > rs1) || (!dirp && rcur + h < rs1)
            h = rs1 - rcur
        end
        err, y5, k7 = dp5_step!(f!, rcur, y, k1, h, tol)
        if !all(isfinite, y5)
            return nothing
        end
        if err <= 1 || abs(h) <= hmin
            rcur += h; y = y5; k1 = k7
            push!(rsm, rcur)
            push!(um, complex(y5[1], y5[2])); push!(dum, complex(y5[3], y5[4]))
            push!(d2m, complex(k7[3], k7[4]))
            h *= clamp(T(9)/10 * err^(-one(T)/5), T(1)/5, T(5))
        else
            h *= max(T(1)/5, T(9)/10 * err^(-one(T)/5))
        end
    end
    if !dirp
        reverse!(rsm); reverse!(um); reverse!(dum); reverse!(d2m)
    end
    return rsm, um, dum, d2m
end

"""
    _gsn_integrate_targets(f!, u0, du0, rs0, targets; tol) -> Vector{(u, du)}

High-order GBS integration (NumericsCore.gbs_integrate_to_targets!) of the complex 2-state ODE
from rs0 to each target r* in sequence — the BigFloat path (GSN_PLAN §3a). `targets` must be
monotone in the direction of integration.
"""
function _gsn_integrate_targets(f!, u0::Complex{T}, du0::Complex{T}, rs0::T, targets::Vector{T}; tol::T) where {T<:Real}
    dir = (isempty(targets) || targets[end] >= rs0) ? 1 : -1
    kmax = max(gbs_kmax(T), 6)
    y0 = T[real(u0), imag(u0), real(du0), imag(du0)]
    Y = gbs_integrate_to_targets!(f!, y0, rs0, targets, dir; tol=tol, kmax=kmax)
    return [(complex(Y[i, 1], Y[i, 2]), complex(Y[i, 3], Y[i, 4])) for i in 1:length(targets)]
end

"""
    _gsn_hermite(rs, u, du, d2u, x) -> (u(x), u′(x))

QUINTIC two-point Hermite dense output on the stored mesh, using the (u, u′, u″) triplet at
both bracketing nodes (u″ comes free from the FSAL stage). Interpolation error is O(h⁶) —
at or below the DP5 step error (order 5) for every accepted step, so interior evaluations
inherit the integration tolerance. (The cubic variant lost ~2 orders on large smooth-region
Riccati steps — measured up to 1.4e-7 on the fixture net — because its O(h⁴) error is BELOW
the integrator order.)
"""
function _gsn_hermite(rs::Vector{T}, u::Vector{Complex{T}}, du::Vector{Complex{T}},
                      d2u::Vector{Complex{T}}, x::T) where {T}
    (x >= rs[1] - 16*eps(T)*(1 + abs(rs[1])) && x <= rs[end] + 16*eps(T)*(1 + abs(rs[end]))) ||
        error("GSN: r*=$x outside the integrated domain [$(rs[1]), $(rs[end])].")
    i = clamp(searchsortedlast(rs, x), 1, length(rs) - 1)
    r0 = rs[i]; r1 = rs[i + 1]; hh = r1 - r0; t = (x - r0)/hh
    f0 = u[i];  g0 = hh*du[i];  c0 = hh^2*d2u[i]
    f1 = u[i+1]; g1 = hh*du[i+1]; c1 = hh^2*d2u[i+1]
    t2 = t*t; t3 = t2*t; t4 = t3*t; t5 = t4*t
    # quintic Hermite basis (value): H0..H5 for f0, g0, c0, f1, g1, c1
    H0 = 1 - 10t3 + 15t4 - 6t5
    H1 = t - 6t3 + 8t4 - 3t5
    H2 = (t2 - 3t3 + 3t4 - t5)/2
    H3 = 10t3 - 15t4 + 6t5
    H4 = -4t3 + 7t4 - 3t5
    H5 = (t3 - 2t4 + t5)/2
    val = H0*f0 + H1*g0 + H2*c0 + H3*f1 + H4*g1 + H5*c1
    # d/dt of the basis, then /hh for d/dx
    D0 = -30t2 + 60t3 - 30t4
    D1 = 1 - 18t2 + 32t3 - 15t4
    D2 = (2t - 9t2 + 12t3 - 5t4)/2
    D3 = 30t2 - 60t3 + 30t4
    D4 = -12t2 + 28t3 - 15t4
    D5 = (3t2 - 8t3 + 5t4)/2
    der = (D0*f0 + D1*g0 + D2*c0 + D3*f1 + D4*g1 + D5*c1)/hh
    return val, der
end

# ------------------------------------------------------------------------------------------------
# Amplitude extraction — ONE 2×2 matching helper replacing the upstream four near-copies
# (Solutions.jl BrefBinc_SN_from_Xin / CrefCinc_SN_from_Xup; GSN_PLAN §3e).
# ------------------------------------------------------------------------------------------------

"""
    _gsn_match_amplitudes(fin, dfin, fout, dfout, k, rsb, r, drdrs, X, dX) -> (c_in, c_out)

Solve the 2×2 boundary match X = c_in·f_in·e^{−ikr*} + c_out·f_out·e^{+ikr*} (and its r*
derivative) at the boundary r* = rsb for the ingoing/outgoing mode amplitudes. `k` is ω at ∞
or p = ω − mΩ_H at the horizon; f/df are the ansatz radial factors at r. [algebra of GSN
Solutions.jl l.522-606]
"""
function _gsn_match_amplitudes(fin, dfin, fout, dfout, k, rsb, drdrs, X, dX)
    A1 = fin * exp(-im*k*rsb)
    A2 = fout * exp(im*k*rsb)
    A3 = (drdrs*dfin - im*k*fin) * exp(-im*k*rsb)
    A4 = (drdrs*dfout + im*k*fout) * exp(im*k*rsb)
    D = A2*A3 - A1*A4
    c_in  = -(A4*X - A2*dX)/D
    c_out = -(-A3*X + A1*dX)/D
    return c_in, c_out
end

"(Bref, Binc) of an In solution from its boundary values at rsout. [GSN BrefBinc_SN_from_Xin]"
function _gsn_BrefBinc(s::Int, m::Int, a, omega, lambda, rsout, X, dX, order::Int, ::Type{T}) where {T}
    rout = _gsn_r_from_rstar(a, rsout)
    cin(i) = _gsn_coeff_in_inf(s, m, a, omega, lambda, i, T)
    cout(i) = _gsn_coeff_out_inf(s, m, a, omega, lambda, i, T)
    fin = _gsn_fansatz(cin, omega, rout; order=order); dfin = _gsn_dfansatz_dr(cin, omega, rout; order=order)
    fout = _gsn_fansatz(cout, omega, rout; order=order); dfout = _gsn_dfansatz_dr(cout, omega, rout; order=order)
    drdrs = _gsn_Delta(a, rout)/(rout^2 + a^2)
    c_in, c_out = _gsn_match_amplitudes(fin, dfin, fout, dfout, omega, rsout, drdrs, X, dX)
    return c_out, c_in                      # Bref = outgoing coefficient, Binc = ingoing
end

"(Cref, Cinc) of an Up solution from its boundary values at rsin. [GSN CrefCinc_SN_from_Xup]"
function _gsn_CrefCinc(s::Int, m::Int, a, omega, lambda, rsin, X, dX, order::Int, ::Type{T}) where {T}
    rin = _gsn_r_from_rstar(a, rsin)
    p = omega - m*_gsn_Omega_h(a)
    gin(i) = _gsn_coeff_in_hor(s, m, a, omega, lambda, i, T)
    gout(i) = _gsn_coeff_out_hor(s, m, a, omega, lambda, i, T)
    fi = _gsn_gansatz(gin, a, rin; order=order); dfi = _gsn_dgansatz_dr(gin, a, rin; order=order)
    fo = _gsn_gansatz(gout, a, rin; order=order); dfo = _gsn_dgansatz_dr(gout, a, rin; order=order)
    drdrs = _gsn_Delta(a, rin)/(rin^2 + a^2)
    c_in, c_out = _gsn_match_amplitudes(fi, dfi, fo, dfo, p, rsin, drdrs, X, dX)
    return c_in, c_out                      # Cref = ingoing coefficient, Cinc = outgoing
end

# ------------------------------------------------------------------------------------------------
# The solve driver (real ω; :auto = Riccati with linear fallback)
# ------------------------------------------------------------------------------------------------

"""
    _gsn_solve_X(s, l, m, a, omega, lambda, bc; rsin, rsout, hor_order, inf_order, tol, form)
        -> _GSNXSolution

Solve the GSN equation for one boundary condition at REAL ω: In (IC at rsin via `_gsn_Xin_bc`,
integrated outward) or Up (IC at rsout via `_gsn_Xup_bc`, integrated inward), then extract the
SN-frame amplitudes at the far boundary. `form=:auto` tries the Riccati form first and falls
back to the linear form on integrator breakdown (upstream `method="auto"` semantics).
Complex ω (QNM rotated-contour machinery, upstream ComplexFrequencies.jl) is NOT ported in
GSN-2 — [HOLE: quarantined to GSN-3 per ruling OQ-4; evidence: upstream solves it on a
ρ-parametrized rotated contour with β± = −arg(ω), −arg(p) segment rotations, which requires
the contour ODE r(ρ) and rsmp optimizer, deferred with the inhomogeneous subsystem].
"""
function _gsn_solve_X(s::Int, l::Int, m::Int, a, omega, lambda, bc::Symbol;
                      rsin::Real, rsout::Real, hor_order::Int=3, inf_order::Int=6,
                      tol::Union{Nothing,Real}=nothing, form::Symbol=:auto)
    isreal(omega) || throw(ArgumentError(
        "GSN backend: complex ω is not ported in GSN-2 (quarantined to GSN-3 with the contour machinery)."))
    T = typeof(float(real(a) + real(omega)))
    aT = T(a); omegaT = T(real(omega)); rsinT = T(rsin); rsoutT = T(rsout)
    lamT = Complex{T}(lambda)
    rsinT < rsoutT || throw(DomainError(rsout, "rsout must be larger than rsin"))
    # Default tolerance: precision-inherited rk_tol for Float64; for BigFloat the radial solve
    # runs the SAME type-generic DP5 mesh at a default 1e-12 truncation [deviation D-GBS,
    # reported with evidence]: the GBS phase-form tier targeted by GSN_PLAN §3a has a MEASURED
    # non-h² (h³-scaling) extrapolation floor near the horizon (probes 2026-07-02/03: tol=1e-8 →
    # 356 RHS calls over [-10,-8] at 256 bits; tol=1e-15 → 214k calls, ladder bottoms then
    # rises; and at tol=1e-10 the 256- and 512-bit amplitudes agreed with each other but sat
    # 5e-6 from the Float64 result — the floor biases accuracy, not just cost). DP5@1e-12 keeps
    # BigFloat runs exact in ARITHMETIC (precision threading, cross-precision determinism) with
    # an algorithmic truncation shared across precisions; deep-tolerance BigFloat integration is
    # deferred to GSN-2b with the contour study.
    tolT = tol === nothing ? (T <: BigFloat ? T(1)/T(10)^12 : T(rk_tol(T))) : T(tol)

    if bc === :In
        X0, dX0 = _gsn_Xin_bc(s, m, aT, omegaT, lamT, rsinT, hor_order, T)
        rs0, rs1 = rsinT, rsoutT
    elseif bc === :Up
        X0, dX0 = _gsn_Xup_bc(s, m, aT, omegaT, lamT, rsoutT, inf_order, T)
        rs0, rs1 = rsoutT, rsinT
    else
        throw(ArgumentError("GSN backend: bc must be :In or :Up (internal Out/Down are amplitude combinations)"))
    end

    fL!(g, rs, y) = _gsn_rhs_linear!(g, rs, y, s, m, aT, omegaT, lamT)
    fR!(g, rs, y) = _gsn_rhs_riccati!(g, rs, y, s, m, aT, omegaT, lamT)

    mesh = nothing; used_form = form
    if form === :auto || form === :riccati
        P0 = -im*log(X0); dP0 = -im*dX0/X0     # [GSN PhiPhiprime_from_XXprime]
        mesh = _gsn_integrate_mesh(fR!, P0, dP0, rs0, rs1; tol=tolT)
        used_form = :riccati
    end
    if mesh === nothing && (form === :auto || form === :linear)
        mesh = _gsn_integrate_mesh(fL!, X0, dX0, rs0, rs1; tol=tolT)
        used_form = :linear
    end
    mesh === nothing && error("GSN: both Riccati and linear integrations broke down for ($s,$l,$m,a=$a,ω=$omega,$bc).")
    rsm, um, dum, d2m = mesh
    # boundary state at the far end for amplitude extraction
    ub, dub = used_form === :riccati ?
        _gsn_XXp_from_Phi(_gsn_hermite(rsm, um, dum, d2m, rs1)...) :
        _gsn_hermite(rsm, um, dum, d2m, rs1)

    if bc === :In
        refA, incA = _gsn_BrefBinc(s, m, aT, omegaT, lamT, rsoutT, ub, dub, inf_order, T)
    else
        refA, incA = _gsn_CrefCinc(s, m, aT, omegaT, lamT, rsinT, ub, dub, hor_order, T)
    end

    return _GSNXSolution{T}(s, l, m, aT, omegaT, lamT, bc, rsinT, rsoutT, hor_order, inf_order,
                            used_form, rsm, um, dum, d2m, Complex{T}(incA), Complex{T}(refA),
                            Dict{T,Tuple{Complex{T},Complex{T}}}())
end

"X = e^{iΦ}, X′ = i·e^{iΦ}·Φ′. [GSN XXprime_from_PhiPhiprime]"
_gsn_XXp_from_Phi(P, dP) = (exp(im*P), im*exp(im*P)*dP)

"""
    _gsn_X_at(sol::_GSNXSolution, rs) -> (X, dX/dr*)

Semianalytic GSN solution evaluation: numerical mesh (Hermite; Riccati meshes are mapped by
X = e^{iΦ}) inside [rsin, rsout]; the asymptotic ansatz continuation outside — In: pure
transmission g·e^{−ipr*} below rsin / Bref·f_out·e^{iωr*} + Binc·f_in·e^{−iωr*} above rsout;
Up: pure transmission f_out·e^{iωr*} above rsout / Cinc·g_out·e^{ipr*} + Cref·g_in·e^{−ipr*}
below rsin. [GSN Solutions.jl semianalytical_Xin/Xup] On the BigFloat path interior values are
GBS-integrated on demand from the IC boundary and cached per r*.
"""
function _gsn_X_at(sol::_GSNXSolution{T}, rs::Real) where {T}
    rsT = T(rs)
    s, m, a, omega, lambda = sol.s, sol.m, sol.a, sol.omega, sol.lambda
    if rsT < sol.rsin
        r = _gsn_r_from_rstar(a, rsT)
        p = omega - m*_gsn_Omega_h(a)
        drdrs = _gsn_Delta(a, r)/(r^2 + a^2)
        gin(i) = _gsn_coeff_in_hor(s, m, a, omega, lambda, i, T)
        gi = _gsn_gansatz(gin, a, r; order=sol.hor_order)
        dgi = _gsn_dgansatz_dr(gin, a, r; order=sol.hor_order)
        ein = exp(-im*p*rsT)
        if sol.bc === :In
            return ein*gi, ein*(drdrs*dgi - im*p*gi)
        else
            gout(i) = _gsn_coeff_out_hor(s, m, a, omega, lambda, i, T)
            go = _gsn_gansatz(gout, a, r; order=sol.hor_order)
            dgo = _gsn_dgansatz_dr(gout, a, r; order=sol.hor_order)
            eout = exp(im*p*rsT)
            X = sol.inc_SN*go*eout + sol.ref_SN*gi*ein
            dX = sol.inc_SN*(drdrs*dgo + im*p*go)*eout + sol.ref_SN*(drdrs*dgi - im*p*gi)*ein
            return X, dX
        end
    elseif rsT > sol.rsout
        r = _gsn_r_from_rstar(a, rsT)
        drdrs = _gsn_Delta(a, r)/(r^2 + a^2)
        cout(i) = _gsn_coeff_out_inf(s, m, a, omega, lambda, i, T)
        fo = _gsn_fansatz(cout, omega, r; order=sol.inf_order)
        dfo = _gsn_dfansatz_dr(cout, omega, r; order=sol.inf_order)
        eout = exp(im*omega*rsT)
        if sol.bc === :Up
            return eout*fo, eout*(drdrs*dfo + im*omega*fo)
        else
            cin(i) = _gsn_coeff_in_inf(s, m, a, omega, lambda, i, T)
            fi = _gsn_fansatz(cin, omega, r; order=sol.inf_order)
            dfi = _gsn_dfansatz_dr(cin, omega, r; order=sol.inf_order)
            ein = exp(-im*omega*rsT)
            X = sol.ref_SN*fo*eout + sol.inc_SN*fi*ein
            dX = sol.ref_SN*(drdrs*dfo + im*omega*fo)*eout + sol.inc_SN*(drdrs*dfi - im*omega*fi)*ein
            return X, dX
        end
    else
        if !isempty(sol.mesh_rs)
            u, du = _gsn_hermite(sol.mesh_rs, sol.mesh_u, sol.mesh_du, sol.mesh_d2u, rsT)
            return sol.form === :riccati ? _gsn_XXp_from_Phi(u, du) : (u, du)
        end
        # BigFloat on-demand path
        haskey(sol.demand, rsT) && return sol.demand[rsT]
        if sol.bc === :In
            X0, dX0 = _gsn_Xin_bc(s, m, a, omega, lambda, sol.rsin, sol.hor_order, T)
            rs0 = sol.rsin
        else
            X0, dX0 = _gsn_Xup_bc(s, m, a, omega, lambda, sol.rsout, sol.inf_order, T)
            rs0 = sol.rsout
        end
        # Riccati (phase) form — the only BigFloat-feasible representation (see _gsn_solve_X)
        fR!(g, x, y) = _gsn_rhs_riccati!(g, x, y, s, m, a, omega, lambda)
        P0 = -im*log(X0); dP0 = -im*dX0/X0
        st = _gsn_integrate_targets(fR!, P0, dP0, rs0, [rsT];
                                    tol=(T <: BigFloat ? T(1)/T(10)^10 : T(rk_tol(T))))  # D-GBS floor
        val = _gsn_XXp_from_Phi(st[1][1], st[1][2])
        sol.demand[rsT] = val
        return val
    end
end

# ------------------------------------------------------------------------------------------------
# Internal gates
# ------------------------------------------------------------------------------------------------

"""
    _gsn_sanity_identity(Xin::_GSNXSolution, Xup::_GSNXSolution) -> (actual, expected)

The In/Up cross-normalization identity Binc/Cinc = p·c₀/(ω·η(r₊)) [arXiv:2306.16469 Eq. (42);
GSN Solutions.jl check_XinXup_sanity]. Used as the auto order-widening criterion and as a test
gate.
"""
function _gsn_sanity_identity(Xin::_GSNXSolution{T}, Xup::_GSNXSolution{T}) where {T}
    s, m, a, omega, lambda = Xin.s, Xin.m, Xin.a, Xin.omega, Xin.lambda
    p = omega - m*_gsn_Omega_h(a)
    c0 = _gsn_eta_coefficient(s, m, a, omega, lambda, 0)
    eta_rp = _gsn_eta(s, m, a, omega, lambda, _gsn_rp(a))
    expected = (p*c0)/(omega*eta_rp)
    actual = Xin.inc_SN/Xup.inc_SN
    return actual, expected
end

"""
    _gsn_wronskian(Xin::_GSNXSolution, Xup::_GSNXSolution, rs) -> Complex

Scaled GSN Wronskian W = (X_in·X_up′ − X_in′·X_up)/η(r) — an r*-independent constant for exact
solutions; its drift across the domain is an integration-quality gate. [GSN Solutions.jl
scaled_Wronskian_GSN; the upstream Phi-variant with undefined variables (defect D4) is DROPPED]
"""
function _gsn_wronskian(Xin::_GSNXSolution{T}, Xup::_GSNXSolution{T}, rs::Real) where {T}
    r = _gsn_r_from_rstar(Xin.a, T(rs))
    Xi, dXi = _gsn_X_at(Xin, rs)
    Xu, dXu = _gsn_X_at(Xup, rs)
    return (Xi*dXu - dXi*Xu)/_gsn_eta(Xin.s, Xin.m, Xin.a, Xin.omega, Xin.lambda, r)
end
