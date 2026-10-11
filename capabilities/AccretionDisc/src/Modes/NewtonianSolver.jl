# Miranda & Rafikov linear disc-mode SOLVER (arXiv:1811.09628 method): integrate the enthalpy master
# equation δh″ = Ψ − C₁δh′ − C₀δh outward and inward from corotation, impose the WKB outgoing-wave
# condition at adaptive boundaries, refine, then reconstruct δΣ, δu_r, δu_φ.
#
# Ported from the legacy `DiscModes/solver_newt.jl` (verified against the authors' `main.py`) together
# with the two solver variants that were FIGURE-LOCAL in the legacy MR19 reproduction scripts
# (`results/emri-disk/figures/1811.09628/fast_solver.jl` — the cached-potential throughput path — and
# `homog_solver.jl` — the homogeneous / passive-wake problem of MR19 §"passive", figures 10-12).
# Both are physics and now live here, once.
#
# The adaptive DP5 grid integrator is `NumericsCore.integrate_grid` (the single vendored ODE family);
# the cached-potential path interpolates with `NumericsCore.local_interp`.

using NumericsCore: integrate_grid, local_interp

"Log-spaced radial grid r_in … r_out with `n` points."
loggrid(d::NewtonianDisc{T}, r_in, r_out, n) where {T} =
    T(10) .^ range(log10(T(r_in)), log10(T(r_out)); length = n)

"Corotation index: argmin |mΩ(r) − ω_p| (the same radius for every m)."
function corotation_index(d::NewtonianDisc, x, m)
    best = 1; bv = abs(m * Omega(d, x[1]) - omega_pat(d, m))
    @inbounds for i in 2:length(x)
        v = abs(m * Omega(d, x[i]) - omega_pat(d, m)); v < bv && (bv = v; best = i)
    end
    best
end

"""
    master_rhs(d, r, y, m, ih) -> NTuple{4}

Master-equation right-hand side as a real 4-vector `y = [Re δh, Im δh, Re δh′, Im δh′]`.
`ih = 1` includes the secondary's forcing potential, `ih = 0` gives the homogeneous problem.
"""
function master_rhs(d::NewtonianDisc, r, y, m, ih)
    Om = Omega(d, r); dOm = dOmegadr(d, r); Sig = Sigma(d, r); dSig = dSigmadr(d, r)
    cs = csadi(d, r); ot = omega_pat(d, m) - m * Om
    Dv = Dfreq(d, r, m); dD = dDdr(d, r, m)
    LSi = LSinv(d, r); dLSi = dLSinvdr(d, r); LTi = LTinv(d, r); dLTi = dLTinvdr(d, r)
    N2 = Nr2(d, r); dN2 = dNr2dr(d, r); t = tc(d, r); dt = dtcdr(d, r)
    if ih == 1
        Ph = Phi_m(d, r, m); dPh = dPhidr_m(d, r, m); d2Ph = d2Phidr2_m(d, r, m)
    else
        Ph = dPh = d2Ph = zero(r)
    end
    C1R, C1I, C0R, C0I, PsiR, PsiI = master_equation_cool(
        r, d.gamma, cs, Om, dOm, Sig, dSig, LTi, dLTi, LSi, dLSi, ot, t, dt, Dv, dD, N2, dN2,
        d.pole_disp, Ph, dPh, d2Ph, m)
    @inbounds (y[3], y[4],
        PsiR - C1R * y[3] + C1I * y[4] - C0R * y[1] + C0I * y[2],
        PsiI - C1R * y[4] - C1I * y[3] - C0R * y[2] - C0I * y[1])
end

"""
    MasterRHS(d, m, ih)

Callable master-equation right-hand side with BOTH an in-place method `(out, r, y)` (zero-allocation,
the form `NumericsCore.integrate_grid` detects and prefers) and the allocating `(r, y)` form.
"""
struct MasterRHS{D,I}
    d::D
    m::I
    ih::I
end
@inline function (g::MasterRHS)(out::AbstractVector, r, y)
    d1, d2, d3, d4 = master_rhs(g.d, r, y, g.m, g.ih)
    @inbounds (out[1] = d1; out[2] = d2; out[3] = d3; out[4] = d4)
    out
end
@inline (g::MasterRHS)(r, y) = collect(master_rhs(g.d, r, y, g.m, g.ih))

"""
    integrate_master(d, x, bc, m, ih; rtol, atol) -> complex (n, 2)

Integrate the master equation on grid `x` outward AND inward from corotation with the state
`bc = (δh, δh′)` there; column 1 is δh, column 2 is δh′.
"""
function integrate_master(d::NewtonianDisc{T}, x, bc, m, ih; rtol = 1e-10, atol = 1e-12) where {T}
    ico = corotation_index(d, x, m)
    y0 = T[real(bc[1]), imag(bc[1]), real(bc[2]), imag(bc[2])]
    frhs = MasterRHS(d, m, ih)
    Yout = integrate_grid(frhs, y0, x, ico, +1; rtol = rtol, atol = atol)
    Yin  = integrate_grid(frhs, y0, x, ico, -1; rtol = rtol, atol = atol)
    n = length(x); f = zeros(Complex{T}, n, 2)
    @inbounds for i in 1:n
        Y = i >= ico ? Yout : Yin
        f[i, 1] = Y[i, 1] + im * Y[i, 2]; f[i, 2] = Y[i, 3] + im * Y[i, 4]
    end
    f[ico, 1] = bc[1]; f[ico, 2] = bc[2]
    f
end

"Coefficients a₁,a₂ making a₁yh₁ + a₂yh₂ + yih satisfy δh′ = C δh at both boundaries."
function homogeneous_coefficients(yh1, yh2, yih, Ci, Co, i_in, i_out)
    m00 = yh1[i_in, 2] - Ci * yh1[i_in, 1];  m01 = yh2[i_in, 2] - Ci * yh2[i_in, 1]
    m10 = yh1[i_out, 2] - Co * yh1[i_out, 1]; m11 = yh2[i_out, 2] - Co * yh2[i_out, 1]
    det = m00 * m11 - m01 * m10
    v0 = Ci * yih[i_in, 1] - yih[i_in, 2]; v1 = Co * yih[i_out, 1] - yih[i_out, 2]
    a1 = (m11 * v0 - m01 * v1) / det; a2 = (-m10 * v0 + m00 * v1) / det
    a1, a2
end

"Adaptive boundary indices: where |δh| first drops below `thresh·|δh(r=1)|` on either side."
function mode_boundaries(f, i1, thresh = 1e-4)
    mx = abs(f[i1, 1]); n = size(f, 1); i_in = 1; i_out = n
    for i in i1:-1:1
        abs(f[i, 1]) < thresh * mx && (i_in = i; break)
    end
    for i in i1:n
        abs(f[i, 1]) < thresh * mx && (i_out = i; break)
    end
    i_in, i_out
end

"WKB outgoing-wave coefficients at the two boundaries (with the m = 1 inner correction of main.py)."
function wkb_boundary_coeffs(d::NewtonianDisc, x, i_in, i_out, m)
    Ci = wkbC(d, x[i_in], m)
    m == 1 && (Ci -= 2im * kwave(d, x[i_in], m))
    Co = wkbC(d, x[i_out], m)
    Ci, Co
end

"""
    solve_m_mode(d, m, x; reint_steps = 2, phgrad_steps = 0, diagnostics = false, rtol = 1e-10, atol = 1e-12) -> complex (n, 2)

Enthalpy perturbation `(δh, δh′)` of azimuthal number `m` on grid `x` by the MR19 method: two
homogeneous plus one inhomogeneous integration from corotation, WKB outgoing-wave conditions at
adaptive boundaries, and `reint_steps` inhomogeneous re-integration refinements. Values outside the
converged boundaries are zeroed. `phgrad_steps=1` applies one MR secant correction; the default
zero preserves the original solver. `diagnostics=true` returns a named tuple with mode, before,
boundary constants, boundary errors, and relative changes. With a phase step, `before` is the
pre-boundary-update baseline, matching the mirrored algorithm. `boundary_coeffs` and
`phase_scale_height` select the published boundary and scale-height convention.
"""
function solve_m_mode(d::NewtonianDisc, m, x; rtol = 1e-10, atol = 1e-12, kwargs...)
    _solve_m_mode(d, m, x; rtol, atol, kwargs...)
end

function _solve_m_mode(d::NewtonianDisc{T}, m, x, cache...; reint_steps = 2,
                       phgrad_steps = 0, diagnostics = false, boundary_coeffs = wkb_boundary_coeffs,
                       phase_scale_height = _phase_scale_height, rtol, atol) where {T}
    phgrad_steps in (0,1) || throw(ArgumentError("phgrad_steps must be 0 or 1 (one MR secant step)"))
    n = length(x); i1 = argmin(abs.(x .- 1))
    fh1 = integrate_master(d, x, (one(Complex{T}), zero(Complex{T})), m, 0, cache...; rtol = rtol, atol = atol)
    fh2 = integrate_master(d, x, (zero(Complex{T}), one(Complex{T})), m, 0, cache...; rtol = rtol, atol = atol)
    fih = integrate_master(d, x, (zero(Complex{T}), zero(Complex{T})), m, 1, cache...; rtol = rtol, atol = atol)
    i_in = 1; i_out = n
    Ci, Co = boundary_coeffs(d, x, i_in, i_out, m)
    a1, a2 = homogeneous_coefficients(fh1, fh2, fih, Ci, Co, i_in, i_out)
    fm = a1 .* fh1 .+ a2 .* fh2 .+ fih
    i_in, i_out = mode_boundaries(fm, i1)
    Ci, Co = boundary_coeffs(d, x, i_in, i_out, m)
    baseline = fm
    for _ in 1:reint_steps
        ico = corotation_index(d, x, m)
        fih = integrate_master(d, x, (fm[ico, 1], fm[ico, 2]), m, 1, cache...; rtol = rtol, atol = atol)
        a1, a2 = homogeneous_coefficients(fh1, fh2, fih, Ci, Co, i_in, i_out)
        fm = a1 .* fh1 .+ a2 .* fh2 .+ fih
        baseline = fm  # main.py:639 uses this pre-boundary-update solution.
        i_in, i_out = mode_boundaries(fm, i1)
        Ci, Co = boundary_coeffs(d, x, i_in, i_out, m)
        a1, a2 = homogeneous_coefficients(fh1, fh2, fih, Ci, Co, i_in, i_out)
        phgrad_steps == 0 && (fm = a1 .* fh1 .+ a2 .* fh2 .+ fih)
    end
    C_before = (Ci,Co)
    before = phgrad_steps == 1 ? baseline : fm
    if phgrad_steps == 1
        step = phase_gradient_refine(d,x,m,fh1,fh2,fih,Ci,Co,i_in,i_out;
                                     baseline,scale_height=phase_scale_height)
        fm = step.mode
        Ci,Co = step.C
    end
    fm[1:i_in-1, :] .= 0; fm[i_out+1:end, :] .= 0
    if diagnostics
        X_before = phase_error_diagnostic(d,before,x,m,i_in,i_out; scale_height=phase_scale_height)
        X_after = phase_error_diagnostic(d,fm,x,m,i_in,i_out; scale_height=phase_scale_height)
        relative_state_change = ntuple(j -> sqrt(sum(abs2,fm[:,j]-before[:,j])/sum(abs2,before[:,j])),2)
        return (; mode=fm,before,C_before,C_after=(Ci,Co),X_before,X_after,
                 indices=(i_in,i_out),relative_C_change=abs.(((Ci,Co).-C_before)./C_before),
                 relative_state_change)
    end
    fm
end

"""
    vars_from_enthalpy(d, f, x, m) -> (δu_r, δu_φ, δΣ)

Reconstruct the velocity and surface-density perturbations from the enthalpy mode `f`
(main.py `vars_from_enthalpy`). Complex vectors of `length(x)`.
"""
function vars_from_enthalpy(d::NewtonianDisc{T}, f, x, m) where {T}
    n = length(x); u_r = zeros(Complex{T}, n); u_phi = zeros(Complex{T}, n); sig = zeros(Complex{T}, n)
    @inbounds for i in 1:n
        r = x[i]; Om = Omega(d, r); kap = kappa(d, r)
        ot = (omega_pat(d, m) - m * Om) + im * d.pole_disp; t = tc(d, r); bt = ot * t
        N2 = Nr2(d, r); LSi = LSinv(d, r); LTi = LTinv(d, r)
        Leff = (LTi - im * d.gamma * bt * LSi) / (1 - im * d.gamma * bt)
        Dv = Dfreq(d, r, m); Dc = Dv - im * d.gamma * bt / (1 - im * d.gamma * bt) * N2
        Ph = Phi_m(d, r, m); dPh = dPhidr_m(d, r, m); cs = csadi(d, r)
        u_r[i] = im / Dc * (ot * (f[i, 2] + dPh) - 2 * m * Om / r * (f[i, 1] + Ph) - Leff * ot * f[i, 1])
        u_phi[i] = 1 / Dc * (kap^2 / (2 * Om) * (f[i, 2] + dPh) -
                             m * (kap^2 - Dc) / (r * ot) * (f[i, 1] + Ph) - Leff * kap^2 / (2 * Om) * f[i, 1])
        sig[i] = Sigma(d, r) / (1 / d.gamma - im * bt) * ((1 - im * bt) * f[i, 1] / cs^2 + LSi * t * u_r[i])
    end
    u_r, u_phi, sig
end

"""
    solve_modes(d, mmax, x; mmin = 1, kwargs...) -> Dict{Int,Matrix}

Solve `m = mmin…mmax` IN PARALLEL across Julia threads (the modes are independent; run with `-t N`
and BLAS pinned to one thread). Returns `m => (n, 2)` complex enthalpy modes.
"""
function solve_modes(d::NewtonianDisc{T}, mmax::Int, x; mmin::Int = 1, kwargs...) where {T}
    return solve_modes(d, mmin:mmax, x; kwargs...)
end

"""
    solve_modes(d, ms::AbstractVector, x; kwargs...) -> Dict{Int,Matrix}

Solve an arbitrary (not necessarily contiguous) list of azimuthal numbers in parallel.
"""
function solve_modes(d::NewtonianDisc{T}, ms::AbstractVector, x; reint_steps = 2, phgrad_steps = 0,
                     rtol = 1e-10, atol = 1e-12) where {T}
    ms = collect(ms)
    out = Vector{Matrix{Complex{T}}}(undef, length(ms))
    Threads.@threads for i in eachindex(ms)
        out[i] = solve_m_mode(d, ms[i], x; phgrad_steps, reint_steps = reint_steps, rtol = rtol, atol = atol)
    end
    Dict(ms[i] => out[i] for i in eachindex(ms))
end

# ---- cached-potential throughput path -------------------------------------------------------------
"""
    PotentialCache(d, m, rin, rout; n = 4000)

Precomputed Φ_m, Φ_m′, Φ_m″ on a log-spaced radial grid. Each evaluation of the softened potential is
an `lap_nphi`-point azimuthal quadrature — the per-right-hand-side bottleneck of the adaptive
integrator — so for parameter sweeps the potential is tabulated once and interpolated
(4-point local Lagrange via `NumericsCore.local_interp` on log₁₀ r; the potential is smooth). The
PHYSICS is unchanged: gate the cached solve against `solve_m_mode` before trusting a sweep.
"""
struct PotentialCache{T}
    lg::Vector{T}
    P::Vector{T}; dP::Vector{T}; d2P::Vector{T}
end
function PotentialCache(d::NewtonianDisc{T}, m, rin, rout; n = 4000) where {T}
    rg = collect(T(10) .^ range(log10(T(rin)), log10(T(rout)); length = n))
    PotentialCache{T}(log10.(rg), [Phi_m(d, r, m) for r in rg], [dPhidr_m(d, r, m) for r in rg],
                      [d2Phidr2_m(d, r, m) for r in rg])
end
@inline function potential_at(c::PotentialCache, r)
    lq = log10(r)
    (local_interp(lq, c.lg, c.P; npts = 4), local_interp(lq, c.lg, c.dP; npts = 4),
     local_interp(lq, c.lg, c.d2P; npts = 4))
end

"Master-equation right-hand side with the tabulated potential (see [`PotentialCache`](@ref))."
struct CachedMasterRHS{D,I,C}
    d::D
    m::I
    ih::I
    cache::C
end
function _cached_rhs(g::CachedMasterRHS, r, y)
    d = g.d
    Om = Omega(d, r); dOm = dOmegadr(d, r); Sig = Sigma(d, r); dSig = dSigmadr(d, r)
    cs = csadi(d, r); ot = omega_pat(d, g.m) - g.m * Om
    Dv = Dfreq(d, r, g.m); dD = dDdr(d, r, g.m)
    LSi = LSinv(d, r); dLSi = dLSinvdr(d, r); LTi = LTinv(d, r); dLTi = dLTinvdr(d, r)
    N2 = Nr2(d, r); dN2 = dNr2dr(d, r); t = tc(d, r); dt = dtcdr(d, r)
    if g.ih == 1
        Ph, dPh, d2Ph = potential_at(g.cache, r)
    else
        Ph = dPh = d2Ph = zero(r)
    end
    C1R, C1I, C0R, C0I, PsiR, PsiI = master_equation_cool(
        r, d.gamma, cs, Om, dOm, Sig, dSig, LTi, dLTi, LSi, dLSi, ot, t, dt, Dv, dD, N2, dN2,
        d.pole_disp, Ph, dPh, d2Ph, g.m)
    @inbounds (y[3], y[4],
        PsiR - C1R * y[3] + C1I * y[4] - C0R * y[1] + C0I * y[2],
        PsiI - C1R * y[4] - C1I * y[3] - C0R * y[2] - C0I * y[1])
end
@inline function (g::CachedMasterRHS)(out::AbstractVector, r, y)
    d1, d2, d3, d4 = _cached_rhs(g, r, y)
    @inbounds (out[1] = d1; out[2] = d2; out[3] = d3; out[4] = d4)
    out
end
@inline (g::CachedMasterRHS)(r, y) = collect(_cached_rhs(g, r, y))

function integrate_master(d::NewtonianDisc{T}, x, bc, m, ih, cache::PotentialCache;
                          rtol = 1e-8, atol = 1e-10) where {T}
    ico = corotation_index(d, x, m)
    y0 = T[real(bc[1]), imag(bc[1]), real(bc[2]), imag(bc[2])]
    frhs = CachedMasterRHS(d, m, ih, cache)
    Yout = integrate_grid(frhs, y0, x, ico, +1; rtol = rtol, atol = atol)
    Yin  = integrate_grid(frhs, y0, x, ico, -1; rtol = rtol, atol = atol)
    n = length(x); f = zeros(Complex{T}, n, 2)
    @inbounds for i in 1:n
        Y = i >= ico ? Yout : Yin
        f[i, 1] = Y[i, 1] + im * Y[i, 2]; f[i, 2] = Y[i, 3] + im * Y[i, 4]
    end
    f[ico, 1] = bc[1]; f[ico, 2] = bc[2]
    f
end

"""
    solve_m_mode(d, m, x, cache::PotentialCache; reint_steps = 2, rtol = 1e-8, atol = 1e-10)

Same algorithm as [`solve_m_mode`](@ref) with the tabulated forcing potential (throughput path for
parameter sweeps over many `m`).
"""
function solve_m_mode(d::NewtonianDisc, m, x, cache::PotentialCache;
                      rtol = 1e-8, atol = 1e-10, kwargs...)
    _solve_m_mode(d, m, x, cache; rtol, atol, kwargs...)
end

# ---- homogeneous (passive-wake) problem, MR19 §passive / Figs. 10-12 ------------------------------
# Pattern frequency ω_p = 0 (no secondary, no resonance inside the domain) and Φ_m ≡ 0: the wave is
# launched at the outer edge and propagates inward. Ported from the legacy figure-local
# `homog_solver.jl`; the ω_p-dependence is the ONLY difference from the forced problem above, so the
# ω_p = 0 forms are given their own thin functions rather than a second copy of the model.
passive_Dfreq(d::NewtonianDisc, s, m) = kappa(d, s)^2 - (m * Omega(d, s))^2
passive_dDdr(d::NewtonianDisc, s, m)  = 2 * kappa(d, s) * dkappadr(d, s) -
                                        2 * m^2 * Omega(d, s) * dOmegadr(d, s)

"WKB radial wavenumber of the passive (ω_p = 0) problem."
function passive_kwave(d::NewtonianDisc, s, m)
    Om = Omega(d, s); ot = -m * Om; t = tc(d, s)
    cseff = ((1 + (d.gamma * ot * t)^2) / (1 + (ot * t)^2))^(one(s) / 4) * csiso(d, s)
    Deff = passive_Dfreq(d, s, m) + (d.gamma * ot * t)^2 / (1 + (d.gamma * ot * t)^2) * Nr2(d, s)
    sqrt(complex(-Deff)) / cseff
end

"Outgoing-wave boundary coefficient of the passive (ω_p = 0) problem."
function passive_wkbC(d::NewtonianDisc, s, m)
    Om = Omega(d, s); ot = -m * Om; t = tc(d, s)
    rek = passive_kwave(d, s, m)
    imk = (d.gamma - 1) / 2 * ot * t / (1 + d.gamma * (ot * t)^2) * rek
    k = rek + im * imk
    w = im * k
    w += passive_dDdr(d, s, m) / (4 * passive_Dfreq(d, s, m)) - dSigmadr(d, s) / (2 * Sigma(d, s)) +
         dcsisodr(d, s) / (2 * csiso(d, s)) - 1 / (2 * s)
    w += (1 + im * d.gamma * ot * t) / 2 / (1 + (d.gamma * ot * t)^2) * LTinv(d, s)
    w
end

struct PassiveMasterRHS{D,I}
    d::D
    m::I
end
function _passive_rhs(g::PassiveMasterRHS, r, y)
    d = g.d; m = g.m
    Om = Omega(d, r); dOm = dOmegadr(d, r); Sig = Sigma(d, r); dSig = dSigmadr(d, r)
    cs = csadi(d, r); ot = -m * Om
    Dv = passive_Dfreq(d, r, m); dD = passive_dDdr(d, r, m)
    LSi = LSinv(d, r); dLSi = dLSinvdr(d, r); LTi = LTinv(d, r); dLTi = dLTinvdr(d, r)
    N2 = Nr2(d, r); dN2 = dNr2dr(d, r); t = tc(d, r); dt = dtcdr(d, r)
    C1R, C1I, C0R, C0I, _, _ = master_equation_cool(
        r, d.gamma, cs, Om, dOm, Sig, dSig, LTi, dLTi, LSi, dLSi, ot, t, dt, Dv, dD, N2, dN2,
        d.pole_disp, zero(r), zero(r), zero(r), m)
    @inbounds (y[3], y[4],
        -C1R * y[3] + C1I * y[4] - C0R * y[1] + C0I * y[2],
        -C1R * y[4] - C1I * y[3] - C0R * y[2] - C0I * y[1])
end
@inline function (g::PassiveMasterRHS)(out::AbstractVector, r, y)
    d1, d2, d3, d4 = _passive_rhs(g, r, y)
    @inbounds (out[1] = d1; out[2] = d2; out[3] = d3; out[4] = d4)
    out
end
@inline (g::PassiveMasterRHS)(r, y) = collect(_passive_rhs(g, r, y))

"""
    passive_mode(d, m, x; rtol = 1e-9, atol = 1e-30) -> complex (n, 2)

Unique outgoing-wave homogeneous ("passive wake") mode, MR19 §passive (Figs. 10-12). The physical
inward-propagating wave GROWS toward small r by angular-momentum-flux conservation, so the
well-conditioned direction is INWARD from r_out: the WKB outgoing-wave relation δh′ = C_o δh is
imposed at `x[end]` and the equation integrated to `x[1]`, then normalised to δh(r_out) = 1 (all
modes launched in phase at the outer edge).
"""
function passive_mode(d::NewtonianDisc{T}, m, x; rtol = 1e-9, atol = 1e-30) where {T}
    n = length(x)
    Co = passive_wkbC(d, x[n], m)
    dh0 = one(Complex{T}); dhp0 = Co * dh0
    y0 = T[real(dh0), imag(dh0), real(dhp0), imag(dhp0)]
    Y = integrate_grid(PassiveMasterRHS(d, m), y0, x, n, -1; rtol = rtol, atol = atol)
    f = zeros(Complex{T}, n, 2)
    @inbounds for i in 1:n
        f[i, 1] = Y[i, 1] + im * Y[i, 2]; f[i, 2] = Y[i, 3] + im * Y[i, 4]
    end
    f[n, 1] = dh0; f[n, 2] = dhp0
    f ./ f[n, 1]
end

"""
    passive_sigma(d, f, x, m) -> δΣ

Surface-density perturbation of a passive (ω_p = 0, Φ = 0) mode `f` — the ω_p = 0 specialisation of
[`vars_from_enthalpy`](@ref).
"""
function passive_sigma(d::NewtonianDisc{T}, f, x, m) where {T}
    n = length(x); sig = zeros(Complex{T}, n)
    @inbounds for i in 1:n
        r = x[i]; Om = Omega(d, r)
        ot = (-m * Om) + im * d.pole_disp; t = tc(d, r); bt = ot * t
        N2 = Nr2(d, r); LSi = LSinv(d, r); LTi = LTinv(d, r)
        Leff = (LTi - im * d.gamma * bt * LSi) / (1 - im * d.gamma * bt)
        Dv = passive_Dfreq(d, r, m); Dc = Dv - im * d.gamma * bt / (1 - im * d.gamma * bt) * N2
        cs = csadi(d, r)
        u_r = im / Dc * (ot * f[i, 2] - 2 * m * Om / r * f[i, 1] - Leff * ot * f[i, 1])
        sig[i] = Sigma(d, r) / (1 / d.gamma - im * bt) * ((1 - im * bt) * f[i, 1] / cs^2 + LSi * t * u_r)
    end
    sig
end
