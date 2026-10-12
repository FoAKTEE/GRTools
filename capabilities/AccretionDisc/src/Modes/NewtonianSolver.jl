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

# Direct forcing remains the circular Laplace harmonic, even when a different
# frequency is supplied for an operator test. Eccentric production supplies a cache.
_forcing_at(d, r, m, ::Nothing) = (Phi_m(d,r,m), dPhidr_m(d,r,m), d2Phidr2_m(d,r,m))
_forcing_value_derivative(d, r, m, ::Nothing) = (Phi_m(d,r,m), dPhidr_m(d,r,m))

function _mode_anchor(d, x, mode, anchor)
    if anchor === :legacy || (anchor === :auto && omega_pat(d,mode) == _azimuthal_number(mode))
        return argmin(abs.(x .- 1))
    elseif anchor === :corotation || anchor === :auto
        return corotation_index(d,x,mode)
    end
    throw(ArgumentError("anchor must be :auto, :legacy or :corotation"))
end

"Log-spaced radial grid r_in … r_out with `n` points."
loggrid(d::NewtonianDisc{T}, r_in, r_out, n) where {T} =
    T(10) .^ range(log10(T(r_in)), log10(T(r_out)); length = n)

"Corotation index: argmin |mΩ(r) − ω| for an integer m or a HarmonicMode."
function corotation_index(d::NewtonianDisc, x, m)
    best = 1; bv = abs(_azimuthal_number(m) * Omega(d, x[1]) - omega_pat(d, m))
    @inbounds for i in 2:length(x)
        v = abs(_azimuthal_number(m) * Omega(d, x[i]) - omega_pat(d, m)); v < bv && (bv = v; best = i)
    end
    best
end

"""
    master_rhs(d, r, y, mode, ih[, cache]) -> NTuple{4}

Master-equation right-hand side as a real 4-vector `y = [Re δh, Im δh, Re δh′, Im δh′]`.
`ih = 1` includes the secondary's forcing potential, `ih = 0` gives the homogeneous problem.
`mode` is an integer circular m or a `HarmonicMode`. The optional cache supplies
Phi, Phi' and Phi''; without it the forcing is the circular Laplace harmonic.
"""
function master_rhs(d::NewtonianDisc, r, y, mode, ih, cache=nothing)
    m = _azimuthal_number(mode)
    Om = Omega(d, r); dOm = dOmegadr(d, r); Sig = Sigma(d, r); dSig = dSigmadr(d, r)
    cs = csadi(d, r); ot = omega_pat(d, mode) - m * Om
    Dv = Dfreq(d, r, mode); dD = dDdr(d, r, mode)
    LSi = LSinv(d, r); dLSi = dLSinvdr(d, r); LTi = LTinv(d, r); dLTi = dLTinvdr(d, r)
    N2 = Nr2(d, r); dN2 = dNr2dr(d, r); t = tc(d, r); dt = dtcdr(d, r)
    if ih == 1
        Ph, dPh, d2Ph = _forcing_at(d, r, m, cache)
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
struct MasterRHS{D,M,I}
    d::D
    m::M
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
function integrate_master(d::NewtonianDisc{T}, x, bc, m, ih; rtol = 1e-10, atol = 1e-12, dense::Bool = true) where {T}
    _integrate_master(d, x, bc, m, MasterRHS(d, m, ih); rtol, atol, dense)
end

function _integrate_master(d::NewtonianDisc{T}, x, bc, m, frhs; rtol, atol, dense = true) where {T}
    ico = corotation_index(d, x, m)
    y0 = T[real(bc[1]), imag(bc[1]), real(bc[2]), imag(bc[2])]
    Yout = integrate_grid(frhs, y0, x, ico, +1; rtol = rtol, atol = atol, dense)
    Yin  = integrate_grid(frhs, y0, x, ico, -1; rtol = rtol, atol = atol, dense)
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

"Adaptive boundary indices: first drop below `thresh*|δh(x[i1])|` on either side."
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
    _azimuthal_number(m) == 1 && (Ci -= 2im * kwave(d, x[i_in], m))
    Co = wkbC(d, x[i_out], m)
    Ci, Co
end

"""
    group_velocity_boundary_coeffs(d, x, i_in, i_out, mode)

Locally isothermal (`beta=0`) radiation coefficients. In a wave region,
`v_g=c_s^2*k_r/(omega-m*Omega)`, so the inner/outer edge needs
`sign(k_r)=(-1/+1)*sign(omega-m*Omega)`. In an evanescent region select
decay into the exterior. Retains the existing WKB amplitude gradient;
corotation and exact turning points at an edge are rejected. Unlike the
legacy m=1 special case, selection depends on the actual local dispersion.
"""
function group_velocity_boundary_coeffs(d::NewtonianDisc, x, i_in, i_out, mode)
    iszero(d.beta) || throw(ArgumentError("group_velocity_bc requires beta=0"))
    map(((i_in,-1),(i_out,1))) do (i,normal)
        r=x[i]; w=omega_pat(d,mode)-_azimuthal_number(mode)*Omega(d,r)
        D=Dfreq(d,r,mode)
        iszero(w) || iszero(D) ? throw(DomainError(r,"move boundary away from corotation/turning point")) : nothing
        branch = D < 0 ? normal*sign(w) : normal
        wkbC(d,r,mode)+(branch-1)*im*kwave(d,r,mode)
    end
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
`phgrad_safeguard=true` applies the descent safeguard of
[`phase_gradient_refine`](@ref) with at most `phgrad_max_backtracks` halvings; the diagnostics
then also carry its audit as `phase_safeguard`.
This direct method forces with the circular Laplace harmonic Φ_m, so a `HarmonicMode` is accepted
here only for the circular planet (`k == 0`, `omega == m`); any other mode needs its own potential
through the cached method or [`solve_lm_mode`](@ref).
"""
function solve_m_mode(d::NewtonianDisc, m, x; rtol = 1e-10, atol = 1e-12, kwargs...)
    m isa HarmonicMode && !(m.k == 0 && m.omega == m.m) &&
        throw(ArgumentError("direct forcing is the circular Phi_m; pass this mode's PotentialCache"))
    _solve_m_mode(d, m, x; rtol, atol, kwargs...)
end

function _solve_m_mode(d::NewtonianDisc{T}, m, x, cache...; reint_steps = 2,
                       phgrad_steps = 0, diagnostics = false, boundary_coeffs = wkb_boundary_coeffs,
                       phase_scale_height = _phase_scale_height, anchor = :auto,
                       phgrad_safeguard::Bool = false, phgrad_max_backtracks::Integer = 6, dense::Bool = true, rtol, atol) where {T}
    phgrad_steps in (0,1) || throw(ArgumentError("phgrad_steps must be 0 or 1 (one MR secant step)"))
    n = length(x); i1 = _mode_anchor(d, x, m, anchor)
    fh1 = integrate_master(d, x, (one(Complex{T}), zero(Complex{T})), m, 0, cache...; rtol = rtol, atol = atol, dense)
    fh2 = integrate_master(d, x, (zero(Complex{T}), one(Complex{T})), m, 0, cache...; rtol = rtol, atol = atol, dense)
    fih = integrate_master(d, x, (zero(Complex{T}), zero(Complex{T})), m, 1, cache...; rtol = rtol, atol = atol, dense)
    i_in = 1; i_out = n
    Ci, Co = boundary_coeffs(d, x, i_in, i_out, m)
    a1, a2 = homogeneous_coefficients(fh1, fh2, fih, Ci, Co, i_in, i_out)
    fm = a1 .* fh1 .+ a2 .* fh2 .+ fih
    i_in, i_out = mode_boundaries(fm, i1)
    Ci, Co = boundary_coeffs(d, x, i_in, i_out, m)
    baseline = fm
    for _ in 1:reint_steps
        ico = corotation_index(d, x, m)
        fih = integrate_master(d, x, (fm[ico, 1], fm[ico, 2]), m, 1, cache...; rtol = rtol, atol = atol, dense)
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
                                     baseline,scale_height=phase_scale_height,
                                     safeguard=phgrad_safeguard,max_backtracks=phgrad_max_backtracks)
        fm = step.mode
        Ci,Co = step.C
    end
    fm[1:i_in-1, :] .= 0; fm[i_out+1:end, :] .= 0
    if diagnostics
        X_before = phase_error_diagnostic(d,before,x,m,i_in,i_out; scale_height=phase_scale_height)
        X_after = phase_error_diagnostic(d,fm,x,m,i_in,i_out; scale_height=phase_scale_height)
        relative_state_change = ntuple(j -> sqrt(sum(abs2,fm[:,j]-before[:,j])/sum(abs2,before[:,j])),2)
        result = (; mode=fm,before,C_before,C_after=(Ci,Co),X_before,X_after,
                 indices=(i_in,i_out),relative_C_change=abs.(((Ci,Co).-C_before)./C_before),
                 relative_state_change)
        return phgrad_safeguard && phgrad_steps == 1 ?
            (; result..., phase_safeguard=step.safeguard) : result
    end
    fm
end

"""
    vars_from_enthalpy(d, f, x, mode[, cache]) -> (δu_r, δu_φ, δΣ)

Reconstruct the velocity and surface-density perturbations from the enthalpy mode `f`
(main.py `vars_from_enthalpy`). Complex vectors of `length(x)`.
An integer mode preserves circular forcing; a `HarmonicMode` supplies the total
frequency. Pass the same potential cache used by the eccentric solve.
For a three-column `lindblad_regular=true` solve, also pass that keyword here:
the third column supplies radial mass flux, including at D=0, and the forcing
is evaluated with the solve's C2 interpolation.
"""
function vars_from_enthalpy(d::NewtonianDisc{T}, f, x, mode, cache=nothing;
                            lindblad_regular::Bool=false) where {T}
    lindblad_regular && return _regular_vars(d, f, x, mode, cache)
    m = _azimuthal_number(mode)
    n = length(x); u_r = zeros(Complex{T}, n); u_phi = zeros(Complex{T}, n); sig = zeros(Complex{T}, n)
    @inbounds for i in 1:n
        r = x[i]; Om = Omega(d, r); kap = kappa(d, r)
        ot = (omega_pat(d, mode) - m * Om) + im * d.pole_disp; t = tc(d, r); bt = ot * t
        N2 = Nr2(d, r); LSi = LSinv(d, r); LTi = LTinv(d, r)
        Leff = (LTi - im * d.gamma * bt * LSi) / (1 - im * d.gamma * bt)
        Dv = Dfreq(d, r, mode); Dc = Dv - im * d.gamma * bt / (1 - im * d.gamma * bt) * N2
        Ph, dPh = _forcing_value_derivative(d, r, m, cache); cs = csadi(d, r)
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

_forcing_at(d, r, m, cache::PotentialCache) = potential_at(cache,r)
_forcing_value_derivative(d, r, m, cache::PotentialCache) = potential_at(cache,r)[1:2]

"""
    solve_lm_mode(d, mode::HarmonicMode, x, cache::PotentialCache;
                  reint_steps=2, phgrad_steps=0, anchor=:auto, kwargs...)

Frequency-aware MR shooting solver. Construct `mode = NewtonianMode(m,l,n)` in
`AccretionDisc.ForcingModes`; its total frequency is `l*n`, in the disc's units
GM=1. Supply the forcing through any `PotentialCache`, including `harmonic_cache`.
The cache and grid must use those same length units; for a Kepler orbit n=a^(-3/2).

By default this corotation shooting method requires m>0 and a positive frequency
whose corotation lies in the grid. `group_velocity_bc=true` opts into signed
frequency radiation boundaries for beta=0 and admits omega<=0, still with m>0.
These modes have no corotation: the existing integration seed (minimum |omega-m
Omega|) is the OUTER grid node, and adaptive edge selection is anchored at the
node nearest r=1 (the reference semimajor axis in the FR units). Both homogeneous
solutions and the particular solution are shot inward from that outer seed.
At omega=0 the same regular Doppler frequency -m*Omega and radiation/decay
conditions apply; there is no division by the total frequency, and mode power
is exactly zero. The global m=1,l=0 response remains boundary dependent.
For omega<=0 the phase-gradient correction is disabled explicitly: its published
corotation-based refinement has no certificate in this sector. `diagnostics`
still exposes the selected boundaries. Positive-frequency phase settings are
unchanged. The opt-in fixes boundary selection, not exterior truncation error.

`:auto` anchors adaptive boundaries at corotation for non-circular frequencies,
and at r=1 when omega=m to retain circular bit identity. `:corotation` explicitly
selects corotation for all modes; `:legacy` explicitly selects r=1. Phase correction
is opt-in, as in `solve_m_mode`; use `phgrad_steps=1` for the FR refinement. Matching
options and a circular Phi_m cache give the identical circular solve. The phase-step descent
safeguard (`phgrad_safeguard=true`) uses outer descent with bounded inner enthalpy
change at m=1, and two-edge descent otherwise. FR22 consumers pass the measured settings
[`FR22_SOLVER_SETTINGS`](@ref) (`solve_lm_mode(d, mode, x, cache; FR22_SOLVER_SETTINGS...)`)
rather than their own. `dense=false` uses the existing grid-capped DP5 integrator;
`dense=true` retains the legacy cubic-Hermite dense output. Neither setting asserts
response convergence: refine tolerances, output grid and forcing cache separately.

`lindblad_regular=true` (locally isothermal `beta=0` only) integrates the regular
first-order enthalpy-plus-potential and radial mass-flux system through D=0. The
forcing uses a global C2 spline of the supplied cache's potential values; its
derivatives are taken from that same spline. This opt-in result has THREE columns:
`(δh, δh′, rΣδu_r)`. Pass `lindblad_regular=true` to `vars_from_enthalpy` to use the
retained flux at resonances without division by D. The legacy two-column result
and all circular defaults are unchanged. Grid and forcing convergence remain the
caller's responsibility; use `dense=false` to avoid cubic dense-output error.
"""
function solve_lm_mode(d::NewtonianDisc, mode::HarmonicMode, x, cache::PotentialCache;
                       lindblad_regular::Bool=false, group_velocity_bc::Bool=false, kwargs...)
    mode.m > 0 && (mode.omega > 0 || group_velocity_bc) ||
        throw(ArgumentError("corotation shooting requires m > 0 and omega > 0"))
    length(x) >= 2 && all(r -> isfinite(r) && r > 0,x) && all(diff(x) .> 0) ||
        throw(ArgumentError("the radial grid must be positive, finite and strictly increasing"))
    (group_velocity_bc && mode.omega <= 0) ||
        mode.m*Omega(d,first(x)) >= mode.omega >= mode.m*Omega(d,last(x)) ||
        throw(ArgumentError("corotation must lie within the radial grid"))
    length(cache.lg) >= 4 && length(cache.lg) == length(cache.P) == length(cache.dP) == length(cache.d2P) &&
        all(diff(cache.lg) .> 0) && first(cache.lg) <= log10(first(x)) && last(cache.lg) >= log10(last(x)) ||
        throw(ArgumentError("the potential cache must have at least four ordered nodes covering the grid"))
    options = (;kwargs...)
    if group_velocity_bc
        iszero(d.beta) || throw(ArgumentError("group_velocity_bc requires beta=0"))
        options = merge(options,(;boundary_coeffs=group_velocity_boundary_coeffs))
        if mode.omega <= 0
            options = merge(options,(;anchor=:legacy,phgrad_steps=0))
        end
    end
    if lindblad_regular
        iszero(d.beta) || throw(ArgumentError("lindblad_regular requires beta=0"))
        return _solve_m_mode(d,mode,x,RegularPotential(cache);
                             rtol=1e-8,atol=1e-10,options...)
    end
    solve_m_mode(d,mode,x,cache; options...)
end

"""
    FR22_SOLVER_SETTINGS

The FR22 keywords of [`solve_lm_mode`](@ref): two re-seeds, one safeguarded phase step anchored at
corotation, `rtol=1e-10, atol=1e-12`, grid-capped DP5 (`dense=false`) and the regular Lindblad
crossing (`lindblad_regular=true`, so `beta=0`). Measured (ECC::fr22-linear-solver, lane w16-g40)
only for the tested scope: 4000 log-spaced radii on [0.1, 5] a (native control [0.05, 5] a), forcing
caches of at least 3997 nodes and the twelve-mode ladder plus the native m = 40 control, where the
torque error is at most 1.8e-5 against a fixed-mesh eighth-order oracle. Other grids, caches or
domains need their own refinement.
"""
const FR22_SOLVER_SETTINGS = (; reint_steps=2, phgrad_steps=1, phgrad_safeguard=true, anchor=:corotation,
                               rtol=1e-10, atol=1e-12, dense=false, lindblad_regular=true)

"Master-equation right-hand side with the tabulated potential (see [`PotentialCache`](@ref))."
struct CachedMasterRHS{D,M,I,C}
    d::D
    m::M
    ih::I
    cache::C
end
_cached_rhs(g::CachedMasterRHS, r, y) = master_rhs(g.d, r, y, g.m, g.ih, g.cache)
@inline function (g::CachedMasterRHS)(out::AbstractVector, r, y)
    d1, d2, d3, d4 = _cached_rhs(g, r, y)
    @inbounds (out[1] = d1; out[2] = d2; out[3] = d3; out[4] = d4)
    out
end
@inline (g::CachedMasterRHS)(r, y) = collect(_cached_rhs(g, r, y))

function integrate_master(d::NewtonianDisc{T}, x, bc, m, ih, cache::PotentialCache;
                          rtol = 1e-8, atol = 1e-10, dense::Bool = true) where {T}
    _integrate_master(d, x, bc, m, CachedMasterRHS(d, m, ih, cache); rtol, atol, dense)
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
