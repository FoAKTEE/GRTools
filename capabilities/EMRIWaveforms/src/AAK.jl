"""
    EMRIWaveforms.AAK

Augmented-analytic-Kerr (AAK) waveform pieces (design §6.8: FEW `summation/aaksummation.jl`; KerrBackwards
`{summation,bessel_ak,parametermap_aak}.jl`):

- `ParameterMapAAK` — the AK PN coordinate-rate polynomials (`drdm/dtdm/dthetadm/dphidm`), the AK azimuthal
  frequency `OmegaPhi`, the mapped fundamental frequencies `_ak_mapped_freqs`, and the Kerr→AK `par_map`
  3×3 frequency-matching solve (NLsolve).
- `ak_bessel_coeffs` — the Peters–Mathews Bessel basis, canonical home HERE (hoisted from Kernels per
  the §5.5/§7.7.1 purity gate: Kernels carries no SpecialFunctions machinery).
- `aak_summation` — thin delegation to `Kernels.aak_summation` injecting `ak_bessel_coeffs` as the
  kernel's `bessel_coeffs` provider, with a `backend::Symbol=:few` no-op selector (the FEW and KB
  legacy forks were byte-identical delegations, now one API).
- `aak_waveform` — top-level trajectory-driven convenience: runs `Trajectory.run_inspiral`, builds the AK
  basis (`λ=arccos Y`, `Φ̄=Φ_r`, `γ=Φ_θ−Φ_r`, `α=Φ_φ−Φ_θ`, `ν=Ω_r/2π`, `γ̇=Ω_θ−Ω_r`) and reduces it
  through `aak_summation`. Generic-inclination (`Y0≠1`) is supported: `Trajectory.pn5_rhs` drives
  the direct `dpdt/dedt/dYdt8H_5PNe10` fluxes for inclined orbits.
"""
module AAK

import Kernels
using SpecialFunctions: besselj
using NumericsCore: MTSUN_SI
using ..Trajectory: EMRIInspiral, PN5FluxODE, run_inspiral, _symmetric_mass_ratio
using Spacetime: kerr_frequencies
using ..Mappings: xI_to_Y
export ak_bessel_coeffs, aak_summation, aak_waveform

"""
    ak_bessel_coeffs(n, e) -> (J0, J1, J2, J3, J4)

The five analytic-kludge Bessel coefficients `(J_{n-2}, J_{n-1}, J_n, J_{n+1}, J_{n+2})` at argument
`ne = n·e` (Peters–Mathews / Barack–Cutler harmonic amplitudes). For `n = 1` the order `n-2 = -1` is taken
via `J_{-1}(x) = -J_1(x)` (DLMF 10.4.1), so `J0 = -J_1(e)`. Upstream: `gpuAAK.cu` AK amplitude loop.

Canonical home (S4-c §5.5 purity hoist, moved verbatim from `Kernels.aak`): Kernels' `aak_summation`
consumes this as its injected `bessel_coeffs` provider and carries no SpecialFunctions dep itself.
"""
function ak_bessel_coeffs(n::Integer, e::Real)
    n ≥ 1 || throw(ArgumentError("ak_bessel_coeffs: harmonic n must be ≥ 1, got $n"))
    ne = n * e
    J0 = n == 1 ? -besselj(1, ne) : besselj(n - 2, ne)
    J1 = besselj(n - 1, ne)
    J2 = besselj(n,     ne)
    J3 = besselj(n + 1, ne)
    J4 = besselj(n + 2, ne)
    return (J0, J1, J2, J3, J4)
end

# =====================================================================================
# AK PN coordinate-rate polynomials + frequency map (KerrBackwards parametermap_aak.jl).
# =====================================================================================
module ParameterMapAAK

import NLsolve

export drdm, dtdm, dthetadm, dphidm, OmegaPhi, par_map

"`drdm(v,e,Y,q)` — AK PN radial-rate polynomial (`ParameterMapAAK.cc:30`)."
function drdm(v, e, Y, q)
    v2 = v*v; v3 = v2*v
    e2 = e*e; e4 = e2*e2
    Y2 = Y*Y; q2 = q*q
    return (16. + 8. * (-3. + e2) * v2 - 16. * (-3. + e2) * q * v3 * Y +
            8. * (33. + 4. * e2 - 3. * e4) * q * v3 * v2 * Y +
            v3 * v3 * (-351. + 132. * q2 + e2 * (-135. + 21. * e2 + 5. * e4 + 2. * (7. + e2) * q2) +
                       2. * (-204. + 13. * e2 * (-3. + e2)) * q2 * Y2) +
            2. * v2 * v2 * (-45. + 3. * e4 + 4. * q2 * (1. - 4. * Y2) + 2. * e2 * q2 * (1. + Y2))) / (16. * v)
end

"`dtdm(v,e,Y,q)` — AK PN coordinate-time-rate polynomial (`ParameterMapAAK.cc:46`)."
function dtdm(v, e, Y, q)
    v2 = v*v; v3 = v2*v; v4 = v2*v2
    e2 = e*e; e4 = e2*e2
    Y2 = Y*Y; q2 = q*q
    return (2. * e4 * (240. + v2 * (-120. + v * (42. * (-27. + 4. * q2) * v + (-6567. + 1996. * q2) * v3 +
                48. * q * (8. + 77. * v2) * Y - 4. * q2 * v * (90. + 1577. * v2) * Y2))) +
            e4 * e2 * (560. + v2 * (-360. + 960. * q * v * Y + 8816. * q * v3 * Y +
                v4 * (-15565. + 24. * q2 * (200. - 629. * Y2)) +
                v2 * (-2742. + 80. * q2 * (5. - 11. * Y2)))) -
            8. * e2 * (-48. + v2 * (8. - 64. * q * v * Y - 688. * q * v3 * Y +
                2. * v2 * (99. + 16. * q2 * (-1. + 2. * Y2)) + v4 * (1233. + 8. * q2 * (-47. + 150. * Y2))))
            + 16. * (16. + v2 * (24. + v2 * (27. * (2. + 5. * v2) - 48. * q * v * Y +
                4. * q2 * (2. - v2 + (-2. + 3. * v2) * Y2))))) / (256. * v4)
end

"`dthetadm(v,e,Y,q)` — AK PN polar-rate polynomial (`ParameterMapAAK.cc:67`)."
function dthetadm(v, e, Y, q)
    v2 = v*v; v3 = v2*v
    e2 = e*e; e4 = e2*e2
    Y2 = Y*Y; q2 = q*q
    return (16. + 8. * (3. + e2) * v2 - 16. * (3. + e2) * q * v3 * Y -
            8. * (3. + e2) * (5. + 3. * e2) * q * v3 * v2 * Y +
            v3 * v3 * (135. - 54. * q2 + e2 * (5. * (27. + 9. * e2 + e4) + 2. * (-38. + e2) * q2) +
                       2. * (57. + 90. * e2 + 13. * e4) * q2 * Y2) +
            2. * v2 * v2 * (27. + 3. * e4 + 2. * q2 * (-1. + 7. * Y2) + 2. * e2 * (9. + q2 * (1. + Y2)))) / (16. * v)
end

"`dphidm(v,e,Y,q)` — AK PN azimuthal-rate polynomial (`ParameterMapAAK.cc:84`)."
function dphidm(v, e, Y, q)
    v2 = v*v; v3 = v2*v
    e2 = e*e; e4 = e2*e2
    Y2 = Y*Y; q2 = q*q
    return (16. + 8. * (3. + e2) * v2 - 16. * q * v3 * (-2. + (3. + e2) * Y) -
            8. * q * v3 * v2 * (-6. + 15. * Y + 3. * e4 * Y + 2. * e2 * (-4. + 7. * Y)) +
            2. * v2 * v2 * (27. + 3. * e4 + 2. * q2 * (-1. + Y) * (1. + 7. * Y) + 2. * e2 * (9. + q2 * (1. + Y2))) +
            v3 * v3 * (5. * e4 * e2 + e4 * (45. + q2 * (2. + 26. * Y2)) +
                       e2 * (135. + 4. * q2 * (-19. + 5. * Y * (-7. + 9. * Y))) +
                       3. * (45. + 2. * q2 * (-9. + Y * (-6. + 19. * Y))))) / (16. * v)
end

"""
    OmegaPhi(v, e, cosiota, s, M) -> Float64

AK azimuthal angular-frequency magnitude `dphidm/dtdm/M` with the prograde (`cosiota>0`) / retrograde
(`cosiota≤0`) sign flip on `(cosiota, s)`. As `v→0` (`p→∞`) it tends to the Newtonian Kepler limit
`v³/M`. (`ParameterMapAAK.cc:100`.)
"""
function OmegaPhi(v, e, cosiota, s, M)
    if cosiota > 0
        return dphidm(v, e, cosiota, s) / dtdm(v, e, cosiota, s) / M
    else
        return dphidm(v, e, -cosiota, -s) / dtdm(v, e, -cosiota, -s) / M
    end
end

"""
    _ak_mapped_freqs(v, M_map, s, e, iota, M) -> (Ω_r, Ω_θ, Ω_φ)

AK fundamental frequencies `d{r,θ,φ}dm / dtdm / (2π)` with the prograde/retrograde sign flip on
`(cos ι, s)` and the extra sign on `Ω_φ` for retrograde. (`ParameterMapAAK.cc::sol_fun`.)
"""
function _ak_mapped_freqs(v, M_map, s, e, iota, M)
    c = cos(iota)
    if c > 0
        Or = drdm(v, e, c, s) / dtdm(v, e, c, s) / (2π)
        Ot = dthetadm(v, e, c, s) / dtdm(v, e, c, s) / (2π)
        Op = dphidm(v, e, c, s) / dtdm(v, e, c, s) / (2π)
    else
        Or = drdm(v, e, -c, -s) / dtdm(v, e, -c, -s) / (2π)
        Ot = dthetadm(v, e, -c, -s) / dtdm(v, e, -c, -s) / (2π)
        Op = -dphidm(v, e, -c, -s) / dtdm(v, e, -c, -s) / (2π)
    end
    return (Or, Ot, Op)
end

"""
    par_map(Omega_r, Omega_theta, Omega_phi, p, M, s, e, iota; ftol=1e-12) -> (; v, M_map, s_map, converged)

Kerr → AK parameter map (`ParameterMapAAK.cc::ParMap`): solve `Ω_i·M_map − Ω_i^AK(v,s_map)·M = 0`,
`i ∈ {r,θ,φ}`, for `(v, M_map, s_map)` from `(1/√p, M, s)` (NLsolve, forward-mode autodiff Jacobian).
"""
function par_map(Omega_r, Omega_theta, Omega_phi, p, M, s, e, iota; ftol::Real = 1e-12)
    function residual!(F, X)
        v, M_map, s_map = X[1], X[2], X[3]
        Or, Ot, Op = _ak_mapped_freqs(v, M_map, s_map, e, iota, M)
        F[1] = Omega_r * M_map - Or * M
        F[2] = Omega_theta * M_map - Ot * M
        F[3] = Omega_phi * M_map - Op * M
        return F
    end
    sol = NLsolve.nlsolve(residual!, [1 / sqrt(p), float(M), float(s)]; ftol = ftol, autodiff = :forward)
    return (v = sol.zero[1], M_map = sol.zero[2], s_map = sol.zero[3], converged = sol.f_converged)
end

end # module ParameterMapAAK

# =====================================================================================
# Unified aak_summation → Kernels.aak_summation (backend selector is a no-op; the FEW and KB legacy
# forks were byte-identical delegations, now one module-scope API).
# =====================================================================================

"""
    aak_summation(t, e, Phivec, gimvec, alpvec, nuvec, gimdotvec, OmegaPhi, lam,
                  M, a, mu, qS, phiS, qK, phiK, dist, nmodes;
                  backend=:few, mich=false, dt=10.0, out_len=nothing, threaded=false)

AAK harmonic summation (`few.summation.aaksummation.AAKSummation.sum` ∪ KerrBackwards
`summation.aak_summation`, one canonical form). Delegates to `Kernels.aak_summation`, injecting
[`ak_bessel_coeffs`](@ref) as the kernel's required `bessel_coeffs` provider (§5.5 purity: the Bessel
basis lives here, not in Kernels); `backend::Symbol` (`:few` default, `:kb`) is a no-op selector
(both legacy forks were byte-identical delegations).
"""
function aak_summation(args...; backend::Symbol = :few, kwargs...)
    backend in (:few, :kb) ||
        throw(ArgumentError("backend must be :few or :kb, got $backend"))
    return Kernels.aak_summation(args...; bessel_coeffs = ak_bessel_coeffs, kwargs...)
end

"""
    _aak_lam_fixup!(lam; fill_val=1e-6) -> lam

AAK polar-angle fixup (`aakwave.py`): nudge `λ = ι` away from `0, π/2, π` by `fill_val` so the
analytic-kludge geometric projection stays regular.
"""
# Phase7-D dedup: canonical body lives in Kernels (core/Kernels/src/aak.jl); private symbol, aliased.
const _aak_lam_fixup! = Kernels._aak_lam_fixup!

"""
    aak_waveform(M, mu, a, p0, e0, Y0, qS, phiS, qK, phiK, dist;
                 T=1.0, dt=10.0, mich=false, nmodes=nothing,
                 Phi_phi0=0.0, Phi_theta0=0.0, Phi_r0=0.0, err=1e-10) -> Vector{ComplexF64}

Augmented-analytic-kludge (AAK) top-level waveform with the 5PN trajectory (`Pn5AAKWaveform`). Runs
[`run_inspiral`](@ref) to a sparse `(t,p,e,xI,Φ_φ,Φ_θ,Φ_r)`, builds the analytic-kludge basis exactly as
`AAKSummation.sum` — `λ = arccos Y` (regularised by [`_aak_lam_fixup!`]), `Φ̄ = Φ_r`, `γ = Φ_θ − Φ_r`,
`α = Φ_φ − Φ_θ`, `ν = Ω_r/2π`, `γ̇ = Ω_θ − Ω_r`, with the `L`-signed polar frequency, frequencies from
[`kerr_frequencies`](@ref) divided by `M·MTSUN_SI` — and reduces it through the grounded
[`aak_summation`](@ref) kernel.

`Y0` is the initial `cos ι` (generic inclination `Y0 ≠ 1` supported via the direct 5PN fluxes in
[`Trajectory`](@ref)). `nmodes` defaults to `max(4, ⌊30 e0⌋)`. `mich=false` returns `h₊ − i h₍×₎`.
"""
function aak_waveform(M, mu, a, p0, e0, Y0, qS, phiS, qK, phiK, dist;
                      T = 1.0, dt = 10.0, mich::Bool = false,
                      nmodes::Union{Integer,Nothing} = nothing,
                      Phi_phi0 = 0.0, Phi_theta0 = 0.0, Phi_r0 = 0.0, err = 1e-10)
    (isfinite(T) && T > 0) || throw(ArgumentError("T must be finite and positive, got $T"))
    insp = EMRIInspiral(PN5FluxODE(float(a)); err = err)
    massratio = _symmetric_mass_ratio(M, mu)
    y0 = [float(p0), float(e0), float(Y0),
          massratio * float(Phi_phi0), massratio * float(Phi_theta0),
          massratio * float(Phi_r0)]
    tr = run_inspiral(insp, float(M), float(mu), float(a), y0;
                      T = T, dt = dt)
    nm = nmodes === nothing ? max(4, floor(Int, 30 * e0)) : Int(nmodes)

    t = tr.t; p = tr.p; e = tr.e; xI = tr.x
    Phi_phi = tr.Phi_phi; Phi_theta = tr.Phi_theta; Phi_r = tr.Phi_r
    n = length(t)
    Msec = M * MTSUN_SI

    OmegaPhi = Vector{Float64}(undef, n)
    OmegaR = Vector{Float64}(undef, n)
    OmegaTheta = Vector{Float64}(undef, n)
    Yv = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        f = kerr_frequencies(a, p[i], e[i], xI[i])
        OmegaPhi[i] = f.Omegaph / Msec
        OmegaTheta[i] = copysign(f.Omegath, xI[i]) / Msec
        OmegaR[i] = f.Omegar / Msec
        Yv[i] = xI_to_Y(a, p[i], e[i], xI[i])
    end

    lam = acos.(Yv)
    _aak_lam_fixup!(lam)
    Phivec = Phi_r
    gimvec = Phi_theta .- Phi_r
    alpvec = Phi_phi .- Phi_theta
    nuvec = OmegaR ./ (2 * π)
    gimdotvec = OmegaTheta .- OmegaR

    return aak_summation(t, e, Phivec, gimvec, alpvec, nuvec, gimdotvec, OmegaPhi, lam,
                         M, a, mu, qS, phiS, qK, phiK, dist, nm; mich = mich, dt = dt)
end

end # module AAK
