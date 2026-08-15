"""
    EMRIWaveforms.Generate

Top-level observed-waveform generation + SSB source/detector frame algebra (design §6.8: FEW
`waveform/generate.jl`; KerrBackwards `waveform.jl`). Two pieces:

- The frame transform — `get_viewing_angles` (detector sky/spin → source viewing angles),
  `ssb_polarization_angle` (the `ψ_LDC` polarization angle), and `to_SSB_frame` (the `2ψ` polarization
  rotation of `(h₊, h_×)`). Closed-form and self-contained.
- `generate_emri_waveform` — the `GenerateEMRIWaveform` dispatcher: it maps `x0 → Y₀` (`Mappings.xI_to_Y`)
  when the generator needs `Y`, calls the per-class generator backend, applies the distance scaling, and —
  for `frame=:detector` — rotates the source-frame polarizations into the SSB frame.

# v0.1 wired paths
Two generators are wired: the source-frame `:FastSchwarzschildEccentricFlux` (composed from
[`fast_schwarzschild_eccentric_flux`] — `run_inspiral`(`SchwarzEccFluxODE`) → `RomanAmplitude` →
`ModeSelect` → `ylms_few` → interpolated mode-sum) and the detector-frame `:Pn5AAKWaveform` (the
trajectory-driven `AAK.aak_waveform`). The AAK backend is an injectable seam (`aak_generator`, default
`AAK.aak_waveform`) so the dispatch / frame composition / return shape stay unit-testable against a stub.
`SlowSchwarzschildEccentricFlux` remains `[FUTURE]`. Generic-inclination Kerr is supported on the AAK
path (the generic-Y 5PN flux is implemented in `Trajectory.pn5_rhs` via the KerrBackwards direct fluxes).
"""
module Generate

using ..Mappings: xI_to_Y
using NumericsCore: Gpc_in_m, MRSUN_SI
import ..Trajectory
using ..Trajectory: EMRIInspiral, run_inspiral
import ..FluxModels
import ..AmplitudeModels
import ..ModeSelect
import ..ModeSum
import ..AAK
import Data
using Harmonics: ylms_few

export get_viewing_angles, ssb_polarization_angle, to_SSB_frame, generate_emri_waveform,
       fast_schwarzschild_eccentric_flux

# =====================================================================================
# Source-frame Schwarzschild eccentric generator (KerrBackwards/fastschwarz_waveform.jl).
# =====================================================================================

"""
    SchwarzEccFluxODE(flux)

`Trajectory.ODEModel` adapter for the Schwarzschild eccentric flux model
(`FluxModels.KerrBackwards.SchwarzEccFlux`): `a = 0`, evolves `x` directly (no `Y` conversion), and its
`evaluate_rhs` returns the per-unit-mass-ratio flux RHS so `run_inspiral`'s mass-ratio adimensionalisation
applies unchanged (`epsilon = 1`).
"""
struct SchwarzEccFluxODE{F} <: Trajectory.ODEModel
    flux::F
end
Trajectory.ode_spin(::SchwarzEccFluxODE) = 0.0
Trajectory.converts_Y(::SchwarzEccFluxODE) = false
Trajectory.evaluate_rhs(m::SchwarzEccFluxODE, p, e, x, p_sep) =
    collect(FluxModels.KerrBackwards.schwarz_ecc_flux_deriv(m.flux, 1.0, 0.0, p, e, 1.0))

"""
    fast_schwarzschild_eccentric_flux(M, mu, p0, e0, theta, phi; T=1.0, dt=10.0, eps=1e-5,
                                      normalize_amps=true, mode_selection=:power,
                                      roman=nothing) -> (t, h)

Time-domain Schwarzschild eccentric EMRI strain `h(t) = h₊ − i h_×` — the `FastSchwarzschildEccentricFlux`
generator (`few/waveform.py`). Composes the already-grounded GrTools pieces: [`run_inspiral`](@ref) with
[`SchwarzEccFluxODE`] (trajectory) → `AmplitudeModels.RomanAmplitude` Teukolsky amplitudes
(`get_amplitudes(...; backend=:few)` carries the amplitude-vector renormalization) → `ModeSelect.KerrBackwards`
power-thresholded mode selection → `Harmonics.ylms_few` angular factors → `ModeSum.KerrBackwards`
interpolated mode-sum. Returns the uniform `dt`-spaced grid `t` and the complex strain `h`.

`mode_selection = :all` keeps every mode; `normalize_amps = false` skips the flux renormalization.
When `roman` is supplied, its `ModeIndexTable` determines the mode grid and takes precedence over the
constructor-only `lmax`/`nmax` keywords.
"""
function fast_schwarzschild_eccentric_flux(M, mu, p0, e0, theta, phi;
        T::Real = 1.0, dt::Real = 10.0, eps::Real = 1e-5,
        normalize_amps::Bool = true, mode_selection::Symbol = :power,
        roman = nothing, lmax::Int = 10, nmax::Int = 30,
        Phi_phi0 = 0.0, Phi_theta0 = 0.0, Phi_r0 = 0.0)

    (isfinite(T) && T > 0) || throw(ArgumentError("T must be finite and positive, got $T"))
    mode_selection in (:power, :all) ||
        throw(ArgumentError("mode_selection must be :power or :all, got $mode_selection"))
    roman === nothing && (roman = AmplitudeModels.RomanAmplitude(; lmax = lmax, nmax = nmax))
    flux = FluxModels.KerrBackwards.SchwarzEccFlux()
    insp = EMRIInspiral(SchwarzEccFluxODE(flux); err = 1e-11, enforce_schwarz_sep = true)
    massratio = Trajectory._symmetric_mass_ratio(M, mu)
    y0 = [float(p0), float(e0), 1.0,
          massratio * float(Phi_phi0), massratio * float(Phi_theta0),
          massratio * float(Phi_r0)]
    traj = run_inspiral(insp, float(M), float(mu), 0.0, y0;
                        T = T, dt = dt)
    t = traj.t; p = traj.p; e = traj.e
    nt = length(t)

    table = roman.table
    modes = NTuple{3,Int}[(table.l_arr[i], table.m_arr[i], table.n_arr[i])
                          for i in 1:table.num_modes]
    l_arr = [mm[1] for mm in modes]; m_arr = [mm[2] for mm in modes]; n_arr = [mm[3] for mm in modes]
    K = length(modes)

    teuk = Matrix{ComplexF64}(undef, nt, K)
    for i in 1:nt
        teuk[i, :] .= AmplitudeModels.get_amplitudes(roman, p[i], e[i]; backend = :few, renormalize = normalize_amps)
    end

    # angular factors over the full ±m table (base m≥0 modes, then m>0 → −m partners).
    mpos = findall(m_arr .!= 0)
    full_l = vcat(l_arr, l_arr[mpos]); full_m = vcat(m_arr, .-m_arr[mpos])
    ylms = ylms_few(full_l, full_m, theta, phi)

    if mode_selection === :all
        ms = collect(m_arr); ns = collect(n_arr)
        keep = collect(1:K)
        temp2 = [k <= (K - length(mpos)) ? k : k + length(mpos) for k in keep]
        teuk_in = teuk; ylms_in = ylms[vcat(keep, temp2)]
    else
        sel = ModeSelect.KerrBackwards.mode_selector(teuk, ylms, l_arr, m_arr, n_arr; eps = eps)
        teuk_in = sel.teuk_modes; ylms_in = sel.ylms; ms = sel.m_arr; ns = sel.n_arr
    end

    h = ModeSum.KerrBackwards.interpolated_mode_sum(t, teuk_in, ylms_in, traj.Phi_phi, traj.Phi_r,
                                                    ms, ns; dt = dt)
    tout = [dt * (i - 1) for i in 1:length(h)]
    return (t = tout, h = h)
end

# =====================================================================================
# Source ↔ detector frame algebra (KerrBackwards/waveform.jl).
# =====================================================================================

"""
    get_viewing_angles(qS, phiS, qK, phiK) -> (theta, phi)

Source-frame viewing angles from the detector-frame sky location `(qS, phiS)` and spin orientation
`(qK, phiK)` (`GenerateEMRIWaveform._get_viewing_angles`). With sky unit vector `R` and spin unit vector
`S`, the azimuth is fixed to `φ = −π/2` and the polar angle is `θ = arccos(−R·S)`.
"""
function get_viewing_angles(qS, phiS, qK, phiK)
    R = (sin(qS) * cos(phiS), sin(qS) * sin(phiS), cos(qS))
    S = (sin(qK) * cos(phiK), sin(qK) * sin(phiK), cos(qK))
    RdotS = R[1] * S[1] + R[2] * S[2] + R[3] * S[3]
    return (acos(clamp(-RdotS, -1, 1)), -π / 2)
end

"""
    ssb_polarization_angle(qS, phiS, qK, phiK) -> psi

Polarization angle `ψ_LDC` rotating the source-frame polarizations into the SSB (detector) frame
(`GenerateEMRIWaveform._to_SSB_frame`). Returns `π/2` on the degenerate `sin(qK) sin(phiS−phiK) = 0` axis.
"""
function ssb_polarization_angle(qS, phiS, qK, phiK)
    up = cos(qS) * sin(qK) * cos(phiS - phiK) - cos(qK) * sin(qS)
    dw = sin(qK) * sin(phiS - phiK)
    return dw != 0.0 ? -atan(up, dw) : 0.5 * π
end

"""
    to_SSB_frame(hp, hc, qS, phiS, qK, phiK) -> (hp_new, hc_new)

Rotate source-frame polarizations `(h₊, h_×)` into the SSB frame by the `2ψ_LDC` polarization rotation:
`h₊' = cos2ψ·h₊ − sin2ψ·h_×`, `h_×' = sin2ψ·h₊ + cos2ψ·h_×`. Scalars or arrays (elementwise).
"""
function to_SSB_frame(hp, hc, qS, phiS, qK, phiK)
    psi = ssb_polarization_angle(qS, phiS, qK, phiK)
    c2 = cos(2 * psi); s2 = sin(2 * psi)
    return (c2 .* hp .- s2 .* hc, s2 .* hp .+ c2 .* hc)
end

# =====================================================================================
# Top-level dispatcher (FEW waveform/generate.jl + KerrBackwards waveform.jl).
# =====================================================================================

# Per-generator metadata (FEW baseclasses): (background, descriptor, frame, needs_Y).
const _GENERATOR_INFO = Dict{Symbol,NamedTuple}(
    :FastSchwarzschildEccentricFlux =>
        (background = :Schwarzschild, descriptor = :eccentric, frame = :source, needs_Y = false),
    :SlowSchwarzschildEccentricFlux =>
        (background = :Schwarzschild, descriptor = :eccentric, frame = :source, needs_Y = false),
    :Pn5AAKWaveform =>
        (background = :Kerr, descriptor = :generic, frame = :detector, needs_Y = true),
)

"""
    generate_emri_waveform(waveform_class, M, mu, a, p0, e0, x0, dist, qS, phiS, qK, phiK;
                           frame=:detector, return_list=false, Phi_phi0=0.0, Phi_theta0=0.0,
                           Phi_r0=0.0, T=1.0, dt=10.0, aak_generator=_aak_hole, kwargs...)

Top-level EMRI waveform dispatcher (`GenerateEMRIWaveform`). `waveform_class` selects the generator;
`x0 → Y₀ = xI_to_Y(a,p0,e0,x0)` is applied when the generator `needs_Y`. The detector-frame
`:Pn5AAKWaveform` generator is produced by `aak_generator` (default = [`AAK.aak_waveform`]), placed on the
`t = dt·(i−1)` grid, distance-scaled (`(dist·Gpc)/(mu·MRSUN_SI)`; identity in the detector frame, inverted
into the source frame), and — for `frame=:detector` — rotated into the SSB frame via
[`to_SSB_frame`](@ref). The source-frame `:FastSchwarzschildEccentricFlux` generator is produced by
[`fast_schwarzschild_eccentric_flux`]. Returns `(t, h)` with `h = h₊ − i h_×` (`return_list=false`) or
`(t, hp, hc)`.
"""
function generate_emri_waveform(waveform_class::Symbol, M, mu, a, p0, e0, x0, dist,
                                qS, phiS, qK, phiK;
                                frame::Symbol = :detector, return_list::Bool = false,
                                Phi_phi0 = 0.0, Phi_theta0 = 0.0, Phi_r0 = 0.0,
                                T::Real = 1.0, dt::Real = 10.0,
                                aak_generator = AAK.aak_waveform, kwargs...)
    haskey(_GENERATOR_INFO, waveform_class) ||
        throw(ArgumentError("unknown waveform_class $waveform_class; choose one of \
                            $(collect(keys(_GENERATOR_INFO)))"))
    frame in (:detector, :source) ||
        throw(ArgumentError("frame must be :detector or :source, got $frame"))
    info = _GENERATOR_INFO[waveform_class]

    if info.frame === :source
        # Source-frame Schwarzschild generator: produced at the source viewing angles, distance-scaled,
        # π-rotated (source-frame convention), then (frame=:detector) rotated into the SSB frame.
        waveform_class === :FastSchwarzschildEccentricFlux ||
            error("[FUTURE] generate_emri_waveform source-frame generator $waveform_class is grid-blocked \
                   (SlowSchwarzschildEccentricFlux not built).")
        theta, phi = get_viewing_angles(qS, phiS, qK, phiK)
        res = fast_schwarzschild_eccentric_flux(M, mu, p0, e0, theta, phi;
                                                T = T, dt = dt,
                                                Phi_phi0 = Phi_phi0,
                                                Phi_theta0 = Phi_theta0,
                                                Phi_r0 = Phi_r0, kwargs...)
        t = res.t
        dist_dimensionless = frame === :detector ? (dist * Gpc_in_m) / (mu * MRSUN_SI) : 1.0
        h = (res.h ./ dist_dimensionless) .* (-1)
    else
        # Detector-frame AAK generator.
        Y0 = info.needs_Y ? xI_to_Y(a, p0, e0, x0) : x0
        hvec = aak_generator(M, mu, a, p0, e0, Y0, qS, phiS, qK, phiK, dist;
                             T = T, dt = dt, Phi_phi0 = Phi_phi0, Phi_theta0 = Phi_theta0,
                             Phi_r0 = Phi_r0, kwargs...)
        t = [dt * (i - 1) for i in 1:length(hvec)]
        dist_dimensionless = frame === :source ? 1.0 / ((dist * Gpc_in_m) / (mu * MRSUN_SI)) : 1.0
        h = hvec ./ dist_dimensionless
    end

    if frame === :detector
        hp, hc = to_SSB_frame(real.(h), -imag.(h), qS, phiS, qK, phiK)
    else
        hp, hc = real.(h), -imag.(h)
    end

    return return_list ? (t = t, hp = hp, hc = hc) : (t = t, h = hp .- im .* hc)
end

end # module Generate
