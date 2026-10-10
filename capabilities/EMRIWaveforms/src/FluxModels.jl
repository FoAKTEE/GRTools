"""
    EMRIWaveforms.FluxModels

Grid-interpolated flux ODE right-hand sides (design §6.8: KerrBackwards `schwarzeccflux.jl`; FEW
`trajectory/kerrflux.jl`). Distinct from the already-built analytic 5PN table (`Trajectory.Pn5Fluxes`):
here the orbit-averaged fluxes are the closed-form PN normalization PLUS a residual interpolated over an
on-disk Teukolsky flux grid loaded through `Data`, from which `(ṗ, ė)` follow algebraically.

# v0.1 TRANSITIONAL two-fork (§4.8 collision survivor `SchwarzEccFlux`)
- `FluxModels.KerrBackwards` — the Schwarzschild eccentric fork: a `GslBicubic` residual over the
  `(y, e)` grid in `FluxNewMinusPNScaled_fixed_y_order.dat`, closed-form `EdotPN`/`LdotPN`, and the
  verbatim `(ṗ, ė)` algebra with the `p ≤ 6 + 2e` separatrix guard.
- `FluxModels.FEW` — the Kerr eccentric-equatorial fork: per-region `TricubicNotAKnot` residuals over
  `KerrEccEqFluxData.h5` (`:A` close / `:B` far) loaded via `Data.read_fluxgrid`, the PN normalizations
  `kerr_PN_alt`/`kerr_pdot_PN`/`kerr_edot_PN`, and the apex-domain `evaluate_rhs`.

The forks are reachable only qualified (`FluxModels.KerrBackwards.*` / `FluxModels.FEW.*`); they are
collapsed only after a per-submodule bit-equivalence gate (§13). The flux-grid loaders are the single
consolidated `Data.read_dat` / `Data.read_fluxgrid` readers (the FEW-vs-KB duplication is de-duped at
the `Data` tier).
"""
module FluxModels

# ==========================================================================================
# KerrBackwards fork — Schwarzschild eccentric flux (KerrBackwards/schwarzeccflux.jl).
# ==========================================================================================
module KerrBackwards

using NumericsCore: GslBicubic
using Spacetime: kerr_frequencies
import Data

export SchwarzEccFlux, schwarz_ecc_flux_deriv, ode_deriv

"Default Schwarzschild flux-residual grid (resolved under the configured `Data` datadir)."
_default_flux_file() = Data.datapath("FluxNewMinusPNScaled_fixed_y_order.dat")

"""
    SchwarzEccFlux(; file=Data.datapath("FluxNewMinusPNScaled_fixed_y_order.dat"))

Schwarzschild eccentric flux-driven inspiral model (`SchwarzEccFlux::deriv_func`,
`src/ode_base_example.cc`). The orbit-averaged `(Ė, L̇z)` are the closed-form PN pieces plus a bicubic
residual over the cached `(y, e)` grid. The residual surfaces `Edot`/`Ldot` are [`GslBicubic`] (GSL
bicubic, as upstream) over `y = log(p − 2e − 2.1)` and `e` node axes loaded with `Data.read_dat`.
"""
struct SchwarzEccFlux
    Edot::GslBicubic{Float64}
    Ldot::GslBicubic{Float64}
end

convert_Y(::SchwarzEccFlux) = false
background(::SchwarzEccFlux) = :Schwarzschild
num_add_args(::SchwarzEccFlux) = 0

function SchwarzEccFlux(; file::AbstractString = _default_flux_file())
    raw = Data.read_dat(file)
    yc = Float64.(raw[:, 1]); ec = Float64.(raw[:, 2])
    Ec = Float64.(raw[:, 3]); Lc = Float64.(raw[:, 4])
    uy = sort(unique(yc)); ue = sort(unique(ec))
    nx = length(uy); ny = length(ue)
    nx * ny == length(yc) ||
        throw(ArgumentError("flux grid $(nx)×$(ny) ≠ $(length(yc)) rows"))
    # File order: row k = (j-1)*nx + i with i the y index, j the e index.
    ZE = Matrix{Float64}(undef, nx, ny)
    ZL = Matrix{Float64}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        k = (j - 1) * nx + i
        ZE[i, j] = Ec[k]
        ZL[i, j] = Lc[k]
    end
    return SchwarzEccFlux(GslBicubic(uy, ue, ZE), GslBicubic(uy, ue, ZL))
end

"`_EdotPN(e, yPN)` — leading PN energy flux `(96+292e²+37e⁴)/(15(1−e²)^3.5)·yPN^5`."
_EdotPN(e, yPN) = (96 + 292 * e^2 + 37 * e^4) / (15 * (1 - e^2)^3.5) * yPN^5

"`_LdotPN(e, yPN)` — leading PN angular-momentum flux `4(8+7e²)/(5(−1+e²)²)·yPN^(7/2)`."
_LdotPN(e, yPN) = (4 * (8 + 7 * e^2)) / (5 * (-1 + e^2)^2) * yPN^(7 / 2)

"""
    _pdot_edot(Edot, Ldot, p, e) -> (pdot, edot)

The verbatim Schwarzschild `(ṗ, ė)` algebra from `(Ė, L̇z)` (`SchwarzEccFlux::deriv_func`). `e == 0`
zeroes `ė`. Pure closed form — no grid access.
"""
function _pdot_edot(Edot, Ldot, p, e)
    pdot = (-2 * (Edot * sqrt((4 * e^2 - (-2 + p)^2) / (3 + e^2 - p)) * (3 + e^2 - p) * p^1.5 +
                  Ldot * (-4 + p)^2 * sqrt(-3 - e^2 + p))) / (4 * e^2 - (-6 + p)^2)
    if e > 0.0
        edot = -((Edot * sqrt((4 * e^2 - (-2 + p)^2) / (3 + e^2 - p)) * p^1.5 *
                  (18 + 2 * e^4 - 3 * e^2 * (-4 + p) - 9 * p + p^2) +
                  (-1 + e^2) * Ldot * sqrt(-3 - e^2 + p) * (12 + 4 * e^2 - 8 * p + p^2)) /
                 (e * (4 * e^2 - (-6 + p)^2) * p))
    else
        edot = 0.0
    end
    return (pdot, edot)
end

"""
    schwarz_ecc_flux_deriv(model::SchwarzEccFlux, epsilon, a, p, e, x)
    ode_deriv(model::SchwarzEccFlux, epsilon, a, p, e, x)

Schwarzschild eccentric flux RHS: `(pdot, edot, xdot, Ω_φ, Ω_θ, Ω_r)`. At/under the separatrix
(`p ≤ 6 + 2e`) every rate is zeroed (upstream guard). `epsilon` is the mass ratio `μ/M`.
"""
function schwarz_ecc_flux_deriv(model::SchwarzEccFlux, epsilon, a, p, e, x)
    if (6.0 + 2.0 * e) >= p
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    f = kerr_frequencies(0.0, p, e, 1.0)      # Schwarzschild geodesic frequencies
    Omega_phi = f.Omegaph
    Omega_r = f.Omegar
    Omega_theta = Omega_phi

    y1 = log(p - 2.0 * e - 2.1)
    yPN = Omega_phi^(2 / 3)

    Edot = -epsilon * (model.Edot(y1, e) * yPN^6 + _EdotPN(e, yPN))
    Ldot = -epsilon * (model.Ldot(y1, e) * yPN^(9 / 2) + _LdotPN(e, yPN))

    pdot, edot = _pdot_edot(Edot, Ldot, p, e)
    return (pdot, edot, 0.0, Omega_phi, Omega_theta, Omega_r)
end

ode_deriv(model::SchwarzEccFlux, epsilon, a, p, e, x) =
    schwarz_ecc_flux_deriv(model, epsilon, a, p, e, x)

end # module KerrBackwards

# ==========================================================================================
# FEW fork — Kerr eccentric-equatorial flux (FEW/trajectory/kerrflux.jl).
# ==========================================================================================
module FEW

using NumericsCore: TricubicNotAKnot
using Spacetime: kerr_frequencies, get_separatrix
using ...Mappings: apex_of_uwyz, apex_of_UWYZ, kerrecceq_forward_map, ELdot_to_PEdot_Jacobian
import Data

export KerrEccEqFlux, kerr_PN_alt, kerr_pdot_PN, kerr_edot_PN,
       interpolate_flux_grids, evaluate_rhs

"Default Kerr-ecc-eq flux grid (resolved under the configured `Data` datadir)."
_default_flux_file() = Data.datapath("KerrEccEqFluxData.h5")

"""
    kerr_PN_alt(p, e) -> (Edot, Ldot)

Leading PN energy/angular-momentum flux normalization, ELQ convention
(`few.trajectory.ode.flux._PN_alt`, arXiv:2201.07044 eq. 91).
"""
function kerr_PN_alt(p, e)
    oneme2 = (1 - e^2)^1.5
    Edot = 32.0 / 5.0 * p^(-5) * oneme2 * (1 + 73 / 24 * e^2 + 37 / 96 * e^4)
    Ldot = 32.0 / 5.0 * p^(-7 / 2) * oneme2 * (1 + 7.0 / 8.0 * e^2)
    return Edot, Ldot
end

"`kerr_pdot_PN(p, e, risco, p_sep)` — PN `pdot` normalization, pex convention (`_pdot_PN`)."
kerr_pdot_PN(p, e, risco, p_sep) =
    (8.0 * (1.0 - e^2)^1.5 * (8.0 + 7.0 * e^2)) /
    (5.0 * p * ((p - risco)^2 - (-risco + p_sep)^2))

"`kerr_edot_PN(p, e, risco, p_sep)` — PN `edot` normalization, pex convention (`_edot_PN`)."
kerr_edot_PN(p, e, risco, p_sep) =
    ((1.0 - e^2)^1.5 * (304.0 + 121.0 * e^2)) /
    (15.0 * p^2 * ((p - risco)^2 - (-risco + p_sep)^2))

"""
    KerrEccEqFlux(; convention="pex", downsample=((1,1,1),(1,1,1)),
                    file=Data.datapath("KerrEccEqFluxData.h5"))

Kerr eccentric-equatorial flux-grid loader (`few.trajectory.ode.flux.KerrEccEqFlux`). Loads region-A
(close) and region-B (far) flux samples via `Data.read_fluxgrid`, reconstructs the apex `(a,p,e,x)` per
`(u,w,z)` node through the grounded `apex_of_uwyz`/`apex_of_UWYZ` maps, divides by the PN normalization,
and builds [`TricubicNotAKnot`] residual interpolants. `convention="ELQ"` stores `Edot`/`Ldot`;
`"pex"` converts to `pdot`/`edot` via `ELdot_to_PEdot_Jacobian`. `downsample = ((su,sw,sz)_A,
(su,sw,sz)_B)` strides the grids in lockstep (grid-convergence / testing).
"""
struct KerrEccEqFlux
    convention::String
    interp::Dict{Symbol,TricubicNotAKnot{Float64}}
    grids::Dict{Symbol,Array{Float64,3}}
end

convert_Y(::KerrEccEqFlux) = false
background(::KerrEccEqFlux) = :Kerr
num_add_args(::KerrEccEqFlux) = 0

# Apex coordinates of one (downsampled) region's (u,w,z) nodes (z fixed via the y=1 equatorial slice).
function _region_apex(u, w, z, is_regionA)
    nu, nw, nz = length(u), length(w), length(z)
    p = Array{Float64}(undef, nu, nw, nz)
    e = Array{Float64}(undef, nu, nw, nz)
    a = Array{Float64}(undef, nu, nw, nz)
    for k in 1:nz, j in 1:nw, i in 1:nu
        ai, pi, ei, _ = is_regionA ? apex_of_uwyz(u[i], w[j], 1.0, z[k]) :
                        apex_of_UWYZ(u[i], w[j], 1.0, z[k], true)
        a[i, j, k] = ai; p[i, j, k] = pi; e[i, j, k] = ei
    end
    return a, p, e
end

function KerrEccEqFlux(;
    convention::String = "pex",
    downsample = ((1, 1, 1), (1, 1, 1)),
    file::AbstractString = _default_flux_file(),
)
    convention in ("pex", "ELQ") || throw(ArgumentError("convention must be \"pex\" or \"ELQ\""))
    interp = Dict{Symbol,TricubicNotAKnot{Float64}}()
    grids = Dict{Symbol,Array{Float64,3}}()

    for (region, ds, is_A, tag) in (
        (:A, downsample[1], true, :A),
        (:B, downsample[2], false, :B),
    )
        schema, ch = Data.read_fluxgrid(file, region; downsample = ds)
        u, w, z = schema.u, schema.w, schema.z
        E, L = ch[:Edot], ch[:Ldot]
        a, p, e = _region_apex(u, w, z, is_A)
        nu, nw, nz = size(E)
        if convention == "ELQ"
            Eg = Array{Float64}(undef, nu, nw, nz)
            Lg = Array{Float64}(undef, nu, nw, nz)
            for k in 1:nz, j in 1:nw, i in 1:nu
                EPN, LPN = kerr_PN_alt(p[i, j, k], e[i, j, k])
                Eg[i, j, k] = E[i, j, k] / EPN
                Lg[i, j, k] = L[i, j, k] / LPN
            end
            interp[Symbol(:Edot, tag)] = TricubicNotAKnot(u, w, z, Eg)
            interp[Symbol(:Ldot, tag)] = TricubicNotAKnot(u, w, z, Lg)
            grids[Symbol(:Edot, tag)] = Eg
            grids[Symbol(:Ldot, tag)] = Lg
        else
            pg = Array{Float64}(undef, nu, nw, nz)
            eg = Array{Float64}(undef, nu, nw, nz)
            for k in 1:nz, j in 1:nw, i in 1:nu
                aa = a[i, j, k]
                xs = aa == 0 ? 1.0 : sign(aa)
                a_abs = abs(aa)
                pdot, edot = ELdot_to_PEdot_Jacobian(
                    a_abs, p[i, j, k], e[i, j, k], xs, E[i, j, k], L[i, j, k] * xs)
                risco = get_separatrix(a_abs, 0.0, 1.0)
                psep = get_separatrix(a_abs, e[i, j, k], 1.0)
                pg[i, j, k] = pdot / kerr_pdot_PN(p[i, j, k], e[i, j, k], risco, psep)
                eg[i, j, k] = edot / kerr_edot_PN(p[i, j, k], e[i, j, k], risco, psep)
            end
            interp[Symbol(:pdot, tag)] = TricubicNotAKnot(u, w, z, pg)
            interp[Symbol(:edot, tag)] = TricubicNotAKnot(u, w, z, eg)
            grids[Symbol(:pdot, tag)] = pg
            grids[Symbol(:edot, tag)] = eg
        end
    end
    return KerrEccEqFlux(convention, interp, grids)
end

"""
    interpolate_flux_grids(flux::KerrEccEqFlux, p, e, x, a; pLSO=nothing) -> (f1, f2)

Total fluxes at apex `(a, p, e, x)`: `(Ė, L̇)` (ELQ) or `(ṗ, ė)` (pex). The apex is mapped to the grid
`(u, w, z)` by `kerrecceq_forward_map`; the region-tagged residual interpolant is multiplied back by the
PN normalization (`KerrEccEqFlux.interpolate_flux_grids`).
"""
function interpolate_flux_grids(flux::KerrEccEqFlux, p::Real, e::Real, x::Real, a::Real;
                                pLSO = nothing)
    a_in = x == -1 ? -a : a
    # FEW flux.py:640-641 takes the separatrix at the apex (a, e, x), not at the grid's signed a_in.
    pl = pLSO === nothing ? get_separatrix(a, e, x) : pLSO
    u, w, _, z, inA = kerrecceq_forward_map(a_in, p, e, 1.0; pLSO = pl,
                                            return_mask = true, kind = "flux")
    tag = inA ? :A : :B
    if flux.convention == "ELQ"
        EPN, LPN = kerr_PN_alt(p, e)
        Edot = -flux.interp[Symbol(:Edot, tag)](u, w, z) * EPN
        Ldot = -flux.interp[Symbol(:Ldot, tag)](u, w, z) * LPN
        a_in < 0 && (Ldot *= -1)
        return Edot, Ldot
    else
        risco = get_separatrix(a_in, 0.0, 1.0)
        pdotPN = kerr_pdot_PN(p, e, risco, pl)
        edotPN = kerr_edot_PN(p, e, risco, pl)
        pdot = -flux.interp[Symbol(:pdot, tag)](u, w, z) * pdotPN
        edot = -flux.interp[Symbol(:edot, tag)](u, w, z) * edotPN
        return pdot, edot
    end
end

"""
    evaluate_rhs(flux::KerrEccEqFlux, a, p, e, x, p_sep) -> Vector{Float64}

Kerr-ecc-eq flux RHS `[ḟ₁, ḟ₂, 0, Ω_φ, Ω_θ, Ω_r]` (`KerrEccEqFlux.evaluate_rhs`): `(ḟ₁, ḟ₂)` from
[`interpolate_flux_grids`] with the cached separatrix `p_sep` as `pLSO`, and `(Ω_φ, Ω_θ, Ω_r)` from the
Kerr geodesic frequencies at `(a, p, e, x)`.
"""
function evaluate_rhs(flux::KerrEccEqFlux, a, p, e, x, p_sep)
    fr = kerr_frequencies(a, p, e, x)
    f1, f2 = interpolate_flux_grids(flux, p, e, x, a; pLSO = p_sep)
    return [f1, f2, 0.0, fr.Omegaph, fr.Omegath, fr.Omegar]
end

end # module FEW

end # module FluxModels
