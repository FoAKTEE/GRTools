# WaSABI-C SF-PN hybrid: interpolant layer + dependent variables + the WF-frequency
# correction δx (`InspiralModels/WaSABI-C.m`, §Dependent variables + §Construct Interpolations).
#
# `WaSABICModel` bundles every interpolant the WaSABI-C closed forms read:
#   * Kerr 1SF energy fluxes FI/FH (2D Hermite over (Log[1-a], √(rISCO[a]/r)), node r^5·F∞ /
#     r^(15/2)·F_H) — identical grids to `0PAKerr`, so reuse `forcing_0PAKerr`.
#   * leading secondary-spin fluxes FIS2/FHS2 (same remap, node r^5·χ2∞ / r^9·χ2_H).
#   * the redshift invariant ΔU (2D Hermite, node ΔU).
#   * the (2,2) 1SF amplitude Re/Im on the cos-spaced "NewGrid".
#   * Schwarzschild 1SF (order-8) and 2SF (Re∘Exp order-8 log-exp) infinity fluxes.
#   * the 3D StopHybridFixedChi critical-frequency table (InterpolationOrder->All).
#
# δx and dδx/dχ are MACHINE-GENERATED from the WL InputForm (see wasabic_codegen.jl); this file
# supplies the heads those generated bodies call: the dependent-variable closed forms (`_wc_x`,
# `_wc_rISCO`, …) and the interpolant partial-derivative shims. `PolyGamma` is the canonical
# `BigSpecial.PolyGamma` complex digamma (AHA — the legacy hand-rolled `_digamma` fork collapsed
# into the one shared primitive); it is re-exported here for the generated bodies.

# --- dependent-variable closed forms (9-arg, matching the WL signature) -----

# x = ω^(2/3) (PN frequency parameter).
_wc_x(ω, φ, ν, m, χ1, χ2, δm, δν, δχ) = ω^(2 / 3)
# rISCO depends only on χ1 (Bardeen–Press–Teukolsky); identical to `kerr_isco`.
_wc_rISCO(ω, φ, ν, m, χ1, χ2, δm, δν, δχ) = kerr_isco(χ1)

"""
    wc_x(ω) = ω^(2/3)

WaSABI-C PN frequency parameter `x` (`WaSABI-C.m` §PN Parameter).
"""
wc_x(ω::Real) = ω^(2 / 3)

"""
    wc_s(ν, χ1, χ2)

Symmetric spin parameter `s = ¼((1+√(1-4ν))² χ1 + (-1+√(1-4ν))² χ2)`
(`WaSABI-C.m` §Spin parameters).
"""
wc_s(ν::Real, χ1::Real, χ2::Real) =
    (1 / 4) * ((1 + sqrt(1 - 4ν))^2 * χ1 + (-1 + sqrt(1 - 4ν))^2 * χ2)

"""
    wc_sigma(ν, χ1, χ2)

Antisymmetric spin parameter `σ = ½(-(1+√(1-4ν)) χ1 - (-1+√(1-4ν)) χ2)`
(`WaSABI-C.m` §Spin parameters).
"""
wc_sigma(ν::Real, χ1::Real, χ2::Real) =
    (1 / 2) * (-((1 + sqrt(1 - 4ν)) * χ1) - (-1 + sqrt(1 - 4ν)) * χ2)

# --- interpolant partial-derivative shims used by the generated bodies ------

partial_0_1(f, x, y) = partial(f, x, y, 0, 1)
partial_1_0(f, x, y) = partial(f, x, y, 1, 0)
partial_1_1(f, x, y) = partial(f, x, y, 1, 1)
partial_0_2(f, x, y) = partial(f, x, y, 0, 2)
partial_2_0(f, x, y) = partial(f, x, y, 2, 0)
partial1d_1(f, x) = derivative(f, x)

# --- the WaSABI-C model: every interpolant the closed forms read ------------

"""
    WaSABICModel

Bundle of every interpolant used by the WaSABI-C SF-PN hybrid closed forms (`WaSABI-C.m`
§Construct Interpolations). Build with [`wasabic_model`](@ref).
"""
struct WaSABICModel{F2}
    FI1SFKerr::Interp2Hermite{Float64}
    FH1SFKerr::Interp2Hermite{Float64}
    FIS2Kerr::Interp2Hermite{Float64}
    FHS2Kerr::Interp2Hermite{Float64}
    DU::Interp2Hermite{Float64}
    Re22::Interp2Hermite{Float64}
    Im22::Interp2Hermite{Float64}
    FI1SFSchw::Interp1{Float64,Float64}
    FI2SFSchw::F2
    ωcrit::Interp3Lagrange{Float64}
end

# Build a 2D Hermite interpolant over the WaSABI-C Kerr grid (Log[1-a], √(rISCO[a]/r)) with the
# per-row column-snap (rows agree to ~1e-15 by construction of the r-grid), matching the WL
# `Interpolation[Table[...], Method->"Hermite", InterpolationOrder->{All,All}]`.
function _kerr_grid2d(a_x::AbstractVector, a_isco::AbstractVector,
                      r::AbstractVector, val::AbstractVector)
    rows = Dict{Float64,Vector{Tuple{Float64,Float64}}}()
    for k in eachindex(a_x)
        c2 = sqrt(kerr_isco(a_isco[k]) / r[k])
        push!(get!(rows, Float64(a_x[k]), Tuple{Float64,Float64}[]), (c2, Float64(val[k])))
    end
    as = sort!(collect(keys(rows)))
    for a in as
        sort!(rows[a], by = first)
    end
    ny = length(rows[as[1]])
    ycanon = [rows[as[1]][j][1] for j in 1:ny]
    g = [((log(1 - as[i]), ycanon[j]), rows[as[i]][j][2])
         for i in eachindex(as) for j in 1:ny]
    return interp2_hermite(g)
end

"""
    wasabic_model(inspiral_datadir, amp_datadir) -> WaSABICModel

Load every WaSABI-C interpolant from `inspiral_datadir` (`InspiralModels/sf_data`) and
`amp_datadir` (`AmplitudeModels/sf_amp_data`). Reuses [`forcing_0PAKerr`](@ref) for the Kerr 1SF
flux grids (identical construction), and `interp2_hermite`/`interp1_order8`/`interp3_lagrange`
for the rest.
"""
function wasabic_model(inspiral_datadir::AbstractString, amp_datadir::AbstractString)
    # Kerr 1SF FI/FH — identical grids to 0PAKerr.
    f0 = forcing_0PAKerr(inspiral_datadir)
    FI1SFKerr, FH1SFKerr = f0.FI, f0.FH

    # secondary-spin fluxes: x-abscissa from the spin table's a, rISCO from the 1SF table's a.
    flux1 = read_wl_data(joinpath(inspiral_datadir, "1SF_Flux", "Kerr_Circ",
                                  "SMRfluxdata2025_36x36.data"))
    spin = read_wl_data(joinpath(inspiral_datadir, "Leading_S2_Flux", "Kerr_Circ",
                                 "SMRfluxdataWithSecondarySpin2025_36x36.m"))
    a1 = [Float64(rec["a"]) for rec in flux1]
    aS = [Float64(rec["a"]) for rec in spin]
    rS = [Float64(rec["r"]) for rec in spin]
    FIS2Kerr = _kerr_grid2d(aS, a1, rS, [Float64(rec["r"])^5 * Float64(rec["chi2inf"]) for rec in spin])
    FHS2Kerr = _kerr_grid2d(aS, a1, rS, [Float64(rec["r"])^9 * Float64(rec["chi2hor"]) for rec in spin])

    # redshift invariant ΔU.
    invar = read_wl_data(joinpath(inspiral_datadir, "1SF_Local_Invariants", "Kerr_Circ",
                                  "KerrCircularEquatorialInvariants2025b_36x36.data"))
    aU = [Float64(rec["a"]) for rec in invar]
    rU = [Float64(rec["r"]) for rec in invar]
    DU = _kerr_grid2d(aU, aU, rU, [Float64(rec["\\[CapitalDelta]U"]) for rec in invar])

    # (2,2) 1SF amplitude on the NewGrid.
    Re22 = interp2_hermite(read_wl_grid2d(joinpath(amp_datadir, "1SF", "Kerr_Circ",
                                                   "ReAmp1SFKerr22NewGrid.m")))
    Im22 = interp2_hermite(read_wl_grid2d(joinpath(amp_datadir, "1SF", "Kerr_Circ",
                                                   "ImAmp1SFKerr22NewGrid.m")))

    # Schwarzschild 1SF (order-8) + 2SF (log-exp order-8) infinity fluxes.
    FI1SFSchw = interp1_order8(read_wl_table(joinpath(inspiral_datadir, "1SF_Flux", "Schwarz_Circ",
                                                      "1SFCircShwarzDotEInf.m")))
    # WaSABI-C's ℱℐint2SFSchw = Function[{logr0}, Re[Exp[interp[logr0]]]] — its argument is
    # already Log[r0] (no inner log, unlike the 1PAT1 `interp_logexp`).
    _f2itp = interp1_order8(read_wl_table_complex(
        joinpath(inspiral_datadir, "2SF_Flux", "Schwarz_Circ", "2SFCircShwarzDotEInf.m")))
    FI2SFSchw = logr0 -> real(exp(_f2itp(logr0)))

    # 3D critical-frequency stop table.
    ωcrit = interp3_lagrange(read_wl_grid3d(joinpath(inspiral_datadir, "stop_conditions",
                                                     "StopHybridFixedChi.m")))

    return WaSABICModel{typeof(FI2SFSchw)}(FI1SFKerr, FH1SFKerr, FIS2Kerr, FHS2Kerr, DU,
                                           Re22, Im22, FI1SFSchw, FI2SFSchw, ωcrit)
end

# --- public WF-frequency correction ----------------------------------------

"""
    wc_δx(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)

WaSABI-C waveform-frequency correction `δx` (`WaSABI-C.m` §WF frequency correction). Machine
-generated closed form over the (2,2) amplitude + Kerr 1SF flux interpolants.
"""
wc_δx(mdl::WaSABICModel, ω, φ, ν, m, χ1, χ2, δm, δν, δχ) =
    _wc_dx(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)

"""
    wc_dδxdχ(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)

`∂δx/∂χ1` (`WaSABI-C.m`), the χ1-derivative of [`wc_δx`](@ref).
"""
wc_dδxdχ(mdl::WaSABICModel, ω, φ, ν, m, χ1, χ2, δm, δν, δχ) =
    _wc_ddxdχ(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)
