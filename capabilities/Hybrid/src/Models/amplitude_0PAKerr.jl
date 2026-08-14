# 0PAKerr amplitude model (`AmplitudeModels/0PAKerr.m`): adiabatic (0PA) waveform strain
# amplitude for circular equatorial Kerr orbits. Only the (2,2) mode is tabulated, as a real
# 1SF amplitude `RealAmp22` on a 2D cos-spaced grid (Hermite, InterpolationOrder->{All,All}),
# indexed by the magic abscissae `(c1, c2)` of `0PAKerr.m`.

"""
    Amp0PAKerr

0PAKerr amplitude model: the real (2,2) 1SF amplitude `RealAmp22` as a 2D global-Lagrange
([`Interp2Hermite`](@ref)) interpolant. Build with [`amplitudes_0PAKerr`](@ref).
"""
struct Amp0PAKerr
    R22::Interp2Hermite{Float64}
end

"""
    amplitudes_0PAKerr(amp_datadir) -> Amp0PAKerr

Load `1SF/Kerr_Circ/1SF_RealAmp_KerrCirc_22_NewGrid.m` from `amp_datadir`
(`AmplitudeModels/sf_amp_data`) and build the 2D Hermite interpolant.
"""
function amplitudes_0PAKerr(amp_datadir::AbstractString)
    g = read_wl_grid2d(joinpath(amp_datadir, "1SF", "Kerr_Circ",
                                "1SF_RealAmp_KerrCirc_22_NewGrid.m"))
    return Amp0PAKerr(interp2_hermite(g))
end

# Magic abscissa constants of the (2,2) amplitude grid remap (`0PAKerr.m`):
#   c1 = OFFSET + SLOPE·Log[1-χ1],   c2 = -1 + 2 √(ω^(2/3) rISCO[χ1] / (1-χ1 ω)^(2/3)).
const _KERR_AMP_C1_OFFSET = 0.7988880369434854
const _KERR_AMP_C1_SLOPE = 0.2901432317650783

"""
    real_amp22(model::Amp0PAKerr, ω, χ1, ν, M) -> Float64

Strain amplitude of the (2,2) mode at orbital frequency `ω`, spin `χ1`, mass ratio `ν`
(M units):
`(1/(2(1-χ1 ω)^(2/3))) M ν ω^(2/3) RealAmp22[c1, c2]` with the `0PAKerr.m` abscissa remap.
Reuses [`kerr_isco`](@ref) for `rISCO[χ1]`.
"""
function real_amp22(model::Amp0PAKerr, ω::Real, χ1::Real, ν::Real, M::Real)
    c1 = _KERR_AMP_C1_OFFSET + _KERR_AMP_C1_SLOPE * log(1 - χ1)
    c2 = -1 + 2 * sqrt(ω^(2 / 3) * kerr_isco(χ1) / (1 - χ1 * ω)^(2 / 3))
    return (1 / (2 * (1 - χ1 * ω)^(2 / 3))) * M * ν * ω^(2 / 3) * model.R22(c1, c2)
end

"""
    mode_keys_0PAKerr() -> Vector{Tuple{Int,Int}}

The `(l, m)` modes of the 0PAKerr model: only `(2, 2)` is tabulated (`0PAKerr.m`).
"""
mode_keys_0PAKerr() = [(2, 2)]
