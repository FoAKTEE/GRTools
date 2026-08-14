# 0PAKerr inspiral forcing (`InspiralModels/0PAKerr.m`): adiabatic (0PA) evolution of a
# circular equatorial Kerr EMRI. The waveform-frequency forcing `Fω[ω,ν,χ1]` is a closed-form
# prefactor times the 1SF horizon + infinity energy fluxes, each tabulated on a 2D 36×36 grid
# (Hermite, InterpolationOrder->{All,All}) in `(Log[1-a], √(rISCO[a]/r))` with node values
# `r^5·F_∞` (infinity) and `r^(15/2)·F_H` (horizon). The state `(ω, φ)` evolves with
# `dω/dt = Fω/M`, `dφ/dt = ω/M`; integration stops when `ω` reaches the critical frequency
# `Ωcrit(χ1) = min(1/((1.05 rISCO[a])^{3/2}+χ1), 1/(6.26^{3/2}+χ1))`. The ODE driver
# (`int_inspiral`) lives in `Drivers`; this leaf owns only the forcing + rhs algebra.

"""
    Forcing0PAKerr

0PAKerr forcing model: the 1SF infinity (`FI`) and horizon (`FH`) energy-flux interpolants as
2D global-Lagrange grids ([`Interp2Hermite`](@ref)) over `(Log[1-a], √(rISCO[a]/r))`. Build with
[`forcing_0PAKerr`](@ref); evaluate the forcing with [`forcing_Fω`](@ref).
"""
struct Forcing0PAKerr
    FI::Interp2Hermite{Float64}   # ℱℐ: node value r^5 · F_∞
    FH::Interp2Hermite{Float64}   # ℱℋ: node value r^(15/2) · F_H
end

"""
    forcing_0PAKerr(sf_datadir) -> Forcing0PAKerr

Load `1SF_Flux/Kerr_Circ/SMRfluxdata2025_36x36.data` from `sf_datadir`
(`InspiralModels/sf_data`) and build the two 2D Hermite flux interpolants. The grid's second
abscissa `√(rISCO[a]/r)` is computed per node and agrees across `a`-rows to ~1e-15; the rows are
snapped onto one canonical column grid so the data forms an exact tensor product.
"""
function forcing_0PAKerr(sf_datadir::AbstractString)
    d = read_wl_data(joinpath(sf_datadir, "1SF_Flux", "Kerr_Circ",
                              "SMRfluxdata2025_36x36.data"))
    # Group records by spin a; each (c2, zI, zH) node carries the remapped flux values.
    rows = Dict{Float64,Vector{NTuple{3,Float64}}}()
    for rec in d
        a = Float64(rec["a"])
        r = Float64(rec["r"])
        c2 = sqrt(kerr_isco(a) / r)
        zI = r^5 * Float64(rec["inf"])
        zH = r^(15 / 2) * Float64(rec["hor"])
        push!(get!(rows, a, NTuple{3,Float64}[]), (c2, zI, zH))
    end
    as = sort!(collect(keys(rows)))
    for a in as
        sort!(rows[a], by = first)   # order each a-row by the c2 abscissa
    end
    ny = length(rows[as[1]])
    ycanon = [rows[as[1]][j][1] for j in 1:ny]   # canonical column grid (rows agree to ~1e-15)
    gI = [((log(1 - as[i]), ycanon[j]), rows[as[i]][j][2]) for i in eachindex(as) for j in 1:ny]
    gH = [((log(1 - as[i]), ycanon[j]), rows[as[i]][j][3]) for i in eachindex(as) for j in 1:ny]
    return Forcing0PAKerr(interp2_hermite(gI), interp2_hermite(gH))
end

"""
    forcing_Fω(model::Forcing0PAKerr, ω, ν, χ1) -> Float64

Waveform-frequency forcing `Fω[ω,ν,χ1]` (`0PAKerr.m`): the closed-form prefactor times the
2D-interpolated horizon + infinity fluxes at `(Log[1-χ1], √(rISCO[χ1]/r))` with
`r = (1/ω - χ1)^(2/3)`. The inspiral's `dω/dt = Fω/M`.
"""
function forcing_Fω(model::Forcing0PAKerr, ω::Real, ν::Real, χ1::Real)
    s = χ1 * ω
    u = 1 - s
    rr = (1 / ω - χ1)^(2 / 3)
    c1 = log(1 - χ1)
    c2 = sqrt(kerr_isco(χ1) / rr)
    fi = model.FI(c1, c2)
    fh = model.FH(c1, c2)
    num = 3 * ν * ((-1 + s) * (-1 - s + 3 * ω^(2 / 3) * u^(1 / 3)))^(3 / 2) *
          (ω^(5 / 3) * u^(1 / 3) * fh + (-1 + s)^2 * fi)
    den = (-χ1 + 1 / ω)^(11 / 3) * (-1 + s) *
          (-1 + 6 * ω^(2 / 3) * u^(1 / 3) - 6 * χ1 * (ω + ω^(5 / 3) * u^(1 / 3)) +
           χ1^2 * (7 * ω^2 + 3 * ω^(4 / 3) * u^(2 / 3)))
    return num / den
end

"""
    omega_crit(model::Forcing0PAKerr, χ1) -> Float64

Critical (terminating) orbital frequency `Ωcrit = min(1/((1.05 rISCO[χ1])^{3/2}+χ1),
1/(6.26^{3/2}+χ1))` (`0PAKerr.m`); integration stops when `ω ≥ Ωcrit`.
"""
function omega_crit(::Forcing0PAKerr, χ1::Real)
    ri = kerr_isco(χ1)
    return min(1 / ((1.05 * ri)^(3 / 2) + χ1), 1 / (6.26^(3 / 2) + χ1))
end

"""
    kerr_inspiral_rhs!(du, u, p, t)

In-place ODE right-hand side for the 0PA Kerr circular-equatorial inspiral, state
`u = (ω, φ)` with parameters `p = (model::Forcing0PAKerr, ν, M, χ1)`:
`dω/dt = Fω(ω; ν, χ1)/M`, `dφ/dt = ω/M` (the `0PAKerr.m` evolution equations, with
[`forcing_Fω`](@ref) the 1SF flux forcing). SciML-signature rhs consumed by the
`Hybrid.Drivers` `Vern9` integrator; `t` is unused (autonomous system).
"""
function kerr_inspiral_rhs!(du, u, p, t)
    m, ν, M, χ1 = p
    ω, φ = u
    du[1] = forcing_Fω(m, ω, ν, χ1) / M
    du[2] = ω / M
    return nothing
end
