"""
    PostNewtonian.SelfForce

PN self-force accessor layer over the loaded series store: the [`PostNewtonianData`](@ref) wrapper +
order metadata, [`postnewtonian_coefficient`](@ref) (the `y^n Log[y]^nL` extraction),
[`ResummedSeriesData`](@ref) (`sf*sd` resummed series), and the [`postnewtonian_expansion`](@ref)
load-by-name entry. Builds on `Series` (`PNSeries`/`Normal`), `Catalog` (discovery), and
`MMAInputForm` (the `.m` loader). `Normal` is EXTENDED here on `ResummedSeriesData` over the single
`Series.Normal` method table. `ResummedSeriesData.sf::Symbolics.Num` is the second Symbolics
quarantine anchor.
"""
module SelfForce

using ..Series: PNSeries, Normal, series_min_order
using ..Catalog: SERIES_DATA_DIR, pn_browse
using ..MMAInputForm: parse_pn_data_file
import ..Series: Normal          # extend the SAME generic on ResummedSeriesData
import Symbolics

export PostNewtonianData, pn_series, pn_start_order, pn_max_order,
    postnewtonian_coefficient, ResummedSeriesData, resummed_series,
    postnewtonian_expansion

# ----------------------------------------------------------------------------
# PostNewtonianData wrapper + order metadata
# ----------------------------------------------------------------------------

"""
    PostNewtonianData(data::Dict{String,Any})

Wrapper holding a loaded PN series together with its metadata. `data` mirrors the upstream
Association, with `data["Series"]::PNSeries` carrying the expansion and the remaining keys
(`"Name"`, `"Description"`, `"Authors"`, `"References"`) the metadata. Index with `pnd[key]`;
`keys(pnd)` lists the available fields.
"""
struct PostNewtonianData
    data::Dict{String,Any}
end

Base.getindex(pnd::PostNewtonianData, key::AbstractString) = pnd.data[key]
Base.keys(pnd::PostNewtonianData) = keys(pnd.data)

"""
    pn_series(pnd) -> PNSeries

The `"Series"` payload of `pnd`.
"""
pn_series(pnd::PostNewtonianData) = pnd.data["Series"]::PNSeries

"""
    pn_start_order(pnd) -> Rational

Starting PN order of the expansion, `nmin // den`. Upstream: `Series[[4]]/den`.
"""
pn_start_order(pnd::PostNewtonianData) = series_min_order(pn_series(pnd))

"""
    pn_max_order(pnd) -> Rational

Maximum (last retained) PN order, `(nmax - 1) // den`. Upstream: `(Series[[5]]-1)/den` — note the
`-1`, which makes this distinct from `Series.series_max_order` (= `nmax // den`).
"""
function pn_max_order(pnd::PostNewtonianData)
    s = pn_series(pnd)
    return (s.nmax - 1) // s.den
end

function Base.show(io::IO, pnd::PostNewtonianData)
    name = get(pnd.data, "Name", "PostNewtonianData")
    print(io, "PostNewtonianData(\"", name, "\"; starting order ", pn_start_order(pnd),
        ", maximum order ", pn_max_order(pnd), ")")
end

# ----------------------------------------------------------------------------
# postnewtonian_coefficient (y^n Log[y]^nL extraction)
# ----------------------------------------------------------------------------

"""
    postnewtonian_coefficient(pnd, n) -> Num

The `y^n` coefficient of the expansion (full, including any `Log[y]` sentinel terms). `n` is a real
PN order; `n*den` must land on the integer grid, else the coefficient is `0`. Upstream:
`SeriesCoefficient[pnd["Series"], n]`.
"""
function postnewtonian_coefficient(pnd::PostNewtonianData, n)
    s = pn_series(pnd)
    num = n * s.den
    isinteger(num) || return Symbolics.Num(0)
    i = Int(num) - s.nmin + 1
    return (1 <= i <= length(s.coeffs)) ? s.coeffs[i] : Symbolics.Num(0)
end

"""
    postnewtonian_coefficient(pnd, n, nL) -> Num

The `y^n Log[y]^nL` coefficient. `nL == 0` strips the log sentinel (`Logvar -> 0`); `nL >= 1` extracts
the coefficient of `Logvar^nL`, where `Logvar = Symbol("Log", series variable)`.
"""
function postnewtonian_coefficient(pnd::PostNewtonianData, n, nL)
    s = pn_series(pnd)
    c = postnewtonian_coefficient(pnd, n)
    logvar = Symbolics.variable(Symbol("Log", s.var))
    nL == 0 && return Symbolics.substitute(c, Dict(logvar => 0))
    return Symbolics.Num(Symbolics.coeff(c, logvar^nL))
end

# ----------------------------------------------------------------------------
# ResummedSeriesData (sf*sd resummed series)
# ----------------------------------------------------------------------------

"""
    ResummedSeriesData(sf::Num, sd::PNSeries)

A resummed power series `sf * sd`: a singular prefactor `sf` (a `Symbolics.Num`, typically a monomial
like `p^-5`) times a truncated series `sd`. Build via [`resummed_series`](@ref), which applies the
upstream collapse rule. Index with `rsd["SingularFactor"]` / `rsd["SeriesData"]`.
"""
struct ResummedSeriesData
    sf::Symbolics.Num
    sd::PNSeries
end

function Base.getindex(rsd::ResummedSeriesData, key::AbstractString)
    key == "SingularFactor" && return rsd.sf
    key == "SeriesData" && return rsd.sd
    throw(KeyError(key))
end

"""
    resummed_series(sf, sd) -> ResummedSeriesData or Num

Construct a resummed series `sf*sd`. If `sd` is a `PNSeries`, return the [`ResummedSeriesData`](@ref)
wrapper; otherwise collapse to `sf*sd` (a `Num`), matching the upstream conditional rewrite
`ResummedSeriesData[sf, sd] /; !MatchQ[sd, _SeriesData] := sf*sd`.
"""
function resummed_series(sf, sd)
    sd isa PNSeries && return ResummedSeriesData(Symbolics.Num(sf), sd)
    return Symbolics.Num(sf) * sd
end

"""
    Normal(rsd::ResummedSeriesData) -> Num

Materialize the resummed series to `sf * Normal(sd)` (the singular prefactor times the Laurent
expression of the series, with its O[] term dropped). Extends `Series.Normal`.
"""
Normal(rsd::ResummedSeriesData) = rsd.sf * Normal(rsd.sd)

# ----------------------------------------------------------------------------
# postnewtonian_expansion (load-by-name entry)
# ----------------------------------------------------------------------------

"""
    postnewtonian_expansion(quantity=""; root=SERIES_DATA_DIR)

Locate and load stored PN series matching `quantity` (a substring filter string, or a vector of them
applied in sequence — leading `!` negates, case-insensitive). On a unique match, load that file and
return a [`PostNewtonianData`](@ref); otherwise return the `Vector{String}` of matching catalog keys
(so an empty/ambiguous query browses). Upstream: `PostNewtonianExpansion[quantity]`.
"""
function postnewtonian_expansion(quantity = ""; root::AbstractString = SERIES_DATA_DIR)
    sublist = pn_browse(quantity, root)
    if length(sublist) == 1
        return PostNewtonianData(parse_pn_data_file(string(root, sublist[1], ".m")))
    end
    return sublist
end

end # module SelfForce
