# 1PAT1ea amplitude model (`AmplitudeModels/1PAT1ea.m`): the 1PA Schwarzschild-circular strain
# amplitudes extended with linear-in-spin channels — secondary spin S2 (`C2σ`, RWZ amplitudes
# with the i-factor for odd l+m and the √((l-1)l(l+1)(l+2))/2 normalization) and primary spin S1
# (`Za`, "slow-χ1" Teukolsky amplitudes) — plus the δm (mass-loss) channel. Built on top of the
# 1PAT1 base ([`Amp1PAT1`](@ref): 0PA `Z_I` + 2SF `Z2_I`).

"""
    Amp1PAT1ea

1PAT1ea amplitude model: the [`Amp1PAT1`](@ref) base (1SF `Z_I`, 2SF `Z2_I`) plus the primary-
spin `Za` (order-8) and secondary-spin `C2σ` (order-8) complex amplitude interpolants,
`l = 2..5`, `m = 1..l`. Build with [`amplitudes_1PAT1ea`](@ref).
"""
struct Amp1PAT1ea
    base::Amp1PAT1
    Za::Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}    # primary spin S1
    C2σ::Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}   # secondary spin S2 (RWZ)
end

"""
    amplitudes_1PAT1ea(amp_datadir) -> Amp1PAT1ea

Load the 1PAT1 base tables plus the primary-spin `1SF/Linear_Kerr_Circ/Slow_chi1_*` and
secondary-spin `Leading_S2/Schwarz_Circ/chi2_2SF_RWZamp*` tables under `amp_datadir`
(`AmplitudeModels/sf_amp_data`), all order-8.
"""
function amplitudes_1PAT1ea(amp_datadir::AbstractString)
    base = amplitudes_1PAT1(amp_datadir)
    ds1 = joinpath(amp_datadir, "1SF", "Linear_Kerr_Circ")
    ds2 = joinpath(amp_datadir, "Leading_S2", "Schwarz_Circ")
    Za = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    C2σ = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    for l in 2:5, m in 1:l
        Za[(l, m)] = interp1_order8(
            read_wl_table_complex(joinpath(ds1, "Slow_chi1_TeukampSchwarzCirc$(l)$(m).m")))
        C2σ[(l, m)] = interp1_order8(
            read_wl_table_complex(joinpath(ds2, "chi2_2SF_RWZampSchwarzCirc$(l)$(m).m")))
    end
    return Amp1PAT1ea(base, Za, C2σ)
end

# Primary-spin amplitude Za_I[l,m] (neg-m: (-1)^l conj), `1PAT1ea.m`.
function _Za_I(model::Amp1PAT1ea, l::Integer, m::Integer, r0::Real)
    m < 0 && return (-1)^l * conj(_Za_I(model, l, -m, r0))
    return model.Za[(l, m)](r0)
end

# Secondary-spin RWZ amplitude C2σ_I[l,m] (neg-m: (-1)^l conj), `1PAT1ea.m`.
function _C2σ_I(model::Amp1PAT1ea, l::Integer, m::Integer, r0::Real)
    m < 0 && return (-1)^l * conj(_C2σ_I(model, l, -m, r0))
    return model.C2σ[(l, m)](r0)
end

# S1 strain channel h1PAampa[l,m] = -2(-Za_I)/(i m √(1/r0³))² (M=1), `1PAT1ea.m`.
function _h1PA_amp_S1(model::Amp1PAT1ea, l::Integer, m::Integer, r0::Real)
    Za = _Za_I(model, l, m, r0)
    return (-2 * (-Za)) / (im * m * sqrt(1 / r0^3))^2
end

# S2 strain channel h1PASpinamp[l,m] = norm·C2σ_I (even l+m) or -i·norm·C2σ_I (odd l+m), with
# norm = √((l-1)l(l+1)(l+2))/2 (`1PAT1ea.m`).
function _h1PA_amp_S2(model::Amp1PAT1ea, l::Integer, m::Integer, r0::Real)
    norm = sqrt((l - 1) * l * (l + 1) * (l + 2)) / 2
    c = _C2σ_I(model, l, m, r0)
    return iseven(l + m) ? norm * c : -im * norm * c
end

# δm channel h1PAampνδm[l,m] = √(1/r0³)(-2 r0^(5/2)/3) h0PAamp'[l,m] (M=1), `1PAT1ea.m`.
_h1PA_amp_δm(model::Amp1PAT1ea, l::Integer, m::Integer, r0::Real) =
    sqrt(1 / r0^3) * (-2 * r0^(5 / 2) / 3) * h0PA_amp_deriv(model.base, l, m, r0)

"""
    h1PA_ampΩ_ea(model::Amp1PAT1ea, l, m, Ω, ν, χt1, χt2, δm) -> ComplexF64

Frequency-domain 1PAT1ea strain amplitude `h1PAAmpΩ[l,m]` (M=1; `1PAT1ea.m`). Common spin part
`c = ν h0 + ν χt1 hS1 + ν χt2 hS2`; even `m`: `M(c + ν²(h1ν + δm hδm))`; odd `m`:
`√(1-4ν) M(c + ν²(2 h0 + h1ν + δm hδm))`, all evaluated at `r0 = Ω^{-2/3}`.
"""
function h1PA_ampΩ_ea(model::Amp1PAT1ea, l::Integer, m::Integer, Ω::Real, ν::Real,
                      χt1::Real, χt2::Real, δm::Real)
    M = 1.0
    r0 = Ω^(-2 / 3)
    h0 = h0PA_amp(model.base, l, m, r0)
    h1ν = h1PA_ampν(model.base, l, m, r0)
    hS1 = _h1PA_amp_S1(model, l, m, r0)
    hS2 = _h1PA_amp_S2(model, l, m, r0)
    hδm = _h1PA_amp_δm(model, l, m, r0)
    c = ν * h0 + ν * χt1 * hS1 + ν * χt2 * hS2
    if iseven(m)
        return M * (c + ν^2 * (h1ν + δm * hδm))
    else
        return sqrt(1 - 4ν) * M * (c + ν^2 * (2 * h0 + h1ν + δm * hδm))
    end
end

"""
    mode_keys_1PAT1ea() -> Vector{Tuple{Int,Int}}

The `(l, m)` modes of the 1PAT1ea model, `l = 2..5`, `m = 1..l`.
"""
mode_keys_1PAT1ea() = [(l, m) for l in 2:5 for m in 1:l]
