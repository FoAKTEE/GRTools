# 1PAT1 amplitude model (`AmplitudeModels/1PAT1.m`): first post-adiabatic (1PA) strain
# amplitudes for circular Schwarzschild orbits. Built from the 1SF Teukolsky amplitudes
# `Z_I[l,m]` (order-8; reused from the 0PA model) and the 2SF amplitudes
# `Z2_I[l,m] = Exp[interp(Log r0)]` (order-3, complex). The frequency-domain amplitude
# `h1PAAmpΩ[l,m]` carries an even/odd-m parity split: odd m gains a `√(1-4ν)` factor and an
# extra `2ν² h0PAamp` term (`1PAT1.m`).

"""
    Amp1PAT1

1PAT1 amplitude model: the 0PA amplitude tables ([`Amp0PASchwarz`](@ref), 1SF `Z_I`) plus the
2SF amplitude interpolants `Z2_I[l,m]` (order-3 over `{Log r0, Log Z}`, complex), `l = 2..5`,
`m = 1..l`. Build with [`amplitudes_1PAT1`](@ref).
"""
struct Amp1PAT1
    base::Amp0PASchwarz
    Z2::Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}
end

"""
    amplitudes_1PAT1(amp_datadir) -> Amp1PAT1

Load the 1SF (`1SF/Schwarz_Circ`, order-8) and 2SF (`2SF/Schwarz_Circ`, order-3) Teukolsky
amplitude tables under `amp_datadir` (`AmplitudeModels/sf_amp_data`) and build the model.
"""
function amplitudes_1PAT1(amp_datadir::AbstractString)
    base = amplitudes_0PASchwarz(joinpath(amp_datadir, "1SF", "Schwarz_Circ"))
    d2 = joinpath(amp_datadir, "2SF", "Schwarz_Circ")
    Z2 = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    for l in 2:5, m in 1:l
        tbl = read_wl_table_complex(joinpath(d2, "2SFTeukampSchwarzCirc$(l)$(m).m"))
        Z2[(l, m)] = interp1_order8(tbl; order = 3)
    end
    return Amp1PAT1(base, Z2)
end

"""
    Z2_I(model::Amp1PAT1, l, m, r0) -> ComplexF64

2SF Teukolsky amplitude `Z2_I[l,m](r0) = Exp[interp(Log r0)]`. Negative `m` uses
`Z2[l,-m] = (-1)^l conj(Z2[l,m])` (`1PAT1.m`).
"""
function Z2_I(model::Amp1PAT1, l::Integer, m::Integer, r0::Real)
    if m < 0
        return (-1)^l * conj(Z2_I(model, l, -m, r0))
    end
    return exp(model.Z2[(l, m)](log(r0)))
end

"""
    h0PA_amp(model::Amp1PAT1, l, m, r0) -> ComplexF64

0PA strain amplitude `h0PAamp[l,m][r0]` (delegates to the 1SF [`Amp0PASchwarz`](@ref) base).
"""
h0PA_amp(model::Amp1PAT1, l::Integer, m::Integer, r0::Real) = h0PA_amp(model.base, l, m, r0)

"""
    h0PA_amp_deriv(model, l, m, r0) -> ComplexF64

`d/dr0` of `h0PAamp[l,m][r0]`. With `h0PAamp = -2 Z_I r0³/m²` (M = 1), this is
`-2/m² (Z_I'(r0) r0³ + 3 r0² Z_I(r0))`. Reuses [`Z_I`](@ref) and [`Z_I_deriv`](@ref).
"""
function h0PA_amp_deriv(base::Amp0PASchwarz, l::Integer, m::Integer, r0::Real)
    Z = Z_I(base, l, m, r0)
    dZ = Z_I_deriv(base, l, m, r0)
    return -2 / m^2 * (dZ * r0^3 + 3 * r0^2 * Z)
end
h0PA_amp_deriv(model::Amp1PAT1, l::Integer, m::Integer, r0::Real) =
    h0PA_amp_deriv(model.base, l, m, r0)

"""
    h1PA_amp(model::Amp1PAT1, l, m, r0) -> ComplexF64

1PA strain amplitude `h1PAamp[l,m][r0] = (-2 Z2_I[l,m][r0])/(i m √(M/r0³))²` (M = 1; `1PAT1.m`).
"""
function h1PA_amp(model::Amp1PAT1, l::Integer, m::Integer, r0::Real)
    M = 1.0
    Z2 = Z2_I(model, l, m, r0)
    return (-2 * Z2) / (im * m * sqrt(M / r0^3))^2
end

"""
    h1PA_ampν(model::Amp1PAT1, l, m, r0) -> ComplexF64

`h1PAampν[l,m][r0] = h0PAamp + h1PAamp + (2/3) r0 (1 - (M/r0)/√(1-3M/r0)) h0PAamp'[r0]`
(M = 1; `1PAT1.m`).
"""
function h1PA_ampν(model::Amp1PAT1, l::Integer, m::Integer, r0::Real)
    M = 1.0
    extra = (2 / 3) * r0 / M * (1 - (M / r0) / sqrt(1 - 3M / r0)) *
            h0PA_amp_deriv(model, l, m, r0)
    return h0PA_amp(model, l, m, r0) + h1PA_amp(model, l, m, r0) + extra
end

"""
    h1PA_ampΩ(model::Amp1PAT1, l, m, Ω, ν) -> ComplexF64

Frequency-domain 1PA strain amplitude `h1PAAmpΩ[l,m]` at orbital frequency `Ω` and mass ratio
`ν` (M = 1; `1PAT1.m`). Even `m`: `M(ν h0 + ν² h1ν)`. Odd `m`:
`√(1-4ν) M(ν h0 + ν²(2 h0 + h1ν))`, where `h0 = h0PAamp[Ω^{-2/3}]`, `h1ν = h1PAampν[Ω^{-2/3}]`.
"""
function h1PA_ampΩ(model::Amp1PAT1, l::Integer, m::Integer, Ω::Real, ν::Real)
    M = 1.0
    r0 = Ω^(-2 / 3)
    h0 = h0PA_amp(model, l, m, r0)
    h1ν = h1PA_ampν(model, l, m, r0)
    if iseven(m)
        return M * (ν * h0 + ν^2 * h1ν)
    else
        return sqrt(1 - 4ν) * M * (ν * h0 + ν^2 * (2 * h0 + h1ν))
    end
end

"""
    mode_keys_1PAT1() -> Vector{Tuple{Int,Int}}

The `(l, m)` modes of the 1PAT1 model, `l = 2..5`, `m = 1..l`.
"""
mode_keys_1PAT1() = [(l, m) for l in 2:5 for m in 1:l]
