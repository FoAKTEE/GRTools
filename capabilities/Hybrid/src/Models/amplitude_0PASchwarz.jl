# 0PASchwarz amplitude model (AmplitudeModels/0PASchwarz.m): adiabatic (0PA) waveform strain
# amplitudes for circular Schwarzschild orbits, built from the Teukolsky amplitude tables
# `Z_I[l,m]` (order-8 interpolated, `l = 2..5`, `m = 1..l`). Negative `m` follows from the
# closed-form symmetry `Z[l,-m] = (-1)^l conj(Z[l,m])`.

"""
    Amp0PASchwarz

0PASchwarz amplitude model: a table of order-8 [`Interp1`](@ref) interpolants for the complex
Teukolsky amplitudes `Z_I[l,m](r0)`, `l = 2..5`, `m = 1..l`. Build with
[`amplitudes_0PASchwarz`](@ref).
"""
struct Amp0PASchwarz
    Z::Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}
end

"""
    amplitudes_0PASchwarz(datadir) -> Amp0PASchwarz

Load the `1SFTeukampSchwarzCirc<l><m>.m` complex amplitude tables from `datadir`
(`sf_amp_data/1SF/Schwarz_Circ`) and build the order-8 `Z_I[l,m]` interpolants for `l = 2..5`,
`m = 1..l`.
"""
function amplitudes_0PASchwarz(datadir::AbstractString)
    Z = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    for l in 2:5, m in 1:l
        tbl = read_wl_table_complex(joinpath(datadir, "1SFTeukampSchwarzCirc$(l)$(m).m"))
        Z[(l, m)] = interp1_order8(tbl)
    end
    return Amp0PASchwarz(Z)
end

"""
    Z_I(model, l, m, r0)

Interpolated Teukolsky amplitude `Z_I[l,m](r0)`. Negative `m` uses `Z[l,-m] = (-1)^l
conj(Z[l,m])` (0PASchwarz.m).
"""
function Z_I(model::Amp0PASchwarz, l::Integer, m::Integer, r0::Real)
    if m < 0
        return (-1)^l * conj(model.Z[(l, -m)](r0))
    end
    return model.Z[(l, m)](r0)
end

"""
    Z_I_deriv(model, l, m, r0)

`d/dr0` of the interpolated Teukolsky amplitude `Z_I[l,m](r0)`. Negative `m` uses
`Z'[l,-m] = (-1)^l conj(Z'[l,m])`. Used by the 1PAT1 amplitude's `h0PAamp'[r0]` term.
"""
function Z_I_deriv(model::Amp0PASchwarz, l::Integer, m::Integer, r0::Real)
    if m < 0
        return (-1)^l * conj(derivative(model.Z[(l, -m)], r0))
    end
    return derivative(model.Z[(l, m)], r0)
end

"""
    h0PA_amp(model, l, m, r0)

Strain amplitude `h0PAamp[l,m][r0] = (-2(-Z_I[l,m][r0]))/(i m √(M/r0³))²` (M = 1).
"""
function h0PA_amp(model::Amp0PASchwarz, l::Integer, m::Integer, r0::Real)
    M = one(float(r0))
    Z = Z_I(model, l, m, r0)
    return (-2 * (-Z)) / (im * m * sqrt(M / r0^3))^2
end

"""
    h0PA_ampΩ(model, l, m, Ω, ν)

Strain amplitude as a function of orbital frequency `Ω`, scaled by mass ratio `ν`:
`h0PAAmpΩ[l,m] = ν M h0PAamp[l,m][Ω^{-2/3}]` (M = 1).
"""
h0PA_ampΩ(model::Amp0PASchwarz, l::Integer, m::Integer, Ω::Real, ν::Real) =
    ν * one(float(Ω)) * h0PA_amp(model, l, m, Ω^(-2 / 3))

"""
    mode_keys_0PASchwarz() -> Vector{Tuple{Int,Int}}

The `(l, m)` modes of the 0PASchwarz model, `l = 2..5`, `m = 1..l`.
"""
mode_keys_0PASchwarz() = [(l, m) for l in 2:5 for m in 1:l]
