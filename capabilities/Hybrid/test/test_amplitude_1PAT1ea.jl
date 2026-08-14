# Characterization of the 1PAT1ea linear-in-spin / δm strain-amplitude model (Hybrid.Models, §4.9).
#
# GROUNDING (no WaSABI live dep — ANALYTIC closed forms + synthetic complex interpolant fixtures):
#   The frequency-domain amplitude h1PAAmpΩ[l,m] is rebuilt from the literal 1PAT1ea.m channels:
#     hS1 = -2(-Za)/(i m √(1/r0³))²,  hS2 = norm·C2σ (even l+m) / -i norm·C2σ (odd l+m) with
#       norm = √((l-1)l(l+1)(l+2))/2,  hδm = √(1/r0³)(-2 r0^{5/2}/3) h0PAamp',
#     c = ν h0 + ν χt1 hS1 + ν χt2 hS2,  even m: c+ν²(h1ν+δm hδm),  odd m: √(1-4ν)(c+ν²(2h0+h1ν+δm hδm)).
#   Za / C2σ are r-linear interpolants (values known in closed form); h0/h1ν/h0PAamp' come from
#   the already-grounded 1PAT1 base.

using Test
using Hybrid.Models

_zampea(l, m, r) = complex(0.01 + 0.001 * r + 0.0001 * (l + m),
                         -0.002 + 0.0005 * r - 0.0003 * (l - m))
_zafun(l, m, r) = complex(0.005 + 0.0007 * r + 0.0002 * (l - m),
                          0.003 - 0.0004 * r + 0.0001 * (l + m))
_csfun(l, m, r) = complex(0.008 - 0.0003 * r + 0.0001 * l,
                          -0.001 + 0.0006 * r - 0.0002 * m)

@testset "Models.amplitude_1PAT1ea" begin
    rs = collect(range(6.2, 32.0; length = 14))
    Z = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    Z2 = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    Za = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    C2σ = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    α, β = 0.1 + 0.2im, -1.5 + 0.1im
    for l in 2:5, m in 1:l
        Z[(l, m)] = interp1_order8(rs, [_zampea(l, m, r) for r in rs])
        Z2[(l, m)] = interp1_order8(log.(rs), [(α + 0.01 * (l + m)) + β * log(r) for r in rs]; order = 3)
        Za[(l, m)] = interp1_order8(rs, [_zafun(l, m, r) for r in rs])
        C2σ[(l, m)] = interp1_order8(rs, [_csfun(l, m, r) for r in rs])
    end
    base1pa = Amp1PAT1(Amp0PASchwarz(Z), Z2)
    model = Amp1PAT1ea(base1pa, Za, C2σ)

    @testset "h1PA_ampΩ_ea literal channel assembly" begin
        ν, χt1, χt2, δm = 1.0e-2, 0.3, -0.2, 0.4
        for (l, m, Ω) in ((2, 2, 0.03), (3, 1, 0.025), (4, 3, 0.02), (5, 2, 0.018))
            r0 = Ω^(-2 / 3)
            h0 = h0PA_amp(base1pa, l, m, r0)
            h1ν = h1PA_ampν(base1pa, l, m, r0)
            Zav = _zafun(l, m, r0); Csv = _csfun(l, m, r0)
            hS1 = (-2 * (-Zav)) / (im * m * sqrt(1 / r0^3))^2
            norm = sqrt((l - 1) * l * (l + 1) * (l + 2)) / 2
            hS2 = iseven(l + m) ? norm * Csv : -im * norm * Csv
            hδm = sqrt(1 / r0^3) * (-2 * r0^(5 / 2) / 3) * h0PA_amp_deriv(base1pa, l, m, r0)
            c = ν * h0 + ν * χt1 * hS1 + ν * χt2 * hS2
            exp_amp = iseven(m) ? (c + ν^2 * (h1ν + δm * hδm)) :
                      sqrt(1 - 4ν) * (c + ν^2 * (2 * h0 + h1ν + δm * hδm))
            @test h1PA_ampΩ_ea(model, l, m, Ω, ν, χt1, χt2, δm) ≈ exp_amp rtol = 1e-9
        end
    end

    @testset "mode_keys" begin
        @test mode_keys_1PAT1ea() == [(l, m) for l in 2:5 for m in 1:l]
    end
end
