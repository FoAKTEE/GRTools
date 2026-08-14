# Characterization of the 1PAT1 first-post-adiabatic strain-amplitude model (Hybrid.Models, §4.9).
#
# GROUNDING (no WaSABI live dep — ANALYTIC closed forms + synthetic complex interpolant fixtures):
#   * Z2_I == Exp[interp(Log r0)]: with the 2SF table set to a complex log-linear L(u)=α+βu the
#     order-3 interpolant is exact, so Z2_I == exp(α+β Log r0) in closed form; neg-m symmetry.
#   * h0PA_amp_deriv == -2/m² (Z'(r0) r0³ + 3 r0² Z) (M=1) on the r-linear base Z table.
#   * h1PA_amp == (-2 Z2)/(i m √(1/r0³))² (M=1).
#   * h1PA_ampν == h0 + h1 + (2/3)r0(1-(1/r0)/√(1-3/r0)) h0PAamp' (M=1).
#   * h1PA_ampΩ even/odd-m parity split (odd m: √(1-4ν) factor + extra 2ν² h0).

using Test
using Hybrid.Models

_zamp1(l, m, r) = complex(0.01 + 0.001 * r + 0.0001 * (l + m),
                         -0.002 + 0.0005 * r - 0.0003 * (l - m))
_dzamp1(l, m, r) = complex(0.001, 0.0005)   # d/dr of _zamp1 (r-linear part)

@testset "Models.amplitude_1PAT1" begin
    rs = collect(range(6.2, 32.0; length = 14))
    Z = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    for l in 2:5, m in 1:l
        Z[(l, m)] = interp1_order8(rs, [_zamp1(l, m, r) for r in rs])
    end
    base = Amp0PASchwarz(Z)

    # 2SF table: complex log-linear L(u) = α + β u over u = Log r0 (order-3 exact).
    α, β = 0.1 + 0.2im, -1.5 + 0.1im
    Z2 = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    for l in 2:5, m in 1:l
        αlm = α + 0.01 * (l + m)        # vary per mode but keep log-linear
        Z2[(l, m)] = interp1_order8(log.(rs), [αlm + β * log(r) for r in rs]; order = 3)
    end
    z2exp(l, m, r) = exp((α + 0.01 * (l + m)) + β * log(r))
    model = Amp1PAT1(base, Z2)

    @testset "Z2_I closed form + neg-m symmetry" begin
        for (l, m, r0) in ((2, 2, 12.0), (3, 1, 9.0), (5, 4, 20.0))
            @test Z2_I(model, l, m, r0) ≈ z2exp(l, m, r0) rtol = 1e-10
            @test Z2_I(model, l, -m, r0) == (-1)^l * conj(Z2_I(model, l, m, r0))
        end
    end

    @testset "h0PA_amp_deriv closed form" begin
        for (l, m, r0) in ((2, 2, 12.0), (3, 2, 15.0), (5, 3, 20.0))
            Z0 = _zamp1(l, m, r0); dZ = _dzamp1(l, m, r0)
            @test h0PA_amp_deriv(model, l, m, r0) ≈ -2 / m^2 * (dZ * r0^3 + 3 * r0^2 * Z0) rtol = 1e-9
        end
    end

    @testset "h1PA_amp / h1PA_ampν closed forms" begin
        for (l, m, r0) in ((2, 2, 12.0), (3, 1, 9.0), (4, 4, 18.0))
            Z2v = z2exp(l, m, r0)
            @test h1PA_amp(model, l, m, r0) ≈ (-2 * Z2v) / (im * m * sqrt(1 / r0^3))^2 rtol = 1e-10
            extra = (2 / 3) * r0 * (1 - (1 / r0) / sqrt(1 - 3 / r0)) * h0PA_amp_deriv(model, l, m, r0)
            @test h1PA_ampν(model, l, m, r0) ≈
                  h0PA_amp(model, l, m, r0) + h1PA_amp(model, l, m, r0) + extra rtol = 1e-10
        end
    end

    @testset "h1PA_ampΩ even/odd-m parity" begin
        ν = 1.0e-2
        for (l, m, Ω) in ((2, 2, 0.03), (3, 3, 0.025), (4, 2, 0.02), (5, 1, 0.018))
            r0 = Ω^(-2 / 3)
            h0 = h0PA_amp(model, l, m, r0); h1ν = h1PA_ampν(model, l, m, r0)
            exp_amp = iseven(m) ? (ν * h0 + ν^2 * h1ν) :
                      sqrt(1 - 4ν) * (ν * h0 + ν^2 * (2 * h0 + h1ν))
            @test h1PA_ampΩ(model, l, m, Ω, ν) ≈ exp_amp rtol = 1e-10
        end
    end

    @testset "mode_keys" begin
        @test mode_keys_1PAT1() == [(l, m) for l in 2:5 for m in 1:l]
    end
end
