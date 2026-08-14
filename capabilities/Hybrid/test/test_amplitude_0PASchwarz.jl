# Characterization of the 0PASchwarz waveform-strain amplitude model (Hybrid.Models, §4.9).
#
# GROUNDING (no WaSABI live dep — ANALYTIC + synthetic complex WL-table fixture):
#   * Z_I / Z_I_deriv negative-m symmetry Z[l,-m]==(-1)^l conj(Z[l,m]): exact identity on any
#     table (and likewise for the derivative).
#   * h0PA_amp vs the literal closed form (-2(-Z))/(i m √(M/r0³))² == -2 Z r0³/m² (M=1).
#   * h0PA_ampΩ == ν M h0PA_amp(·, Ω^{-2/3}).
#   * mode_keys_0PASchwarz == [(l,m) for l in 2:5 for m in 1:l] literal.
#   * Data-load path: synthetic complex `.m` tables round-tripped via Data.read_wl_table_complex.

using Test
using Hybrid.Models

# Synthetic complex amplitude Z_lm(r0) — degree-1 complex polynomial (exact under order-8).
_zfun(l, m, r) = complex(0.01 + 0.001 * r + 0.0001 * (l + m),
                         -0.002 + 0.0005 * r - 0.0003 * (l - m))

function _build_amp(rs)
    Z = Dict{Tuple{Int,Int},Interp1{Float64,ComplexF64}}()
    for l in 2:5, m in 1:l
        Z[(l, m)] = interp1_order8(rs, [_zfun(l, m, r) for r in rs])
    end
    return Amp0PASchwarz(Z)
end

# Write a flat WL complex `{ {x, re ± im*I}, ... }` table.
function _write_wl_complex(path, rs, vals)
    parts = String[]
    for (r, z) in zip(rs, vals)
        re, im = real(z), imag(z)
        imstr = im >= 0 ? "+ $(im)*I" : "- $(-im)*I"
        push!(parts, "{$(r), $(re) $(imstr)}")
    end
    write(path, "{" * join(parts, ", ") * "}")
    return path
end

@testset "Models.amplitude_0PASchwarz" begin
    rs = collect(range(6.2, 32.0; length = 14))
    model = _build_amp(rs)

    @testset "mode_keys" begin
        @test mode_keys_0PASchwarz() == [(l, m) for l in 2:5 for m in 1:l]
    end

    @testset "Z_I value + negative-m symmetry" begin
        for (l, m, r0) in ((2, 2, 12.0), (3, 1, 9.0), (4, 3, 20.0), (5, 2, 15.0))
            @test Z_I(model, l, m, r0) ≈ _zfun(l, m, r0) rtol=1e-10
            # exact symmetry identity
            @test Z_I(model, l, -m, r0) == (-1)^l * conj(Z_I(model, l, m, r0))
            @test Z_I_deriv(model, l, -m, r0) == (-1)^l * conj(Z_I_deriv(model, l, m, r0))
        end
    end

    @testset "h0PA_amp closed form" begin
        for (l, m, r0) in ((2, 2, 12.0), (3, 1, 9.0), (5, 4, 18.0))
            Z = Z_I(model, l, m, r0)
            literal = (-2 * (-Z)) / (im * m * sqrt(1.0 / r0^3))^2
            @test h0PA_amp(model, l, m, r0) ≈ literal rtol=1e-12
            @test h0PA_amp(model, l, m, r0) ≈ -2 * Z * r0^3 / m^2 rtol=1e-12  # simplified
        end
    end

    @testset "h0PA_ampΩ frequency form" begin
        ν = 1.0e-3
        for (l, m, Ω) in ((2, 2, 0.03), (3, 2, 0.025), (4, 4, 0.02))
            r0 = Ω^(-2 / 3)
            @test h0PA_ampΩ(model, l, m, Ω, ν) ≈ ν * 1.0 * h0PA_amp(model, l, m, r0) rtol=1e-12
        end
    end

    @testset "Data-load round-trip (synthetic complex WL tables)" begin
        dir = mktempdir()
        for l in 2:5, m in 1:l
            vals = [_zfun(l, m, r) for r in rs]
            _write_wl_complex(joinpath(dir, "1SFTeukampSchwarzCirc$(l)$(m).m"), rs, vals)
        end
        loaded = amplitudes_0PASchwarz(dir)
        for (l, m, r0) in ((2, 2, 12.0), (3, 3, 15.0), (5, 5, 25.0))
            @test Z_I(loaded, l, m, r0) ≈ _zfun(l, m, r0) rtol=1e-9
        end
    end
end
