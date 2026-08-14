# Characterization of the 0PASchwarz adiabatic forcing model (Hybrid.Models, §4.9).
#
# GROUNDING (no WaSABI live dep — ANALYTIC + tiny synthetic WL-table fixture):
#   * forcing prefactor / forcing_F0 / omega_dot vs the literal closed form
#     3(1-3/r0)^{3/2}√(1/r0)/(1-6/r0)·F1(r0) (0PASchwarz.m).
#   * flux_1SF == FEI(r0)+FEH(r0): additive identity on Interp1s built from known polynomial
#     tables (each interpolant reproduces its degree-≤8 polynomial exactly).
#   * inspiral_rhs!: du[1]==omega_dot, du[2]==Ω/M exactly; out-of-place == in-place.
#   * stop_r0 == 6.01 literal.
#   * Data-load path: write synthetic WL `.m` flux tables, load via forcing_0PASchwarz
#     (Data.read_wl_table), assert the interpolants round-trip the table node values.

using Test
using Hybrid.Models

# Write a flat WL `{ {x,y}, ... }` table to `path`.
function _write_wl_table(path, xs, ys)
    body = join(("{$(x), $(y)}" for (x, y) in zip(xs, ys)), ", ")
    write(path, "{" * body * "}")
    return path
end

@testset "Models.inspiral_0PASchwarz" begin
    # Synthetic degree-≤8 polynomial flux tables (exactly reproduced by order-8 Interp1).
    fI(r) = 0.02 + 1.0e-3 * r - 5.0e-5 * r^2 + 1.0e-6 * r^3
    fH(r) = 0.005 - 2.0e-4 * r + 3.0e-6 * r^2
    rs = collect(range(6.2, 32.0; length = 14))
    mI = interp1_order8(rs, fI.(rs))
    mH = interp1_order8(rs, fH.(rs))
    model = Forcing0PASchwarz(mH, mI)

    @testset "flux_1SF additive identity" begin
        for r0 in (8.0, 12.5, 20.0, 30.0)
            @test flux_1SF(model, r0) ≈ mI(r0) + mH(r0) rtol=0
            @test flux_1SF(model, r0) ≈ fI(r0) + fH(r0) rtol=1e-10
        end
    end

    @testset "forcing prefactor / F0 / omega_dot closed form" begin
        for r0 in (8.0, 12.5, 20.0)
            F1 = fI(r0) + fH(r0)
            pref = 3 * (1 - 3 / r0)^(3 / 2) * sqrt(1 / r0) / (1 - 6 / r0)
            @test forcing_F0(model, r0) ≈ pref * F1 rtol=1e-10
            ν, M = 1.0e-4, 1.0
            @test omega_dot(model, r0, ν, M) ≈ (ν / M) * pref * F1 rtol=1e-10
        end
    end

    @testset "stop_r0" begin
        @test stop_r0(model, 1.0e-3) == 6.01
        @test stop_r0(model, 1.0e-6) == 6.01
    end

    @testset "inspiral_rhs (in-place == out-of-place == defn)" begin
        ν, M = 2.0e-4, 1.0
        Ω, φ = 0.03, 1.2
        r0 = Ω^(-2 / 3)
        du = zeros(2)
        inspiral_rhs!(du, [Ω, φ], (model, ν, M), 0.0)
        @test du[1] == omega_dot(model, r0, ν, M)
        @test du[2] == Ω / M
        oop = inspiral_rhs(model, Ω, φ, ν, M)
        @test oop[1] == du[1]
        @test oop[2] == du[2]
    end

    @testset "Data-load round-trip (synthetic WL tables)" begin
        dir = mktempdir()
        _write_wl_table(joinpath(dir, "1SFCircShwarzDotEInf.m"), rs, fI.(rs))
        _write_wl_table(joinpath(dir, "1SFCircShwarzDotEHorizon.m"), rs, fH.(rs))
        loaded = forcing_0PASchwarz(dir)
        # Interpolants reproduce the table node values to machine precision.
        for (k, r) in enumerate(rs)
            @test loaded.FEI(r) ≈ fI(r) rtol=1e-10
            @test loaded.FEH(r) ≈ fH(r) rtol=1e-10
        end
        @test flux_1SF(loaded, 15.0) ≈ fI(15.0) + fH(15.0) rtol=1e-10
    end
end
