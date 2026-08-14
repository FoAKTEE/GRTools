# Characterization of the 1PAT1 first-post-adiabatic forcing model (Hybrid.Models, §4.9).
#
# GROUNDING (no WaSABI live dep — ANALYTIC closed forms + synthetic WL-table fixtures, with the
# already-grounded Hybrid.BindingEnergy.E_FLx_deriv reused for the EFLx' factor):
#   * flux_1SF / flux_1SF_deriv: additive identity FEI(±)+FEH(±) on degree-≤8 polynomial tables
#     (order-8 Interp1 reproduces the polynomials and their derivatives exactly).
#   * flux_2SF == Re[Exp[F2(Log r0)]]: with the 2SF table set to {Log r0, log g(r0)} for a
#     log-linear g(r0)=e^a r0^b, F2 is linear and flux_2SF == g(r0) in closed form.
#   * forcing_F0 == pref·F1 (the literal 0PA prefactor, shared with 0PASchwarz).
#   * forcing_F1 == pref·(F2I + 2(1-3/r0)^{3/2}/(1-6/r0)·F1·E_FLx_deriv(z,1/r0)) — the literal
#     1PA closed form, with E_FLx_deriv taken from the already-grounded BindingEnergy leaf.
#   * omega_dot == (νF0+ν²F1)/M and stop_r0 == max(6.25,6.05+4.5ν-1.5ν²) literals.

using Test
using Hybrid.Models
using Hybrid.BindingEnergy: E_FLx_deriv

_write_wl(path, xs, ys) =
    (write(path, "{" * join(("{$(x), $(y)}" for (x, y) in zip(xs, ys)), ", ") * "}"); path)

@testset "Models.inspiral_1PAT1" begin
    # Synthetic degree-≤8 polynomial 1SF flux tables over r0.
    fI(r) = 0.02 + 1.0e-3 * r - 5.0e-5 * r^2 + 1.0e-6 * r^3
    fH(r) = 0.005 - 2.0e-4 * r + 3.0e-6 * r^2
    dfI(r) = 1.0e-3 - 1.0e-4 * r + 3.0e-6 * r^2
    dfH(r) = -2.0e-4 + 6.0e-6 * r
    rs = collect(range(6.2, 32.0; length = 14))
    FEH = interp1_order8(rs, fH.(rs))
    FEI = interp1_order8(rs, fI.(rs))

    # 2SF complex table {Log r0, log g(r0)} with g(r0) = e^a r0^b (log-linear → exact).
    a2, b2 = -3.0, -3.0
    g2(r) = exp(a2) * r^b2
    F2 = interp1_order8(log.(rs), [complex(a2 + b2 * log(r), 0.0) for r in rs])

    # Redshift interpolant z over x = 1/r0 (degree-≤8 poly; E_FLx_deriv evaluates z(x)).
    zfun(x) = 0.5 + 0.3 * x - 0.1 * x^2 + 0.02 * x^3
    xs = collect(range(0.02, 0.18; length = 12))
    z = interp1_order8(xs, zfun.(xs))

    model = Forcing1PAT1(FEH, FEI, F2, z)

    @testset "flux identities" begin
        for r0 in (8.0, 12.5, 20.0, 30.0)
            @test flux_1SF(model, r0) ≈ fI(r0) + fH(r0) rtol = 1e-10
            @test flux_1SF_deriv(model, r0) ≈ dfI(r0) + dfH(r0) rtol = 1e-8
            @test flux_2SF(model, r0) ≈ g2(r0) rtol = 1e-10
        end
    end

    @testset "forcing F0 / F1 closed forms" begin
        for r0 in (8.0, 12.5, 20.0)
            F1 = fI(r0) + fH(r0)
            pref = 3 * (1 - 3 / r0)^(3 / 2) * sqrt(1 / r0) / (1 - 6 / r0)
            @test forcing_F0(model, r0) ≈ pref * F1 rtol = 1e-10
            coupling = (2 * (1 - 3 / r0)^(3 / 2)) / (1 - 6 / r0)
            expF1 = pref * (g2(r0) + coupling * F1 * E_FLx_deriv(z, 1 / r0))
            @test forcing_F1(model, r0) ≈ expF1 rtol = 1e-10
            ν, M = 1.0e-3, 1.0
            @test omega_dot(model, r0, ν, M) ≈
                  (ν * forcing_F0(model, r0) + ν^2 * forcing_F1(model, r0)) / M rtol = 0
        end
    end

    @testset "stop_r0 ν-dependent" begin
        for ν in (1.0e-4, 1.0e-2, 0.1, 0.2)
            @test stop_r0(model, ν) == max(6.25, 6.05 + 4.5 * ν - 1.5 * ν^2)
        end
    end

    @testset "Data-load round-trip" begin
        dir = mktempdir()
        d1 = mkpath(joinpath(dir, "1SF_Flux", "Schwarz_Circ"))
        d2 = mkpath(joinpath(dir, "2SF_Flux", "Schwarz_Circ"))
        dz = mkpath(joinpath(dir, "1SF_Local_Invariants", "Schwarz_Circ"))
        _write_wl(joinpath(d1, "1SFCircShwarzDotEInf.m"), rs, fI.(rs))
        _write_wl(joinpath(d1, "1SFCircShwarzDotEHorizon.m"), rs, fH.(rs))
        # 2SF complex table: {Log r0, re + im*I}
        open(joinpath(d2, "2SFCircShwarzDotEInf.m"), "w") do io
            print(io, "{" * join(("{$(log(r)), $(a2 + b2 * log(r)) + 0.0*I}" for r in rs), ", ") * "}")
        end
        _write_wl(joinpath(dz, "1SFCircShwarzRedshfit.m"), xs, zfun.(xs))
        loaded = forcing_1PAT1(dir)
        @test flux_1SF(loaded, 15.0) ≈ fI(15.0) + fH(15.0) rtol = 1e-10
        @test flux_2SF(loaded, 15.0) ≈ g2(15.0) rtol = 1e-10
        @test forcing_F1(loaded, 15.0) ≈ forcing_F1(model, 15.0) rtol = 1e-9
    end
end
