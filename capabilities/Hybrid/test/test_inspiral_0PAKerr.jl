# Characterization of the 0PAKerr adiabatic (equatorial Kerr) forcing model (Hybrid.Models, §4.9).
#
# GROUNDING (no WaSABI live dep — ANALYTIC closed forms + synthetic 2D-grid fixtures + the light
# live Spacetime.kerr_isco for the Kerr abscissa remap):
#   * forcing_Fω == the literal num/den closed form of 0PAKerr.m, evaluated with the (already
#     grounded) Interp2Hermite fluxes at (Log[1-χ1], √(rISCO[χ1]/r)), r=(1/ω-χ1)^{2/3}.
#   * omega_crit == min(1/((1.05 rISCO)^{3/2}+χ1), 1/(6.26^{3/2}+χ1)) vs Spacetime.kerr_isco.
#   * kerr_inspiral_rhs! == (Fω/M, ω/M).
#   * forcing_0PAKerr loader: a synthetic .data file with constant-c2 columns (r = rISCO/c2²)
#     forms an exact tensor product; the remapped node values r^5·inf / r^{15/2}·hor reproduce a
#     separable linear surface, which the loaded Interp2Hermite returns exactly.

using Test
using Hybrid.Models
using Spacetime: kerr_isco

@testset "Models.inspiral_0PAKerr" begin
    # --- synthetic separable 2D Hermite fluxes (struct-level forcing_Fω algebra) ----------
    Xn = log.(1 .- [0.1, 0.3, 0.5, 0.7]); Yn = [0.2, 0.4, 0.6, 0.8, 1.0]
    PI(X, Y) = 1.0 + 0.5X + 0.2Y
    PH(X, Y) = 0.7 - 0.3X + 0.4Y
    gI = [((X, Y), PI(X, Y)) for X in Xn for Y in Yn]
    gH = [((X, Y), PH(X, Y)) for X in Xn for Y in Yn]
    model = Forcing0PAKerr(interp2_hermite(gI), interp2_hermite(gH))

    @testset "forcing_Fω literal num/den" begin
        ν = 1.0e-3
        for (ω, χ1) in ((0.04, 0.1), (0.05, 0.3), (0.03, 0.5))
            s = χ1 * ω; u = 1 - s
            rr = (1 / ω - χ1)^(2 / 3)
            c1 = log(1 - χ1); c2 = sqrt(kerr_isco(χ1) / rr)
            fi = model.FI(c1, c2); fh = model.FH(c1, c2)
            num = 3 * ν * ((-1 + s) * (-1 - s + 3 * ω^(2 / 3) * u^(1 / 3)))^(3 / 2) *
                  (ω^(5 / 3) * u^(1 / 3) * fh + (-1 + s)^2 * fi)
            den = (-χ1 + 1 / ω)^(11 / 3) * (-1 + s) *
                  (-1 + 6 * ω^(2 / 3) * u^(1 / 3) - 6 * χ1 * (ω + ω^(5 / 3) * u^(1 / 3)) +
                   χ1^2 * (7 * ω^2 + 3 * ω^(4 / 3) * u^(2 / 3)))
            @test forcing_Fω(model, ω, ν, χ1) ≈ num / den rtol = 1e-12
        end
    end

    @testset "omega_crit vs kerr_isco" begin
        for χ1 in (0.0, 0.3, 0.7, 0.9)
            ri = kerr_isco(χ1)
            @test omega_crit(model, χ1) ≈
                  min(1 / ((1.05 * ri)^(3 / 2) + χ1), 1 / (6.26^(3 / 2) + χ1)) rtol = 1e-12
        end
    end

    @testset "kerr_inspiral_rhs!" begin
        ν, M, χ1 = 1.0e-3, 1.0, 0.3; ω, φ = 0.05, 0.7
        du = zeros(2)
        kerr_inspiral_rhs!(du, [ω, φ], (model, ν, M, χ1), 0.0)
        @test du[1] == forcing_Fω(model, ω, ν, χ1) / M
        @test du[2] == ω / M
    end

    @testset "forcing_0PAKerr loader (synthetic tensor-product .data)" begin
        as = [0.1, 0.3, 0.5]; ts = [0.3, 0.6, 0.9]
        recs = String[]
        for a in as, t in ts
            r = kerr_isco(a) / t^2
            X = log(1 - a)
            inf = PI(X, t) / r^5
            hor = PH(X, t) / r^(15 / 2)
            push!(recs, "<|\"a\" -> $(repr(a)), \"r\" -> $(repr(r)), " *
                        "\"inf\" -> $(repr(inf)), \"hor\" -> $(repr(hor))|>")
        end
        dir = mkpath(joinpath(mktempdir(), "1SF_Flux", "Kerr_Circ"))
        write(joinpath(dir, "SMRfluxdata2025_36x36.data"), "{" * join(recs, ", ") * "}")
        loaded = forcing_0PAKerr(dirname(dirname(dir)))
        for a in as, t in ts
            X = log(1 - a)
            @test loaded.FI(X, t) ≈ PI(X, t) rtol = 1e-9
            @test loaded.FH(X, t) ≈ PH(X, t) rtol = 1e-9
        end
        # linear surface ⇒ exact also off-node
        @test loaded.FI(log(1 - 0.2), 0.5) ≈ PI(log(1 - 0.2), 0.5) rtol = 1e-9
    end
end
