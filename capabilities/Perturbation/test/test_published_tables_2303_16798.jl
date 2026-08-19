# PUBLISHED-TABLE regression: arXiv:2303.16798 (Skoupý, Lukes-Gerakopoulos, Drummond, Hughes),
# "Asymptotic gravitational-wave fluxes from a spinning test body on generic orbits around a Kerr
# black hole" — Tables I, II and III of the paper, computed here by Perturbation.SpinFluxes.
#
# These three tables are the headline numeric reproduction of the paper; until now they were locked
# ONLY in the legacy TeukolskySpinFluxes test suite. The goldens below are the PUBLISHED numbers
# (paper.tex tab:error_generic L613-619, tab:spherical L850-865, tab:generic L878-903), not a
# self-captured fixture.
#
# Conventions carried over from the legacy lock:
#   * Table III sign: C_Julia = −C_paper (the flux |C|² is convention-independent).
#   * Table II I=150°, m=2 has a PAPER TYPO: F^Jz prints 1.658e-5 but is 1.658e-4 — that row's
#     F^Jz is asserted against the corrected value and the typo is recorded here, not silently fixed.
#   * The tolerances are the reproduction levels established in the legacy lock (Table III 0.28 %
#     worst, Table II ~2.5 % worst on the retrograde rows, Table I ~0.3 % for the truncated sum).
#
# Cost note: these are inhomogeneous-Teukolsky mode sums; the file is the slowest in the suite
# (a few minutes at -t4). The subsets locked here are the legacy ones; the FULL 24-row Table III
# and the FULL 12-row Table II run as reproduction drivers under reproduce/arXiv-2303.16798/.
using Test
using Perturbation
using Perturbation: KerrSpinOrbitCorrection, KerrSpinOrbitCorrectionSpherical,
                    _orbit_derivatives_spherical

@testset "published tables — arXiv:2303.16798" begin

    @testset "Table III (tab:generic) — linear-in-spin amplitudes C^±_{S,lmnk}" begin
        # a=0.9, p=12, e=0.2, I=30° (x=cos30°). C^+ = AmplitudesCorrection.I, C^- = .H, sign-flipped.
        oc = KerrSpinOrbitCorrection(0.9, 12.0, 0.2, cos(pi / 6), 8, 8)
        rows = [(2, 1, 0, 1,  4.8962e-6 - 1.6020e-6im, -5.2716e-6 - 2.7823e-7im),
                (2, 1, 1, 1,  9.9514e-6 - 2.7846e-6im,  3.8592e-6 + 8.7391e-7im),
                (2, 2, 0, 0, -1.9890e-5 + 5.8986e-6im, -7.5636e-6 - 3.6977e-6im),
                (2, 2, 1, 0, -3.6535e-5 + 1.0473e-5im, -2.8727e-5 - 9.5114e-6im),
                (2, 2, 2, 0, -2.8239e-5 + 8.6430e-6im, -2.1354e-5 - 7.5839e-6im),
                (3, 2, 1, 1,  3.8727e-6 + 7.7467e-6im, -1.2893e-7 - 2.4943e-7im),
                (3, 3, 1, 0, -8.9746e-6 - 1.7781e-5im,  5.0351e-7 + 3.8648e-6im),
                (4, 3, 2, 1, -5.6428e-6 + 4.2331e-6im, -8.4756e-8 + 9.4155e-8im)]
        for (l, m, n, k, Cp, Cm) in rows
            md = TeukolskySpinModeCorrection(l, m, n, k, oc)
            @test isapprox(-md.AmplitudesCorrection.I, Cp; rtol = 2e-3)
            @test isapprox(-md.AmplitudesCorrection.H, Cm; rtol = 2e-3)
        end
    end

    @testset "Table II (tab:spherical) — spherical-orbit fluxes F^E_{S,m}, F^{Jz}_{S,m}" begin
        # a=0.9, p=10, e=0, I=30°; the m=2 row, summed over l=2..5 and k=-4..4.
        oc  = KerrSpinOrbitCorrectionSpherical(0.9, 10.0, 0.0, cos(pi / 6), 8)
        odv = _orbit_derivatives_spherical(0.9, 10.0, 0.0, cos(pi / 6), 8; h = 1e-5)
        FE = 0.0 + 0im; FJ = 0.0 + 0im
        for l in 2:5, k in -4:4
            md = TeukolskySpinModeSphericalCorrection(l, 2, k, oc, odv)
            FE += md.FluxesCorrection.Energy.I + md.FluxesCorrection.Energy.H
            FJ += md.FluxesCorrection.AngularMomentum.I + md.FluxesCorrection.AngularMomentum.H
        end
        @test isapprox(real(FE), -2.702e-6; rtol = 1e-3)     # published
        @test isapprox(real(FJ), -6.431e-5; rtol = 1e-3)     # published
        @test real(FE) < 0 && real(FJ) < 0                    # prograde signs
    end

    @testset "Table I (tab:error_generic) — generic-orbit summed flux F^E_{S,m}, row 1" begin
        # a=0.9, p=10, e=0.1, I=15°, m=2 → published −2.8259e-6. The n-sum converges by |n|≤3 and
        # the l-sum by l≤3; the truncated l=2..3, |n|,|k|≤2 partial sits ~0.3 % from the published
        # value (the fully converged Julia value is −2.819e-6, ~0.24 %, the MST-analytic level).
        oc = KerrSpinOrbitCorrection(0.9, 10.0, 0.1, cos(15 * pi / 180), 8, 8)
        FE = 0.0 + 0im
        for l in 2:3, n in -2:2, k in -2:2
            FE += let md = TeukolskySpinModeCorrection(l, 2, n, k, oc)
                md.FluxesCorrection.Energy.I + md.FluxesCorrection.Energy.H
            end
        end
        @test isapprox(real(FE), -2.8259e-6; rtol = 5e-3)
    end
end
