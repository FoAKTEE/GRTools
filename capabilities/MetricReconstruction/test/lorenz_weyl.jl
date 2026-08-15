# Grounding for Lorenz.LorenzWeyl — the gauge-invariant ψ4 of the a≠0 Lorenz reconstruction. ANALYTIC /
# characterization gates (no heavy oracle: SWSH S from Harmonics, Teukolsky R from Perturbation — both built):
#
# (a) θ-FACTORIZATION identity: weyl_psi4_lorenz(rec,r,θ;l) == ρ⁴·weyl_psi4_lorenz_amp(rec,r)·{}_{-2}S_{lm}(θ)
#     (assembly identity — the θ-dependence factors through ρ⁴ and the spheroidal harmonic).
# (b) BRANCH selection + module-1 tie: weyl_psi4_lorenz_amp(rec,r;bc=:Up)==valK(rec.Pm2,r) for r>r0 and
#     (bc=:In)==valK(rec.Pm2,r) for r<r0 (the HomFieldK amplitude×radial).
# (c) ψ4 RECONSTRUCTION identity: the amp radial part R_{-2}(r) == Aup·Teukolsky s=−2 R^up(r) — i.e. P_{-2}
#     IS the Teukolsky Weyl mode (reldiff vs GrTools Perturbation TeukolskyRadial < 1e-9, the App-H ψ4 check).

using Test
import MetricReconstruction.Lorenz as L
import Perturbation: TeukolskyRadial
import Harmonics: SpinWeightedSpheroidalHarmonicS

@testset "Lorenz.LorenzWeyl (gauge-invariant ψ4 cross-check)" begin
    a = 0.5; r0 = 6.0; l = 2; m = 2
    sol, aux = L.closure_kerr(a, r0, m; lmax=max(2l, l+8), return_aux=true)
    rec = L.reconstruct_kerr(a, r0, l; m=m, aux=aux, sol=sol)

    @testset "(a) θ-factorization ψ4 == ρ⁴·amp·S(θ)" begin
        for r in (r0+2.0, r0-1.5), θ in (0.4, 1.1, 2.3)
            bc = r > r0 ? :Up : :In
            amp = L.weyl_psi4_lorenz_amp(rec, r; bc=bc)
            Sm2 = real(SpinWeightedSpheroidalHarmonicS(-2, l, m, Float64(a*rec.jd.ω))(θ, 0.0))
            ρ = -1/(r - im*a*cos(θ))
            lhs = L.weyl_psi4_lorenz(rec, r, θ; bc=bc, l=l)
            @test abs(lhs - ρ^4*amp*Sm2) / (abs(ρ^4*amp*Sm2) + 1e-300) < 1e-12
        end
    end

    @testset "(b) branch selection == HomFieldK valK(rec.Pm2)" begin
        for r in (r0+2.0, r0+4.0)
            @test abs(L.weyl_psi4_lorenz_amp(rec, r; bc=:Up) - L.valK(rec.Pm2, r)) /
                  (abs(L.valK(rec.Pm2, r)) + 1e-300) < 1e-12
        end
        for r in (r0-1.0, r0-2.0)
            @test abs(L.weyl_psi4_lorenz_amp(rec, r; bc=:In) - L.valK(rec.Pm2, r)) /
                  (abs(L.valK(rec.Pm2, r)) + 1e-300) < 1e-12
        end
    end

    @testset "(c) ψ4 reconstruction: amp/Aup == Teukolsky s=−2 R^up (P_{-2}≡R_{-2})" begin
        Rteuk = TeukolskyRadial(-2, l, m, a, rec.jd.ω; lambda=rec.jd.Λ)
        for r in (r0+1.0, r0+3.0)
            ampR = L.weyl_psi4_lorenz_amp(rec, r; bc=:Up) / rec.Pm2.Aup
            @test abs(ampR - Rteuk.Up(r)) / (abs(Rteuk.Up(r)) + 1e-300) < 1e-9
        end
    end
end
