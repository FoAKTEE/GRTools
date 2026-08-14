# Grounding for Lorenz.KerrSMatrices — the a≠0 closure S-matrices (smatrices_kerr) + the radial-operator chains
# and the lplp/mpmp projected tetrad components.
#
# (1) smatrices_kerr CHARACTERIZATION vs the legacy KerrLorenzCirc.smatrices_kerr oracle (the port reuses the
#     SAME App-H banded formulas already inside closure_kerr — validated App-H a=0.99 lower 6.85e-8). The expected
#     entries are a ONE-SHOT capture (corpus-free):
#       cd packages/KerrLorenzCirc && julia --project=. -e '
#         using KerrLorenzCirc; const K=KerrLorenzCirc; S=K.smatrices_kerr(0.5,6.0,2;lmax=12)
#         for (l,lp) in ((2,2),(2,4),(4,2),(3,3)); println(S.s2lplp(l,lp), …); end'
# (2) a=0 REDUCTION (analytic): at a=0 the b-matrices = I and the aω-mixing vanishes ⇒ S.s2mpmp/s0lplp are the
#     identity, S.s2lplp(l,l)=Λ̂(l), S.s1mpmp(l,l)=−λ̂₂(l), S.s0mpmp_h≡0.
# (3) tetrad_lplp_kerr/tetrad_mpmp_kerr CHARACTERIZATION vs legacy (same capture script).
# (4) a→0 REGRESSION (G-K11): tetrad_lplp/mpmp_kerr(0) collapse to the a=0 SchwarzschildComponents.tetrad_modes
#     lplp/mpmp (the validated a=0 diagonal).
# (5) sm1_third / skappa_second analytic vs central finite difference of the module-1 accessors.

using Test
import MetricReconstruction.Lorenz as L

@testset "Lorenz.KerrSMatrices (a≠0 closure S-matrices + lplp/mpmp components)" begin

    @testset "(1) smatrices_kerr(0.5,6,2) vs legacy KLC oracle" begin
        S = L.smatrices_kerr(0.5, 6.0, 2; lmax=max(2*2, 2+8))
        ref = Dict(
            (:s2lplp,2,2)=>9.375882640855467,    (:s2lplp,2,4)=>0.0012396923462869188,
            (:s2lplp,4,2)=>-0.0018852439786235476,
            (:s2mpmp,2,2)=>0.9999233789476186,   (:s2mpmp,2,4)=>0.00012059721080922477,
            (:s2mpmp,3,3)=>0.9998453823315746,
            (:s1lplp,2,2)=>-4.791828620503162,   (:s1lplp,2,4)=>-0.00030044619504083353,
            (:s1mpmp,2,2)=>-1.955977398665349,   (:s1mpmp,3,3)=>-3.1270051190554855,
            (:s0lplp,2,2)=>0.9999999986983268,   (:s0lplp,2,4)=>5.102299788760325e-5,
            (:s0mpmp_κ,2,2)=>4.686088931622465,  (:s0mpmp_κ,3,3)=>10.714173572672388,
            (:s0mpmp_h,2,2)=>-1.5720931279593808,(:s0mpmp_h,2,4)=>-0.08500805226710709,
        )
        fns = (s2lplp=S.s2lplp, s2mpmp=S.s2mpmp, s1lplp=S.s1lplp, s1mpmp=S.s1mpmp,
               s0lplp=S.s0lplp, s0mpmp_κ=S.s0mpmp_κ, s0mpmp_h=S.s0mpmp_h)
        for ((sym,l,lp), v) in ref
            got = getfield(fns, sym)(l, lp)
            @test abs(got - v) / (abs(v) + 1e-12) < 1e-9
        end
    end

    @testset "(2) a=0 reduction (b=I, analytic)" begin
        S0 = L.smatrices_kerr(0.0, 6.0, 2; lmax=12)
        for l in (2,3,4)
            σ = (-1)^(l+2)   # σof(l) at m=2
            @test abs(S0.s2mpmp(l,l) - 1.0) < 1e-12
            @test abs(S0.s0lplp(l,l) - 1.0) < 1e-12
            # s2lplp b=I collapse: Λ̂(l) from bget(2,l) + σ·Λ̂(l) from bget(−2,l) ⇒ (1+σ)Λ̂(l)
            # (even-l → 2Λ̂, odd-l → 0 by the σ=−1 cancellation).
            @test abs(S0.s2lplp(l,l) - (1+σ)*L.Λ̂(l)) < 1e-10
            @test abs(S0.s1mpmp(l,l) - (-L.λ̂2(l))) < 1e-10
            @test abs(S0.s0mpmp_h(l,l)) < 1e-12
            for lp in (l-2, l+2)
                lp < 2 && continue
                @test abs(S0.s2mpmp(l,lp)) < 1e-12      # off-diagonal vanishes at a=0
                @test abs(S0.s0lplp(l,lp)) < 1e-12
            end
        end
    end

    @testset "(3) tetrad_lplp/mpmp_kerr(0.5,6,2) vs legacy KLC oracle" begin
        ref = Dict(
            (:lplp,7.0)=>0.49180152849016756 + 0.1994847469131672im,
            (:mpmp,7.0)=>-5.171223616613378 - 9.213644641927916im,
            (:lplp,9.0)=>0.17295754166041094 + 0.1289531925444104im,
            (:mpmp,9.0)=>-2.9657468223897503 - 14.014205614882798im,
        )
        sol, aux = L.closure_kerr(0.5, 6.0, 2; lmax=12, return_aux=true)
        recs = Dict{Int,Any}()
        for l in aux.ls
            recs[l] = L._attach_recmeta(L.reconstruct_kerr(0.5, 6.0, l; m=2, aux=aux, sol=sol), l, aux)
        end
        S = L.smatrices_kerr(0.5, 6.0, 2; lmax=12)
        for r in (7.0, 9.0)
            lk = L.tetrad_lplp_kerr(0.5, 6.0, 2, 2, r; recs=recs, S=S, ls=aux.ls)
            mk = L.tetrad_mpmp_kerr(0.5, 6.0, 2, 2, r; recs=recs, S=S, ls=aux.ls)
            @test abs(lk - ref[(:lplp,r)]) / abs(ref[(:lplp,r)]) < 1e-8
            @test abs(mk - ref[(:mpmp,r)]) / abs(ref[(:mpmp,r)]) < 1e-8
        end
    end

    @testset "(4) a→0 regression: tetrad_*_kerr(0) == a=0 tetrad_modes (G-K11)" begin
        R3 = L.build_recon(6.0, 2; m=2)
        for r in (8.0, 9.0)
            M3 = L.tetrad_modes(R3, r)
            lk = L.tetrad_lplp_kerr(0.0, 6.0, 2, 2, r)
            mk = L.tetrad_mpmp_kerr(0.0, 6.0, 2, 2, r)
            @test abs(lk - M3.lplp) / abs(M3.lplp) < 3.2e-3
            @test abs(mk - M3.mpmp) / abs(M3.mpmp) < 3.2e-3
        end
    end

    @testset "(5) sm1_third / skappa_second analytic vs FD" begin
        sol, aux = L.closure_kerr(0.5, 6.0, 2; lmax=12, return_aux=true)
        rec = L.reconstruct_kerr(0.5, 6.0, 2; m=2, aux=aux, sol=sol)
        r = 8.5; h = 1e-5
        P = L.valK(rec.Pm1, r); Pp = L.derK(rec.Pm1, r); Ppp = L.der2K(rec.Pm1, r)
        sm1 = L.sm1_third(rec, r, P, Pp, Ppp)
        fd_sm1 = (L.der2K(rec.Pm1, r+h) - L.der2K(rec.Pm1, r-h)) / (2h)
        @test abs(sm1 - fd_sm1) / abs(fd_sm1) < 1e-4
        κ = L.κvalK(rec.κ, r); κp = L.κderK(rec.κ, r)
        sk = L.skappa_second(rec, r, κ, κp)
        fd_sk = (L.κderK(rec.κ, r+h) - L.κderK(rec.κ, r-h)) / (2h)
        @test abs(sk - fd_sk) / abs(fd_sk) < 1e-4
    end
end
