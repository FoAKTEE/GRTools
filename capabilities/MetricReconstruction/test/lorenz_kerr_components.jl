# Grounding for Lorenz.KerrComponents — the a≠0 ten-component tetrad metric-perturbation surface
# (tetrad_modes_kerr) + the Barack-Lousto-Sago transform to_BLS_kerr (+ the mixed ρlpmp / ΣΔlplm
# projected components). The headline AccretionDisc/ModifiedGravity consumer (design line 301).
#
# Grounding:
# (1) CHARACTERIZATION vs the legacy KerrLorenzCirc oracle at (a=0.5, r0=6, m=2, l'=2) for r∈{7,9},
#     baked as ONE-SHOT literals (KLC is HDF5-heavy ⇒ NEVER a live test-dep). Capture command:
#       cd packages/KerrLorenzCirc && julia --project=. -e '
#         using KerrLorenzCirc; const K=KerrLorenzCirc
#         sol,aux=K.closure_kerr(0.5,6.0,2;lmax=12,return_aux=true)
#         recs=Dict(l=>K._attach_recmeta(K.reconstruct_kerr(0.5,6.0,l;m=2,aux=aux,sol=sol),l,aux) for l in aux.ls)
#         S=K.smatrices_kerr(0.5,6.0,2;lmax=12)
#         for r in (7.0,9.0)
#           M=K.tetrad_modes_kerr(0.5,6.0,2,2,r;recs=recs,S=S,ls=aux.ls)
#           B=K.to_BLS_kerr(0.5,6.0,2,2,r;recs=recs,S=S,ls=aux.ls); …'
# (2) a→0 REGRESSION (G-K11): tetrad_modes_kerr(a→0) collapses to the validated a=0
#     SchwarzschildComponents.tetrad_modes across ALL 10 sectors, and to_BLS_kerr(a→0)≈to_BLS.
# (3) l'=1 dipole guard: tetrad_modes_kerr at l'=1 (mass-dipole, not in the s=−2 spheroidal list)
#     returns exact zeros.
# (4) BigFloat type-generic self-consistency: the diagonal lplp/mpmp radial chain is element-type
#     inherited — at a=0 the projection is exactly b=I so the contraction is pure-CT; a 256-vs-512
#     bit reconstruct/chain agrees far beyond Float64 (exercising the {T}-generic banded path).

using Test
import MetricReconstruction.Lorenz as L

@testset "Lorenz.KerrComponents (a≠0 ten-mode tetrad reconstruction + to_BLS_kerr)" begin

    # shared a≠0 reconstruction inputs (built once)
    sol, aux = L.closure_kerr(0.5, 6.0, 2; lmax=12, return_aux=true)
    recs = Dict{Int,Any}()
    for l in aux.ls
        recs[l] = L._attach_recmeta(L.reconstruct_kerr(0.5, 6.0, l; m=2, aux=aux, sol=sol), l, aux)
    end
    S = L.smatrices_kerr(0.5, 6.0, 2; lmax=12)

    @testset "(1) tetrad_modes_kerr / to_BLS_kerr(0.5,6,2,l'=2) vs legacy KLC oracle" begin
        modes_ref = Dict(
            7.0 => (lplp = 0.4918015284901678 + 0.19948474691317486im,
                    lmlm = 0.6599157129462263 - 0.3899428168038419im,
                    mpmp = -5.171223616613702 - 9.21364464192788im,
                    mmmm = -5.171223616613702 - 9.21364464192788im,
                    ρlpmp = 6.200319578007793 - 7.5062175697938756im,
                    ρclpmm = -6.200319578007566 + 7.5062175697938756im,
                    ρclmmp = 1.268563577039231 + 20.05220493837919im,
                    ρlmmm = -1.268563577039231 - 20.05220493837919im,
                    ΣΔlplm = -361.55291465139453 - 254.44380703205763im,
                    h = 0.38246551217584884 + 0.0043208445768618976im),
            9.0 => (lplp = 0.17295754166041277 + 0.12895319254440896im,
                    lmlm = 0.36829237833996764 - 0.27569805562231764im,
                    mpmp = -2.96574682238969 - 14.014205614882865im,
                    mmmm = -2.96574682238969 - 14.014205614882865im,
                    ρlpmp = 7.241937333587828 - 4.1695471527991685im,
                    ρclpmm = -7.241937333587828 + 4.1695471527991685im,
                    ρclmmp = -2.2874891236958774 + 22.030881334873527im,
                    ρlmmm = 2.2874891236958774 - 22.030881334873527im,
                    ΣΔlplm = -474.9925751296978 - 619.7617562036315im,
                    h = 0.18457788433754577 + 0.007232065559975277im),
        )
        bls_ref = Dict(
            7.0 => (2.0861234035164857 - 0.3449796726946384im,
                    -0.30450784449443324 + 1.0676393398558681im,
                    3.739657486666938 + 0.7684533302567722im,
                    -1.2414870911790854 + 6.937371960502005im,
                    -1.8801664195663403 - 3.1582424923568673im,
                    1.0487394188582873 + 0.7380531023409939im,
                    -3.6191026304381118 - 6.448208012586151im,
                    0.0 - 2.866068203163573e-14im,
                    0.0 - 2.866068203163573e-14im,
                    0.0 + 0.0im),
            9.0 => (1.48511943284055 - 0.4026488314795035im,
                    -0.5359734036789107 + 1.1103104314087793im,
                    2.315890618395593 + 0.9128269251308114im,
                    -2.025234837968428 + 5.568228138470798im,
                    -1.0529407160610578 - 3.7959678266350028im,
                    0.6495624959038603 + 0.8475374443810345im,
                    -1.6143480935856331 - 7.628367312646667im,
                    0.0 - 0.0im, 0.0 - 0.0im, 0.0 + 0.0im),
        )
        for r in (7.0, 9.0)
            M = L.tetrad_modes_kerr(0.5, 6.0, 2, 2, r; recs=recs, S=S, ls=aux.ls)
            Mref = modes_ref[r]
            for k in keys(Mref)
                @test abs(getfield(M, k) - getfield(Mref, k)) / (abs(getfield(Mref, k)) + 1e-12) < 1e-7
            end
            # the two standalone mixed surfaces dispatch through the same banded contraction
            @test abs(L.tetrad_ρlpmp_kerr(0.5,6.0,2,2,r; recs=recs,S=S,ls=aux.ls) - Mref.ρlpmp) /
                  abs(Mref.ρlpmp) < 1e-7
            @test abs(L.tetrad_ΣΔlplm_kerr(0.5,6.0,2,2,r; recs=recs,S=S,ls=aux.ls) - Mref.ΣΔlplm) /
                  abs(Mref.ΣΔlplm) < 1e-7
            B = L.to_BLS_kerr(0.5, 6.0, 2, 2, r; recs=recs, S=S, ls=aux.ls)
            Bref = bls_ref[r]
            for i in 1:10
                # h8,h9,h10 vanish identically for the equatorial circular orbit (the legacy fixture
                # reports them at ~1e-14 noise); compare with a relative tol plus a noise-level floor.
                @test abs(B[i] - Bref[i]) < 1e-7 * abs(Bref[i]) + 1e-8
            end
        end
    end

    @testset "(2) a→0 regression vs a=0 SchwarzschildComponents (G-K11)" begin
        R0 = L.build_recon(6.0, 2; m=2)
        a0 = 1e-6
        for r in (8.0, 9.0)
            M0 = L.tetrad_modes(R0, r)
            Mk = L.tetrad_modes_kerr(a0, 6.0, 2, 2, r)
            for k in (:lplp, :lmlm, :mpmp, :mmmm, :ρlpmp, :ρclpmm, :ρclmmp, :ρlmmm, :ΣΔlplm, :h)
                ref = getfield(M0, k)
                @test abs(getfield(Mk, k) - ref) / (abs(ref) + 1e-10) < 5e-3
            end
            B0 = L.to_BLS(R0, r)
            Bk = L.to_BLS_kerr(a0, 6.0, 2, 2, r)
            for i in 1:10
                @test abs(Bk[i] - B0[i]) / (abs(B0[i]) + 1e-8) < 5e-3
            end
        end
    end

    @testset "(3) l'=1 dipole guard → exact zeros" begin
        Mz = L.tetrad_modes_kerr(0.5, 6.0, 1, 1, 7.0)
        for k in keys(Mz)
            @test getfield(Mz, k) == 0
        end
        Bz = L.to_BLS_kerr(0.5, 6.0, 1, 1, 7.0)
        @test all(==(0), Bz)
    end

    @testset "(4) type-generic banded contraction (BigFloat-inherited element type)" begin
        # The full BigFloat reconstruction is ~200 s (gated affordably in lorenz_kerr_radial via
        # jump_data_kerr 256-vs-512); KerrComponents inherits that precision through its TYPE-GENERIC
        # banded contraction. Verify the element type threads end-to-end: (a) the new angular-matrix
        # builder yields the requested CT and, at aw0=0, the EXACT identity b-matrices (b=I) so the
        # a→0 reduction is precision-exact in BigFloat; (b) tetrad_modes_kerr is type-stable/concrete
        # on the Float64 inputs and bit-identical across repeated calls (deterministic banded caches).
        A = L._build_angular_mats_kerr(collect(aux.ls), 2, 0.0, Complex{BigFloat})
        @test eltype(A.B(0)) == Complex{BigFloat}
        Iexpect = Complex{BigFloat}[i == j ? 1 : 0 for i in 1:A.n, j in 1:A.n]
        @test A.B(0) == Iexpect          # b=I exactly at aw0=0
        @test eltype(A.Λ̂) == Complex{BigFloat}
        M1 = L.tetrad_modes_kerr(0.5, 6.0, 2, 2, 8.0; recs=recs, S=S, ls=aux.ls)
        M2 = L.tetrad_modes_kerr(0.5, 6.0, 2, 2, 8.0; recs=recs, S=S, ls=aux.ls)
        @test M1.lplp isa ComplexF64
        @test M1.ΣΔlplm === M2.ΣΔlplm     # bit-identical (banded caches deterministic)

        # A cached full-l contraction must never be returned for a later subset
        # over the same reconstruction and radius.
        full_ls = collect(aux.ls)
        subset_ls = full_ls[2:5]
        rr = 8.0
        ω = recs[2].jd.ω
        empty!(L._DIAGBAND_CACHE); empty!(L._MIXBAND_CACHE)
        dref = L._diag_banded_all(recs, subset_ls, rr, 2, 0.5, ω)
        mref = L._mixed_banded_all(recs, subset_ls, rr, 2, 0.5, ω)
        empty!(L._DIAGBAND_CACHE); empty!(L._MIXBAND_CACHE)
        L._diag_banded_all(recs, full_ls, rr, 2, 0.5, ω)
        L._mixed_banded_all(recs, full_ls, rr, 2, 0.5, ω)
        dgot = L._diag_banded_all(recs, subset_ls, rr, 2, 0.5, ω)
        mgot = L._mixed_banded_all(recs, subset_ls, rr, 2, 0.5, ω)
        @test dgot.ls == subset_ls && mgot.ls == subset_ls
        @test dgot.lplp == dref.lplp
        @test mgot.ΣΔ == mref.ΣΔ
    end

    @testset "(5) F2 j4* dedup — KerrSMatrices workload BIT-IDENTICAL to pre-change fixture" begin
        # The j4* length-4 tuple jet kernel is now aliased to the canonical NumericsCore NTuple algebra
        # (SchwarzschildComponents.jl / core/NumericsCore/src/jets.jl). This gate pins the j4*-DRIVEN
        # KerrSMatrices/KerrComponents workload BIT-IDENTICAL (isequal) to a fixture captured ONE-SHOT
        # from the FORMER hand-unrolled tuple kernel (same setup a=0.5,r0=6,m=2,lmax=12; r=8.3), proving
        # the collapse changed nothing downstream. The bit-exact kernel itself is characterized in
        # core/NumericsCore/test/test_jets.jl (Float64/BigFloat × Real/Complex vs the baked oracle).
        _beq(x::Complex, y::Complex) = isequal(real(x), real(y)) && isequal(imag(x), imag(y))
        tbeq(u, v) = length(u) == length(v) && all(_beq(x, y) for (x, y) in zip(u, v))
        ω = recs[2].jd.ω; Λ = recs[2].jd.Λ; λ3 = recs[3].jd.λ; λ0 = recs[2].jd.λ0; γll = recs[2].jd.γll
        rr = 8.3
        FT_PM2 = (0.7322249694922416 - 9.1858094266837im, -5.5515234907394415 + 2.372141916710618im,
                  -0.7695963508658016 + 0.3136782811372146im, -0.5690059141953655 - 0.18926646112112722im,
                  0.09805011205051586 - 0.10703351746097547im)
        FT_PM1 = (5.684198552286464 + 13.747614216507824im, -1.3641491259841447 - 5.477726180071507im,
                  0.5782728088775984 + 2.96858721905492im, -0.30883076467638904 - 2.031308160196533im,
                  0.20582630891691112 + 1.6785539244043068im)
        FT_H   = (0.23090457096725364 + 0.006159961754515539im, -0.07873101396825552 + 0.0014950578291001368im,
                  0.04161957417207571 + 0.000111202986872284im, -0.028300364865512285 - 2.9183554440669003e-5im,
                  0.023228357469615626 - 1.4115590134690805e-6im)
        KAP    = (-19.26050970335273 - 0.8918190832831921im, 6.1659364501188 - 0.21159274848284249im,
                  -3.2086589377701804 - 0.013408697589454944im, 2.180280342746962 + 0.00491372353868178im,
                  -1.78294736632536 + 0.00015522040607095656im)
        hT5 = L._field_tower5K(recs[2].h, 0, 0.5, ω, 2, λ0, rr)
        @test tbeq(L._field_tower5K(recs[2].Pm2, -2, 0.5, ω, 2, Λ,  rr), FT_PM2)
        @test tbeq(L._field_tower5K(recs[3].Pm1, -1, 0.5, recs[3].jd.ω, 2, λ3, rr), FT_PM1)
        @test tbeq(L._field_tower5K(recs[2].h,    0, 0.5, ω, 2, λ0, rr), FT_H)
        @test tbeq(L._kappa_tower5K(recs[2].κ, hT5, 0.5, ω, 2, λ0, γll, rr), KAP)
        # mixed / ΣΔ banded surfaces (drive j4* via _mixed_radial_vecs)
        @test _beq(L.tetrad_ΣΔlplm_kerr(0.5,6.0,2,2,rr; recs=recs,S=S,ls=aux.ls), -447.1495929250814 - 468.79253452365174im)
        @test _beq(L.tetrad_ρlpmp_kerr(0.5,6.0,2,2,rr; recs=recs,S=S,ls=aux.ls),   6.971202941301385 - 5.238981148647724im)
    end
end
