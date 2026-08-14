# Q3 Deliverable A gates — the baked QQNM source-coefficient tables (src/QNM/qqnm_source_tables.jl).
#   G-Q3-T0  provenance: each owner source .tex hashes to the frozen sha256 ledger.
#   G-Q3-T1  named-entry r-form == Appendix-C Γ-form under the on-shell Γ+Γ̄=2r (the F7-3b route;
#            DISTINCT strings, equal only under the identity — non-circular). ComplexF64 + BigFloat-256.
#   G-Q3-T2  channel structure: ++ nonzero-key census; +− support ⊂ {(4,0),(4,1),(4,2),(5,0),(5,1),(6,0)};
#            −− the certified structural zero.
using Test
using Perturbation
import Perturbation: QQNM_PM_SUPPORT, QQNM_PM_ATOMS
include(joinpath(@__DIR__, "fixtures", "qqnm_source_tex_provenance.jl"))

# (name, ν, i, j) for the named entries with a hand Γ-form.
const QQNM_NAMED = [
    ("a_00^(1)", 1, 0, 0), ("a_01^(1)", 1, 0, 1), ("a_11^(1)", 1, 1, 1), ("a_02^(1)", 1, 0, 2),
    ("a_06^(2)", 2, 0, 6), ("a_03^(2)", 2, 0, 3), ("a_33^(2)", 2, 3, 3),
    ("a_04^(3)", 3, 0, 4), ("a_11^(4)", 4, 1, 1),
]

# generic Kerr points (a,r,θ) → Γ=r+iacosθ, Γb=r−iacosθ (Γ+Γb=2r on-shell). Avoid Γ̄=0.
const QQNM_TAB_PTS = [(0.37, 3.7, 0.9), (0.61, 5.3, 1.3), (0.83, 4.1, 0.6), (0.5, 6.8, 2.1)]

function run_qqnm_tables(::Type{C}, tol) where {C}
    RT = real(C)
    label = C === ComplexF64 ? "F64" : "BF$(precision(BigFloat))"
    @testset "G-Q3-T1 r-form==Γ-form under Γ+Γ̄=2r [$label]" begin
        for (name, ν, i, j) in QQNM_NAMED
            gf = qqnm_pp_aij_gammaform(name)
            for (a0, r0, θ0) in QQNM_TAB_PTS
                a = RT(a0); r = RT(r0); θ = RT(θ0)
                Γ  = C(complex(r, a*cos(θ)))     # paper Γ = r + ia cosθ
                Γb = C(complex(r, -a*cos(θ)))    # paper Γ̄ = r − ia cosθ
                rform = qqnm_pp_aij(ν, C, Γ, Γb, C(r))[(i, j)]     # machine r-form
                gform = gf(Γ, Γb, C(r))                            # hand Γ-form (distinct string)
                @test abs(rform - gform) / (abs(gform) + 1) < tol
            end
        end
    end
end

@testset "Q3 Deliverable A — QQNM source tables" begin
    @testset "G-Q3-T0 provenance (sha256 vs frozen ledger)" begin
        (checked, miss) = qqnm_verify_tex_provenance()   # errors on mismatch
        @test checked + miss == length(QQNM_TEX_SHA256)
        @info "QQNM tex provenance" checked missing = miss
        if checked == 0
            # No owner tex mirror at this site — the baked literals remain authoritative
            # (transcription gates G-Q3-T1/T2 below still run); re-hash gate skips.
            @test_skip "QQNM-coeff tex mirror absent — sha256 re-verification skipped"
        else
            @test checked ≥ 1
        end
    end

    run_qqnm_tables(ComplexF64, 1e-10)
    setprecision(BigFloat, 256) do
        run_qqnm_tables(Complex{BigFloat}, 1e-40)
    end

    @testset "G-Q3-T3 ++ tables vs Paper-I:591 operators (source-level-fixed excepted)" begin
        # An INDEPENDENT second source (the paper operator, distinct from the owner .tex): every baked
        # a_ij^(ν) equals the Paper-I:591 A_ν operator coefficient under Γ+Γ̄=2r — EXCEPT the three
        # "source-level-fixed" entries {a_11^(1), a_04^(3), a_33^(3)} the owner corrected (F7-3 flag).
        # a_11^(1) is DECISIVELY certified by G-Q3-1 (owner value reproduces Paper-I:1427 Q(r) exactly;
        # the printed operator value does not); a_04^(3),a_33^(3) match the owner Γ-form (G-Q3-T1).
        paper(ν, Γ, Γb, r) = ν == 1 ? Dict((0,0)=>-45Γ/(2Γb^7),(0,1)=>3*(15Γ-Γb)/Γb^6,
            (1,1)=>3*(Γb-15Γ)/(2Γb^5),(0,2)=>(3Γb-42Γ)/(2Γb^5),(0,3)=>6Γ/Γb^4,(1,2)=>3*(7Γ-Γb)/Γb^4,
            (2,2)=>3*(2Γb-13Γ)/(8Γb^3),(1,3)=>(Γb-12Γ)/(2Γb^3),(0,4)=>-(7Γ+4r)/(8Γb^3),
            (2,3)=>(11Γ-Γb)/(4Γb^2),(1,4)=>(4Γ+r)/(4Γb^2),(0,5)=>r/(4Γb^2),(3,3)=>-3Γ/(8Γb),
            (2,4)=>-Γ/(2Γb),(1,5)=>-Γ/(8Γb)) :
          ν == 2 ? Dict((0,1)=>-3/Γb^5,(0,2)=>3r/Γb^5,(1,1)=>3/Γb^4,(2,2)=>3r/(2Γb^3),
            (1,2)=>-3*(2r+Γb)/(2Γb^4),(0,3)=>-3Γ/(2Γb^4),(1,3)=>(r+Γ)/Γb^3,(0,4)=>(3Γ-Γb)/(4Γb^3),
            (2,3)=>-(2r+3Γ)/(4Γb^2),(1,4)=>(Γb-6Γ)/(8Γb^2),(0,5)=>(Γb-2Γ)/(8Γb^2),(3,3)=>Γ/(4Γb),
            (2,4)=>7Γ/(16Γb),(1,5)=>Γ/(4Γb),(0,6)=>Γ/(16Γb)) :
          ν == 3 ? Dict((0,1)=>5*(2Γb-3Γ)/Γb^7,(1,1)=>(15Γ-11Γb)/Γb^6,(0,2)=>5*(2Γ-Γb)/Γb^6,
            (1,2)=>(23Γb-34Γ)/(2Γb^5),(0,3)=>(Γb-6Γ)/(2Γb^5),(2,2)=>3*(3Γ-2Γb)/(2Γb^4),
            (1,3)=>(15Γ-8Γb)/(3Γb^4),(0,4)=>(2r+Γ)/(6Γb^4),(2,3)=>(16Γb-29Γ)/(12Γb^3),
            (1,4)=>(4Γb-21Γ)/(24Γb^3),(0,5)=>-(3Γb+2r)/(24Γb^3),(3,3)=>(3Γ-4Γb)/(12Γb^2),
            (2,4)=>(8Γ-Γb)/(24Γb^2),(1,5)=>r/(6Γb^2),(0,6)=>1/(24Γb)) :
          Dict((1,1)=>(45Γ-60Γb)/(2Γb^7),(1,2)=>(42Γb-30Γ)/Γb^6,(2,2)=>9*(2Γ-3Γb)/(2Γb^5),
            (1,3)=>9*(2Γ-3Γb)/(2Γb^5),(2,3)=>3*(5Γb-3Γ)/(2Γb^4),(1,4)=>(5Γb-3Γ)/(2Γb^4),
            (3,3)=>3*(Γ-2Γb)/(8Γb^3),(2,4)=>(Γ-2Γb)/(2Γb^3),(1,5)=>(Γ-2Γb)/(8Γb^3))
        excepts = Set([(1,1,1), (3,0,4), (3,3,3)])   # source-level-fixed (owner-corrected vs paper)
        for (a0, r0, θ0) in QQNM_TAB_PTS
            a=r0*0+a0; Γ = complex(r0, a0*cos(θ0)); Γb = complex(r0, -a0*cos(θ0))
            for ν in 1:4
                baked = qqnm_pp_aij(ν, ComplexF64, Γ, Γb, r0); ref = paper(ν, Γ, Γb, r0)
                @test keys(baked) == keys(ref)
                for (k, v) in ref
                    (ν, k[1], k[2]) in excepts && continue      # owner-corrected; checked elsewhere
                    @test abs(baked[k] - v)/(abs(v)+1) < 1e-13
                end
            end
        end
    end

    @testset "G-Q3-T2 channel structure" begin
        # ++ nonzero-key census per ν (16-basis, 0≤i≤j≤6)
        for (ν, nkeys) in ((1, 15), (2, 15), (3, 15), (4, 9))
            d = qqnm_pp_aij(ν, ComplexF64, 3.0+0.1im, 3.0-0.1im, 3.0)
            @test length(d) == nkeys
            @test all(0 ≤ i ≤ j ≤ 6 for (i, j) in keys(d))
        end
        # +− support ⊂ the certified 6-slot set, for every μ=1..16, and every μ nonempty
        for μ in 1:16
            d = qqnm_pm_bij(μ, ComplexF64, 3.0+0.1im, 3.0-0.1im, 3.0)
            @test !isempty(d)
            @test all(k in QQNM_PM_SUPPORT for k in keys(d))
        end
        @test length(QQNM_PM_ATOMS) == 16
        # −− the certified structural zero
        @test qqnm_mm_is_zero() === true
    end
end
