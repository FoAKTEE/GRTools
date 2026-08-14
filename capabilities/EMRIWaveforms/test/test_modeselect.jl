# ModeSelect characterization.
#
# Grounding:
#  (A) ANALYTIC hand-computed mode power — single-point K=3 sets whose per-mode |A·y|² are exact integers,
#      so the cumulative-(1−eps) keep/drop ordering is hand-derivable and asserted.
#  (B) CAPTURED FIXTURE — legacy KerrBackwards.mode_selector run ONCE (heavy pkg, never a live dep) on a
#      fixed nt=2,K=4 slice; kept indices baked verbatim.
#  (C) select_modes_by_list set-membership identity + FEW-fork "all"/list/nothing branches.

using Test
import EMRIWaveforms.ModeSelect
import EMRIWaveforms.ModeSelect: select_modes_by_list
const KBsel = ModeSelect.KerrBackwards
const FEWsel = ModeSelect.FEW

@testset "ModeSelect" begin
    # ---- (A) hand-computed power threshold ----
    @testset "hand-computed power threshold" begin
        # nt=1, K=3, all m=0 (no partners): |A·y|² = [1, 4, 0.25], total 5.25.
        teuk = ComplexF64[1.0 2.0 0.5]
        ylms = ComplexF64[1, 1, 1]
        l = [2, 3, 4]; m = [0, 0, 0]; n = [0, 0, 0]
        # eps=0.05 ⇒ thresh 4.9875: keep strongest (mode2,P=4) then mode1 (cum 4<thresh); drop mode3.
        @test KBsel.mode_selector(teuk, ylms, l, m, n; eps = 0.05).keep_modes == [1, 2]
        # eps=0.5 ⇒ thresh 2.625: keep only mode2 (cum 4 ≥ thresh after first).
        r05 = KBsel.mode_selector(teuk, ylms, l, m, n; eps = 0.5)
        @test r05.keep_modes == [2]
        @test length(r05.ylms) == 2 * length(r05.keep_modes)   # ylmkeep = [keep, temp2] layout
        # eps→0 keeps every mode.
        @test KBsel.mode_selector(teuk, ylms, l, m, n; eps = 1e-12).keep_modes == [1, 2, 3]
    end

    # ---- (A) partner (m<0) keeps the shared base mode ----
    @testset "partner pairing keeps base mode" begin
        # K=2: mode1 m=0, mode2 m=2. powers base=[0.1, 0.5], partner of mode2 = 9 (dominant).
        teuk = ComplexF64[0.1 0.5]              # |A|·|y|=|A| with ylms below
        # base ylms: [1,1]; partner ylm for the single m>0 mode chosen so |conj(A2)·yp|² = 9.
        ylms = ComplexF64[1, 1, 18.0]           # conj(A2)=0.5 ⇒ |0.5·18|²=81? recompute below
        l = [2, 2]; m = [0, 2]; n = [0, 0]
        # |A2|=0.5 ⇒ partner power = |conj(0.5)·yp|² = (0.5·yp)². Want 9 ⇒ yp=6.
        ylms = ComplexF64[1, 1, 6.0]
        # row = [|0.1·1|², |0.5·1|², (0.5·6)²] = [0.01, 0.25, 9.0]; total 9.26.
        # eps=0.1 ⇒ thresh 8.334: first kept element is the partner (P=9) ⇒ base mode2 kept; cum 9 ≥ thresh.
        r = KBsel.mode_selector(teuk, ylms, l, m, n; eps = 0.1)
        @test r.keep_modes == [2]
    end

    @testset "partner mapping is invariant under arbitrary base-mode permutation" begin
        # Canonical FEW order has m=0 first. The dominant contribution is the negative-m partner of
        # physical mode id=30 (base index 3), so the selector must map that partner back to id=30.
        teuk = ComplexF64[0.1 0.2 0.5 0.4]
        l = [10, 20, 30, 40]; m = [0, 0, 2, 1]; n = [1, 2, 3, 4]
        base_y = ComplexF64[1, 2, 3, 4]
        partner_y = ComplexF64[6, 0.1] # partners of base indices 3 and 4
        ylms = vcat(base_y, partner_y)
        ref = KBsel.mode_selector(teuk, ylms, l, m, n; eps = 0.4)
        @test ref.keep_modes == [3]
        @test ref.l_arr == [30]
        @test ref.ylms == ComplexF64[3, 6]

        # Interleave positive and zero m. The partner block follows positive modes in their new
        # actual order, independently of where those positive modes occur in the base block.
        perm = [3, 1, 4, 2]
        positive_old = perm[findall(>(0), m[perm])]
        old_partner_col = Dict(3 => 1, 4 => 2)
        ylms_perm = vcat(base_y[perm], partner_y[[old_partner_col[j] for j in positive_old]])
        got = KBsel.mode_selector(teuk[:, perm], ylms_perm, l[perm], m[perm], n[perm]; eps = 0.4)
        @test got.keep_modes == [1]
        @test got.l_arr == ref.l_arr
        @test got.m_arr == ref.m_arr
        @test got.n_arr == ref.n_arr
        @test got.ylms == ref.ylms
    end

    # ---- (B) captured legacy KerrBackwards.mode_selector ----
    @testset "captured KB fixture" begin
        teuk = ComplexF64[1.0 2.0 0.5 3.0;
                          0.2 4.0 1.5 0.1]
        l = [2, 2, 3, 4]; m = [0, 2, 1, 0]; n = [0, 0, 0, 0]
        ylms = ComplexF64[1, 1, 1, 1, 0.7, 0.3]
        @test KBsel.mode_selector(teuk, ylms, l, m, n; eps = 1e-2).keep_modes == [1, 2, 3, 4]
        @test KBsel.mode_selector(teuk, ylms, l, m, n; eps = 1e-5).keep_modes == [1, 2, 3, 4]
        @test KBsel.mode_selector(teuk, ylms, l, m, n; eps = 0.5).keep_modes  == [2, 4]
    end

    # ---- (C) select_modes_by_list membership identity ----
    @testset "select_modes_by_list" begin
        l = [2, 2, 3, 4]; m = [0, 2, 1, 0]; k = [0, 0, 0, 0]; n = [0, 0, 0, 0]
        @test select_modes_by_list(l, m, k, n, [(2, 2, 0, 0), (4, 0, 0, 0)]) == [2, 4]
        @test select_modes_by_list(l, m, k, n, [(3, 1, 0, 0)]) == [3]
        @test_throws ArgumentError select_modes_by_list(l, m, k, n, [(9, 9, 0, 0)])
    end

    # ---- (C) FEW fork branches ----
    @testset "FEW fork" begin
        teuk = ComplexF64[1.0 2.0 0.5 3.0;
                          0.2 4.0 1.5 0.1]
        l = [2, 2, 3, 4]; m = [0, 2, 1, 0]; k = [0, 0, 0, 0]; n = [0, 0, 0, 0]
        ylms = ComplexF64[1, 1, 1, 1, 0.7, 0.3]
        # nothing-branch delegates to KB (same keep_modes), threads k through.
        rnothing = FEWsel.mode_selector(teuk, ylms, l, m, k, n; eps = 0.5)
        @test rnothing.keep_modes == [2, 4]
        @test rnothing.k_arr == [0, 0]
        # "all" keeps every mode.
        rall = FEWsel.mode_selector(teuk, ylms, l, m, k, n; mode_selection = "all")
        @test rall.keep_modes == [1, 2, 3, 4]
        @test rall.ylms == ComplexF64[1, 1, 1, 1, 1, 0.7, 0.3, 1]
        # explicit list ≡ select_modes_by_list.
        rlist = FEWsel.mode_selector(teuk, ylms, l, m, k, n; mode_selection = [(2, 2, 0, 0), (4, 0, 0, 0)])
        @test rlist.keep_modes == [2, 4]
        @test rlist.ylms == ComplexF64[1, 1, 0.7, 1]
        @test_throws ArgumentError FEWsel.mode_selector(teuk, ylms, l, m, k, n; mode_selection = 3)
        @test_throws ArgumentError FEWsel.mode_selector(teuk, ylms, l, m, k, n; eps = -0.1)
        @test_throws ArgumentError KBsel.mode_selector(teuk, ylms, l, m, n; eps = Inf)
    end
end
