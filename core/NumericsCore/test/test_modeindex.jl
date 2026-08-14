# Characterization of the mode-index primitives.
#
# Oracle 1 (gamma<->lm): H1Lorenz test-dep DROPPED (pulls HDF5); the flat work-counter is the closed
#   form gamma = l(l+1)/2 + m, m in 0:l (H1Lorenz driver.jl:17-26, read not run) — grounded by that
#   closed form + decode round-trip + literal reference rows.
# Oracle 2 (ModeIndexTable): heavy KerrBackwards test-dep DROPPED; re-grounded WITHOUT a live dep by
#   DETERMINISTIC ENUMERATION of the documented FEW SchwarzschildEccentric construction
#   (ref-code .../few/utils/baseclasses.py:223-323, ported in packages/KerrBackwards/src/modeindex.jl):
#   - closed-form counts: l in 2:10, m in 0:l, n in -30:30 (61 vals) =>
#       num_modes = 61*sum_{l=2}^{10}(l+1) = 61*63 = 3843; num_m0 = 9*61 = 549;
#       num_m_1_up = 3843-549 = 3294; full length = 3843+3294 = 7137.
#   - m0sort + mirror structure invariants (m==0 first, then m>0, then the m>0 block mirrored with m<0).
#   - literal reference rows traced from the construction (e.g. row 1 = (2,0,-30), row 550 = (2,1,-30)).
# Oracle 3 (a_il, ladder): documented closed forms — SecondOrderRicci `a_il` (base 2^(-1/2), /sqrt(l(l+1))
#   for i in {4,5,8,9}, /sqrt((l-1)l(l+1)(l+2)) for i in {7,10}); ladder sqrt(l(l+1)),
#   sqrt((l-1)l(l+1)(l+2)).
# Every expected value is a live LIGHT-oracle value or an analytic relation — never the code's own output.

using NumericsCore

@testset "modeindex" begin
    @testset "gamma<->lm closed form (gamma=l(l+1)/2+m) + round-trip" begin
        # literal reference rows (H1Lorenz driver.jl: lm_to_gamma(l,m)=l*(l+1)÷2+m; m in 0:l):
        @test NumericsCore.lm_to_gamma(0, 0) == 0
        @test NumericsCore.lm_to_gamma(1, 0) == 1
        @test NumericsCore.lm_to_gamma(2, 0) == 3
        @test NumericsCore.lm_to_gamma(2, 2) == 5
        @test NumericsCore.gamma_to_lm(0) == (0, 0)
        @test NumericsCore.gamma_to_lm(5) == (2, 2)
        for g in 0:200
            l, m = NumericsCore.gamma_to_lm(g)
            @test 0 <= m <= l                                  # selection rule m in 0:l
            @test l*(l+1) ÷ 2 + m == g                         # decode matches the closed form
            @test NumericsCore.lm_to_gamma(l, m) == g          # round-trip gamma->lm->gamma
        end
        for l in 0:15, m in 0:l
            @test NumericsCore.lm_to_gamma(l, m) == l*(l+1) ÷ 2 + m   # encode closed form
            @test NumericsCore.gamma_to_lm(NumericsCore.lm_to_gamma(l, m)) == (l, m)
        end
        @test NumericsCore.mode_list(6) == [(l, m) for l in 0:6 for m in 0:l]   # H1Lorenz mode_list form (m in 0:l)
    end

    @testset "ModeIndexTable analytic enumeration (lmax=10, nmax=30)" begin
        t = NumericsCore.mode_index_table(lmax=10, nmax=30)

        # --- closed-form counts (61 n-values; sum_{l=2}^{10}(l+1)=63) ---
        @test t.num_modes     == 61 * 63          # = 3843
        @test t.num_m0        == 9 * 61           # = 549  (one m==0 column per l in 2:10)
        @test t.num_m_1_up    == 3843 - 549       # = 3294
        @test t.num_m_zero_up == 3843
        @test length(t.l_arr) == 3843 + 3294      # = 7137
        @test length(t.m_arr) == 7137
        @test length(t.n_arr) == 7137

        # --- selection rules: l in 2:10, |m|<=l, |n|<=30 ---
        @test all(2 .<= t.l_arr .<= 10)
        @test all(abs.(t.m_arr) .<= t.l_arr)
        @test all(abs.(t.n_arr) .<= 30)

        # --- m0sort + mirror block structure ---
        nm  = t.num_modes; n0 = t.num_m0
        @test all(t.m_arr[1:n0] .== 0)                         # m==0 block first
        @test all(t.m_arr[n0+1:nm] .> 0)                       # then the m>0 block
        @test t.l_arr[nm+1:end] == t.l_arr[n0+1:nm]            # appended block mirrors m>0 block: l same
        @test t.m_arr[nm+1:end] == -t.m_arr[n0+1:nm]           #   m negated
        @test t.n_arr[nm+1:end] == t.n_arr[n0+1:nm]            #   n same

        # --- literal reference rows traced from the documented (l outer, m middle, n inner) enumeration ---
        @test (t.l_arr[1],    t.m_arr[1],    t.n_arr[1])    == (2,  0, -30)   # first m==0 entry
        @test (t.l_arr[549],  t.m_arr[549],  t.n_arr[549])  == (10, 0,  30)   # last  m==0 entry
        @test (t.l_arr[550],  t.m_arr[550],  t.n_arr[550])  == (2,  1, -30)   # first m>0  entry
        @test (t.l_arr[3843], t.m_arr[3843], t.n_arr[3843]) == (10, 10, 30)   # last  m>0  entry
        @test (t.l_arr[3844], t.m_arr[3844], t.n_arr[3844]) == (2, -1, -30)   # first mirrored (m<0)
        @test (t.l_arr[7137], t.m_arr[7137], t.n_arr[7137]) == (10,-10, 30)   # last  mirrored (m<0)

        # --- index_map round-trips every row; keys are unique (1-based positions) ---
        @test length(t.index_map) == 7137
        @test all(t.index_map[(t.l_arr[i], t.m_arr[i], t.n_arr[i])] == i for i in 1:7137)
        @test t.index_map[(2, 0, -30)] == 1
        @test t.index_map[(2, 1, -30)] == 550

        # --- special_index_map: m>=0 -> self; mirrored m<0 -> its m>0 partner (i - num_m_1_up) ---
        @test t.special_index_map[(2, 0, -30)] == 1
        @test t.special_index_map[(2, 1, -30)] == 550        # m>0 maps to itself
        @test t.special_index_map[(2,-1, -30)] == 550        # m<0 maps to its m>0 partner
        @test all(t.special_index_map[(t.l_arr[i], t.m_arr[i], t.n_arr[i])] == i for i in 1:nm)
        @test all(t.special_index_map[(t.l_arr[i], t.m_arr[i], t.n_arr[i])] == i - t.num_m_1_up
                  for i in nm+1:7137)
    end

    @testset "a_il closed form (SecondOrderRicci)" begin
        base = sqrt(0.5)
        for l in 2:8
            for i in (1, 2, 3, 6)
                @test NumericsCore.a_il(i, l) == base
            end
            for i in (4, 5, 8, 9)
                @test NumericsCore.a_il(i, l) ≈ base/sqrt(l*(l+1)) rtol=1e-14
            end
            for i in (7, 10)
                @test NumericsCore.a_il(i, l) ≈ base/sqrt((l-1)*l*(l+1)*(l+2)) rtol=1e-14
            end
        end
        # BigFloat precision inheritance
        ab = NumericsCore.a_il(BigFloat, 4, 5)
        @test ab isa BigFloat
        @test ab ≈ sqrt(big(0.5))/sqrt(big(5*6)) rtol=1e-50
    end

    @testset "ladder factors closed form" begin
        for l in 2:10
            @test NumericsCore.lambdahat(l) ≈ sqrt(l*(l+1.0)) rtol=1e-14
            @test NumericsCore.Lambdahat(l) ≈ sqrt((l-1.0)*l*(l+1.0)*(l+2.0)) rtol=1e-14
            # a_il consistency: i=4 ladder == base/lambdahat
            @test NumericsCore.a_il(4, l) ≈ sqrt(0.5)/NumericsCore.lambdahat(l) rtol=1e-13
            @test NumericsCore.a_il(7, l) ≈ sqrt(0.5)/NumericsCore.Lambdahat(l) rtol=1e-13
        end
        @test NumericsCore.lambdahat(BigFloat, 5) isa BigFloat
    end
end
