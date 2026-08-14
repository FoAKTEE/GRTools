# CoupledSets — leaf solve-plan layer for the Lorenz field solver: the per-parity CoupledSet
# descriptor + the CSETS table + parity selection + the γ↔(l,m) mode enumeration.
#
# Grounding (no live legacy dep — pure wiring/enumeration vs hand-enumerated tables):
#   • select_coupled_set vs the hand-enumerated parity table that mirrors H1Lorenz
#     coupled_sets.jl select_correct_coupled_set (even→even / odd→odd / l==1 even→even_dipole /
#     l==1 odd→odd_dipole / l==0 even→monopole; the (l+m)even,m==0 even-static swap is the caller's
#     responsibility, exposed via select_coupled_set(l,m; static=true)).
#   • ode!-field OBJECT IDENTITY: each CoupledSet.ode! === the corresponding consumed
#     MetricReconstruction.Lorenz.LorenzRadial RHS (gates the monopole_ode! arity-4 canonical form).
#   • γ↔(l,m) round-trip gamma_to_lm∘lm_to_gamma == identity, and mode_list(lmax) length ==
#     (lmax+1)(lmax+2)/2 enumerated by hand.

using Test
using SelfForce
import MetricReconstruction
const LR = MetricReconstruction.Lorenz.LorenzRadial
const LF = SelfForce.LorenzField

@testset "LorenzField.CoupledSets" begin

    @testset "select_coupled_set parity table (hand-enumerated)" begin
        # (l+m) even ⇒ even sector; l==0 monopole, l==1 even_dipole, else even.
        @test LF.select_coupled_set(0, 0).kind === :monopole
        @test LF.select_coupled_set(1, 1).kind === :even_dipole   # l+m=2 even, l==1
        @test LF.select_coupled_set(2, 0).kind === :even          # l+m even, l>=2
        @test LF.select_coupled_set(2, 2).kind === :even
        @test LF.select_coupled_set(3, 1).kind === :even
        @test LF.select_coupled_set(4, 2).kind === :even
        # (l+m) odd ⇒ odd sector; l==1 odd_dipole, else odd.
        @test LF.select_coupled_set(1, 0).kind === :odd_dipole    # l+m=1 odd, l==1
        @test LF.select_coupled_set(2, 1).kind === :odd           # l+m odd, l>=2
        @test LF.select_coupled_set(3, 0).kind === :odd
        @test LF.select_coupled_set(3, 2).kind === :odd
        # the even-static swap (l+m even, m==0, n==0 ⇒ ω==0): exposed via static=true
        @test LF.select_coupled_set(2, 0; static=true).kind === :even_static
        @test LF.select_coupled_set(4, 0; static=true).kind === :even_static
        # static flag must NOT swap when m != 0 or l<2 (caller only swaps (l+m)even & m==0 & l>=2)
        @test LF.select_coupled_set(2, 2; static=true).kind === :even
        @test LF.select_coupled_set(0, 0; static=true).kind === :monopole
        @test LF.select_coupled_set(1, 1; static=true).kind === :even_dipole
        @test LF.select_coupled_set(3, 0; static=true).kind === :odd   # l+m odd ⇒ never static-even
    end

    @testset "ode! field === consumed LorenzRadial RHS (object identity)" begin
        @test LF.CSETS.even.ode!        === LR.even_ode!
        @test LF.CSETS.even_dipole.ode! === LR.even_dipole_ode!
        @test LF.CSETS.even_static.ode! === LR.even_static_ode!
        @test LF.CSETS.monopole.ode!    === LR.monopole_ode!     # arity-4 canonical generic fn
        @test LF.CSETS.odd.ode!         === LR.odd_ode!
        @test LF.CSETS.odd_dipole.ode!  === LR.odd_dipole_ode!
        # select_coupled_set returns the SAME descriptor objects (===) as the CSETS table.
        @test LF.select_coupled_set(2, 0).ode! === LR.even_ode!
        @test LF.select_coupled_set(2, 1).ode! === LR.odd_ode!
        @test LF.select_coupled_set(0, 0).ode! === LR.monopole_ode!
        @test LF.select_coupled_set(1, 0).ode! === LR.odd_dipole_ode!
        @test LF.select_coupled_set(1, 1).ode! === LR.even_dipole_ode!
    end

    @testset "CoupledSet field counts (num_coupled / num_gauge)" begin
        # the BLS coupled-field counts per sector (H1Lorenz coupled_sets.jl literals).
        @test (LF.CSETS.even.num_coupled,        LF.CSETS.even.num_gauge)        == (5, 2)
        @test (LF.CSETS.odd.num_coupled,         LF.CSETS.odd.num_gauge)         == (2, 1)
        @test (LF.CSETS.odd_dipole.num_coupled,  LF.CSETS.odd_dipole.num_gauge)  == (1, 1)
        @test (LF.CSETS.monopole.num_coupled,    LF.CSETS.monopole.num_gauge)    == (3, 1)
        @test (LF.CSETS.even_dipole.num_coupled, LF.CSETS.even_dipole.num_gauge) == (4, 2)
        @test (LF.CSETS.even_static.num_coupled, LF.CSETS.even_static.num_gauge) == (3, 4)
        # sources vector length matches num_coupled
        @test length(LF.CSETS.even.sources)        == 5
        @test length(LF.CSETS.odd.sources)         == 2
        @test length(LF.CSETS.odd_dipole.sources)  == 1
        @test length(LF.CSETS.monopole.sources)    == 3
        @test length(LF.CSETS.even_dipole.sources) == 4
        @test length(LF.CSETS.even_static.sources) == 3
    end

    @testset "γ ↔ (l,m) enumeration round-trip + mode_list length" begin
        for l in 0:8, m in 0:l
            γ = LF.lm_to_gamma(l, m)
            @test LF.gamma_to_lm(γ) == (l, m)
        end
        # γ enumerates 0,1,2,... contiguously over mode_list ordering
        modes = LF.mode_list(8)
        for (i, (l, m)) in enumerate(modes)
            @test LF.lm_to_gamma(l, m) == i - 1
        end
        for lmax in 0:8
            @test length(LF.mode_list(lmax)) == (lmax + 1) * (lmax + 2) ÷ 2
        end
    end

    @testset "out-BC wrappers route to the right BC series (uniform signature)" begin
        # the _out_* closures present a uniform (a0s,l,ω,p; rs_out) signature; each routes to the
        # right BC series. Smoke: ω=0 static-friendly call returns (R, dR, rs) with R-length ==
        # num_coupled for monopole/odd_dipole (the simple ω-independent sectors).
        T = Float64
        Rm, dRm, _ = LF.CSETS.monopole.out_bc!(T[1,0,0], 0, T(0), T(10))
        @test length(Rm) == 3
        Rod, dRod, _ = LF.CSETS.odd_dipole.out_bc!(T[1,0], 1, T(0), T(10))
        @test length(Rod) == 1
    end
end
