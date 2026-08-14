# Grounding for psi4_from_hertz (Radpsi4fromHzpotmode, Sourceincoordinates.nb SC l.408/436).
#
# Physics oracle (Shah & van de Meent 2015 reconstruction identity, notebook IN173-176/IN188):
# the ψ4 reconstructed from the (l,-m,-n,-k) Hertz potential equals the point-particle ψ4
# amplitude for (l,m,n,k) divided by r^4, in the vacuum region r>p. The reference values are the
# BHPToolkit Teukolsky paclet `TeukolskyPointParticleMode[-2,l,m,orbit][r]/r^4` captured by
# oracle/psi4_from_hertz_oracle.wl (which also derives the reduced closed-form kernel xAct-free,
# Psi4check + DerHzpotsubstf, and checks recon vs toolkit: reldiff 6.7e-6 @r=20, 7e-8 @r=50,
# 4e-9 @r=100 — the reconstruction converges to the toolkit ψ4 away from the particle). The Hertz
# front-end runs via the Perturbation GRTOOLS live dep; paclet values are frozen literals. Values +
# tolerances transcribed verbatim from the legacy weyl_from_hertz_test.jl. No value from new code.
#
# Circular Schwarzschild a=0, p=10, l=2, m=2, n=k=0. Toolkit ψ4/r^4 reference (paclet):
#   r=20  1.4806277372676876e-5 + 1.5593229011404426e-5 i
#   r=50  1.436133987182408e-5  + 1.5130981988029799e-5 i
#   r=100 -9.428525966798886e-6 - 6.27623614400317e-6  i

using ModifiedGravity, Test

@testset "psi4_from_hertz" begin
    f = psi4_from_hertz(2, -2, 0, 0, 0.0, 10.0, 0.0, 1.0)   # reconstruct (l,-m,-n,-k)
    @testset "reconstruction == Teukolsky ψ4 / r^4 (vacuum)" begin
        # tolerance follows the reconstruction convergence (tighter farther from the particle)
        @test f(20.0)  ≈ (1.4806277372676876e-5 + 1.5593229011404426e-5im) rtol=1e-5
        @test f(50.0)  ≈ (1.436133987182408e-5  + 1.5130981988029799e-5im) rtol=1e-6
        @test f(100.0) ≈ (-9.428525966798886e-6 - 6.27623614400317e-6im)   rtol=1e-6
    end

    # h44_from_hertz (Radh44fromHzpotmode, SC l.368/430/432): reconstructed (4,4)=m̄m̄ tetrad metric,
    # native GHP assembly (thorn_prime + schw_spincoeff, ghp_operators.jl). Two grounds land together:
    #
    #  (a) NECESSARY — h44 → 2·ψ4/ω² as r→∞ (Radh44asympt, SC l.438): the m̄m̄ metric mode asymptotes to
    #      twice the reconstructed ψ4 over ω² (ψ4=−½∂_t²h44 at leading order). Checked vs
    #      psi4_from_hertz (BHPToolkit-grounded) at large r — validates signcoef=(−1)^m, the ξ⁴=r⁴
    #      Radf† profile, the GHP weight-sum SΨ=4, and the leading þ′þ′ operator structure.
    #
    #  (b) SUFFICIENT (finite r) — vs the independent wolframscript oracle
    #      oracle/h44_from_hertz_oracle.wl: the SAME recovered operator (SC l.368) composed SYMBOLICALLY
    #      (Mathematica D[] chain/product rule) and evaluated numerically with the BHPToolkit Teukolsky
    #      paclet RadHZ† (Up radial) — a fully independent implementation of the assembly (symbolic D[]
    #      vs Julia jet towers + the s=+2 radial-Teukolsky ODE for RadHZ†″). The þ′ operator itself is
    #      grounded standalone in test_ghp_operators.jl (GHP commutator + Ricci identity), so this
    #      certifies the finite-r ASSEMBLY. Reference literals (l=2,m=−2, a=0, p=10) transcribed from
    #      the oracle print (double-precision paclet). No expected value taken from the Julia code.
    @testset "h44_from_hertz: native GHP m̄m̄ assembly" begin
        fh = h44_from_hertz(2, -2, 0, 0, 0.0, 10.0, 0.0, 1.0)
        _, _, ω = ModifiedGravity.hertz_potential_mode(2, -2, 0, 0, 0.0, 10.0, 0.0, 1.0)

        @testset "sufficient: finite-r vs wolframscript operator oracle" begin
            # oracle/h44_from_hertz_oracle.wl printout, mode (l,m)=(2,−2):
            @test fh(12.0)  ≈ (0.08855530028315524  - 0.07645267029093028im)  rtol=1e-6
            @test fh(20.0)  ≈ (0.03765356624708456  - 0.024807065873447047im) rtol=1e-6
            @test fh(50.0)  ≈ (0.009251822349186273 + 0.008714717919847483im) rtol=1e-6
            @test fh(100.0) ≈ (-0.0049849139989386415 - 0.003245793126587966im) rtol=1e-6
        end

        @testset "necessary: h44 → 2·ψ4/ω² asymptotic (Radh44asympt)" begin
            fp = psi4_from_hertz(2, -2, 0, 0, 0.0, 10.0, 0.0, 1.0)
            @test fh(2000.0)  ≈ 2*fp(2000.0)/ω^2  rtol=1e-3
            @test fh(40000.0) ≈ 2*fp(40000.0)/ω^2 rtol=1e-5
        end
    end
end
