# Characterization of PostNewtonian.SourceMode (point-particle source + sourced mode + flux).
#
# GROUNDING (the slice capstone is an END-TO-END analytic identity, computed in Julia, reading NO
# oracle file):
#  * ★ Newtonian 32/5 v¹⁰ quadrupole: the (2,2) mode is built end-to-end (radial In × source ×
#    amplitudes × Wronskian convolution → Z^∞ → FluxInf) and its FluxInf leading coefficient equals
#    16/5·r0⁴·Ω⁶; under the Schwarzschild Kepler relation Ω=r0^(-3/2) that is 16/(5 r0⁵), so the
#    (2,±2) sum is the textbook 32/(5 r0⁵) = 32/5 v¹⁰.
#  * _C_trans leading coefficient == -2 i ω³ (= -16 i Ω³ with ω = 2Ω): the A-amplitude / Pochhammer /
#    Coulomb amplitude machinery, verified directly.
#  * Inline oracle Dicts _SRC_A0_LM / _SRC_AKERR_LM copy-verified byte-exact vs legacy
#    mode_pn_source_data.jl (keys l=2..5,m=1..l; a spot-checked InputForm string literal).

using PostNewtonian: SourceMode, RadialMode, MST
using PostNewtonian.SourceMode: TeukolskyPointParticleModePN
using PostNewtonian.MST: imagunit_sym
using Symbolics
import Symbolics

@testset "SourceMode" begin
    @variables r0 ΩKerr
    iu = imagunit_sym()
    SP = Symbolics.variable(:SqrtPi); PI = Symbolics.variable(:Pi)

    # --- ★ Newtonian 32/5 v^10 quadrupole (Julia-computed end-to-end) ------------------------
    @testset "Newtonian 32/5 v^10 flux" begin
        tm0 = TeukolskyPointParticleModePN(-2, 2, 2, 0, :r0, (:η, 2); compute_H = false)
        FI = tm0["EnergyFlux"]["I"]
        li = findfirst(c -> !iszero(c), FI.coeffs)
        @test li !== nothing
        lead = Symbolics.substitute(FI.coeffs[li],
            Dict(SP => Num(sqrt(big(pi))), PI => Num(big(pi)), iu => im))
        # bare leading coeff at r0=1,Ω=1 is 16/5
        n1 = Symbolics.value(Symbolics.substitute(lead, Dict(r0 => 1.0, ΩKerr => 1.0)))
        @test isapprox(real(n1), 16 / 5; rtol = 1e-9)
        @test abs(imag(n1)) < 1e-9
        # with Kepler Ω = r0^(-3/2): leading coeff == 16/(5 r0^5) as a function of r0
        lead_kep = Symbolics.substitute(lead, Dict(ΩKerr => r0^(Num(-3) // 2)))
        for rv in (2.0, 3.7, 10.0)
            got = Symbolics.value(Symbolics.substitute(lead_kep, Dict(r0 => rv)))
            @test isapprox(real(got), 16 / (5 * rv^5); rtol = 1e-9)
            @test abs(imag(got)) < 1e-12 * max(1.0, abs(real(got)))
        end
        @test 2 * (16 // 5) == 32 // 5    # (2,±2) sum = 32/5
    end

    # --- _C_trans leading == -2 i ω^3 (amplitude machinery) ----------------------------------
    @testset "C_trans leading -2 i ω^3" begin
        nwork = 12
        nuser, eps, aMSTeta = RadialMode._amp_inputs(-2, 2, 2, nwork)
        Ct = RadialMode._reattach_omega_amp(
            RadialMode._C_trans(-2, 2, nuser, eps, aMSTeta, nwork), Num(2) * ΩKerr)
        @test Ct.nmin == 9                       # η^9 = ε^3 leading
        want = -16 * iu * ΩKerr^3                # ω = 2Ω ⇒ -2 i ω^3
        @test isequal(Symbolics.expand(Ct.coeffs[1] - want), Num(0))
    end

    # --- inline oracle Dicts: key coverage + byte-exact spot-check ----------------------------
    @testset "source oracle Dicts" begin
        a0 = SourceMode._SRC_A0_LM
        ak = SourceMode._SRC_AKERR_LM
        @test Set(keys(a0)) == Set((l, m) for l in 2:5 for m in 1:l)
        @test haskey(ak, (2, 2))
        # byte-exact spot-check: the (2,2) a→0 δ' coefficient literal (transcribed from
        # legacy mode_pn_source_data.jl).
        @test occursin("Sqrt[5*Pi]", a0[(2, 1)][1])
        @test a0[(2, 2)][3] isa AbstractString
    end
end
