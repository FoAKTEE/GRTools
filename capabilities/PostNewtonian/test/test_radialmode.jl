# Characterization of PostNewtonian.RadialMode (η-graded symbolic radial + amplitude assembly).
#
# GROUNDING:
#  * Symbolic special-function CLOSED VALUES (analytic): the digamma at positive integers and the
#    Hurwitz/Riemann ζ — ψ(1)=-γ, ψ(2)=1-γ, ζ(2)=π²/6. These are exact and independent of the port.
#  * Captured one-shot fixtures from LEGACY TeukolskyPN.TeukolskyRadialPN (baked numeric literals,
#    NOT a live test dep): the In/Up radial pair structure (nmin/nmax) and leading coefficients for
#    Schwarzschild s=-2,l=2,m=2 — In leading 4/5·r⁴ω⁴ at η⁴, Up leading -3𝕚/(rω) at η⁻¹.

using PostNewtonian: RadialMode, MST
using PostNewtonian.RadialMode: TeukolskyRadialPN
using PostNewtonian.MST: imagunit_sym
using Symbolics
import Symbolics

@testset "RadialMode" begin
    EG = Symbolics.variable(:EulerGamma)
    PI = Symbolics.variable(:Pi)

    # --- symbolic special-fn closed values (analytic) ----------------------------------------
    @testset "special-fn closed values" begin
        @test isequal(RadialMode._digamma_int(1), -EG)              # ψ(1) = -γ
        @test isequal(Symbolics.expand(RadialMode._digamma_int(2) - (1 - EG)), Num(0))  # ψ(2)=1-γ
        @test isequal(Symbolics.expand(RadialMode._hurwitz_zeta(2, 1) - (1 // 6) * PI^2), Num(0))  # ζ(2)=π²/6
    end

    # --- TeukolskyRadialPN In/Up pair (captured legacy fixture) -------------------------------
    @testset "TeukolskyRadialPN Schwarzschild l2m2" begin
        R = TeukolskyRadialPN(-2, 2, 2, 0 // 1, :ω, (:η, 4))
        @test R isa Dict{String,PNSeries}
        rin = R["In"]
        rup = R["Up"]
        @test (rin.nmin, rin.nmax) == (4, 8)
        @test (rup.nmin, rup.nmax) == (-1, 3)

        iu = imagunit_sym()
        r = Symbolics.variable(:r); ω = Symbolics.variable(:ω)
        sub = Dict(r => 1.0, ω => 1.0, iu => im)
        # In leading: 4/5 r⁴ ω⁴  -> 0.8
        cin = Symbolics.value(Symbolics.substitute(rin.coeffs[1], sub))
        @test isapprox(real(ComplexF64(cin)), 0.8; atol = 1e-25)
        # Up leading: -3𝕚/(rω) -> -3im
        cup = Symbolics.value(Symbolics.substitute(rup.coeffs[1], sub))
        @test isapprox(ComplexF64(cup), -3im; atol = 1e-12)
    end
end
