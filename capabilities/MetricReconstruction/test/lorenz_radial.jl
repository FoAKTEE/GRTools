# Grounding for Lorenz.LorenzRadial — the a=0 first-order Lorenz-gauge radial field-equation RHS.
#
# Expected values are LITERAL fixtures captured one-shot from the legacy oracles (heavy pkgs are
# NEVER live test-deps): the self-contained pure-arithmetic ODE source of H1Lorenz (orbit.jl +
# ode_systems.jl) and DiscModes (lorenz_metric.jl), included verbatim into bare modules and
# evaluated at (rs=12.3, ω=0.35) on a deterministic state vector y[k]=sin(0.7k)+0.3cos(1.1k).
# Capture script: scratchpad/capture_lorenzradial.jl (run with bare `julia`, prints these arrays).
#
# Two grounding axes:
#   (1) ARITY GATE — H1Lorenz monopole_ode!(g,rs,y,l,ω) is bit-identical across l∈{0,1,2,5,17}
#       (l is a DEAD argument) AND bit-equals DiscModes' arity-4 monopole_ode!(g,rs,y,ω). The
#       captured MONO_H1_l* arrays were verified identical to MONO_DM ⇒ the fork collapse to a
#       canonical arity-4 + trivial arity-5 method is sound.
#   (2) RHS VALUES — even/even_dipole/even_static/odd/odd_dipole/monopole RHS match the captured
#       legacy arrays. The GrTool tortoise inverse (Spacetime.tortoise_schw_inverse) was checked
#       bit-identical to the legacy r_of_rs at rs=12.3, so equality holds to full precision.

using Test
using MetricReconstruction.Lorenz.LorenzRadial: even_ode!, even_dipole_ode!, even_static_ode!,
                                                 monopole_ode!, odd_ode!, odd_dipole_ode!

@testset "Lorenz.LorenzRadial (a=0 Lorenz-gauge radial RHS)" begin
    rs = 12.3
    ω  = 0.35
    mkstate(n) = [sin(0.7k) + 0.3cos(1.1k) for k in 1:n]
    rtol = 1e-12

    # ---- captured legacy fixtures (literals) ----
    EVEN_L2 = [0.8088993948118564, -0.06420789882432898, 0.24278828916237932, -0.06318100691815223, -0.5865059948260289, -0.006711417022450969, -0.8745945420908184, 0.056519178922452874, 0.6583143081152043, -0.05563936111383816, 1.0963640953804163, -0.28191729522193987, -0.6523650043180811, -0.03507248609826565, -0.8844546026755444, 0.06720538595808578, 0.2080195907654674, 0.044914257403586934, 0.6906191077764792, -0.131900327124877]
    EVEN_DIPOLE_L1 = [0.8088993948118564, -0.09090443025850833, 0.24278828916237932, -0.08257877445468229, -0.5865059948260289, -0.001983742502610285, -0.8745945420908184, 0.08855791196634868, 0.6583143081152043, -0.002727538930031005, 1.0963640953804163, -0.17772937283943474, -0.6523650043180811, -0.044325957435229435, -0.8844546026755444, 0.1045120945787063]
    EVEN_STATIC_L2 = [0.8088993948118564, -0.02074027603022273, 0.24278828916237932, -0.03606911197877775, -0.5865059948260289, 0.07846383302148201, -0.8745945420908184, 0.024753630547624805, 0.6583143081152043, -0.19742622762456075, 1.0963640953804163, -0.2311552629833652]
    ODD_L2 = [0.8088993948118564, -0.038314130498883076, 0.24278828916237932, -0.017995627121819853, -0.5865059948260289, -0.041439031428473175, -0.8745945420908184, 0.04221512259669343]
    ODD_DIPOLE_L1 = [0.8088993948118564, -0.06663762351183755, 0.24278828916237932, -0.048419066471478076]
    MONO = [0.8088993948118564, -0.10409574338988456, 0.24278828916237932, -0.09169460244986237, -0.5865059948260289, -0.000546450398322313, -0.8745945420908184, 0.10654627730698177, 0.6583143081152043, 0.012945630008524067, 1.0963640953804163, -0.15814621912704638]

    @testset "RHS values vs legacy fixture" begin
        let y = mkstate(20), g = zeros(20)
            even_ode!(g, rs, y, 2, ω); @test g ≈ EVEN_L2 rtol=rtol
        end
        let y = mkstate(16), g = zeros(16)
            even_dipole_ode!(g, rs, y, 1, ω); @test g ≈ EVEN_DIPOLE_L1 rtol=rtol
        end
        let y = mkstate(12), g = zeros(12)
            even_static_ode!(g, rs, y, 2, 0.0); @test g ≈ EVEN_STATIC_L2 rtol=rtol
        end
        let y = mkstate(8), g = zeros(8)
            odd_ode!(g, rs, y, 2, ω); @test g ≈ ODD_L2 rtol=rtol
        end
        let y = mkstate(4), g = zeros(4)
            odd_dipole_ode!(g, rs, y, 1, ω); @test g ≈ ODD_DIPOLE_L1 rtol=rtol
        end
        let y = mkstate(12), g = zeros(12)
            monopole_ode!(g, rs, y, ω); @test g ≈ MONO rtol=rtol
        end
    end

    @testset "monopole arity gate (l-invariance + arity-4≡arity-5)" begin
        y = mkstate(12)
        g4 = zeros(12); monopole_ode!(g4, rs, y, ω)
        @test g4 ≈ MONO rtol=rtol
        for l in (0, 1, 2, 5, 17)
            g5 = zeros(12); monopole_ode!(g5, rs, y, l, ω)
            @test g5 == g4                      # arity-5 method delegates exactly; l is dead
        end
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        setprecision(BigFloat, 256) do
            y = [sin(big"0.7" * k) + big"0.3" * cos(big"1.1" * k) for k in 1:8]
            g = zeros(BigFloat, 8)
            odd_ode!(g, big"12.3", y, 2, big"0.35")
            @test eltype(g) === BigFloat
            @test Float64.(g) ≈ ODD_L2 rtol=1e-12   # BigFloat path agrees with Float64 fixture
        end
    end
end
