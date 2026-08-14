# Grounding for the static (ω=0) TeukolskyRadial path's §9 BigFloat inheritance (design
# GRTOOLS_DESIGN_CAPABILITY §9 chain: "Perturbation (MST/CF/renorm-ν/static, all type-generic
# {T})"). Two gates:
#   (i)  Float64 REGRESSION — literals baked by running the PRE-{T} COMMITTED code (tree d817c05,
#        teukolsky_radial_static.jl with its ComplexF64 caps) at 3 (s,l,m,a,r) points chosen to
#        cover all static branches: τ≠0 (s=-2,l=2,m=2,a=0.9), τ=0 & s=0 (Up-amplitude special
#        case), τ=0 & s>0 (the Γ·Pochhammer norm + the c≤0-integer regularized-2F1 branch).
#        The {T} refactor must keep the Float64 path BIT-IDENTICAL (==), matching the
#        leg_radial_static baked-fixture discipline.
#   (ii) BigFloat INHERITANCE — BigFloat(256) inputs must give Complex{BigFloat} values/derivs/
#        norms, and the 256-bit vs 512-bit evaluations must self-converge far beyond Float64.
#        Measured 256v512 agreement (post-{T}, all 3 points × value/deriv × In/Up): worst
#        = 1.29e-76; tolerance set 1e-75 (~8x headroom over the worst measurement, still ~60
#        decades beyond Float64 — NOT tuned to a marginal pass). BigFloat-vs-Float64 baked
#        consistency: worst measured = 1.29e-15; rtol set 1e-14 (~8x).

# Baked 2026-07-02 from the committed pre-{T} static path (git tree d817c05); repr-exact Float64.
const REF_STATIC_PREC = [
    (s=-2, l=2, m=2, a=0.9, r=6.0,                                      # τ ≠ 0 branch
     vin=369.1925636560661 - 492.5256758223582im, vup=0.19343450885492436 + 0.047526813547457476im,
     din=333.3493546796276 - 370.2525958027509im, dup=-0.03606314012075073 - 0.01882060846846364im,
     normin=-0.028661351625646325 - 0.08256523812478475im, normup=1.147078669352809 + 0.0im),
    (s=0, l=2, m=0, a=0.5, r=8.0,                                       # τ = 0, s = 0 branch
     vin=97.50000000000003 + 0.0im, vup=0.0029541955596372473 + 0.0im,
     din=28.000000000000007 + 0.0im, dup=-0.001277296516313022 + 0.0im,
     normin=1.0 + 0.0im, normup=0.1924500897298753 + 0.0im),
    (s=2, l=3, m=0, a=0.3, r=5.0,                                       # τ = 0, s > 0 branch
     vin=4.193139346887674 + 0.0im, vup=0.00027981656970150717 + 0.0im,
     din=1.0482848367219184 + 0.0im, dup=-0.0004393419327558083 + 0.0im,
     normin=0.016666666666666666 + 0.0im, normup=0.020734608948592765 + 0.0im),
]

@testset "Radial: static {T} — Float64 regression bit-identical to pre-refactor bake" begin
    for p in REF_STATIC_PREC
        R = Perturbation.TeukolskyRadial(p.s, p.l, p.m, p.a, 0.0)
        @test R.In isa Perturbation.TeukolskyRadialFunction{Float64}
        @test R.In.norm == p.normin
        @test R.Up.norm == p.normup
        vin = R.In(p.r); vup = R.Up(p.r)
        @test vin isa ComplexF64 && vup isa ComplexF64
        @test vin == p.vin
        @test vup == p.vup
        @test Perturbation.teukolsky_radial_deriv(R.In, p.r) == p.din
        @test Perturbation.teukolsky_radial_deriv(R.Up, p.r) == p.dup
    end
end

@testset "Radial: static {T} — BigFloat(256) inheritance + 256v512 self-convergence" begin
    RTOL_CONV = 1e-75    # measured worst 256v512 relerr = 1.29e-76 (see header); NOT tuned
    RTOL_F64X = 1e-14    # measured worst BigFloat-vs-baked-Float64 relerr = 1.29e-15
    for p in REF_STATIC_PREC
        v = Dict{Int,NTuple{4,Complex{BigFloat}}}()
        for prec in (256, 512)
            setprecision(BigFloat, prec) do
                aB = BigFloat(p.a); rB = BigFloat(p.r)   # exact Float64 promotions: same inputs
                R = Perturbation.TeukolskyRadial(p.s, p.l, p.m, aB, zero(BigFloat))
                @test R.In isa Perturbation.TeukolskyRadialFunction{BigFloat}
                @test R.In.a isa BigFloat && R.In.omega isa BigFloat
                @test R.In.norm isa Complex{BigFloat} && R.Up.norm isa Complex{BigFloat}
                @test R.In.lambda isa Complex{BigFloat}
                @test all(x -> x isa Complex{BigFloat}, values(R.In.amplitudes))
                @test all(x -> x isa Complex{BigFloat}, values(R.Up.amplitudes))
                vin = R.In(rB); vup = R.Up(rB)
                din = Perturbation.teukolsky_radial_deriv(R.In, rB)
                dup = Perturbation.teukolsky_radial_deriv(R.Up, rB)
                @test vin isa Complex{BigFloat} && vup isa Complex{BigFloat}
                @test din isa Complex{BigFloat} && dup isa Complex{BigFloat}
                v[prec] = (vin, vup, din, dup)
            end
        end
        worst = 0.0
        for (v256, v512, f64) in zip(v[256], v[512], (p.vin, p.vup, p.din, p.dup))
            relerr = Float64(abs(v256 - v512) / abs(v512))
            worst = max(worst, relerr)
            @test relerr < RTOL_CONV                      # 256-bit self-converged vs 512-bit
            @test isapprox(ComplexF64(v512), f64; rtol=RTOL_F64X)  # consistency vs Float64 bake
        end
        @info "static {T} BigFloat self-convergence (s=$(p.s),l=$(p.l),m=$(p.m),a=$(p.a)): worst 256v512 relerr=$worst"
    end
end
