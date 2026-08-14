# Gates for the uniform-grid FD sweep (fd_uniform.jl) + type-generic Gauss-Legendre rule
# (gauss_legendre.jl) — Move D / D1.
#
# Grounding:
#  - set_d1!/set_d2!/ko_filter!: port of the ORIGINAL TeukEvolution.jl/src/Radial.jl:27-162.
#    INTERIOR stencils are coefficient-identical to fd1c/fd2c (interp.jl:569,572) — asserted below
#    (DEDUP). Whole-function oracle: FD_UNIFORM_REF, baked from the upstream Radial run into
#    fixtures/numericscore_fd_uniform_ref.jl (bake script minrun/bake_numericscore_fd.jl).
#  - gl_nodes_weights: mirrors the self-contained _lgl_nodes_weights (dg.jl); oracle = GL24_REF,
#    baked from FastGaussQuadrature.gausslegendre(24). BigFloat exactness is a NEW capability
#    (upstream gausslegendre was Float64-only).
#
# Exact polynomial degrees (STATED and tested per stencil row, incl. the negative "one higher fails"):
#   set_d1! interior + every edge row : exact ≤ 4   set_d2! interior + every edge row : exact ≤ 5
#   KO dissipation term               : annihilates ≤ 5 interior

using Test
using NumericsCore
using NumericsCore: set_d1!, set_d2!, ko_filter!, KO6_KERNEL, gl_nodes_weights, fd1c, fd2c,
    _legendre_PdP

include("fixtures/numericscore_fd_uniform_ref.jl")   # FD_UNIFORM_REF, GL24_REF (no live FGQ dep)

# --- helpers ------------------------------------------------------------------------------------
# extract the 5 interior weights of set_d1!/set_d2! by impulse (dr=1, read an interior node)
function interior_weights(setd!, offsets)
    nx, i0, dr = 12, 6, 1.0
    ws = Float64[]
    for off in offsets
        v = zeros(nx, 1); v[i0+off, 1] = 1.0
        dv = similar(v); setd!(dv, v, dr)
        push!(ws, dv[i0, 1])
    end
    return ws
end
# fd1c/fd2c weights over the same offsets, extracted by indicator probes (r=0, h=1)
fdc_weights(fdc, offsets) = [fdc(x -> (x == Float64(off) ? 1.0 : 0.0), 0.0, 1.0) for off in offsets]

# per-row polynomial errors for a monomial x^deg on a centered uniform grid
function poly_errs(setd!, deg, order; nx=14, dr=0.1)
    x = [(i - (nx + 1) / 2) * dr for i in 1:nx]
    v = reshape(x .^ deg, nx, 1)
    dv = similar(v); setd!(dv, v, dr)
    if order == 1
        tru = deg == 0 ? zeros(nx) : deg .* x .^ (deg - 1)
    else
        tru = deg <= 1 ? zeros(nx) : deg * (deg - 1) .* x .^ (deg - 2)
    end
    intr = maximum(abs.(dv[3:nx-2, 1] .- tru[3:nx-2]))
    edges = [abs(dv[r, 1] - tru[r]) for r in (1, 2, nx - 1, nx)]
    return intr, edges
end

exact_int(x) = iseven(x) ? 2 / (x + 1) : 0.0        # ∫_{-1}^{1} t^x dt

@testset "fd_uniform + Gauss-Legendre (Move D / D1)" begin

    # ---- GATE 1: interior stencil == fd1c/fd2c (DEDUP coefficient identity) ---------------------
    @testset "interior == fd1c/fd2c coefficient identity" begin
        offs = -2:2
        w_d1 = interior_weights(set_d1!, offs)
        w_d2 = interior_weights(set_d2!, offs)
        @test w_d1 == fdc_weights(fd1c, offs)              # exact Float64 identity
        @test w_d2 == fdc_weights(fd2c, offs)
        # explicit expected rationals (documentation of the shared coefficients)
        @test w_d1 == Float64.([1//12, -2//3, 0, 2//3, -1//12])
        @test w_d2 == Float64.([-1//12, 4//3, -5//2, 4//3, -1//12])
        println("DEDUP interior weights  d1=$(w_d1)  d2=$(w_d2)  (== fd1c/fd2c: ",
                w_d1 == fdc_weights(fd1c, offs) && w_d2 == fdc_weights(fd2c, offs), ")")
    end

    # ---- GATE 2: polynomial exactness — exact degree per stencil incl. every edge row ----------
    @testset "polynomial exactness (exact degrees stated)" begin
        pass, fail = 1e-9, 1e-6
        # set_d1!: exact for degree ≤ 4, fails at 5 (interior AND each of the 4 edge rows)
        for d in 0:4
            intr, edges = poly_errs(set_d1!, d, 1)
            @test intr ≤ pass
            @test all(e ≤ pass for e in edges)
        end
        intr5, edges5 = poly_errs(set_d1!, 5, 1)
        @test intr5 > fail                                # interior NOT exact at degree 5
        @test all(e > fail for e in edges5)               # every edge row NOT exact at degree 5
        println("set_d1! exact ≤4: deg5 interior err=$(round(intr5,sigdigits=3)), ",
                "edge errs=$(round.(edges5,sigdigits=3))")
        # set_d2!: exact for degree ≤ 5, fails at 6 (interior AND each edge row)
        for d in 0:5
            intr, edges = poly_errs(set_d2!, d, 2)
            @test intr ≤ pass
            @test all(e ≤ pass for e in edges)
        end
        intr6, edges6 = poly_errs(set_d2!, 6, 2)
        @test intr6 > fail
        @test all(e > fail for e in edges6)
        println("set_d2! exact ≤5: deg6 interior err=$(round(intr6,sigdigits=3)), ",
                "edge errs=$(round.(edges6,sigdigits=3))")
    end

    # ---- GATE 3: Kreiss-Oliger kernel sum-zero + annihilates degree ≤ 5 interior ---------------
    @testset "KO sum-zero + annihilates degree ≤5" begin
        @test sum(KO6_KERNEL) == 0
        @test collect(KO6_KERNEL) == [1, -6, 15, -20, 15, -6, 1]
        nx = 16
        idx = Float64.(1:nx)
        for d in 0:5                                       # degree ≤5 → dissipation term vanishes
            v = reshape(idx .^ d, nx, 1); v0 = copy(v)
            ko_filter!(v, similar(v), 1.0)
            @test maximum(abs.(v[4:nx-3, 1] .- v0[4:nx-3, 1])) ≤ 1e-6 * maximum(abs.(v0[4:nx-3, 1]) .+ 1)
        end
        v6 = reshape(idx .^ 6, nx, 1); v60 = copy(v6)       # degree 6 → genuinely modified
        ko_filter!(v6, similar(v6), 1.0)
        @test maximum(abs.(v6[4:nx-3, 1] .- v60[4:nx-3, 1])) > 1.0
        # sum-zero manifestation: constant field left unchanged interior (already the d=0 case)
        println("KO6 kernel sum = ", sum(KO6_KERNEL), "; degree-6 interior |Δ|max > 1 confirmed")
    end

    # ---- GATE 4: upstream Radial char fixture reproduced (Float64 == / ≤1e-15) ------------------
    @testset "upstream Radial char fixture (== / ≤1e-15)" begin
        R = FD_UNIFORM_REF
        v = copy(R.v)
        dv1 = similar(v); set_d1!(dv1, v, R.dr)
        dv2 = similar(v); set_d2!(dv2, v, R.dr)
        vf = copy(v); ko_filter!(vf, similar(vf), R.eps_KO)
        e1 = maximum(abs.(dv1 .- R.d1)); e2 = maximum(abs.(dv2 .- R.d2)); ef = maximum(abs.(vf .- R.filt))
        @test dv1 == R.d1                                  # bit-exact (same rational coefficients)
        @test dv2 == R.d2
        @test vf == R.filt
        @test e1 ≤ 1e-15 && e2 ≤ 1e-15 && ef ≤ 1e-15
        println("upstream char max|Δ|: set_d1!=$(e1)  set_d2!=$(e2)  ko_filter!=$(ef)")
    end

    # ---- GATE 5: Gauss-Legendre vs FastGaussQuadrature fixture (≤1e-15) -------------------------
    @testset "GL(24) vs FastGaussQuadrature fixture (≤1e-15)" begin
        x, w = gl_nodes_weights(Float64, 24)
        en = maximum(abs.(x .- GL24_REF.nodes)); ew = maximum(abs.(w .- GL24_REF.weights))
        @test en ≤ 1e-15
        @test ew ≤ 1e-15
        println("GL(24) vs FGQ: max|Δnodes|=$(en)  max|Δweights|=$(ew)")
    end

    # ---- GATE 6: GL Newton residual ≤ 4 eps(T) --------------------------------------------------
    @testset "GL Newton residual ≤ 4 eps(T)" begin
        for n in (5, 12, 24)
            x, _ = gl_nodes_weights(Float64, n)
            res = maximum(abs(begin P, dP = _legendre_PdP(n, xi); P / dP end) for xi in x)
            @test res ≤ 4 * eps(Float64)
        end
        x24, _ = gl_nodes_weights(Float64, 24)
        r24 = maximum(abs(begin P, dP = _legendre_PdP(24, xi); P / dP end) for xi in x24)
        println("GL Newton residual max |P_n/P'_n| (n=24) = ", r24, "  (4eps=$(4eps(Float64)))")
    end

    # ---- GATE 7: quadrature exactness on degree 2n-1 (and NOT 2n) -------------------------------
    @testset "GL quadrature exact to degree 2n-1" begin
        for n in (5, 8)
            x, w = gl_nodes_weights(Float64, n)
            for d in 0:(2n - 1)                            # exact ≤ 2n-1
                q = sum(w[i] * x[i]^d for i in 1:n)
                @test isapprox(q, exact_int(d); atol=1e-12)
            end
            q2n = sum(w[i] * x[i]^(2n) for i in 1:n)        # degree 2n NOT exact
            @test abs(q2n - exact_int(2n)) > 1e-6
        end
        x5, w5 = gl_nodes_weights(Float64, 5)
        println("GL(5): ∫x^9 (exact 0) = ", sum(w5[i]*x5[i]^9 for i in 1:5),
                "; ∫x^10 err (deg 2n) = ", abs(sum(w5[i]*x5[i]^10 for i in 1:5) - exact_int(10)))
    end

    # ---- GATE 8: BigFloat-256 self-consistency + eltype preservation ----------------------------
    @testset "BigFloat GL 256 vs 512 (≤1e-70) + FD/KO eltype-preserving" begin
        # GL nodes/weights self-consistent across 256 vs 512 bits
        n = 12
        x256, w256 = setprecision(BigFloat, 256) do; gl_nodes_weights(BigFloat, n); end
        x512, w512 = setprecision(BigFloat, 512) do; gl_nodes_weights(BigFloat, n); end
        dnodes = maximum(abs.(x256 .- x512)); dweights = maximum(abs.(w256 .- w512))
        @test dnodes ≤ 1e-70
        @test dweights ≤ 1e-70
        @test eltype(x256) == BigFloat
        # BigFloat Newton residual to ~eps(BigFloat 256) and degree 2n-1 exactness at 256 bits
        setprecision(BigFloat, 256) do
            xb, wb = gl_nodes_weights(BigFloat, n)
            resb = maximum(abs(begin P, dP = _legendre_PdP(n, xi); P / dP end) for xi in xb)
            @test resb ≤ 4 * eps(BigFloat)
            q = sum(wb[i] * xb[i]^(2n - 1) for i in 1:n)   # odd degree → exact 0
            @test abs(q) ≤ 1e-70
        end
        println("GL BigFloat 256v512: max|Δnodes|=$(Float64(dnodes))  max|Δweights|=$(Float64(dweights))")

        # FD/KO eltype-preserving on Complex{BigFloat}; genuine-precision path matches Float64 oracle
        setprecision(BigFloat, 256) do
            R = FD_UNIFORM_REF
            vC = Complex{BigFloat}.(R.v)                    # real part = fixture field, imag 0
            drB = BigFloat(R.dr)
            dvC = similar(vC); set_d1!(dvC, vC, drB)
            dv2C = similar(vC); set_d2!(dv2C, vC, drB)
            vfC = copy(vC); ko_filter!(vfC, similar(vfC), BigFloat(R.eps_KO))
            @test eltype(dvC) == Complex{BigFloat}
            @test eltype(dv2C) == Complex{BigFloat}
            @test eltype(vfC) == Complex{BigFloat}
            # real part reproduces the Float64 Radial oracle (rel ≤1e-14) — genuine compute, not stub
            rel1 = maximum(abs.(Float64.(real.(dvC)) .- R.d1)) / maximum(abs.(R.d1))
            relf = maximum(abs.(Float64.(real.(vfC)) .- R.filt)) / maximum(abs.(R.filt))
            @test rel1 ≤ 1e-14
            @test relf ≤ 1e-14
            println("Complex{BigFloat} eltype preserved; real-part vs Float64 oracle rel: ",
                    "set_d1!=$(rel1)  ko_filter!=$(relf)")
        end
    end
end
