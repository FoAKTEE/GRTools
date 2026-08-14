# Module — SphereCollocation: GL-collocation swal + ð/ð′ ladder / Laplacian / filter / killtop
# matrices (TeukEvolution.jl/src/Sphere.jl port, Move D / D2).
#
# Grounding:
#  (A) ADAPTER identity — swal(s,m,l,y) at the ny=24 GL nodes equals the upstream Sphere.swal_vals
#      columns (baked fixture) ≤1e-13. Simultaneously validates the RUN-PINNED convention adapter
#      (D_dedup.md §3: √(2π)(−1)^(max(m,−s)−m)·Re ₛY) AND the type-generic GL grid.
#  (B) CHARACTERIZATION — each collocation matrix (vals/laplacian/raising/lowering/killtop/filter)
#      equals the baked upstream Sphere matrix ≤1e-12 (Float64), for the production (s,m).
#  (C) ANALYTIC ð-consistency — the raising/lowering matrices act on a spin-s swal column as the
#      exact ladder factor × the spin-(s±1) column; the Laplacian matrix is diagonal with the
#      −(l−s)(l+s+1) eigenvalue (the spectral-factor unit checks) ≤1e-10.
#  (D) {T} — every builder runs at BigFloat-256 with the eltype preserved and agrees with the
#      Float64 path ≤1e-13; orthonormality Vᵀ·diag(w)·V ≈ I on the retained num_l block.

using Test
using Harmonics
using LinearAlgebra
# Baked upstream oracle: TeukEvolution.jl/src/Sphere.jl matrices + GL grids at ny=24 (captured once
# via scratchpad/npghp_moveD/TeukEvolution_run/minrun/bake_harmonics_sphere.jl).
include(joinpath(@__DIR__, "fixtures", "sphere_collocation_ref.jl"))  # SPHERE_NY/YV/COS/SIN/W/MATS

_maxabsdiff(A, B) = maximum(abs.(A .- B))
# scale-relative error: abs for O(1) matrices, relative for the large-entry Laplacian
# (spectral factors up to (l−s)(l+s+1) ≈ 90 ⇒ Float64 roundoff ~1e-12 abs at rel ~5e-14).
_reldiff(A, B) = _maxabsdiff(A, B) / max(1.0, maximum(abs.(B)))

@testset "SphereCollocation" begin
    ny = SPHERE_NY                       # 24
    SM = collect(keys(SPHERE_MATS))      # (-2,2),(-2,-2),(0,2),(2,2)

    # ---- num_l truncation buffer -------------------------------------------------------------
    @testset "num_l" begin
        @test num_l(24) == 8                                   # 24 − 2·2 − 2·6
        @test num_l(24; max_s = 2, max_m = 6) == 8
        @test num_l(30) == 14
        @test_throws AssertionError num_l(20)                  # needs ny > 2(max_s+max_m)+4 = 20
    end

    # ---- (A) adapter identity: swal at fixture GL nodes == upstream swal_vals columns ---------
    @testset "adapter identity vs upstream swal_vals (≤1e-13)" begin
        worst = 0.0
        for (s, m) in SM
            vals = SPHERE_MATS[(s, m)].vals
            nl = size(vals, 2)
            lmin = max(abs(s), abs(m))
            for k in 1:nl
                l = lmin + k - 1
                for j in 1:ny
                    d = abs(swal(s, m, l, SPHERE_YV[j]) - vals[j, k])
                    worst = max(worst, d)
                end
            end
        end
        println("(A) adapter-vs-swal_vals max|Δ| = ", worst)
        @test worst < 1e-13
    end

    # independent human-readable anchors (D-R2 §3 full-precision spot values, Jacobi-poly oracle)
    @testset "adapter spot anchors vs Jacobi-poly oracle" begin
        @test swal(-2, 2, 2, -0.7) ≈ 1.142372804735827 rtol = 1e-13
        @test swal(0, 2, 2, -0.3) ≈ 0.8811037112621873 rtol = 1e-13
        @test swal(0, 2, 2, -0.7) ≈ 0.4938052667284015 rtol = 1e-6   # D-R2 printed 0.493805
        # n<0 guard: l below max(|m|,|s|) ⇒ 0
        @test swal(2, 0, 1, 0.3) == 0.0                              # l=1 < |s|=2
    end

    # ---- grids Y_vals / cos_vals / sin_vals + GL weights -------------------------------------
    @testset "GL grids vs fixture (≤1e-13)" begin
        @test _maxabsdiff(Y_vals(ny), SPHERE_YV) < 1e-13
        @test _maxabsdiff(cos_vals(ny), SPHERE_COS) < 1e-13
        @test _maxabsdiff(sin_vals(ny), SPHERE_SIN) < 1e-13
        # inner_product uses the same GL weights → cross-check against baked upstream weights
        _, wport = Harmonics.NumericsCore.gl_nodes_weights(Float64, ny)
        @test _maxabsdiff(wport, SPHERE_W) < 1e-13
    end

    # ---- (B) collocation matrices vs baked upstream fixture (Float64, scale-relative ≤1e-12) ---
    @testset "collocation matrices vs upstream fixture (≤1e-12)" begin
        for (s, m) in SM
            ref = SPHERE_MATS[(s, m)]
            dvals = _reldiff(swal_vals(ny, s, m), ref.vals)
            dlap  = _reldiff(swal_laplacian_matrix(ny, s, m), ref.lap)
            drai  = _reldiff(swal_raising_matrix(ny, s, m), ref.raise)
            dlow  = _reldiff(swal_lowering_matrix(ny, s, m), ref.lower)
            dkil  = _reldiff(swal_killtop_matrix(ny, s, m, SPHERE_KILLTOP_N), ref.killtop2)
            dfil  = _reldiff(swal_filter_matrix(ny, s, m), ref.filt)
            println("(B) (s=$s,m=$m) rel|Δ|: vals=$dvals lap=$dlap raise=$drai lower=$dlow killtop=$dkil filt=$dfil")
            @test dvals < 1e-12
            @test dlap < 1e-12
            @test drai < 1e-12
            @test dlow < 1e-12
            @test dkil < 1e-12
            @test dfil < 1e-12
        end
    end

    # ---- (C) analytic ð-consistency + spectral-factor unit checks (≤1e-10) --------------------
    @testset "spectral factors: Laplacian eigenvalue + ð/ð′ ladders" begin
        yv = Y_vals(ny)
        l0 = 4                          # in every operator's l-window for these (s,m); l0 ≥ max(|m|,|s±1|)
        for (s, m) in SM
            col = [swal(s, m, l0, y) for y in yv]                 # spin-s harmonic sampled

            # Laplacian: transpose(lap)·col == −(l0−s)(l0+s+1)·col
            lap = swal_laplacian_matrix(ny, s, m)
            eig = -(l0 - s) * (l0 + s + 1)
            @test _maxabsdiff(transpose(lap) * col, eig .* col) < 1e-10

            # raising ð: transpose(raise)·col == +√((l0−s)(l0+s+1))·swal(s+1,l0)
            raise = swal_raising_matrix(ny, s, m)
            up = [swal(s + 1, m, l0, y) for y in yv]
            fac_up = sqrt((l0 - s) * (l0 + s + 1))
            @test _maxabsdiff(transpose(raise) * col, fac_up .* up) < 1e-10

            # lowering ð′: transpose(lower)·col == −√((l0+s)(l0−s+1))·swal(s−1,l0)
            lower = swal_lowering_matrix(ny, s, m)
            dn = [swal(s - 1, m, l0, y) for y in yv]
            fac_dn = -sqrt((l0 + s) * (l0 - s + 1))
            @test _maxabsdiff(transpose(lower) * col, fac_dn .* dn) < 1e-10
        end
        # explicit factor literals for (s,m)=(-2,2), l0=4 (unit check of the formulas):
        @test -(4 - (-2)) * (4 + (-2) + 1) == -18                 # Laplacian eigenvalue
        @test sqrt((4 - (-2)) * (4 + (-2) + 1)) == sqrt(18.0)     # raise factor
        @test -sqrt((4 + (-2)) * (4 - (-2) + 1)) == -sqrt(14.0)   # lower factor
    end

    # ---- orthonormality: Vᵀ·diag(w)·V ≈ I on the retained num_l block -------------------------
    @testset "orthonormality Vᵀ diag(w) V ≈ I" begin
        _, wv = Harmonics.NumericsCore.gl_nodes_weights(Float64, ny)
        for (s, m) in SM
            V = swal_vals(ny, s, m)
            G = transpose(V) * Diagonal(wv) * V
            off = _maxabsdiff(G, Matrix{Float64}(I, size(G)...))
            println("orthonormality (s=$s,m=$m) max|G−I| = ", off)
            @test off < 1e-9
        end
        # inner_product API matches the same Gram entries (δ on the retained block)
        V = swal_vals(ny, -2, 2)
        @test inner_product(V[:, 1], V[:, 1]) ≈ 1.0 atol = 1e-9
        @test abs(inner_product(V[:, 1], V[:, 3])) < 1e-9
    end

    # ---- angular_matrix_mult! apply convention -----------------------------------------------
    @testset "angular_matrix_mult! (v[i]·M[i,j] → v[j])" begin
        s, m = -2, 2
        mat = swal_filter_matrix(ny, s, m)
        nx = 3
        f = reshape(collect(1.0:(nx*ny)), nx, ny)
        f_m = similar(f)
        angular_matrix_mult!(f_m, f, mat)
        # per row i: f_m[i,:] == transpose(mat) * f[i,:]
        for i in 1:nx
            @test _maxabsdiff(f_m[i, :], transpose(mat) * f[i, :]) < 1e-10
        end
        @test size(f_m) == (nx, ny)
    end

    # ---- (D) {T}: BigFloat-256 run + eltype preserved + Float64 agreement (≤1e-13) ------------
    @testset "BigFloat-256 genericity + Float64 agreement" begin
        setprecision(BigFloat, 256) do
            s, m = -2, 2
            N = SPHERE_KILLTOP_N

            # grids
            yb = Y_vals(BigFloat, ny)
            @test eltype(yb) == BigFloat
            @test _maxabsdiff(Float64.(yb), SPHERE_YV) < 1e-13
            @test eltype(cos_vals(BigFloat, ny)) == BigFloat
            @test _maxabsdiff(Float64.(sin_vals(BigFloat, ny)), SPHERE_SIN) < 1e-13

            # swal value at BigFloat node
            vb = swal(-2, 2, 2, big"-0.7")
            @test vb isa BigFloat
            @test Float64(vb) ≈ 1.142372804735827 rtol = 1e-13

            # every matrix builder at BigFloat, eltype preserved + Float64-path agreement
            for (name, fbig, ff64) in (
                ("vals", swal_vals(BigFloat, ny, s, m), swal_vals(ny, s, m)),
                ("lap", swal_laplacian_matrix(BigFloat, ny, s, m), swal_laplacian_matrix(ny, s, m)),
                ("raise", swal_raising_matrix(BigFloat, ny, s, m), swal_raising_matrix(ny, s, m)),
                ("lower", swal_lowering_matrix(BigFloat, ny, s, m), swal_lowering_matrix(ny, s, m)),
                ("killtop", swal_killtop_matrix(BigFloat, ny, s, m, N), swal_killtop_matrix(ny, s, m, N)),
                ("filt", swal_filter_matrix(BigFloat, ny, s, m), swal_filter_matrix(ny, s, m)),
            )
                @test eltype(fbig) == BigFloat
                d = _reldiff(Float64.(fbig), ff64)
                println("(D) BigFloat vs Float64 (s=$s,m=$m) $name rel|Δ| = ", d)
                @test d < 1e-13
            end

            # inner_product BigFloat eltype + orthonormality at high precision
            Vb = swal_vals(BigFloat, ny, s, m)
            ip = inner_product(Vb[:, 1], Vb[:, 1])
            @test ip isa BigFloat
            @test Float64(ip) ≈ 1.0 atol = 1e-13
        end
    end
end
