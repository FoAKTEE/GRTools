# Characterization + validation of the 1D nodal-DG spectral-element core (src/dg.jl).
#
# Grounding:
#  - LGL nodes/weights + D/M operators: oracle = the ORIGINAL DG-Teulolsky DGCore.RefElement
#    (FastGaussQuadrature.gausslobatto based), BAKED into fixtures/dg_lgl_ref.jl (captured once by
#    scratchpad/capture_dg_lgl_ref.jl in that repo's own env). NumericsCore's self-contained
#    Newton-on-P'_N rule must reproduce it to <=1e-14 rel for N=2..16 (same mathematical objects).
#  - "DGCore reference element" / "Lagrange interpolation" / "Mesh with forced interface"
#    testsets: ported essentially verbatim from DG-Teulolsky/test/runtests.jl (same tolerances).
#  - Advection superconvergence: DG-Teulolsky/src/Advection.jl ported below as TEST-SUPPORT code
#    (not runtime src) — the Cao-Yang model problem; DG-upwind gives the textbook N+1 global-L2
#    rate and 2N+1 downwind-node superconvergence (same rate tolerances as the original repo).
#  - BigFloat LGL exactness: NEW capability (the original was Float64-only via gausslobatto) —
#    at BigFloat(256) the N+1-point rule integrates degree 2N-1 exactly to 1e-70-class residual.

using NumericsCore
using NumericsCore: _legendre_PdP, _lgl_nodes_weights, rk4_step!
using LinearAlgebra: mul!, norm

# Baked from the ORIGINAL DG-Teulolsky DGCore (captured once) — no live FGQ/DG-Teulolsky dep.
include("fixtures/dg_lgl_ref.jl")

# ---------------------------------------------------------------------------------------------
# Test-support: DG-Teulolsky/src/Advection.jl ported verbatim (module wrapper dropped; DGCore/
# TimeStepping names now live flat in NumericsCore). Linear advection u_t + c u_x = 0, periodic,
# upwind flux — the canonical Cao-Yang superconvergence model validating the DG core machinery.
# ---------------------------------------------------------------------------------------------
function advection_errors(N::Int, K::Int; c=1.0, tfinal=1.0, dt=2e-4, u0=x -> sinpi(2x))
    el = RefElement(N)
    Np = el.Np
    h = 1.0 / K
    x = [(k - 1) * h + (el.r[i] + 1) / 2 * h for i in 1:Np, k in 1:K]
    invMc1 = el.invM[:, 1]
    invMcN = el.invM[:, el.Np]
    Du = Vector{Float64}(undef, Np)
    rhs! = function (du, u, t)
        U = reshape(u, Np, K)
        dU = reshape(du, Np, K)
        @inbounds for k in 1:K
            mul!(Du, el.D, view(U, :, k))
            km1 = k == 1 ? K : k - 1
            kp1 = k == K ? 1 : k + 1
            if c >= 0                       # upwind from the left
                jL = U[1, k] - U[Np, km1]   # (F-F*) on left face
                for i in 1:Np
                    dU[i, k] = -c * (2 / h) * Du[i] - c * (2 / h) * invMc1[i] * jL
                end
            else                            # upwind from the right
                jR = U[Np, k] - U[1, kp1]
                for i in 1:Np
                    dU[i, k] = -c * (2 / h) * Du[i] + c * (2 / h) * invMcN[i] * jR
                end
            end
        end
        return du
    end
    uv = vec([u0(xi) for xi in x])
    rk4!(rhs!, uv, (0.0, tfinal), dt)
    U = reshape(uv, Np, K)
    ue = [u0(xi) for xi in x]              # exact at t=period
    num = 0.0
    den = 0.0
    @inbounds for k in 1:K, i in 1:Np
        w = el.w[i] * h / 2
        num += w * (U[i, k] - ue[i, k])^2
        den += w * ue[i, k]^2
    end
    km = cld(K, 2)
    downwind = abs(U[Np, km] - ue[Np, km])
    return sqrt(num / den), downwind
end

function advection_rates(N::Int, Ks; kw...)
    l2 = Float64[]
    dn = Float64[]
    for K in Ks
        e = advection_errors(N, K; kw...)
        push!(l2, e[1])
        push!(dn, e[2])
    end
    fit(es) = begin
        x = Float64[]
        y = Float64[]
        for (K, e) in zip(Ks, es)
            (1e-13 < e < 0.5) && (push!(x, log(float(K))); push!(y, log(e)))
        end
        n = length(x)
        n < 2 && return NaN
        sx = sum(x); sy = sum(y)
        -(n * sum(x .* y) - sx * sy) / (n * sum(abs2, x) - sx^2)
    end
    return fit(l2), fit(dn), l2, dn
end

@testset "dg (DG spectral-element core)" begin

    @testset "LGL characterization vs baked FGQ original (N=2..16)" begin
        # ours vs the captured FastGaussQuadrature-based DGCore.RefElement — same mathematical
        # objects, required to agree to <=1e-14 (rel, with abs floor 1 on the O(1) nodes).
        devmax_rw = 0.0
        devmax_DM = 0.0
        for (N, ref) in sort(collect(DG_LGL_REF); by=first)
            r, w = _lgl_nodes_weights(Float64, N)
            @test length(r) == N + 1
            dr = maximum(abs.(r .- ref.r) ./ max.(abs.(ref.r), 1.0))
            dw = maximum(abs.(w .- ref.w) ./ abs.(ref.w))
            @test dr <= 1e-14
            @test dw <= 1e-14
            devmax_rw = max(devmax_rw, dr, dw)
            # whole-matrix Frobenius hashes of D and M for every N
            el = RefElement(N)
            fD, fM = DG_DM_FRO[N]
            @test abs(norm(el.D) - fD) / fD <= 1e-14
            @test abs(norm(el.M) - fM) / fM <= 1e-14
        end
        # full D and M matrices on the representative subset
        for (N, ref) in sort(collect(DG_DM_REF); by=first)
            el = RefElement(N)
            dD = norm(el.D - ref.D) / norm(ref.D)
            dM = norm(el.M - ref.M) / norm(ref.M)
            @test dD <= 1e-14
            @test dM <= 1e-14
            devmax_DM = max(devmax_DM, dD, dM)
        end
        println("dg LGL vs FGQ-original agreement: max node/weight rel dev = ",
                round(devmax_rw, sigdigits=3), ", max D/M Frobenius rel dev = ",
                round(devmax_DM, sigdigits=3))
    end

    @testset "LGL BigFloat(256) exactness (NEW capability)" begin
        setprecision(BigFloat, 256) do
            for N in (4, 8, 12, 16)
                r, w = _lgl_nodes_weights(BigFloat, N)
                @test eltype(r) === BigFloat && precision(r[2]) == 256
                # Sum(w) = ∫1 dx = 2 to full precision
                @test abs(sum(w) - 2) < big"1e-70"
                # the N+1-point LGL rule is exact through degree 2N-1:
                #   ∫_{-1}^{1} x^(2N-1) dx = 0   and   ∫_{-1}^{1} x^(2N-2) dx = 2/(2N-1)
                @test abs(sum(w .* r .^ (2N - 1))) < big"1e-70"
                @test abs(sum(w .* r .^ (2N - 2)) - 2 / big(2N - 1)) < big"1e-70"
                # ... and genuinely NOT exact at degree 2N (we sit at the exactness boundary)
                @test abs(sum(w .* r .^ (2N)) - 2 / big(2N + 1)) > big"1e-30"
            end
            # a full BigFloat reference element: spectral operators inherit the precision
            el = RefElement(6; T=BigFloat)
            @test eltype(el.D) === BigFloat
            f = el.r .^ 5
            df = 5 .* el.r .^ 4
            @test maximum(abs.(el.D * f .- df)) < big"1e-70"
            @test maximum(abs.(el.M * ones(BigFloat, 7) .- el.w)) < big"1e-70"
        end
    end

    # --- DG-Teulolsky/test/runtests.jl "DGCore reference element" — ported verbatim ---
    @testset "DGCore reference element" begin
        for N in 2:10
            el = RefElement(N)
            @test el.Np == N + 1
            @test el.r[1] ≈ -1 atol = 1e-13
            @test el.r[end] ≈ 1 atol = 1e-13
            @test sum(el.w) ≈ 2 atol = 1e-12
            # exact differentiation of polynomials up to degree N
            for p in 0:N
                f = el.r .^ p
                df = p == 0 ? zeros(N + 1) : p .* el.r .^ (p - 1)
                @test maximum(abs.(el.D * f .- df)) < 1e-9
            end
            # mass-matrix row sums equal LGL weights (∫ℓ_i dx = w_i)
            @test maximum(abs.(el.M * ones(N + 1) .- el.w)) < 1e-11
            # symmetry of mass matrix
            @test maximum(abs.(el.M - el.M')) < 1e-12
        end
    end

    # --- DG-Teulolsky/test/runtests.jl "Lagrange interpolation" — ported verbatim ---
    @testset "Lagrange interpolation" begin
        el = RefElement(6)
        # interpolating a degree-≤N polynomial is exact
        coef = [0.3, -1.1, 0.7, 0.2, -0.05]      # degree 4 ≤ 6
        poly(x) = sum(coef[k] * x^(k - 1) for k in eachindex(coef))
        u = poly.(el.r)
        for ξ in -0.9:0.37:0.9
            @test lagrange_interp(el, u, ξ) ≈ poly(ξ) atol = 1e-11
        end
        # interpolation at a node returns the node value
        @test lagrange_interp(el, u, el.r[3]) ≈ u[3] atol = 1e-12
    end

    # --- DG-Teulolsky/test/runtests.jl "Mesh with forced interface" — ported verbatim ---
    @testset "Mesh with forced interface" begin
        m = build_mesh(1.0, 50.0, 16; forced=[30.0])
        @test 30.0 in m.edges
        @test m.edges[1] ≈ 1.0 && m.edges[end] ≈ 50.0
        @test length(m.h) == 16
        # node_coords spans the mesh and respects element edges
        el = RefElement(4)
        ρ = node_coords(m, el)
        @test size(ρ) == (16, 5)
        @test ρ[1, 1] ≈ 1.0 && ρ[end, end] ≈ 50.0
        @test all(ρ[k, end] ≈ m.edges[k+1] for k in 1:m.K)
        # element_of locates interior points and interfaces
        @test element_of(m, 1.0) == 1
        @test element_of(m, 49.999) == 16
        @test element_of(m, prevfloat(30.0)) == element_of(m, 30.0)   # forced edge belongs left

        # Endpoints are immutable even when the nearest uniform edge would have been an endpoint.
        edge_forced = build_mesh(0.0, 1.0, 4; forced=[0.01, 0.99])
        @test edge_forced.edges[1] == 0.0
        @test edge_forced.edges[end] == 1.0
        @test all(f -> f in edge_forced.edges, (0.01, 0.99))
        @test all(>(0), edge_forced.h)

        # Nearby requests receive distinct slots; impossible requests fail explicitly.
        injective = build_mesh(0.0, 1.0, 4; forced=[0.24, 0.26])
        @test 0.24 in injective.edges && 0.26 in injective.edges
        @test length(injective.edges) == 5 && allunique(injective.edges)
        @test_throws ArgumentError build_mesh(0.0, 1.0, 2; forced=[0.2, 0.8])
        @test_throws ArgumentError build_mesh(0.0, 1.0, 1; forced=[0.5])
        @test build_mesh(0.0, 1.0, 1; forced=[0.0, 1.0]).edges == [0.0, 1.0]
    end

    @testset "rk4! and cfl_dt" begin
        # u' = -u, u(0)=1 → exp(-1): classical 4th-order convergence of the ported stepper
        rhs!(du, u, t) = (du .= .-u; du)
        err(dt) = begin
            u = [1.0]
            _, tend = rk4!(rhs!, u, (0.0, 1.0), dt)
            @test tend ≈ 1.0 atol = 1e-12     # lands exactly on the end time
            abs(u[1] - exp(-1.0))
        end
        e1, e2 = err(0.1), err(0.05)
        @test e1 / e2 > 12                    # ≈ 2^4 = 16 for order 4 (measured 16.68)
        @test e2 < 1e-7                       # measured 1.998e-8
        # callback sees every step including step 0
        nsteps = Ref(0)
        rk4!(rhs!, [1.0], (0.0, 1.0), 0.25; callback=(u, t, s) -> (nsteps[] = s))
        @test nsteps[] == 4
        # rk4_step! advances a single step in place (matches one rk4! step)
        u = [1.0]; work = ntuple(_ -> similar(u), 5)
        rk4_step!(rhs!, u, 0.0, 0.1, work)
        @test u[1] ≈ exp(-0.1) atol = 1e-7

        # A positive span smaller than the old absolute 1e-14 guard still takes one shortened step.
        unit_rhs!(du, u, t) = (du .= 1; du)
        tiny = [0.0]
        _, tiny_end = rk4!(unit_rhs!, tiny, (0.0, 1e-16), 1.0)
        @test tiny_end === 1e-16
        @test tiny[1] ≈ 1e-16 rtol = 8eps()

        # Non-divisible spans land exactly on the endpoint with the expected shortened last step.
        times = Float64[]
        final = [0.0]
        _, final_t = rk4!(unit_rhs!, final, (0.0, 0.25), 0.1;
                          callback=(u, t, s) -> push!(times, t))
        @test final_t === 0.25
        @test times == [0.0, 0.1, 0.2, 0.25]
        @test final[1] ≈ 0.25 atol = 2eps()

        for bad_dt in (0.0, -0.1, Inf, NaN)
            @test_throws ArgumentError rk4!(rhs!, [1.0], (0.0, 1.0), bad_dt)
        end
        @test_throws ArgumentError rk4!(rhs!, [1.0], (1.0, 0.0), 0.1)
        # CFL helper is the plain formula
        @test cfl_dt(0.5, 0.1) == 0.05
        @test cfl_dt(0.5, 0.1, 2.0) == 0.025
    end

    # --- DG-Teulolsky/test/runtests.jl "Advection superconvergence" — ported verbatim ---
    @testset "Advection superconvergence (DG core validation)" begin
        # Textbook Cao–Yang model: DG-upwind gives N+1 in L2 and 2N+1 at downwind nodes.
        for N in (2, 3)
            l2r, dnr, _, _ = advection_rates(N, (8, 16, 32); dt=2e-4)
            @test isapprox(l2r, N + 1; atol=0.35)        # global L2 rate ≈ N+1
            @test dnr > l2r + 1.0                        # downwind node superconverges
        end
    end
end
