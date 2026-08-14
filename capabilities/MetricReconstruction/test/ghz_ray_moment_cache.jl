using Test
using MetricReconstruction
using MetricReconstruction.GHZ

const GHZ_RMC = MetricReconstruction.GHZ
if !isdefined(GHZ_RMC, :RayMomentCache)
    Base.include(GHZ_RMC, joinpath(@__DIR__, "..", "src", "GHZ", "RayMomentCache.jl"))
end

const RMC_PRIMARY_FP = UInt128(0x2c98fcc9f6a6478aa937bb1783ab7c1a)
const RMC_PARTNER_FP = UInt128(0xcabf3cfe89d84914afe86b22e44a17ea)
const RMC_PAIR_FP = UInt128(0x77165819c780423f8aa11eaa89ffb642)
const RMC_BOUNDARY_FP = UInt128(0xe548388c6551442a81167d0586908f23)

function rmc_middle(::Val{N}, ::Type{T}; l=2, m=1,
        omega=Complex{T}(T(0.2), T(-0.1)), source_fingerprint,
        Sll=q -> zero(q), Slm=q -> zero(q), Sln=q -> zero(q),
        rb=T(2), head_x0=zero(Complex{T}), head_dx0=zero(Complex{T}),
        middle_x0=zero(Complex{T}), middle_dx0=zero(Complex{T})) where
        {N,T<:AbstractFloat}
    src = TetradSourceMode{N}(Sll; Slm, Sln, M=one(T), l, m, omega,
        source_fingerprint)
    head = solve_irg_corrector_head(src,
        InnerRayIVP(src; rb, x_mmb=head_x0, dr_x_mmb=head_dx0))
    solve_irg_corrector_middle(head,
        MiddleRayIVP(head; rb, x_nm=middle_x0, dr_x_nm=middle_dx0))
end

function rmc_mode(primary, partner; x0=0, pair_fingerprint=RMC_PAIR_FP,
        boundary_fingerprint=RMC_BOUNDARY_FP)
    pair = PhysicalCorrectorMiddlePair(primary, partner; pair_fingerprint)
    boundary = TailRayIVP(pair; x_nn=x0, boundary_fingerprint,
        policy=FiniteRayFreeIVP())
    solve_irg_corrector_tail(pair, boundary)
end

function rmc_zero_mode(::Type{T}=Float64; l=2, m=1,
        omega=Complex{T}(T(0.2), T(-0.1)), x0=zero(Complex{T}),
        primary_kwargs=NamedTuple(), partner_kwargs=NamedTuple()) where {T<:AbstractFloat}
    primary = rmc_middle(Val(2), T; l, m, omega,
        source_fingerprint=RMC_PRIMARY_FP, primary_kwargs...)
    partner = rmc_middle(Val(2), T; l, m=-m, omega=-conj(omega),
        source_fingerprint=RMC_PARTNER_FP, partner_kwargs...)
    rmc_mode(primary, partner; x0)
end

function rmc_manufactured(::Type{T}=Float64) where {T<:AbstractFloat}
    root3 = sqrt(T(3))
    omega = Complex{T}(T(0.2), zero(T))
    primary = rmc_middle(Val(2), T; l=2, m=2, omega,
        Sll=r -> -T(6) + T(8) / r,
        Slm=r -> root3 * (r - T(2)) / r - T(2) / (root3 * r^2),
        Sln=r -> one(r), middle_dx0=inv(root3),
        source_fingerprint=RMC_PRIMARY_FP + UInt128(10))
    partner = rmc_middle(Val(2), T; l=2, m=-2, omega=-conj(omega),
        Slm=r -> -T(4) / (root3 * r^2), middle_dx0=T(2) / root3,
        source_fingerprint=RMC_PARTNER_FP + UInt128(10))
    rmc_mode(primary, partner; x0=Complex{T}(T(1.5), T(-0.25)),
        pair_fingerprint=RMC_PAIR_FP + UInt128(10),
        boundary_fingerprint=RMC_BOUNDARY_FP + UInt128(10))
end

function rmc_compare_scalar(cache, mode, grid; rtol, atol)
    for r in grid
        @test GHZ_RMC.x_mmb(cache, r) ≈ x_mmb(mode.pair.primary.head, r; rtol, atol)
        @test GHZ_RMC.dr_x_mmb(cache, r) ≈
            dr_x_mmb(mode.pair.primary.head, r; rtol, atol)
        @test GHZ_RMC.x_nm(cache, r) ≈ x_nm(mode.pair.primary, r; rtol, atol)
        @test GHZ_RMC.dr_x_nm(cache, r) ≈ dr_x_nm(mode.pair.primary, r; rtol, atol)
        @test GHZ_RMC.x_nmb(cache, r) ≈ x_nmb(mode.pair, r; rtol, atol)
        @test GHZ_RMC.dr_x_nmb(cache, r) ≈ dr_x_nmb(mode.pair, r; rtol, atol)
        @test GHZ_RMC.x_nn(cache, r) ≈ x_nn(mode, r; rtol, atol)
    end
end

function rmc_raw_boundary(bc::TailRayIVP{N,T};
        primary_source_fingerprint=bc.primary_source_fingerprint,
        boundary_fingerprint=bc.boundary_fingerprint) where {N,T}
    TailRayIVP{N,T}(bc.M, bc.l, bc.m, bc.omega, primary_source_fingerprint,
        bc.partner_source_fingerprint, bc.pair_fingerprint, bc.primary_head_x0,
        bc.primary_head_dx0, bc.primary_middle_x0, bc.primary_middle_dx0,
        bc.partner_head_x0, bc.partner_head_dx0, bc.partner_middle_x0,
        bc.partner_middle_dx0, bc.rb, bc.x0, boundary_fingerprint, bc.policy)
end

function rmc_raw_cache(cache::GHZ_RMC.RayMomentCache{N,T};
        provenance=cache.provenance, radii=copy(cache.radii), A=copy(cache.A),
        dA=copy(cache.dA), B=copy(cache.B), dB=copy(cache.dB), C=copy(cache.C),
        dC=copy(cache.dC), Nvalues=copy(cache.N),
        diagnostics=GHZ_RMC.RayMomentQuadratureDiagnostics{T}(
            copy(cache.quadrature_diagnostics.primary_sll),
            copy(cache.quadrature_diagnostics.primary_slm),
            copy(cache.quadrature_diagnostics.partner_sll),
            copy(cache.quadrature_diagnostics.partner_slm),
            copy(cache.quadrature_diagnostics.primary_sln)), atol=cache.atol,
        rtol=cache.rtol) where {N,T}
    GHZ_RMC.RayMomentCache{N,T}(provenance, radii, A, dA, B, dB, C, dC,
        Nvalues, diagnostics, atol, rtol, Val(:validated))
end

@testset "GHZ — forward ray-moment cache" begin
    @testset "full-chain manufactured oracle and exact lookup" begin
        mode = rmc_manufactured()
        grid = [2.0, 2.25, 3.0, 3.75, 5.0]
        cache = GHZ_RMC.materialize_grid(mode, grid; rtol=2e-12)
        @test cache isa GHZ_RMC.RayMomentCache{2,Float64}
        @test perturbative_order(cache) == 2
        @test cache.radii == grid
        rmc_compare_scalar(cache, mode, grid; rtol=3e-10, atol=2e-12)

        index = findfirst(==(3.0), grid)
        @test cache.A[index] ≈ 1.0 atol=4e-13
        @test cache.B[index] ≈ inv(sqrt(3.0)) atol=4e-13
        @test cache.C[index] ≈ -2 / sqrt(3.0) atol=5e-13
        @test cache.N[index] ≈ 55 / 9 + 13im / 30 rtol=2e-10

        components = GHZ_RMC.ray_moment_components(cache, 3.0)
        @test components isa GHZ_RMC.RayMomentComponents{2,Float64}
        @test perturbative_order(components) == 2
        @test GHZ_RMC.corrector_tetrad_components(cache, 3.0) == components
        @test (components.A, components.dA, components.B, components.dB,
            components.C, components.dC, components.N) ==
            (cache.A[index], cache.dA[index], cache.B[index], cache.dB[index],
             cache.C[index], cache.dC[index], cache.N[index])
        @test components.grid_index == index
        @test components.r == 3.0
        @test components.quadrature_diagnostic ==
            GHZ_RMC.ray_quadrature_diagnostics(cache, 3.0)
        @test (components.atol, components.rtol) == (cache.atol, cache.rtol)

        provenance = components.provenance
        @test perturbative_order(provenance) == 2
        @test (provenance.M, provenance.l, provenance.m, provenance.omega,
            provenance.rb, provenance.tail_x0) ==
            (mode.boundary.M, mode.boundary.l, mode.boundary.m, mode.boundary.omega,
             mode.boundary.rb, mode.boundary.x0)
        @test (provenance.primary_source_fingerprint,
            provenance.partner_source_fingerprint, provenance.pair_fingerprint,
            provenance.tail_boundary_fingerprint) ==
            (mode.boundary.primary_source_fingerprint,
             mode.boundary.partner_source_fingerprint, mode.boundary.pair_fingerprint,
             mode.boundary.boundary_fingerprint)
        @test provenance.policy isa FiniteRayFreeIVP
        @test provenance.normalization isa EinsteinRHS
        @test provenance.gauge isa GHZIRG
        @test provenance.tetrad isa KinnersleyMinusPlus
        @test provenance.slicing isa SchwarzschildOutgoing
        @test provenance.phase isa ExpMinusIomegaUPlusImPhi

        diagnostics = cache.quadrature_diagnostics
        streams = (diagnostics.primary_sll, diagnostics.primary_slm,
            diagnostics.partner_sll, diagnostics.partner_slm,
            diagnostics.primary_sln)
        @test all(values -> values[1] == 0.0, streams)
        @test all(issorted, streams)
        @test !hasproperty(diagnostics, :total)
    end

    @testset "one-cell diagnostics equal one local indicator" begin
        mode = rmc_manufactured()
        a, b = 2.0, 3.0
        atol, rtol = 0.0, 2e-12
        cache = GHZ_RMC.materialize_grid(mode, [a, b]; atol, rtol)
        primary, partner = mode.pair.primary, mode.pair.partner
        phead, qhead = primary.head, partner.head
        _, ep_sll = GHZ_RMC._ray_sll_moments(
            phead, a, b, atol, rtol, "direct primary S_ll")
        _, ep_slm = GHZ_RMC._ray_slm_moments(
            primary, a, b, atol, rtol, "direct primary S_lm")
        _, eq_sll = GHZ_RMC._ray_sll_moments(
            qhead, a, b, atol, rtol, "direct partner S_ll")
        _, eq_slm = GHZ_RMC._ray_slm_moments(
            partner, a, b, atol, rtol, "direct partner S_lm")
        _, ep_sln = GHZ_RMC._ray_sln_I2_over_b(
            phead.source, a, b, atol, rtol)
        point = GHZ_RMC.ray_quadrature_diagnostics(cache, b)
        @test point.primary_sll == ep_sll
        @test point.primary_slm == ep_slm
        @test point.partner_sll == eq_sll
        @test point.partner_slm == eq_slm
        @test point.primary_sln == ep_sln
    end

    @testset "BigFloat-256 and short-cell cancellation" begin
        scale = floatmax(Float64)
        scale_mode = rmc_zero_mode(Float64; x0=0.75scale)
        scale_cache = GHZ_RMC.materialize_grid(scale_mode, [2.0, 4.0])
        @test scale_cache.N[2] == ComplexF64(0.375scale, 0.0)

        # Raw q^2/q^3 moments overflow Float16 on this ray even though every
        # normalized recurrence value and the scalar correctness oracle are finite.
        long_mode = rmc_zero_mode(Float16; x0=one(Float16))
        long_grid = Float16[2, 400]
        long_cache = GHZ_RMC.materialize_grid(long_mode, long_grid)
        @test long_cache.N[2] == x_nn(long_mode, Float16(400))
        @test long_cache.N[2] == Complex{Float16}(Float16(2) / Float16(400))
        @test all(iszero, (long_cache.A[2], long_cache.dA[2], long_cache.B[2],
            long_cache.dB[2], long_cache.C[2], long_cache.dC[2]))

        # Multiplying H*b before the short-cell factor used to overflow here,
        # although the exact homogeneous recurrence remains representable.
        a = 1.0e200
        b = nextfloat(a)
        dlarge = ComplexF64(1.0e120)
        zero_sll = zeros(ComplexF64, 6)
        zero_slm = zeros(ComplexF64, 3)
        head_value, head_derivative = GHZ_RMC._ray_head_step(
            0.0im, dlarge, a, b, zero_sll)
        middle_value, middle_derivative = GHZ_RMC._ray_middle_step(
            0.0im, dlarge, 0.0im, a, b, zero_sll, zero_slm, 1.0)
        middle_reference = dlarge * ((b - a) *
            (1 + a / b + (a / b)^2) / 3)
        @test head_value == dlarge * ((b - a) * (a / b))
        @test head_derivative == (a / b)^2 * dlarge
        @test middle_value ≈ middle_reference rtol=4eps(Float64)
        @test isfinite(real(middle_derivative))
        @test isfinite(real(GHZ_RMC._ray_integral_A_over_b(
            0.0im, dlarge, a, b, zero_sll)))
        @test isfinite(real(GHZ_RMC._ray_integral_B_over_b(
            0.0im, dlarge, 0.0im, a, b, zero_sll, zero_slm, 1.0)))

        setprecision(BigFloat, 256) do
            mode = rmc_manufactured(BigFloat)
            grid = BigFloat[big"2", big"2.000000000000000000000000000001",
                big"2.0000000000000001", big"2.75", big"3.5"]
            cache = GHZ_RMC.materialize_grid(mode, grid; rtol=big"1e-55")
            @test eltype(cache.A) == Complex{BigFloat}
            rmc_compare_scalar(cache, mode, grid; rtol=big"2e-48", atol=big"1e-58")
            @test cache.quadrature_diagnostics.primary_sll isa AbstractVector{BigFloat}
        end

        setprecision(BigFloat, 2048) do
            y = BigFloat(1) / BigFloat(8)
            phi_reference = y - log1p(y)
            F_reference = (one(BigFloat) + y)^2 / BigFloat(4) +
                BigFloat(3) / BigFloat(4) - inv(one(BigFloat) + y) -
                BigFloat(3) * log1p(y) / BigFloat(2)
            phi_series = GHZ_RMC._ray_delta_minus_log1p(y)
            F_series = GHZ_RMC._ray_F_kernel(y)
            relative_tolerance = sqrt(eps(BigFloat))
            @test GHZ_RMC._ray_series_cap(y) > 256
            @test abs(phi_series - phi_reference) / abs(phi_reference) <=
                relative_tolerance
            @test abs(F_series - F_reference) / abs(F_reference) <=
                relative_tolerance
            @test precision(phi_series) == 2048
            @test precision(F_series) == 2048

            high_precision_mode = rmc_zero_mode(BigFloat; x0=one(BigFloat))
            high_precision_grid = BigFloat[BigFloat(2), BigFloat(9) / BigFloat(4)]
            high_precision_cache = GHZ_RMC.materialize_grid(
                high_precision_mode, high_precision_grid; rtol=relative_tolerance)
            @test high_precision_cache.N[2] == Complex{BigFloat}(
                high_precision_grid[1] / high_precision_grid[2])
            @test precision(real(high_precision_cache.N[2])) == 2048
        end

        mode = rmc_zero_mode(Float64; x0=0.3 - 0.2im,
            primary_kwargs=(; head_x0=0.4 + 0.1im, head_dx0=-0.2 + 0.05im,
                middle_x0=0.1 - 0.3im, middle_dx0=0.25 + 0.15im,
                Sll=q -> 0.2q - 0.1, Slm=q -> -0.3 + 0.05q,
                Sln=q -> 0.4 - 0.2im),
            partner_kwargs=(; head_x0=-0.2 + 0.3im, head_dx0=0.1 - 0.2im,
                middle_x0=-0.4 + 0.2im, middle_dx0=0.15 - 0.05im,
                Sll=q -> -0.1q, Slm=q -> 0.2 + 0.03q))
        grid = [2.0, nextfloat(2.0), 2.0 + 1e-12, 2.000001, 2.2]
        cache = GHZ_RMC.materialize_grid(mode, grid; rtol=2e-12)
        rmc_compare_scalar(cache, mode, grid; rtol=2e-9, atol=2e-13)

        for y in (eps(Float64), 1e-12, 1e-8)
            phi_reference = Float64(BigFloat(y) - log1p(BigFloat(y)))
            F_reference = Float64((one(BigFloat) + BigFloat(y))^2 / 4 + big"0.75" -
                inv(one(BigFloat) + BigFloat(y)) - big"1.5" * log1p(BigFloat(y)))
            @test GHZ_RMC._ray_delta_minus_log1p(y) ≈ phi_reference rtol=5e-15
            @test GHZ_RMC._ray_F_kernel(y) ≈ F_reference rtol=8e-15
            @test GHZ_RMC._ray_tF_kernel(y) ≈
                F_reference / (1 + y) rtol=8e-15
        end
        for y in (0.5, 199.0, 1.0e200)
            ybig = BigFloat(y)
            scaled_F_reference = Float64(((one(BigFloat) + ybig)^2 / 4 +
                big"0.75" - inv(one(BigFloat) + ybig) -
                big"1.5" * log1p(ybig)) / (one(BigFloat) + ybig))
            @test GHZ_RMC._ray_tF_kernel(y) ≈ scaled_F_reference rtol=8e-15
        end
        setprecision(BigFloat, 256) do
            tolerance16 = BigFloat(sqrt(eps(Float16)))
            for y16 in (nextfloat(Float16(0.125)), Float16(0.13), Float16(0.2),
                    Float16(0.625), nextfloat(Float16(0.625)))
                ybig = BigFloat(y16)
                phi_reference = ybig - log1p(ybig)
                F_reference = (one(BigFloat) + ybig)^2 / BigFloat(4) +
                    BigFloat(3) / BigFloat(4) - inv(one(BigFloat) + ybig) -
                    BigFloat(3) * log1p(ybig) / BigFloat(2)
                tF_reference = F_reference / (one(BigFloat) + ybig)
                phi16 = BigFloat(GHZ_RMC._ray_delta_minus_log1p(y16))
                F16 = BigFloat(GHZ_RMC._ray_F_kernel(y16))
                tF16 = BigFloat(GHZ_RMC._ray_tF_kernel(y16))
                @test abs(phi16 - phi_reference) / abs(phi_reference) <= tolerance16
                @test abs(F16 - F_reference) / abs(F_reference) <= tolerance16
                @test abs(tF16 - tF_reference) / abs(tF_reference) <= tolerance16
            end
        end
    end

    @testset "piecewise sources on an explicit cell boundary" begin
        jump = 3.0
        primary = rmc_middle(Val(2), Float64; l=3, m=1,
            omega=0.13 - 0.04im, source_fingerprint=RMC_PRIMARY_FP + 20,
            Sll=q -> q < jump ? 0.2q - 0.1im : -0.3q + 0.4im,
            Slm=q -> q < jump ? -0.15 + 0.07q : 0.25 - 0.03q,
            Sln=q -> q < jump ? 0.5 - 0.2im : -0.1 + 0.3im,
            head_x0=0.2, head_dx0=-0.1, middle_x0=0.3im,
            middle_dx0=0.05)
        partner = rmc_middle(Val(2), Float64; l=3, m=-1,
            omega=-conj(0.13 - 0.04im), source_fingerprint=RMC_PARTNER_FP + 20,
            Sll=q -> q < jump ? -0.1q : 0.2q + 0.1im,
            Slm=q -> q < jump ? 0.08q : -0.12q + 0.2im,
            head_x0=-0.15, head_dx0=0.04im, middle_x0=0.1 - 0.2im,
            middle_dx0=-0.06)
        mode = rmc_mode(primary, partner; x0=0.2 - 0.3im,
            pair_fingerprint=RMC_PAIR_FP + 20,
            boundary_fingerprint=RMC_BOUNDARY_FP + 20)
        grid = [2.0, jump, 4.0, 5.0]
        cache = GHZ_RMC.materialize_grid(mode, grid; rtol=2e-11)
        # The prefix oracle has no breakpoint API.  Compare it only through the jump, then
        # restart the same scalar IVPs at that explicit boundary for the post-jump cells.
        rmc_compare_scalar(cache, mode, grid[1:2]; rtol=2e-9, atol=2e-11)
        jump_index = 2
        restarted_head = solve_irg_corrector_head(primary.head.source,
            InnerRayIVP(primary.head.source; rb=jump, x_mmb=cache.A[jump_index],
                dr_x_mmb=cache.dA[jump_index]))
        restarted_primary = solve_irg_corrector_middle(restarted_head,
            MiddleRayIVP(restarted_head; rb=jump, x_nm=cache.B[jump_index],
                dr_x_nm=cache.dB[jump_index]))
        parity = GHZ_RMC._barred_middle_parity(Float64, primary.head.source.m)
        restarted_partner_head = solve_irg_corrector_head(partner.head.source,
            InnerRayIVP(partner.head.source; rb=jump,
                x_mmb=x_mmb(partner.head, jump),
                dr_x_mmb=dr_x_mmb(partner.head, jump)))
        restarted_partner = solve_irg_corrector_middle(restarted_partner_head,
            MiddleRayIVP(restarted_partner_head; rb=jump,
                x_nm=parity * conj(cache.C[jump_index]),
                dr_x_nm=parity * conj(cache.dC[jump_index])))
        restarted = rmc_mode(restarted_primary, restarted_partner;
            x0=cache.N[jump_index], pair_fingerprint=RMC_PAIR_FP + 21,
            boundary_fingerprint=RMC_BOUNDARY_FP + 21)
        rmc_compare_scalar(cache, restarted, grid[2:end]; rtol=3e-9, atol=3e-11)
    end

    @testset "deterministic linear source sampling and source-free lookup" begin
        counters = [Ref(0) for _ in 1:5]
        counted(index, value) = q -> begin
            counters[index][] += 1
            value + zero(q)
        end
        primary = rmc_middle(Val(2), Float64; source_fingerprint=RMC_PRIMARY_FP + 30,
            Sll=counted(1, 0.2), Slm=counted(2, -0.1),
            Sln=counted(3, 0.05 - 0.02im))
        partner = rmc_middle(Val(2), Float64; m=-1, omega=-conj(0.2 - 0.1im),
            source_fingerprint=RMC_PARTNER_FP + 30,
            Sll=counted(4, -0.15), Slm=counted(5, 0.08))
        mode = rmc_mode(primary, partner; pair_fingerprint=RMC_PAIR_FP + 30,
            boundary_fingerprint=RMC_BOUNDARY_FP + 30)
        function calls_for(grid)
            foreach(counter -> (counter[] = 0), counters)
            cache = GHZ_RMC.materialize_grid(mode, grid; rtol=1e-11)
            sum(counter[] for counter in counters), cache
        end
        base_calls, _ = calls_for([2.0])
        four_calls, cache4 = calls_for(collect(2.0:0.5:3.5))
        seven_calls, _ = calls_for(collect(2.0:0.5:5.0))
        @test four_calls > base_calls
        @test seven_calls - base_calls == 2 * (four_calls - base_calls)

        before = sum(counter[] for counter in counters)
        GHZ_RMC.ray_moment_components(cache4, 3.0)
        GHZ_RMC.x_mmb(cache4, 3.0)
        GHZ_RMC.x_nn(cache4, 3.0)
        @test sum(counter[] for counter in counters) == before
    end

    @testset "grid, source, tolerance, provenance, and raw rejection" begin
        mode = rmc_manufactured()
        for grid in (Float64[], [2.1, 3.0], [2.0, 2.0], [2.0, 1.9],
                [2.0, Inf], [2.0, NaN], [2.0, true], (2.0, 3.0),
                [2.0, 2.5 + 0im])
            @test_throws ArgumentError GHZ_RMC.materialize_grid(mode, grid)
        end
        for kwargs in ((; atol=0.0, rtol=0.0), (; atol=-1.0), (; rtol=Inf),
                (; atol=true), (; rtol=true))
            @test_throws ArgumentError GHZ_RMC.materialize_grid(mode, [2.0, 3.0];
                kwargs...)
        end

        late_bad = rmc_zero_mode(Float64;
            primary_kwargs=(; Sln=q -> q > 2.5 ? Inf : 0.0))
        @test_throws ArgumentError GHZ_RMC.materialize_grid(late_bad, [2.0, 3.0])
        late_bool = rmc_zero_mode(Float64;
            primary_kwargs=(; Sll=q -> q > 2.5 ? true : 0.0))
        @test_throws ArgumentError GHZ_RMC.materialize_grid(late_bool, [2.0, 3.0])

        bad_boundary = rmc_raw_boundary(mode.boundary;
            primary_source_fingerprint=mode.boundary.primary_source_fingerprint + 1)
        bad_mode = GHZ_RMC.IRGCorrectorTailMode{2,Float64,typeof(mode.pair)}(
            mode.pair, bad_boundary, GHZIRG())
        @test_throws ArgumentError GHZ_RMC.materialize_grid(bad_mode, [2.0, 3.0])

        cache = GHZ_RMC.materialize_grid(mode, [2.0, 3.0, 4.0])
        for radius in (2.5, 1.0, Inf, NaN, true, Complex{Bool}(true, false))
            @test_throws ArgumentError GHZ_RMC.x_nn(cache, radius)
        end
        @test GHZ_RMC.x_nn(cache, 2) == cache.N[1]

        snapshots = (cache.radii, cache.A, cache.dA, cache.B, cache.dB,
            cache.C, cache.dC, cache.N,
            cache.quadrature_diagnostics.primary_sll,
            cache.quadrature_diagnostics.primary_slm,
            cache.quadrature_diagnostics.partner_sll,
            cache.quadrature_diagnostics.partner_slm,
            cache.quadrature_diagnostics.primary_sln)
        for values in snapshots
            @test_throws ArgumentError setindex!(values, values[1], 1)
        end
        @test_throws ArgumentError cache.N._storage
        detached_N = copy(cache.N)
        sealed_cache = rmc_raw_cache(cache; Nvalues=detached_N)
        detached_N[2] += 100
        @test sealed_cache.N[2] == cache.N[2]

        bad_A = copy(cache.A)
        bad_A[2] = ComplexF64(Inf, 0)
        @test_throws ArgumentError GHZ_RMC.x_mmb(rmc_raw_cache(cache; A=bad_A), 3.0)
        bad_grid = copy(cache.radii)
        bad_grid[2] = bad_grid[1]
        @test_throws ArgumentError GHZ_RMC.x_nn(rmc_raw_cache(cache; radii=bad_grid), 3.0)
        bad_primary_diagnostic = copy(cache.quadrature_diagnostics.primary_sll)
        bad_primary_diagnostic[2] = -1.0
        @test_throws ArgumentError GHZ_RMC.RayMomentQuadratureDiagnostics{Float64}(
            bad_primary_diagnostic,
            copy(cache.quadrature_diagnostics.primary_slm),
            copy(cache.quadrature_diagnostics.partner_sll),
            copy(cache.quadrature_diagnostics.partner_slm),
            copy(cache.quadrature_diagnostics.primary_sln))
    end

    @testset "method-table ambiguity" begin
        ambiguities = Test.detect_ambiguities(GHZ_RMC; recursive=true)
        ours = filter(pair -> occursin("RayMoment", sprint(show, pair)), ambiguities)
        @test isempty(ours)
    end
end
