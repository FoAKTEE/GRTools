using Test
using AccretionDisc.Modes
using AccretionDisc.Modes: _kepler_newton
using NumericsCore: bisect

@testset "FR22 Kepler orbit" begin
    @test EccentricOrbit(1, 0, 1) isa EccentricOrbit{Float64}
    for (a, e, n) in ((0, 0, 1), (1, -0.1, 1), (1, 1, 1), (1, 0, 0), (Inf, 0, 1))
        @test_throws DomainError EccentricOrbit(a, e, n)
    end
    @test_throws DomainError kepler_E(Inf, 0.2)
    @test_throws DomainError kepler_E(0.3, NaN)
    @test_throws ArgumentError kepler_E(0.3, 0.2; maxiter=-1)

    function check_orbits(T)
        tol = eps(T)
        for e in T.((0, 1//10, 1//4, 1//2, 9//10, 99//100)),
            M in range(-T(pi), T(pi); length=101)
            E = kepler_E(M, e)
            @test abs(E - e*sin(E) - M) <= 4tol
            @test kepler_E(-M, e) == -E
            @test abs(kepler_E(M, e; maxiter=0) - E) <= 16tol / (1-e)
            # Independent bisection at ≥256 bits, without range reduction.
            oracle = setprecision(BigFloat, max(256, 2precision(T))) do
                mm, ee = BigFloat(M), BigFloat(e)
                bisect(x -> x-ee*sin(x)-mm, mm-ee, mm+ee; atol=eps(BigFloat))
            end
            @test abs(E - oracle) <= 8tol / (1-e)
            p = planet_orbit(M, e; a=T(2), n=T(3)/2)
            @test p.r_p ≈ 2(1-e*cos(E)) rtol=8tol
            @test p.r_p ≈ 2(1-e^2)/(1+e*cos(p.phi_p)) rtol=32tol/(1-e)
            @test sin(p.phi_p) ≈ sqrt(1-e^2)*sin(E)/(1-e*cos(E)) atol=32tol/(1-e)
            @test p.dE_dt * (1-e*cos(E)) ≈ T(3)/2 rtol=8tol
            @test p.r_p^2*p.dphi_dt ≈ 6sqrt(1-e^2) rtol=8tol
            @test planet_orbit(-M, e; a=2, n=T(3)/2).r_p == p.r_p
            @test planet_orbit(-M, e).phi_p == -p.phi_p
        end
        for e in T.((0, 1//2, 99//100))
            Ms = range(-10T(pi), 10T(pi); length=1001)
            Es = kepler_E.(Ms, e)
            @test all(diff(Es) .> 0)
            @test all(abs(E-e*sin(E)-M) <= 8tol*max(1,abs(M)) for (E,M) in zip(Es,Ms))
            @test all(diff([planet_orbit(M,e).phi_p for M in Ms]) .> 0)
        end
        for M in T.((-100, -1, 0, 1, 100))
            p = planet_orbit(M, zero(T); a=T(2), n=T(3))
            @test (p.r_p, p.phi_p, p.E, p.dr_dt, p.dphi_dt) == (2, M, M, 0, 3)
            q = planet_orbit(M, tol; a=T(2))
            @test abs(q.r_p-2) <= 4tol
            @test abs(q.phi_p-M) <= 4tol*max(1,abs(M))
        end
        @test kepler_E(T(1), T(1)-sqrt(tol)) isa T
        @test planet_orbit(T(1), T(1)/2).r_p isa T
        # Newton itself, not the bracket, converges (quadratically) for e ≤ 1/2: the
        # fallback must not mask a broken Newton step.
        for e in T.((1//10, 1//4, 1//2)), m in range(T(pi)/50, 49T(pi)/50; length=49)
            E, steps = _kepler_newton(m, e, 32)
            @test steps !== nothing && steps <= 8
            @test abs(E - e*sin(E) - m) <= 4tol
        end
    end
    check_orbits(Float64)
    setprecision(BigFloat, 256) do
        check_orbits(BigFloat)
        @test EccentricOrbit(1, big"0.25", 1) isa EccentricOrbit{BigFloat}
    end
end

@testset "ZZ22 wavelet owner" begin
    using NumericsCore: gl_nodes_weights
    for T in (Float64, BigFloat)
        setprecision(BigFloat, 256) do
            nodes, weights = gl_nodes_weights(T, 128)
            tol = T === Float64 ? T(2e-12) : T(10)^(-65)
            for e in T.((0, 1//10, 1//4, 1//2)), t1 in T.((0, 1, 3)), dir in (-1, 1)
                orb = EccentricOrbit(one(T), e, one(T))
                cs, age = T(1)/10, T(2)
                p = planet_orbit(t1, orb)
                w = zz22_wavelet(t1, t1+age, orb, dir; cs)
                @test w.R ≈ p.r_p+dir*cs*age rtol=8eps(T)
                # Independent quadrature of ZZ22 (7), split at the corotation kink.
                # Trailing-wave convention gives dtheta/dR=-|Omega_pat-Omega(R)|/cs.
                rc = p.dphi_dt^(-T(2)/3)
                lo, hi = minmax(p.r_p, w.R)
                bounds = lo < rc < hi ? (lo, rc, hi) : (lo, hi)
                integral = zero(T)
                for (a, b) in zip(bounds[1:end-1], bounds[2:end])
                    radii = (a+b)/2 .+ (b-a)/2 .* nodes
                    integral += (b-a)/2 * sum(weights .* abs.(p.dphi_dt .- radii.^(-T(3)/2))) / cs
                end
                theta = p.phi_p+p.dphi_dt*age-sign(w.R-p.r_p)*integral
                @test abs(w.theta-theta) <= tol
                launch = zz22_wavelet(t1, t1, orb, dir; cs)
                @test launch.R == p.r_p
                @test abs(launch.theta-p.phi_p) <= 8eps(T)*max(1,abs(p.phi_p))
                @test w.R isa T && w.theta isa T
                # Physical length/time scaling: a -> 2a, n -> 3n, cs -> 6cs.
                scaled = zz22_wavelet(t1/3, (t1+age)/3, EccentricOrbit(T(2),e,T(3)), dir; cs=6cs)
                @test abs(scaled.R/2-w.R) <= 64eps(T)
                @test abs(scaled.theta-w.theta) <= 256eps(T)*max(1,abs(w.theta))
            end
            orb = EccentricOrbit(one(T), T(1)/2, one(T))
            @test zz22_wavelet(zero(T), T(10), orb, -1; cs=T(1)/10) === nothing
            out = zz22_track(orb, T(3), 1; cs=T(1)/10, dt=T(1)/10)
            @test length(out) == 31
            @test first(out)[1] == 0 && last(out)[1] == 30
            @test all(0 <= row[3] <= 2 for row in out)
            @test eltype(out) == Tuple{Int,T,T}
            @test length(zz22_track(orb,T(3),1;cs=T(1)/10,dt=T(2))) == 2
        end
    end
    orb = EccentricOrbit(1, 0.25, 1)
    @test_throws ArgumentError zz22_wavelet(0, 1, orb, 0; cs=0.1)
    @test_throws DomainError zz22_wavelet(1, 0, orb, 1; cs=0.1)
    @test_throws DomainError zz22_wavelet(0, 1, orb, 1; cs=0)
    @test_throws DomainError zz22_track(orb,1,1;cs=0.1,dt=0)
end
