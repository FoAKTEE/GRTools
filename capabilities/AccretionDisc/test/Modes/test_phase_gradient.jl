using Test
using AccretionDisc

@testset "MR phase-gradient refinement" begin
    P = AccretionDisc.Modes
    @testset "analytic phase and boundary windows, inherited precision" begin
        for T in (Float64, BigFloat)
            d = NewtonianDisc{T}(h0=T(1)/10, q=one(T), p=one(T), gamma=T(7)/5,
                                beta=T(80), soft=T(6)/100, pole_disp=T(1)/10^6)
            x = loggrid(d, T(1)/20, T(5), 1200)
            h = exp.(im .* (2 .* x .+ x.^2 ./ 2))
            f = hcat(h, im .* (2 .+ x) .* h)
            err = phase_gradient_error(d, f, x, 10)
            @test eltype(err) == T
            @test err ≈ [2+r-real(P.kwave(d,r,10)) for r in x]
            @test phase_gradient_error(d, (2+3im).*f, x, 10) ≈ err
            ii = argmin(abs.(x .- (x[1]+P.csiso(d,x[1])/P.Omega(d,x[1]))))
            io = argmin(abs.(x .- (x[end]-P.csiso(d,x[end])/P.Omega(d,x[end]))))
            expected = (maximum(err[1:ii-1])-minimum(err[1:ii-1]),
                        maximum(err[io:end-1])-minimum(err[io:end-1]))
            @test collect(phase_error_diagnostic(d,f,x,10,1,length(x))) ≈ collect(expected)
        end
    end
    @testset "mirrored main.py one-step oracle and m=1 outer-only branch" begin
        # Captured by make_author_fixture.py, executing the mirrored function AST and secant loop.
        # Source SHA256: 36987b08433dff679dfece8513817efab269b2b3dc1ddf2520614a7f59e98e1a.
        d = NewtonianDisc(beta=80.0)
        x = loggrid(d, .05, 5., 1200)
        h = exp.(im .* (2 .* x .+ x.^2 ./ 2))
        h1 = hcat(h, im .* (2 .+ x) .* h)
        h2 = conj.(h1)
        ih = hcat(.03 .* x.^2 .+ .02im .* x, .06 .* x .+ .02im)
        Ci, Co = -.4+2.1im, .2+7.1im
        a1,a2 = P.homogeneous_coefficients(h1,h2,ih,Ci,Co,1,length(x))
        baseline = a1.*h1 .+ a2.*h2 .+ ih
        for (m, X, C) in (
            (1, [0.01303871548256108,4.356569738865822],
                [Ci,0.2601960456465893+9.236959620518665im]),
            (10, [147.52826621115764,11.912153824075546],
                [2793.6614914377988-14666.722829903036im,-1075.231346819605-38170.7128132527im]))
            @test collect(phase_error_diagnostic(d,baseline,x,m,1,length(x))) ≈ X rtol=1e-10
            got = phase_gradient_refine(d,x,m,h1,h2,ih,Ci,Co,1,length(x); baseline)
            @test collect(got.C) ≈ C rtol=1e-7
            m == 1 && (@test got.C[1] == Ci)
            if m == 1
                db=NewtonianDisc{BigFloat}(h0=big(1)/10,q=big(1),p=big(1),
                    gamma=big(7)/5,beta=big(80),soft=big(6)/100,pole_disp=big(1)/10^6)
                resultb=phase_gradient_refine(db,BigFloat.(x),m,Complex{BigFloat}.(h1),
                    Complex{BigFloat}.(h2),Complex{BigFloat}.(ih),Complex{BigFloat}(Ci),
                    Complex{BigFloat}(Co),1,length(x);baseline=Complex{BigFloat}.(baseline))
                @test eltype(resultb.mode) == Complex{BigFloat}
                @test collect(resultb.C) ≈ C rtol=1e-7
            end
            # The caller's baseline is used, even if it precedes the updated boundary combination.
            other = phase_gradient_refine(d,x,m,h1,h2,ih,Ci,Co,1,length(x); baseline=h1)
            @test other.C != got.C
        end
    end
    @testset "MR19 paper conventions" begin
        d = NewtonianDisc(h0=.1/sqrt(1.4), gamma=1.4, beta=1e8)
        @test P.csadi(d,1.) ≈ .1
        r=.3; m=10
        # Independent finite difference of the paper's logarithmic amplitude, Eq. (A3).
        amplitude_log(s) = log(abs((P.Dfreq(d,s,m)+P.Nr2(d,s))/(s*P.Sigma(d,s)*P.kwave(d,s,m))))
        dr=1e-5*r
        expected=im*P.kwave(d,r,m)+(amplitude_log(r+dr)-amplitude_log(r-dr))/(4dr)
        Ci,_=P.mr19_boundary_coeffs(d,[r,5.],1,2,m)
        @test Ci ≈ expected rtol=1e-10
    end
    @testset "solve_m_mode phase-step wiring (lander review R1)" begin
        # The cached path at beta=80 on 2000 points: boundaries stay at the grid ends, so the
        # pre-step baseline is the frozen default output, and the diagnostics must be the
        # stand-alone diagnostic of `before` and `mode`. Characterised values (review log
        # default-compat): m=10 X ratios (0.043, 0.239), dC/C (7.4e-3, 7.6e-5); m=1 C_in fixed.
        d = NewtonianDisc(beta=80.0, lap_nphi=800)
        x = loggrid(d, .05, 5., 2000)
        for m in (1, 10)
            cache = PotentialCache(d, m, .05, 5.; n=2000)
            off = solve_m_mode(d, m, x, cache)
            s = solve_m_mode(d, m, x, cache; phgrad_steps=1, diagnostics=true)
            @test s.indices == (1, length(x))
            @test isequal(s.before, off)
            @test collect(s.C_before) == collect(P.wkb_boundary_coeffs(d, x, 1, length(x), m))
            @test collect(s.X_before) == collect(phase_error_diagnostic(d, off, x, m, 1, length(x)))
            @test collect(s.X_after) == collect(phase_error_diagnostic(d, s.mode, x, m, 1, length(x)))
            @test s.mode != off
            if m == 1
                @test s.C_after[1] == s.C_before[1] && s.C_after[2] != s.C_before[2]
            else
                @test all(s.C_after .!= s.C_before)
                @test all(s.relative_C_change .< 0.01)       # MR19 App. A: a change of <~ 1 %
                @test s.X_after[1] < s.X_before[1] / 10      # characterised 0.043 at this grid
            end
        end
        @test_throws ArgumentError solve_m_mode(d, 10, x; phgrad_steps=2)
    end
end
