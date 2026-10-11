# C-66: direct include until registration by the boundary/context integration lane.
using Test
using LinearAlgebra
using Printf
using AccretionDisc
using AccretionDisc.SpiralWaves
const BoundarySW = AccretionDisc.SpiralWaves
isdefined(BoundarySW, :radiation_impedance_row) ||
    Base.include(BoundarySW, joinpath(@__DIR__, "../src/SpiralWaves/boundary.jl"))

module BoundaryEmission
    include(joinpath(@__DIR__, "../src/SpiralWaves/generated/eccentric_adiabatic.jl"))
end

# Independent local acoustic model with an arbitrary real amplitude trace. A=1
# here; a positive GR redshift/radial metric factor changes speed, never its sign.
function boundary_acoustic(σ, κ, z, α)
    C = [α inv(z); -(σ^2-κ^2)*z α]
    dC = [zero(σ) zero(σ); -2σ*z zero(σ)]
    C, dC
end

# Exact manufactured outgoing state exp((α+iq)x)[1,iq exp(ηx)]. Its basis
# varies, so midpoint errors create a measurable incoming component in the
# interior. S=y'-Cy is analytic, not a discrete manufactured source.
function boundary_reflection(n, σ, κ, outward; α=0.2, η=0.3, T=Float64)
    σ, κ, α, η = T.((σ, κ, α, η))
    q = (outward ? 1 : -1) * sign(σ) * sqrt(σ^2-κ^2)
    x = collect(range(zero(T), one(T); length=n+1))
    exact(r) = exp((α+im*q)*r) .* Complex{T}[1, im*q*exp(η*r)]
    coeff(r) = boundary_acoustic(σ, κ, exp(η*r), α)
    rm = (x[1:end-1] + x[2:end])/2
    Cs = [complex.(first(coeff(r))) for r in rm]
    Ss = [Complex{T}[0, im*q*η*exp(η*r)*exp((α+im*q)*r)] for r in rm]
    edge = outward ? last(x) : first(x)
    C, dC = coeff(edge)
    _, _, row = outgoing_branch(C, dC; outward)
    y = outward ? solve_midpoint(x, Cs, Ss, (1, 0), exact(first(x))[1], row, zero(Complex{T})) :
                   solve_midpoint(x, Cs, Ss, row, zero(Complex{T}), (1, 0), exact(last(x))[1])
    # Never evaluate R on the imposed boundary (where R=0 is tautological).
    i = n÷2+1
    C, dC = coeff(x[i])
    amp = wkb_amplitudes(C, dC, y[:, i]; outward)
    abs(amp.R), maximum(abs, y[:, i]-exact(x[i]))
end

function boundary_emitted(bg, m, r, ω)
    z = zero(r)
    hm = (htt=z, htt_r=z, htphi=z, htphi_r=z, hphiphi=z, hphiphi_r=z, hrr=z, hzz=z)
    inp, omt, vr = r2_context(bg, max(m, 1), r; hmodes=hm, landau=false)
    # Adapt existing circular on-shell background inputs, not its frequency.
    p = (; (s => getproperty(inp, s) for s in fieldnames(typeof(inp)))...)
    p = merge(p, (enu=exp(inp.nu), epsi=exp(inp.psi), emu=exp(inp.mu),
                  omt=omt, v_r=vr, m=m, omega_d=ω-m*Omega(bg, r)))
    E = BoundaryEmission
    [E._eccentric_adiabatic_C00(p) E._eccentric_adiabatic_C01(p);
     E._eccentric_adiabatic_C10(p) E._eccentric_adiabatic_C11(p)]
end

@testset "C-66 arbitrary-frequency boundary" begin
    @testset "physical edge normal, negative/small omega, m=0, turning limits" begin
        for T in (Float64, BigFloat), σ in T.((-3, -1//100000, 1//100000, 3)), outward in (false, true)
            κ = abs(σ)/2
            C, dC = boundary_acoustic(σ, κ, T(2), T(20))
            ka, Ya, kt, Yt = wkb_split(C, dC; outward)
            vg = imag(ka)/σ
            @test (outward ? 1 : -1)*vg > 0
            @test C*collect(Ya) ≈ ka*collect(Ya)
            _, _, row = outgoing_branch(C, dC; outward)
            @test abs(row[1]*Ya[1]+row[2]*Ya[2]) ≤ 10eps(T)
            y = collect(Ya) + T(1//5)*collect(Yt)
            a = wkb_amplitudes(C, dC, y; outward)
            @test a.R ≈ T(1//5) rtol=100eps(T)
            # Large common trace must not select the evanescent branch.
            @test abs(imag(ka)) < abs(real(ka))
        end
        for T in (Float64, BigFloat), δ in T.((1//100, 1//1000000)), outward in (false, true)
            C, dC = boundary_acoustic(one(T), sqrt(one(T)-δ), one(T), T(1//10))
            ka, _, _, _ = wkb_split(C, dC; outward)
            @test (outward ? 1 : -1)*imag(ka) > 0
            @test abs(imag(ka)) ≈ sqrt(δ) rtol=T(1//100000000)
            # Evanescent side, with roots that really decay/grow.
            C, dC = boundary_acoustic(one(T), sqrt(one(T)+δ), one(T), zero(T))
            ka, _, _, _ = wkb_split(C, dC; outward)
            @test iszero(imag(ka))
            @test (outward ? 1 : -1)*real(ka) < 0
        end
        # An omega-dependent amplitude trace moves Re(lambda) only: k_r(omega), hence v_g, is
        # unchanged. Here alpha*alpha_omega/q outweighs q_omega, so a group slope without
        # tr(C)_omega would flip the selection (lander review R1, mutant M10).
        for T in (Float64, BigFloat), σ in T.((-1, 1)), outward in (false, true)
            C, dC = boundary_acoustic(σ, T(1//2), one(T), T(4))
            ka, _, _, _ = wkb_split(C, dC - sign(σ)*I; outward)
            @test (outward ? 1 : -1)*imag(ka)/σ > 0
        end
        C, dC = boundary_acoustic(1.0, 1.0, 1.0, 0.0)
        @test_throws DomainError wkb_split(C, dC; outward=true)
        C, dC = boundary_acoustic(1.0, 2.0, 1.0, 10.0)
        @test_throws DomainError wkb_split(C, dC; outward=true) # no exterior decay
        C, dC = boundary_acoustic(0.0, 1.0, 1.0, 0.0)
        @test real(first(wkb_split(C, dC; outward=true))) < 0 # regular static evanescence
        C, dC = boundary_acoustic(2.0, 1.0, 1.0, 0.0)
        @test_throws DomainError wkb_split(C, zero(dC); outward=true)
        @test_throws ArgumentError wkb_split(C .+ 0.01im, dC; outward=true)
        @test_throws ArgumentError wkb_split(C, dC; outward=true, bc=:bad)
        @test_throws DimensionMismatch wkb_split(C, zeros(3,3); outward=true)
        C, dC = boundary_acoustic(-1.0,sqrt(1-1e-6),1.0,0.1)
        old = first(outgoing_branch(complex.(C),-1.0; outward=true))
        corrected = first(outgoing_branch(C,dC; outward=true))
        @test imag(old) > 0 > imag(corrected) # scalar classifier's near-turning counterexample
        for bits in (128,256)
            setprecision(BigFloat,bits) do
                C, dC = boundary_acoustic(big"2",big"1",big"1",big"0")
                k, Y, _ = outgoing_branch(C,dC; outward=true)
                @test precision(real(k)) == bits
                @test precision(real(Y[1])) == bits
            end
        end
        println("PASS numerical: edge-normal group velocity, evanescence, turning-point guards, amplitude split")
    end

    @testset "printed sign is conditional on edge and Doppler frequency" begin
        for σ in (-2.0, 2.0), outward in (false, true)
            C, dC = boundary_acoustic(σ, 1.0, 1.0, 0.0)
            kp = first(outgoing_branch(C, dC; outward, bc=:paper))
            ko = first(outgoing_branch(C, dC; outward))
            @test imag(kp) < 0
            @test (kp == ko) == ((outward ? 1 : -1)*σ < 0)
        end
        println("PASS numerical: printed Im(lambda)<0 outgoing iff n*omega_d<0 in acoustic WKB")
    end

    @testset "manufactured reflection is globally second order" begin
        cases = (("above", 3.0, 2, 0.5, 0.5), ("below", 0.5, 2, 1.0, 0.5),
                 ("negative", -1.0, 2, 0.5, 0.5), ("axisymmetric_k1", 2.0, 0, 1.0, 0.5),
                 ("small", 1e-5, 2, 0.5, 0.5),
                 ("near_turn", sqrt(1+0.05^2), 0, 1.0, 1.0))
        for (name, ω, m, Ω, κ) in cases, outward in (false, true)
            σ = ω-m*Ω
            errors = [boundary_reflection(n, σ, κ, outward)[1] for n in (32,64,128,256)]
            ratios = errors[1:end-1]./errors[2:end]
            @test all(3.8 .< ratios .< 4.2)
            @test last(errors) < 1e-5
            @printf("REFLECTION %-16s %-5s %.9e %.9e %.9e %.9e order=%.6f\n",
                    name, outward ? "outer" : "inner", errors..., log2(last(ratios)))
        end
        r64 = boundary_reflection(64, -2, 1//2, false)[1]
        rb = boundary_reflection(64, -2, 1//2, false; T=BigFloat)[1]
        @test rb isa BigFloat
        @test abs(rb-r64) < 1e-12
        println("PASS numerical: interior incoming/outgoing amplitude tends to zero at order 2; BigFloat banded parity")
    end

    @testset "impedance gap, gauge covariance and exact evanescent transport" begin
        for T in (Float64, BigFloat)
            ka, kt, η = T(-2), T(3), T(1//5)
            Ya, Yt, dYa = (one(T), T(1//4)), (zero(T), one(T)), (zero(T), η)
            row = BoundarySW.radiation_impedance_row(ka, Ya, kt, Yt, dYa)
            β = η/(kt-ka)
            y = Ya .+ β .* Yt
            @test abs(row[1]*y[1]+row[2]*y[2]) < 20eps(T)
            oldrow, _ = BoundarySW._impedance_row(Ya,Yt,kt,dYa,(zero(T),zero(T)),(zero(T),zero(T)))
            @test abs(oldrow[1]*y[1]+oldrow[2]*y[2]) > η/100
            ga, gt, dga = T(3), T(7), T(2)
            row2 = BoundarySW.radiation_impedance_row(ka, ga.*Ya, kt, gt.*Yt, dga.*Ya .+ ga.*dYa)
            @test collect(row2) ≈ ga.*collect(row) rtol=20eps(T)
            r0 = BoundarySW.radiation_impedance_row(ka, Ya, kt, Yt, (zero(T),zero(T)))
            @test r0 == (-Ya[2], Ya[1])
            @test_throws DomainError BoundarySW.radiation_impedance_row(ka, Ya, ka, Ya, dYa)
            # Real triangular basis P=[1 0; ηx 1] has N21=η, others zero.
            # P[1,β] exp(ka*x) is an EXACT homogeneous outgoing/decaying state.
            errs = T[]
            for n in (32,64,128)
                x = collect(range(zero(T), one(T); length=n+1))
                rm = (x[1:end-1]+x[2:end])/2
                Cs = [Complex{T}[ka 0; η*r*(ka-kt) kt] for r in rm]
                Ss = [zeros(Complex{T},2) for r in rm]
                row = BoundarySW.radiation_impedance_row(ka, (one(T),η), kt, Yt, dYa)
                y = solve_midpoint(x,Cs,Ss,(1,0),one(Complex{T}),row,zero(Complex{T}))
                i=n÷2+1
                push!(errs, abs(y[2,i]/y[1,i] - η*x[i] - β))
            end
            @test all(T(3.8) .< errs[1:2]./errs[2:3] .< T(4.2))
        end
        println("PASS numerical: homogeneous impedance uses eigenvalue gap; gauge covariance; second-order evanescent solve")
    end

    @testset "R3 general adiabatic emission at arbitrary frequency" begin
        for T in (Float64, BigFloat)
            bg = DiscBackground{T}(; χ=T(1//2), aEMRI=T(15), h0=T(1//20))
            for (m, ω) in ((2, 3Omega(bg,T(8))), (2, Omega(bg,T(35))/2),
                           (2, -Omega(bg,T(8))), (2,T(1//100000000)), (0,Omega(bg,T(8))))
                for (r,outward) in ((T(8),false),(T(35),true))
                    C = boundary_emitted(bg,m,r,ω)
                    h = max(abs(ω),Omega(bg,r))*eps(T)^(3//4)
                    dC = imag.(boundary_emitted(bg,m,r,ω+im*h))/h
                    ka, Ya, _, _ = wkb_split(C,dC; outward)
                    @test maximum(abs, C*collect(Ya)-ka*collect(Ya)) <
                          500eps(T)*maximum(abs,C)
                    if !iszero(imag(ka))
                        # Independent centred differentiation of the selected root.
                        hfd = max(abs(ω),Omega(bg,r))*cbrt(eps(T))
                        kp = first(wkb_split(boundary_emitted(bg,m,r,ω+hfd),dC; outward))
                        km = first(wkb_split(boundary_emitted(bg,m,r,ω-hfd),dC; outward))
                        slope = imag(kp-km)/(2hfd)
                        @test (outward ? 1 : -1)*slope > 0
                    else
                        @test (outward ? 1 : -1)*real(ka) < 0
                    end
                end
            end
        end
        println("PASS numerical: R3 emitted C, fixed-disc omega derivative, both edges, negative/small omega and m=0")
    end
end
