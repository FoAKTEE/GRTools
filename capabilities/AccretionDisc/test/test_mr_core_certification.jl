# ECC::mr19-core-certification: thermal endpoints supplement the frozen MR fixture gates
# in Modes/test_newtonian_planet.jl. This is numerical characterization, not a new
# FR master-equation certificate or a claim to reproduce all sixteen MR19 figures.
using Test
using AccretionDisc

@testset "MR core thermal endpoints" begin
    M = AccretionDisc.Modes
    for T in (Float64, BigFloat), q in (0//1, 3//4, 1//1)
        @testset "$T q=$q" begin
            disc(beta; gamma = T(7) / 5) = NewtonianDisc{T}(
                h0 = T(1) / 10, q = T(q), p = one(T), gamma = gamma,
                beta = T(beta), soft = T(3) / 50, pole_disp = T(1) / 10^6)
            d0 = disc(0)
            x = T[3//10, 7//10, 13//10, 3]
            # Prescribed response jets isolate reconstruction from BVP integration error.
            f = hcat(complex.(x, x ./ 2), complex.(one.(x), -x))
            for m in (2, 5)
                _, _, sig0 = vars_from_enthalpy(d0, f, x, m)
                expected0 = [M.Sigma(d0, r) * f[i, 1] / M.csiso(d0, r)^2
                             for (i, r) in enumerate(x)]
                @test eltype(sig0) == Complex{T}
                @test isapprox(sig0, expected0; rtol = 100eps(T))
                for r in x
                    @test M.tc(d0, r) == M.dtcdr(d0, r) == 0
                    k0 = sqrt(complex(-M.Dfreq(d0, r, m))) / M.csiso(d0, r)
                    @test isapprox(M.kwave(d0, r, m), k0; rtol = 100eps(T))
                    # Instantaneous cooling removes gamma from the runtime operator,
                    # including nonzero temperature gradients and the forcing term.
                    y = T[1//3, -1//5, 2//7, 1//11]
                    rhs0 = collect(M.master_rhs(d0, r, y, m, 1))
                    rhs1 = collect(M.master_rhs(disc(0; gamma = T(5) / 3), r, y, m, 1))
                    @test eltype(rhs0) == T
                    @test isapprox(rhs0, rhs1; rtol = 1000eps(T))
                end
                # Finite sequences approach the two distinct endpoints. No literal Inf
                # is passed to the generated rational coefficients; beta=80 remains
                # the captured finite-beta MR19 proxy, not an exact infinity assertion.
                for (label, betas) in ((:isothermal, (1//10^6, 1//10^8)),
                                       (:adiabatic, (10^6, 10^8)))
                    errors = T[]
                    for beta in betas
                        d = disc(beta)
                        ur, _, sig = vars_from_enthalpy(d, f, x, m)
                        err = zero(T)
                        for (i, r) in enumerate(x)
                            ot = M.omega_pat(d, m) - m * M.Omega(d, r) + im * d.pole_disp
                            if label === :isothermal
                                expected_sig = expected0[i]
                                expected_k = M.kwave(d0, r, m)
                            else
                                expected_sig = M.Sigma(d, r) *
                                    (f[i, 1] / M.csadi(d, r)^2 + im * M.LSinv(d, r) * ur[i] / ot)
                                expected_k = sqrt(complex(-(M.Dfreq(d, r, m) + M.Nr2(d, r)))) / M.csadi(d, r)
                            end
                            @test isapprox(M.kwave(d, r, m), expected_k; rtol = T(1) / 10^7)
                            err = max(err, abs(sig[i] - expected_sig) / abs(expected_sig))
                        end
                        push!(errors, err)
                    end
                    @test errors[2] < errors[1] / 50
                    @test errors[2] < T(1) / 10^5
                    println("thermal T=$T q=$q m=$m endpoint=$label reconstruction_rel=$errors")
                end
            end
        end
    end
end
