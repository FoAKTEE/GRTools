# SpinningSecondary/action_parts.jl — linear-in-spin ACTION parts + the J/pex gauge transforms.
#
# Grounded with BAKED fixtures captured once from the legacy `SpinningOrbit.KerrSpinActionParts`
# (packages/SpinningOrbit/src/linear_parts.jl, itself oracle-verified against the KerrSpinningFluxes
# `tests.nb` golden at δJ 7.6e-15 / J1pex 2.6e-12 / Ω1J 1.2e-14), plus structural identities that
# hold independently of the fixture: δJφ|_C ≡ 0, the Jacobian shapes, and the ∂Jφ/∂C = (0,1,0) row.
using Test
using Spacetime
using Spacetime.SpinningSecondary: spin_linear_actions, spin_actions_pex,
                                   spin_frequencies_action_gauge, spin_action_jacobian,
                                   spin_frequency_jacobian, spin_orbital_params_jacobian,
                                   spin_linear_orbital_params, spin_linear_frequencies,
                                   spin_linear_constants,
                                   spin_frequencies_constants_gauge,
                                   spin_frequencies_turning_point_gauge

# BAKED legacy fixture: (a, p, e, x) => (J1C, J1pex, Omega1J, pex1CTilde)
const _ACTION_REF = (
    ((0.99, 8.0, 0.1, 0.8),
     (0.8613657052897491, 0.22709577224478045, 0.0),
     [0.0005538635642057255, 0.3370416108887221, 0.4441396662390246],
     [0.022239654067965416, 0.031084566005145302, 0.034391910229075645],
     (-2.921464597969925, -0.0021622341777338128, 0.00541216941230111)),
    ((0.95, 5.0, 0.5, 0.5),
     (1.3177706144109236, 0.11863622431936217, 0.0),
     [0.07214609671547123, 0.6918218915591069, -0.08370609218662751],
     [0.01402796383879944, 0.07815265299028773, 0.10511191011606318],
     (-1.657159763332185, -0.009574565906878357, 0.032221406201730154)),
    ((0.9, 12.0, 0.2, 0.5),
     (1.0180653719760326, 0.0946702289053723, 0.0),
     [0.0018729832528654455, 0.6406377129922547, 0.19697548460284428],
     [0.011970133432548326, 0.016697083304758893, 0.01772689891580126],
     (-3.2539205578088453, -0.0009436845689590467, 0.0046583947222814765)),
    ((0.1, 8.0, 0.3, -0.6),
     (1.7006431639644715, -0.016680159384048036, 0.0),
     [0.033445505625128824, 0.23626590872331626, -0.4072714714987146],
     [0.000999197200244573, 0.07305453699659781, -0.07177421019381859],
     (-2.0679831349528026, 0.000518723478803107, 0.0010980629566186553)),
)

@testset "SpinningSecondary action parts (J / pex gauges)" begin
    @testset "baked legacy fixture: J1C, J1pex, Omega1J, pex1CTilde" begin
        for (orb, J1C, J1pex, Om1J, pex1) in _ACTION_REF
            a, p, e, x = orb
            @test collect(spin_linear_actions(a, p, e, x)) ≈ collect(J1C) rtol = 1e-11
            @test spin_actions_pex(a, p, e, x) ≈ J1pex rtol = 1e-10
            @test spin_frequencies_action_gauge(a, p, e, x) ≈ Om1J rtol = 1e-10
            @test collect(spin_linear_orbital_params(a, p, e, x)) ≈ collect(pex1) rtol = 1e-11
        end
    end

    @testset "structural identities" begin
        for (orb, _, _, _, _) in _ACTION_REF
            a, p, e, x = orb
            @test spin_linear_actions(a, p, e, x)[3] == 0        # Jφ ≡ Lz is held fixed with C
            dJ = spin_action_jacobian(a, p, e, x)
            @test size(dJ) == (3, 3)
            @test dJ[3, :] == [0.0, 1.0, 0.0]                    # ∂Jφ/∂(En,Lz,K) = (0,1,0)
            @test size(spin_frequency_jacobian(a, p, e, x)) == (3, 3)
            @test size(spin_orbital_params_jacobian(a, p, e, x)) == (3, 3)
            @test all(isfinite, dJ)
            @test all(isfinite, spin_frequency_jacobian(a, p, e, x))
            @test all(isfinite, spin_orbital_params_jacobian(a, p, e, x))
        end
    end

    @testset "generic pex1CTilde reduces to the equatorial branch at x = ±1" begin
        # the equatorial closed form (z1 = 0 ⟹ intΣr = 1/r², δx = 0) restated here as the oracle
        function eq_ref(a, p, e)
            En, Lz, Q = constants_of_motion(a, p, e, 1.0)
            K = Q + (Lz - a * En)^2; sK = sqrt(K)
            r1 = p / (1 - e); r2 = p / (1 + e)
            d1 = -((r1^2 + a^2) * En - a * Lz) / (r1 * sK)
            d2 = -((r2^2 + a^2) * En - a * Lz) / (r2 * sK)
            ((2r2^2 * d1 + 2r1^2 * d2) / (r1 + r2)^2, (2r2 * d1 - 2r1 * d2) / (r1 + r2)^2, 0.0)
        end
        for (a, p, e) in ((0.9, 10.0, 0.2), (0.5, 8.0, 0.4), (0.99, 6.0, 0.1))
            @test collect(spin_linear_orbital_params(a, p, e, 1.0)) ≈ collect(eq_ref(a, p, e)) rtol = 1e-12
            @test spin_linear_orbital_params(a, p, e, 1.0)[3] == 0
        end
    end

    @testset "C / pex frequency gauges agree with a finite-difference Jacobian" begin
        # the analytic ∂Ω/∂C and ∂(p,e,x)/∂C must reproduce a central-difference Jacobian of
        # kerr_frequencies / constants_of_motion (converted from (E,Lz,Q) to (E,Lz,K)) — the
        # construction the legacy arXiv:2603.13482 figure script used, to FD accuracy (~1e-6).
        freqvec(a, p, e, x) = (f = kerr_frequencies(a, p, e, x); [f.Omegar, f.Omegath, f.Omegaph])
        constvec(a, p, e, x) = collect(constants_of_motion(a, p, e, x))
        function fdjac(f, a, p, e, x; h = 1e-6)
            J = zeros(3, 3)
            J[:, 1] = (f(a, p + h, e, x) .- f(a, p - h, e, x)) ./ (2h)
            J[:, 2] = (f(a, p, e + h, x) .- f(a, p, e - h, x)) ./ (2h)
            J[:, 3] = (f(a, p, e, x + h) .- f(a, p, e, x - h)) ./ (2h)
            J
        end
        for (a, p, e, x) in ((0.95, 6.0, 0.5, 0.5), (0.9, 10.0, 0.2, 0.7))
            JO = fdjac(freqvec, a, p, e, x)
            JC = fdjac(constvec, a, p, e, x)
            En, Lz, _ = constants_of_motion(a, p, e, x)
            JK = JC[3, :] .+ 2 * (Lz - a * En) * (JC[2, :] .- a * JC[1, :])       # Q → K
            JCK = vcat(transpose(JC[1, :]), transpose(JC[2, :]), transpose(JK))
            @test isapprox(spin_frequency_jacobian(a, p, e, x), JO / JCK; rtol = 1e-6, atol = 1e-8)
            @test isapprox(spin_orbital_params_jacobian(a, p, e, x), inv(JCK); rtol = 1e-6, atol = 1e-8)
            Om1 = collect(spin_linear_frequencies(a, p, e, x))
            @test isapprox(spin_frequencies_constants_gauge(a, p, e, x),
                           Om1 .- (JO / JCK) * collect(spin_linear_constants(a, p, e, x));
                           rtol = 1e-6, atol = 1e-10)
            @test isapprox(spin_frequencies_turning_point_gauge(a, p, e, x),
                           Om1 .- JO * collect(spin_linear_orbital_params(a, p, e, x));
                           rtol = 1e-6, atol = 1e-10)
        end
    end

    @testset "near the separatrix CTilde stays finite while C / J / pex grow" begin
        # the message of arXiv:2603.13482 fig:delta_Omega, at a=0.95, e=0.5, x=0.5
        a, e, x = 0.95, 0.5, 0.5
        ps = separatrix(a, e, x)
        far  = ps + 10.0
        near = ps + 1e-3
        ctf = abs.(collect(spin_linear_frequencies(a, far, e, x)))
        ctn = abs.(collect(spin_linear_frequencies(a, near, e, x)))
        @test all(isfinite, ctn) && maximum(ctn) < 10 * maximum(ctf) + 1     # CTilde stays bounded
        for f in (spin_frequencies_constants_gauge, spin_frequencies_action_gauge)
            @test maximum(abs.(f(a, near, e, x))) > 1e3 * maximum(abs.(f(a, far, e, x)))
        end
        # the turning-point gauge diverges only mildly — much less than C / J
        @test maximum(abs.(spin_frequencies_turning_point_gauge(a, near, e, x))) <
              maximum(abs.(spin_frequencies_constants_gauge(a, near, e, x)))
    end

    @testset "δJ_r|_{p,e,x} is O(e²) near circularity (arXiv:2506.20726 fig:deltaJr)" begin
        # the published message: the radial-action shift at fixed (p,e,x) vanishes like e²
        for (a, p, x) in ((0.99, 8.0, 0.8), (0.99, 3.0, 0.8), (0.1, 8.0, 0.8))
            y1 = abs(spin_actions_pex(a, p, 1e-3, x)[1])
            y2 = abs(spin_actions_pex(a, p, 1e-2, x)[1])
            slope = log(y2 / y1) / log(10.0)
            @test isapprox(slope, 2.0; atol = 0.05)
        end
    end
end
