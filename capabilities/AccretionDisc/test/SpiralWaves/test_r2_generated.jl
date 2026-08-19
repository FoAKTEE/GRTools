# test_r2_generated.jl — EGD::sw-mode-solver / EGD::sw-diagnostics: the ROUTE-R2 generated diagnostics
# of arXiv:2608.17003 (src/SpiralWaves/r2_generated.jl, emitted by the emri-disk generator
# codes/cas/R2_to_julia.wls through the SHARED tree walker codes/cas/cas_to_julia.wl — the same walker
# that emits master_generated.jl, so nothing is transcribed twice).
#
# Oracle: the 50-point R1-vs-R2 gate fixture (EGD::R1-vs-R2-gate, 30 digits) baked by
# fixtures/capture/bake_r2_diag_fixture.py into fixtures/r2_diag_fixture.jl. It carries, per point, the
# FULL route-R2 input set at the field point, so the generated bodies are checked DIRECTLY against
# Mathematica (independent of background.jl), and separately END-TO-END through DiscBackground /
# master_inputs / dv_dr (which additionally gates the new `dv_dr` and `rotation_derivs` of master.jl).
using Test
using AccretionDisc
using AccretionDisc.SpiralWaves
using AccretionDisc.SpiralWaves: MasterInputs, _r2_delta_u_phi, _r2_amf, _r2_torque_density,
                                 rotation_derivs, dv_dr, delta_u_phi
using Spacetime: kerr_equatorial_potentials

include(joinpath(@__DIR__, "fixtures", "r2_diag_fixture.jl"))

_rp(::Type{T}, s::AbstractString) where {T<:Real} = parse(T, s)
_cp(::Type{T}, p::Tuple{<:AbstractString,<:AbstractString}) where {T<:Real} =
    Complex{T}(parse(T, p[1]), parse(T, p[2]))
_reldev(a, b) = abs(a - b) / max(abs(b), eps(real(typeof(float(abs(b))))))

# MasterInputs assembled from the fixture's own background columns (no GrTools background involved).
# Only the fields the route-R2 bodies consume are filled from the fixture; the route-R1-only fields
# (W, Omega, Omega_r, mu_r, omt_r, nu_r and the h-derivatives) get their fixture values too so the
# struct is fully physical, but they are never read by _r2_*.
function _inputs_r2(::Type{T}, pt) where {T<:Real}
    g(k) = _rp(T, getproperty(pt, k))
    W = 1 / sqrt(1 - g(:v)^2)
    MasterInputs{T,Complex{T}}(pt.m, g(:nu), g(:nup), g(:psi), g(:psip), g(:mu), g(:mup), g(:omtp),
        W, g(:v), g(:Omega), g(:Omp), Complex{T}(g(:omd), zero(T)),
        g(:p), g(:pp), g(:eps), g(:epsp), g(:Gam1),
        g(:htt), g(:http), g(:htp), g(:htpp), g(:hpp), g(:hppp), g(:hrr), g(:hzz))
end

@testset "SpiralWaves r2_generated (arXiv:2608.17003 routes R2: δû_φ, AMF (72a), torque density (72b))" begin

    @testset "(a) generated bodies vs the Mathematica R2 fixture — Float64" begin
        worst_duphi = 0.0; worst_amf = 0.0; worst_tq = 0.0
        for pt in R2_DIAG_FIXTURE
            inp = _inputs_r2(Float64, pt)
            omt = _rp(Float64, pt.omt); vp = _rp(Float64, pt.vp)
            y0 = _cp(Float64, pt.sample_y0); y1 = _cp(Float64, pt.sample_y1)
            # δû_φ is AFFINE in (y0, y1): the fixture stores the three coefficients
            c0 = _cp(Float64, pt.deltaUphi_coeffs[1]); c1 = _cp(Float64, pt.deltaUphi_coeffs[2])
            ch = _cp(Float64, pt.deltaUphi_coeffs[3])
            got = _r2_delta_u_phi(inp, omt, vp, y0, y1)
            worst_duphi = max(worst_duphi, _reldev(got, c0 * y0 + c1 * y1 + ch))
            worst_amf = max(worst_amf, _reldev(_r2_amf(inp, omt, vp, y0, y1), _cp(Float64, pt.R2_AMF_at_sample)))
            worst_tq = max(worst_tq, _reldev(_r2_torque_density(inp, omt, vp, y0, y1),
                                             _cp(Float64, pt.R2_TorqueDensity_at_sample)))
        end
        @info "r2_generated fixture (Float64, 50 points): worst rel dev δû_φ $worst_duphi, AMF $worst_amf, torque density $worst_tq"
        @test worst_duphi < 1e-11
        @test worst_amf < 1e-11
        @test worst_tq < 1e-9
        pt = R2_DIAG_FIXTURE[1]
        @test (@inferred _r2_amf(_inputs_r2(Float64, pt), _rp(Float64, pt.omt), _rp(Float64, pt.vp),
                                 _cp(Float64, pt.sample_y0), _cp(Float64, pt.sample_y1))) isa ComplexF64
    end

    @testset "(a') generated bodies vs the fixture — BigFloat(256)" begin
        setprecision(BigFloat, 256) do
            worst = big(0.0)
            for pt in R2_DIAG_FIXTURE[1:8]
                inp = _inputs_r2(BigFloat, pt)
                omt = _rp(BigFloat, pt.omt); vp = _rp(BigFloat, pt.vp)
                y0 = _cp(BigFloat, pt.sample_y0); y1 = _cp(BigFloat, pt.sample_y1)
                c0 = _cp(BigFloat, pt.deltaUphi_coeffs[1]); c1 = _cp(BigFloat, pt.deltaUphi_coeffs[2])
                ch = _cp(BigFloat, pt.deltaUphi_coeffs[3])
                worst = max(worst, _reldev(_r2_delta_u_phi(inp, omt, vp, y0, y1), c0 * y0 + c1 * y1 + ch))
                worst = max(worst, _reldev(_r2_amf(inp, omt, vp, y0, y1), _cp(BigFloat, pt.R2_AMF_at_sample)))
                worst = max(worst, _reldev(_r2_torque_density(inp, omt, vp, y0, y1),
                                           _cp(BigFloat, pt.R2_TorqueDensity_at_sample)))
            end
            @info "r2_generated fixture (BigFloat 256, 8 points): worst rel dev = $(Float64(worst))"
            @test worst < 1e-28          # the fixture itself carries 30 significant digits
        end
    end

    @testset "(b) dv_dr / rotation_derivs vs the fixture background (end-to-end)" begin
        worst_v = 0.0; worst_dv = 0.0; worst_dΩ = 0.0
        for pt in R2_DIAG_FIXTURE
            bg = DiscBackground(; χ = _rp(Float64, pt.chi), aEMRI = _rp(Float64, pt.aEMRI),
                                h0 = _rp(Float64, pt.h0), 𝔭 = _rp(Float64, pt.pdisk), 𝔮 = _rp(Float64, pt.qdisk))
            r = _rp(Float64, pt.r)
            d = rotation_derivs(bg, r)
            worst_v = max(worst_v, _reldev(d.v, _rp(Float64, pt.v)))
            worst_dv = max(worst_dv, _reldev(d.dv, _rp(Float64, pt.vp)))
            worst_dΩ = max(worst_dΩ, _reldev(d.dΩ, _rp(Float64, pt.Omp)))
            @test dv_dr(bg, r) === d.dv
            @test dOmega_dr(bg, r) === d.dΩ
        end
        @info "rotation_derivs vs fixture (Float64, 50 points): worst rel dev v $worst_v, dv/dr $worst_dv, dΩ/dr $worst_dΩ"
        @test worst_v < 1e-13
        @test worst_dv < 1e-9
        @test worst_dΩ < 1e-13
    end

    @testset "(c) delta_u_phi end-to-end through DiscBackground + master_inputs" begin
        worst = 0.0
        for pt in R2_DIAG_FIXTURE
            bg = DiscBackground(; χ = _rp(Float64, pt.chi), aEMRI = _rp(Float64, pt.aEMRI),
                                h0 = _rp(Float64, pt.h0), 𝔭 = _rp(Float64, pt.pdisk), 𝔮 = _rp(Float64, pt.qdisk))
            r = _rp(Float64, pt.r)
            hm = (htt = _rp(Float64, pt.htt), htt_r = _rp(Float64, pt.http),
                  htphi = _rp(Float64, pt.htp), htphi_r = _rp(Float64, pt.htpp),
                  hphiphi = _rp(Float64, pt.hpp), hphiphi_r = _rp(Float64, pt.hppp),
                  hrr = _rp(Float64, pt.hrr), hzz = _rp(Float64, pt.hzz))
            y = (_cp(Float64, pt.sample_y0), _cp(Float64, pt.sample_y1))
            c0 = _cp(Float64, pt.deltaUphi_coeffs[1]); c1 = _cp(Float64, pt.deltaUphi_coeffs[2])
            ch = _cp(Float64, pt.deltaUphi_coeffs[3])
            got = delta_u_phi(bg, pt.m, r, y; hmodes = hm, landau = false)
            worst = max(worst, _reldev(got, c0 * y[1] + c1 * y[2] + ch))
        end
        @info "delta_u_phi end-to-end vs fixture (Float64, 50 points): worst rel dev = $worst"
        @test worst < 1e-8
    end

    @testset "(d) Newtonian limit: δû_φ → r·(50a) (ERRATUM E-M1)" begin
        # (50a) is the ORTHONORMAL component; the covariant one is r times it. In the Newtonian limit
        # (M/a• → 0) the R2 covariant δû_φ must approach r · delta_u_phi_newtonian with Φ̂ = −ĥ_tt/2.
        devs = Float64[]
        for a in (1e4, 1e5)
            bg = DiscBackground(; χ = 0.0, aEMRI = a, h0 = 0.05, 𝔭 = 0.0, 𝔮 = 0.0)
            m = 3
            r = 0.7a
            hm = singular_field_mode(bg, m, r)
            y = (0.3 - 0.1im, -0.2 + 0.4im)
            rel = delta_u_phi(bg, m, r, y; hmodes = hm, landau = false)
            Φ = -hm.htt / 2
            newt = r * delta_u_phi_newtonian(bg, m, r, y, Φ; landau = false)
            push!(devs, abs(rel - newt) / abs(newt))
        end
        @info "δû_φ Newtonian limit: rel dev a•=1e4 $(devs[1]), a•=1e5 $(devs[2]) (ratio $(devs[1]/devs[2]))"
        @test devs[1] < 1e-2
        @test devs[2] < devs[1]                    # → 0 as M/a• → 0
        @test 5 < devs[1] / devs[2] < 20           # first order in M/a•
    end
end
