# Grounding for Lorenz.SchwarzschildRadial — the a=0 full radial reconstruction (legacy radial.jl).
#
# (1) a=0 ANALYTIC self-consistency gate `radial_gates(r0)` (no oracle): the reconstructed P_{-2},P_{-1},h
#     reproduce the closure_a0 jumps (<1e-4), h satisfies its s=0 ODE (<1e-8), κ reproduces its jumps (<1e-4)
#     and satisfies the inhomogeneous ODE L[κ]=½r²h (<1e-6). The App-H a=0 jumps enter through closure_a0.
# (2) a→0-LIMIT cross-check vs module 1 (KerrRadial reconstruct_kerr at a=0): the a=0 Schwarzschild field
#     val(reconstruct(r0,l)) agrees with the spin-general valK(reconstruct_kerr(0,r0,l)) — two independent
#     reconstructions of the SAME a=0 metric, grounded at the closure G-K4 regression level (<1e-4). The
#     P_{-2} jumps are computed from the SAME jumps_pm/orbit_background ⇒ agree to machine precision.

using Test
import MetricReconstruction.Lorenz as L

@testset "Lorenz.SchwarzschildRadial (a=0 radial reconstruction)" begin

    @testset "a=0 self-consistency: radial_gates(6) ls=2:4" begin
        maxjump, maxres_h, maxjump_k, maxres_k = L.radial_gates(6.0; ls=2:4, m=2)
        @test maxjump   < 1e-4
        @test maxres_h  < 1e-8
        @test maxjump_k < 1e-4
        @test maxres_k  < 1e-6
    end

    # (3) S3-b Float64 REGRESSION (pre-{T} fixture): literals captured ONE-SHOT from the committed
    #     Float64-capped code @ d420e97 (capture: jobs/160c1aaf/tmp/capture_ref.jl — jump_data_radial(6,2,2)
    #     + reconstruct(6,2) val/der/der2 at r=4.5,7.5 + κval/κder). Guards that removing the Float64 cap
    #     leaves the Float64 path numerically unchanged. Hom/jump chain deterministic ⇒ 1e-12; the κ VoP
    #     values go through adaptive quadgk ⇒ 1e-10.
    @testset "Float64 regression vs pre-{T} fixture (d817c05-era code)" begin
        relc(x, ref) = abs(x - ref)/abs(ref)
        jd = L.jump_data_radial(6.0, 2, 2)
        @test relc(jd.ω,  1.36082763487954339e-01) < 1e-12
        @test relc(jd.Δ0, 2.40000000000000000e+01) < 1e-12
        @test relc(jd.Pm2[1], ComplexF64(-9.34165202732989974e-01, -3.20351931486590189e+01)) < 1e-12
        @test relc(jd.Pm2[2], ComplexF64(-1.40124780400316169e+01, -1.01063406956798154e+01)) < 1e-12
        @test relc(jd.Pm1[1], ComplexF64(-2.69039578395170167e+02,  1.09834947940785170e+02)) < 1e-12
        @test relc(jd.Pm1[2], ComplexF64( 6.72598945920681643e+01, -3.66116493181650782e+01)) < 1e-12
        @test relc(jd.h[2],   ComplexF64(-5.72057020539855743e-01,  0.0)) < 1e-12
        @test relc(jd.κ[1],   ComplexF64(-3.08910791073022324e+01,  5.83819215194125718e-10)) < 1e-12
        @test relc(jd.κ[2],   ComplexF64( 3.21782074060650203e+01, -5.32178745515921037e-11)) < 1e-12
        rec = L.reconstruct(6.0, 2; m=2)
        # Float64 in ⇒ Float64 fields out (the {T} parametrization must NOT change the Float64 contract).
        @test rec.Pm2 isa L.HomField{Float64}
        @test rec.κ   isa L.KappaField{Float64}
        @test rec.Pm2.Aup isa ComplexF64 && rec.κ.W̃ isa ComplexF64
        @test relc(rec.Pm2.Aup, ComplexF64(3.50533463285308614e-02, -1.17960805942203401e-02)) < 1e-12
        @test relc(rec.Pm2.Ain, ComplexF64(3.15532588227842789e-02,  2.11074367615612865e-02)) < 1e-12
        fix = Dict(  # (field, r) => (val, der, der2) — pre-change literals
            (:Pm2,4.5) => (ComplexF64(3.01078591617468616e+00, 3.04380230788600992e+00),
                           ComplexF64(3.46781427435639333e+00, 3.86559915297893752e+00),
                           ComplexF64(2.64996493977241299e+00, 3.69820180232934348e+00)),
            (:Pm2,7.5) => (ComplexF64(2.10156318636485206e+00, -1.36635545042305289e+01),
                           ComplexF64(-5.88271747837567283e+00, 2.65196708832127559e+00),
                           ComplexF64(-6.16122743522963812e-01, 1.51502582048645140e-01)),
            (:Pm1,4.5) => (ComplexF64(1.42908968435033898e+01, 5.61471648213265340e+00),
                           ComplexF64(1.23020544850018236e+01, 5.32600204901324226e+00),
                           ComplexF64(6.39760609441508343e+00, 3.59110994109836312e+00)),
            (:Pm1,7.5) => (ComplexF64(-1.36821781852024316e+02, 1.04518910920535419e+02),
                           ComplexF64(4.06283042968029804e+01, -9.73784110638545286e+00),
                           ComplexF64(-1.94215587052256211e+01, 6.06405179802333461e+00)),
            (:h,  4.5) => (ComplexF64(3.16591686526208449e-01, 1.98210857381253547e-03),
                           ComplexF64(1.74370580300438671e-01, 1.00691731615557001e-03),
                           ComplexF64(4.13561483243546149e-02, 3.11671728218908262e-04)),
            (:h,  7.5) => (ComplexF64(3.09990932202069858e-01, 6.09258585085384330e-03),
                           ComplexF64(-1.21643613060624600e-01, 1.66071988102026863e-03),
                           ComplexF64(7.27511402255616169e-02, 1.53016130580416881e-04)))
        for ((nm, r), (v0, v1, v2)) in fix
            f = getproperty(rec, nm)
            @test relc(L.val(f, r),  v0) < 1e-12
            @test relc(L.der(f, r),  v1) < 1e-12
            @test relc(L.der2(f, r), v2) < 1e-12
        end
        κfix = Dict(4.5 => (ComplexF64(-9.36319682136186771e+00, -2.72982699700711096e-01),
                            ComplexF64(-5.00165677614668702e+00, -1.46658977147316516e-01)),
                    7.5 => (ComplexF64(-2.50408548402809679e+01, -8.47360811634009670e-01),
                            ComplexF64(9.31165820047700343e+00, -2.27484173155455383e-01)))
        for (r, (kv, kd)) in κfix
            @test relc(L.κval(rec.κ, r), kv) < 1e-10
            @test relc(L.κder(rec.κ, r), kd) < 1e-10
        end
    end

    # (4) S3-b {T}-GENERIC surface: BigFloat inherits (§9 contract, GROTOOLS_DESIGN_CAPABILITY.md:518-525).
    #     Grounded THREE ways, never on the new code's own single output:
    #       (a) BigFloat(256) result == the PRE-change Float64 fixture at Float64 level (<1e-12);
    #       (b) 256↔512-bit SELF-CONVERGENCE far beyond Float64 (<1e-60; the unchanged upstream chain
    #           closure_a0/jumps_pm/TeukolskyRadial measured 1e-75..1e-78 at these args, conv_probe.jl);
    #       (c) hard type assertions (BigFloat in ⇒ Complex{BigFloat} out, no silent demotion).
    #     The hom surface is probed through the same branch/jump_amps/HomField chain reconstruct uses;
    #     the κ VoP FIELD (Abel-tail quadgk over ~7500M in s) is not rebuilt in BigFloat here — its
    #     quadrature tolerances (rtol=1e-9) bound it at Float64 level anyway; the κ JUMPS (the physics
    #     closure) ARE covered in full BigFloat below.
    @testset "{T}-generic: BigFloat inheritance + 256↔512 self-convergence" begin
        l = 2; m = 2
        probe = function (bits)
            setprecision(BigFloat, bits) do
                r0 = BigFloat(6)
                jd = L.jump_data_radial(r0, l, m)
                brm2 = L.branch(-2, l, m, jd.ω, jd.Λ)
                Pm2 = L.HomField{BigFloat}(r0, L.jump_amps(brm2, r0, jd.Pm2...)..., brm2)
                br0 = L.branch(0, l, m, jd.ω, jd.λ0)
                h = L.HomField{BigFloat}(r0, L.jump_amps(br0, r0, jd.h...)..., br0)
                v = (L.val(Pm2, BigFloat(75)/10), L.der(Pm2, BigFloat(75)/10), L.val(Pm2, BigFloat(45)/10),
                     L.val(h, BigFloat(75)/10),  L.der(h, BigFloat(45)/10))
                (jd=jd, Pm2=Pm2, v=v)
            end
        end
        p256 = probe(256); p512 = probe(512)
        # (c) types: BigFloat in ⇒ BigFloat through the whole surface
        @test p256.jd.ω isa BigFloat && p256.jd.Λ isa BigFloat
        @test p256.jd.Pm1[1] isa Complex{BigFloat} && p256.jd.κ[1] isa Complex{BigFloat}
        @test p256.Pm2 isa L.HomField{BigFloat}
        @test p256.v[1] isa Complex{BigFloat} && precision(real(p256.v[1])) == 256
        # (a) Float64-level agreement with the PRE-change fixture (same literals as testset (3))
        relc(x, ref) = Float64(abs(x - ref)/abs(ref))
        @test relc(p256.jd.Pm1[1], ComplexF64(-2.69039578395170167e+02, 1.09834947940785170e+02)) < 1e-12
        @test relc(p256.jd.κ[1],   ComplexF64(-3.08910791073022324e+01, 5.83819215194125718e-10)) < 1e-12
        @test relc(p256.v[1], ComplexF64(2.10156318636485206e+00, -1.36635545042305289e+01)) < 1e-12
        @test relc(p256.v[3], ComplexF64(3.01078591617468616e+00,  3.04380230788600992e+00)) < 1e-12
        @test relc(p256.v[4], ComplexF64(3.09990932202069858e-01,  6.09258585085384330e-03)) < 1e-12
        # (b) 256↔512 self-convergence: BigFloat genuinely inherited (way beyond Float64 eps)
        for k in eachindex(p256.v)
            @test relc(p256.v[k], p512.v[k]) < 1e-60
        end
        for fld in (:Pm2, :Pm1, :κ), j in 1:2
            @test relc(getproperty(p256.jd, fld)[j], getproperty(p512.jd, fld)[j]) < 1e-60
        end
        @test relc(p256.jd.h[2], p512.jd.h[2]) < 1e-60
    end

    @testset "a→0 limit: reconstruct == reconstruct_kerr(0) (G-K4-grounded)" begin
        r0 = 6.0; m = 2
        for l in 2:3
            rec  = L.reconstruct(r0, l; m=m)
            recK = L.reconstruct_kerr(0.0, r0, l; m=m)
            # P_{-2} jumps: identical jumps_pm/orbit_background ⇒ machine-precision match.
            jd  = L.jump_data_radial(r0, l, m)
            jdK = L.jump_data_kerr(0.0, r0, l, m)
            @test abs(jd.Pm2[1] - jdK.Pm2[1]) / (abs(jdK.Pm2[1]) + 1e-300) < 1e-8
            @test abs(jd.Pm2[2] - jdK.Pm2[2]) / (abs(jdK.Pm2[2]) + 1e-300) < 1e-8
            # homogeneous field VALUE away from orbit: two reconstructions of the same a=0 field (G-K4 level).
            for r in (r0+2.0, r0-1.5)
                vP2  = L.val(rec.Pm2, r);  vP2K  = L.valK(recK.Pm2, r)
                @test abs(vP2 - vP2K) / (abs(vP2K) + 1e-300) < 1e-4
            end
        end
    end
end
