# Characterization of the AAK kludge harmonic mode-sum.
# Grounding (NO live heavy dep — KerrBackwards/FEW forbidden as test deps; SpecialFunctions forbidden
# too, design §5.5 purity: Kernels deps == {NumericsCore, KernelAbstractions}):
#   - Peters–Mathews Bessel coefficients enter ONLY through the injected `bessel_coeffs` provider
#     (canonical implementation: EMRIWaveforms.AAK.ak_bessel_coeffs, grounded THERE vs
#     SpecialFunctions.besselj). Here the provider is the BAKED fixture
#     test/fixtures/aak_bessel_ref.jl — the pre-hoist committed Kernels.ak_bessel_coeffs values at
#     exactly the (n,e) pairs these tests drive (bit-exact keys; a KeyError = input drift).
#   - _aak_rot_coeff: HAND-COMPUTED from the n/L/S cross-product algebra:
#       * iota=0 (L∥S) ⇒ cosrot=1, sinrot=0 ⇒ (1,0,0,1).
#       * iota=π/4, θ_S=0 (n=ẑ), θ_K=π/2,φ_K=0 (S=x̂), α=π/2 ⇒ L=(√2/2,-√2/2,0),
#         nxL=(√2/2,√2/2,0), nxS=(0,1,0), cosrot=√2/2, sinrot=-√2 ⇒ rot0=0, rot1=-1 ⇒ (0,-1,1,0).
#     (Both confirmed against the legacy KerrBackwards._aak_rot_coeff oracle during fixture capture.)
#   - _aak_point / aak_summation: CAPTURED FIXTURE — one-shot run of legacy KerrBackwards._aak_point /
#     .aak_summation (heavy pkg, run only to capture, NEVER a live dep) with the exact inputs below; the
#     ComplexF64 outputs are baked verbatim. Capture script: scratchpad/cap.jl via
#     `julia --project=packages/KerrBackwards`.
#   - PARALLEL == SERIAL bit-identity (===) over dense nodes.

using Test
using Kernels
const KA_aak = Kernels   # internal (unexported) helpers accessed via the module

include(joinpath(@__DIR__, "fixtures", "aak_bessel_ref.jl"))   # const AAK_BESSEL_REF
# Baked provider: KeyError on any (n,e) the pre-hoist kernel did not request ⇒ drift guard.
fix_bessel(n, e) = AAK_BESSEL_REF[(Int(n), Float64(e))]

@testset "aak" begin
    @testset "_aak_rot_coeff: hand-computed cross-product algebra" begin
        r1 = KA_aak._aak_rot_coeff(0.0, 1.0, 0.3, 0.9, 1.1, 0.4)        # iota=0 ⇒ (1,0,0,1)
        @test all(isapprox.(r1, (1.0, 0.0, 0.0, 1.0); atol = 1e-12))
        r2 = KA_aak._aak_rot_coeff(pi/4, 0.0, 0.0, pi/2, 0.0, pi/2)     # ⇒ (0,-1,1,0)
        @test all(isapprox.(r2, (0.0, -1.0, 1.0, 0.0); atol = 1e-12))
    end

    @testset "_aak_point: captured KerrBackwards fixture" begin
        args = (0.3, 0.5, 0.2, 0.4, 1e-3, 1e-4, 2e-3, 0.6, 1000.0)
        kw = (M_phys = 1e6, S_phys = 0.5, mu = 10.0, qS = 0.7, phiS = 0.8,
              qK = 0.9, phiK = 1.1, dist = 1.0, nmodes = 3)
        pf = KA_aak._aak_point(fix_bessel, args...; kw..., mich = false)
        pt = KA_aak._aak_point(fix_bessel, args...; kw..., mich = true)
        @test pf ≈ (-3.896516254482147e-23 - 8.627431394501252e-23im) rtol = 1e-12
        @test pt ≈ (-3.6781674568316085e-24 + 4.487738531574858e-24im) rtol = 1e-12
    end

    @testset "aak_summation: captured fixture + parallel === serial" begin
        N = 5
        t = collect(0.0:10.0:10.0 * (N - 1))
        e = [0.3 + 0.001 * i for i in 1:N]
        Phivec = [0.1 * i for i in 1:N]
        gimvec = [0.2 + 0.01 * i for i in 1:N]
        alpvec = [0.4 + 0.02 * i for i in 1:N]
        nuvec = [1e-3 + 1e-5 * i for i in 1:N]
        gimdotvec = [1e-4 * i for i in 1:N]
        OmegaPhi = [2e-3 + 1e-5 * i for i in 1:N]
        lam = [0.6 + 0.01 * i for i in 1:N]
        h = aak_summation(t, e, Phivec, gimvec, alpvec, nuvec, gimdotvec, OmegaPhi, lam,
                          1e6, 0.5, 10.0, 0.7, 0.8, 0.9, 1.1, 1.0, 3;
                          bessel_coeffs = fix_bessel, mich = false, dt = 10.0, out_len = 4)
        # captured legacy KerrBackwards.aak_summation output (verbatim):
        ref = ComplexF64[
            -9.094314566842473e-23 - 1.363760452319579e-23im,
            -8.218650921610496e-23 - 4.2340802446905075e-23im,
            -6.485260605810542e-23 - 6.682712817814172e-23im,
            -4.07364080939573e-23  - 8.460927210783241e-23im,
        ]
        @test h ≈ ref rtol = 1e-12

        ht = aak_summation(t, e, Phivec, gimvec, alpvec, nuvec, gimdotvec, OmegaPhi, lam,
                           1e6, 0.5, 10.0, 0.7, 0.8, 0.9, 1.1, 1.0, 3;
                           bessel_coeffs = fix_bessel, mich = false, dt = 10.0, out_len = 4,
                           threaded = true)
        @test all(ht[i] === h[i] for i in eachindex(h))

        # bessel_coeffs is a REQUIRED keyword (Kernels carries no Bessel machinery, §5.5 purity):
        @test_throws UndefKeywordError aak_summation(
            t, e, Phivec, gimvec, alpvec, nuvec, gimdotvec, OmegaPhi, lam,
            1e6, 0.5, 10.0, 0.7, 0.8, 0.9, 1.1, 1.0, 3; mich = false, dt = 10.0, out_len = 4)
    end
end
