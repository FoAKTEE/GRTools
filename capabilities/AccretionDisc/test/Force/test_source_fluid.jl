# Force.Source + Force.Fluid — effective force f_μ → 𝒬-source V^μ → 𝔉(r), fluid reconstruction.
# GROUNDING (Phase-7B: committed tests reach ONLY test/fixtures/):
#   • FIXTURES test/fixtures/h1-l{2..6}m2.h5 = the a=0 Lorenz BLS corpus (iter94,
#     ref-code/h1Lorenz/run_r8/fields) BYTE-COPIED on the window [r0=8, 17.6] by the one-shot
#     capture scratchpad/capture_s55_accretiondisc_force.jl (2026-07-02); the capture VERIFIED all
#     per-point pipeline quantities on the window are bit-identical to the full corpus.
#   • ORACLE literals = the CURRENT legacy packages/DiscEMRIForce output captured in the same run
#     (force_mode/Vvec/source_F at r≈11, a=0, m=2, r_p=8, lmax=6). NOTE: the legacy repo's OWN
#     cov_source_fluid literals are STALE (written 2026-06-02, before the Jun-15 source.jl
#     iter287/wing changes) — the capture re-ran the legacy oracle as it stands.
#   • tolerances: f_r rel<1e-8 and V_r rel<1e-6 are the legacy suite's frozen cov tolerances
#     (measured port deltas 4.7e-10 / 3.2e-9 — the analytic-Spacetime.christoffel-vs-legacy-FD
#     swap); source_F rel<1e-5 (new literal; measured 1.6e-7, ~60× margin); δp=(ε+p)𝔥 exact
#     identity <1e-12 (legacy); midplane θ-symmetry <1e-6 (legacy).
using Test
using AccretionDisc
using AccretionDisc.Force

const FIXDIR_SF = joinpath(@__DIR__, "..", "fixtures")

@testset "Force.Source/Fluid (a=0 corpus fixtures)" begin
    d = disc_a0(); m = 2; r_p = 8.0
    ω0 = m * disc_OmegaK(d, r_p)
    modes = load_modes(2:6, m; dir = FIXDIR_SF)
    g = modes[2].r_right
    k = argmin(abs.(g .- 11.0))

    @testset "G-D1 loader anchor (corpus fixture)" begin
        v = modes[2].val_right[1][1]                       # h̄¹(r0), l=2, m=2
        @test abs(v - (3.16711 - 0.23710im)) < 1e-4        # the legacy G-D1 anchor tolerance
        @test v == 3.1671076118154446 - 0.2371030070352107im  # byte-copied corpus value (capture)
        @test modes[2].r0 ≈ 8.0 atol = 1e-10
        @test modes[2].parity === :even && modes[3].parity === :odd
    end

    @testset "force_mode (single l=2 mode) vs captured legacy oracle" begin
        f = force_mode(d, modes[2], k, 2, 2, ω0, 0.0)
        @test length(f) == 4 && all(isfinite, f)
        fnorm = sqrt(sum(abs2, f))
        @test fnorm > 0
        @test abs(f[3]) / fnorm < 1e-6                     # midplane reflection symmetry ⇒ f_θ ≈ 0
        f_r_ref = -5.4945916149896774e-5 - 2.42039678673161e-6im   # captured 2026-07-02
        @test abs(f[2] - f_r_ref) / abs(f_r_ref) < 1e-8    # legacy frozen cov tolerance
    end

    @testset "Vvec (𝒬-projected source, l=2..6) vs captured legacy oracle" begin
        V = Vvec(d, modes, k, m, r_p, ω0, 0.0; lmax = 6)
        Vnorm = sqrt(sum(abs2, V))
        @test Vnorm > 0
        @test abs(V[3]) / Vnorm < 1e-6                     # midplane V_θ ≈ 0
        V_r_ref = -0.0086494387869375782 + 0.13142423478335336im   # captured 2026-07-02
        @test abs(V[2] - V_r_ref) / abs(V_r_ref) < 1e-6    # legacy frozen cov tolerance
    end

    @testset "source_F vs captured legacy oracle" begin
        sF = source_F(d, modes, k, m, r_p, ω0; lmax = 6)
        sF_ref = 1.694190742488599e-6 - 4.2784496376152511e-5im    # captured 2026-07-02
        @test abs(sF - sF_ref) / abs(sF_ref) < 1e-5        # measured 1.6e-7 (christoffel swap), ~60× margin
    end

    @testset "reconstruct / FluidVars / Jadv" begin
        hval = 1.0e-3 + 2.0e-3im; hp = -3.0e-4 + 1.0e-4im
        R = reconstruct(d, modes, k, m, r_p, ω0, hval, hp; lmax = 6)
        @test R isa FluidVars
        r = g[k]
        εp = eden0(d, r) * (1 + csound(d, r)^2)
        @test abs(R.δp - εp * hval) / abs(εp * hval) < 1e-12   # EXACT identity δp = (ε+p)𝔥 (legacy)
        @test R.δU_wave != R.δU                            # wave part drops the direct metric force
        @test all(isfinite, R.δU) && all(isfinite, R.δU_wave) && isfinite(R.δp) && isfinite(R.δε)
        Jf = Jadv(d, modes, k, m, r_p, ω0, hval, hp; lmax = 6, wave = false)
        Jw = Jadv(d, modes, k, m, r_p, ω0, hval, hp; lmax = 6, wave = true)
        @test isfinite(Jf) && isfinite(Jw) && Jf != Jw
    end

    @testset "disc(a) background constructor" begin
        @test disc(0.0) == disc_a0()                       # byte-identical fields (same T)
        @test disc_OmegaK(disc(0.0), 8.0) == disc_OmegaK(disc_a0(), 8.0)
        @test lindblad_locations(disc(0.0), 8.0, 2) == lindblad_locations(disc_a0(), 8.0, 2)
        da = disc(0.5)
        @test da.a == 0.5 && disc_OmegaK(da, 8.0) != disc_OmegaK(disc_a0(), 8.0)
        db = disc(big"0.5"; T = BigFloat)
        @test db isa RelativisticDisc{BigFloat}
    end
end
