# test_survey.jl — EGD::sw-survey-parallel: the threaded per-m / per-configuration torque survey of
# arXiv:2608.17003 Sec. III (src/SpiralWaves/survey.jl). Gates:
#   (a) threaded == serial BITWISE over m = 1…8 (per-task workspaces + fixed-order serial reduction);
#   (b) each surveyed mode reproduces the stand-alone solve → amf/torque_density/corotation_split path;
#   (c) the fixed-order sums equal the serial sums of the per-mode values, and feed `m_cut`;
#   (d) `survey_configs` over several discs; argument validation;
#   (e) the compute-usage self-check (mission policy: ≤ 24 CPUs per job, no GPU) is reported.
using Test
using AccretionDisc
using AccretionDisc.SpiralWaves

_bg_ref() = DiscBackground(; χ = 0.0, aEMRI = 15.0, h0 = 0.05, 𝔭 = 0.0, 𝔮 = 0.0)
const GK = (nL = 8, nR = 8, nC = 100, n_anc = 800)      # small but genuinely adaptive grids

@testset "SpiralWaves survey: threaded m-sweep and configuration sweep (Sec. III)" begin

    @testset "(e) compute-usage self-check" begin
        b = compute_budget()
        @info "compute budget: $(b.note)"
        @test b.used == min(Threads.nthreads(), 24)
        @test b.cap == 24
        @test compute_budget(nthreads = 64).ok == false
        @test compute_budget(nthreads = 64).used == 24
        @test compute_budget(nthreads = 4).ok == true
    end

    bg = _bg_ref()
    ms = 1:8
    s_thread = torque_survey(bg, ms; threads = true, grid_kwargs = GK, verbose = true)
    s_serial = torque_survey(bg, ms; threads = false, grid_kwargs = GK)

    @testset "(a) threaded == serial, BITWISE" begin
        @test s_thread.ms == s_serial.ms == collect(ms)
        for (a, b) in zip(s_thread.modes, s_serial.modes)
            @test a.m == b.m
            @test a.T === b.T                       # bitwise
            @test a.T_L === b.T_L
            @test a.T_c === b.T_c
            @test a.F_in === b.F_in && a.F_out === b.F_out
            @test a.nodes == b.nodes
        end
        @test s_thread.T_total === s_serial.T_total
        @test s_thread.T_L_total === s_serial.T_L_total
        @info "survey m = $(collect(ms)) on $(Threads.nthreads()) threads: $(round(s_thread.seconds; digits = 2)) s threaded vs $(round(s_serial.seconds; digits = 2)) s serial; T = $(s_thread.T_total), T_L = $(s_thread.T_L_total); per-mode T_L = $(mode_lindblad_torques(s_thread))"
    end

    @testset "(b) a surveyed mode == the stand-alone pipeline" begin
        m = 4
        g = build_grid(bg, m; GK...)
        sol = solve_mode(bg, m; grid = g, threads = false)
        F = amf(bg, sol); D = torque_density(bg, sol)
        sp = corotation_split(bg, sol; flux = F, density = D)
        mt = s_serial.modes[findfirst(==(m), s_serial.ms)]
        @test mt.T === sp.T
        @test mt.T_L === sp.T_L
        @test mt.T_c === sp.T_c
        @test mt.F_in === F[1] && mt.F_out === F[end]
        @test mt.nodes == length(g.r)
        @test mt.residual === sol.residual
    end

    @testset "(c) fixed-order reduction, profiles and the m-cutoffs" begin
        @test s_serial.T_total == sum(mode_torques(s_serial))
        @test s_serial.T_L_total == sum(mode_lindblad_torques(s_serial))
        @test s_serial.profiles === nothing
        sp = torque_survey(bg, 1:3; threads = false, keep_profiles = true, grid_kwargs = GK)
        @test length(sp.profiles) == 3
        for (i, pr) in enumerate(sp.profiles)
            @test length(pr.r) == length(pr.F_J) == length(pr.density) == sp.modes[i].nodes
            @test pr.F_J[1] === sp.modes[i].F_in
        end
        TL = mode_lindblad_torques(s_serial)
        @test all(>(0), TL)                              # every mode carries AM outward
        @test 1 ≤ m_cut(TL, 0.95) ≤ length(TL)
        @test m_cut(TL, 0.45) ≤ m_cut(TL, 0.95)
    end

    @testset "(d) survey_configs and argument validation" begin
        cfgs = [(χ = 0.0, aEMRI = 15.0, h0 = 0.05), (χ = 0.9, aEMRI = 15.0, h0 = 0.05)]
        out = survey_configs(cfgs, 2:3; threads = false, grid_kwargs = GK)
        @test length(out) == 2
        @test all(s -> s.ms == [2, 3], out)
        @test out[1].T_L_total != out[2].T_L_total       # spin changes the torque
        out2 = survey_configs([bg], 2:3; threads = false, grid_kwargs = GK)
        @test out2[1].T_L_total === out[1].T_L_total     # DiscBackground and NamedTuple agree bitwise
        @test_throws ArgumentError torque_survey(bg, Int[]; grid_kwargs = GK)
        @test_throws ArgumentError torque_survey(bg, 0:2; grid_kwargs = GK)
    end
end
