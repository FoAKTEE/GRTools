# test_selfforce_levels.jl — the LEVEL-SEQUENCE seam between the m-mode fan-out (row 47) and
# the self-force extraction (row 44), landed by ESSF node `fig-convergence` (GLOBAL_DAG row 50).
#
# WHY IT EXISTS. `mode_self_force_history` / `mode_force_convergence` were written against an
# `AmrHistory`, which owns a live `ModeSolution`. The fan-out deliberately does NOT ship that
# object back (`Solve/parallel.jl`: 597 KB of process-local solver state per mode); a worker
# returns a `ModeResult` carrying the per-level `AmrLevel` diagnostics and `x_p` and nothing
# else. So the per-AMR-level `F_r^m` sequence — the ONE quantity `fig:ssf_convergence` plots —
# was unreachable from a fanned-out sweep. These methods close that seam by taking the pieces
# the payload actually carries:
#
#     mode_self_force(orbit, config, m, level::AmrLevel, x_p; phi_p)
#     mode_self_force_history(orbit, config, m, levels; x_p, phi_p)
#     mode_force_convergence(orbit, config, m, levels; x_p, component, phi_p)
#
# and the `AmrHistory` methods become one-line delegations to them, so the two routes are the
# SAME code and cannot drift. The gates below are exactly that statement:
#   (l1) the level route reproduces the `AmrHistory` route BITWISE, field by field;
#   (l2) the fan-out payload (`ModeResult.levels`, `.x_p`) reproduces it BITWISE — the seam
#        this node needs;
#   (l3) the angular chart is read from the CONFIG (`block_lattice_bounds`), not assumed, and
#        the `∂_z = 2z ∂_{z²}` factor is actually applied;
#   (l4) every refusal of the `AmrHistory` methods survives the delegation (short table, bad
#        component, non-penetrating chart, a config whose mode number is not `m`).
#
# The sweeps here are deliberately SHORT (`iterations = 4`, five levels — one more than
# `mode_force_convergence` needs): this file gates the seam, not the physics. The physics
# oracles of `F_r^m` are `test_selfforce.jl` (row 44) and the fan-out's own bitwise parity is
# `test_parallel.jl` (row 47).

using Test
using ScalarSelfForce

fcv_cfg(; spin = 0.5, r0 = 10.0, m = 1, L = 0, P = 4) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, P = P, L = L)

# a copy of one config group with one field replaced (test_selfforce.jl has the same helper
# under its own `sf_` prefix; the two are deliberately NOT shared through `Main`, which is how
# the package suite acquired its `cdiff` method-overwrite warning).
function fcv_setfield(cfg::RunConfig{T}, group::Symbol, field::Symbol, value) where {T}
    g = getfield(cfg, group)
    G = typeof(g)
    g2 = G((f === field ? value : getfield(g, f) for f in fieldnames(G))...)
    return RunConfig{T}((f === group ? g2 : getfield(cfg, f)
                         for f in fieldnames(RunConfig{T}))...)
end

const FCV_ORB = CircularOrbit(0.5, 10.0)
const FCV_M = 1
const FCV_CFG = fcv_cfg(; m = FCV_M)
const FCV_HIST = solve_mode_amr(FCV_ORB, FCV_CFG, FCV_M; iterations = 4)   # 5 levels

# every field of `mode_self_force`'s result, compared bitwise (`===` on the reals, `==` on the
# element list) — a tolerance here would hide exactly the drift the delegation must prevent.
function fcv_same(a, b)
    keys(a) == keys(b) || return false
    for k in keys(a)
        va, vb = getproperty(a, k), getproperty(b, k)
        if va isa AbstractVector
            va == vb || return false
        else
            isequal(va, vb) || return false
        end
    end
    return true
end

@testset "(l1) the level route reproduces the AmrHistory route bitwise" begin
    @test length(FCV_HIST.levels) == 5
    ref = mode_self_force_history(FCV_HIST)
    got = mode_self_force_history(FCV_ORB, FCV_CFG, FCV_M, FCV_HIST.levels; x_p = FCV_HIST.x_p)
    @test length(got) == length(ref)
    @test all(fcv_same(r, g) for (r, g) in zip(ref, got))

    # the single-level method too, on every level
    for k in eachindex(FCV_HIST.levels)
        @test fcv_same(mode_self_force(FCV_HIST; level = k),
                       mode_self_force(FCV_ORB, FCV_CFG, FCV_M, FCV_HIST.levels[k],
                                       FCV_HIST.x_p))
    end

    cref = mode_force_convergence(FCV_HIST)
    cgot = mode_force_convergence(FCV_ORB, FCV_CFG, FCV_M, FCV_HIST.levels; x_p = FCV_HIST.x_p)
    @test keys(cref) == keys(cgot)
    for k in keys(cref)
        @test isequal(getproperty(cref, k), getproperty(cgot, k))
    end
    # and the table is the one fig:ssf_convergence plots: ndofs is the abscissa, values the
    # per-level F_r^m, diffs their successive differences, resolution_error the last of them.
    @test cgot.ndofs == [l.mesh.ndofs for l in FCV_HIST.levels]
    @test cgot.values == [f.Fr for f in mode_self_force_history(FCV_HIST)]
    @test cgot.resolution_error == cgot.diffs[end]

    # a non-default component travels through the delegation as well
    ct = mode_force_convergence(FCV_HIST; component = :Ft)
    cg = mode_force_convergence(FCV_ORB, FCV_CFG, FCV_M, FCV_HIST.levels;
                                x_p = FCV_HIST.x_p, component = :Ft)
    @test ct.values == cg.values
end

@testset "(l2) the fan-out payload reproduces the same table bitwise" begin
    # ONE process, the sealed serial driver (the fan-out's own bitwise oracle) — what comes
    # back is exactly what a worker ships: AmrLevels + x_p, no ModeSolution.
    set = solve_modes_serial(FCV_ORB, fcv_cfg(; m = 0), [FCV_M]; iterations = 4)
    @test isempty(set.failures)
    res = set[FCV_M]
    @test res isa ModeResult
    @test res.m == FCV_M
    @test res.x_p == FCV_HIST.x_p
    @test length(res.levels) == 5

    mcfg = mode_config(fcv_cfg(; m = 0), FCV_M)     # what the worker itself solved with
    cgot = mode_force_convergence(FCV_ORB, mcfg, FCV_M, res.levels; x_p = res.x_p)
    cref = mode_force_convergence(FCV_HIST)
    @test cgot.ndofs == cref.ndofs
    @test cgot.values == cref.values
    @test cgot.diffs == cref.diffs
    @test cgot.resolution_error == cref.resolution_error
end

@testset "(l3) the angular chart comes from the config, not from an assumption" begin
    lvl = FCV_HIST.levels[end]
    # production chart is z² with the puncture ON the equator, so ∂_z Ψ = 2z ∂_{z²}Ψ = 0 and
    # F_θ is a bitwise zero — the sealed statement of row 44.
    @test FCV_HIST.solution.dom.chart === :z2
    @test FCV_HIST.x_p[2] == 0.0
    @test mode_self_force(FCV_ORB, FCV_CFG, FCV_M, lvl, FCV_HIST.x_p).Ftheta == 0
    # move the evaluation point off the equator IN THE SAME CHART: the 2√(z²) factor must
    # appear, so F_θ becomes the raw-ingredient answer at dpsi_z = 2√(0.25)·dpsi_w.
    # `lvl.psi/dpsi_r/dpsi_w` are the MEAN over the containing elements (row 43 seq 6), so the
    # raw-ingredient comparison asks the level route for the SAME reduction
    off = mode_self_force(FCV_ORB, FCV_CFG, FCV_M, lvl, (FCV_HIST.x_p[1], 0.25);
                          reduction = :mean)
    raw = mode_self_force(FCV_ORB, FCV_M, lvl.psi, lvl.dpsi_r, 2 * sqrt(0.25) * lvl.dpsi_w)
    @test off.Ftheta == raw.Ftheta
    @test off.Ftheta != 0
end

@testset "(l4) the refusals survive the delegation" begin
    short = FCV_HIST.levels[1:3]
    @test_throws ArgumentError mode_force_convergence(FCV_ORB, FCV_CFG, FCV_M, short;
                                                      x_p = FCV_HIST.x_p)
    @test length(mode_self_force_history(FCV_ORB, FCV_CFG, FCV_M, short;
                                         x_p = FCV_HIST.x_p)) == 3
    @test_throws ArgumentError mode_force_convergence(FCV_ORB, FCV_CFG, FCV_M, FCV_HIST.levels;
                                                      x_p = FCV_HIST.x_p, component = :nope)
    # the config must be the one that solved this mode: a mode number mismatch is a caller
    # error the level route CAN see (the AmrLevels carry no mode number of their own).
    @test_throws ArgumentError mode_force_convergence(FCV_ORB, FCV_CFG, FCV_M + 1,
                                                      FCV_HIST.levels; x_p = FCV_HIST.x_p)
    # and the radial-chart contract of row 44 is checked here too
    bad = fcv_setfield(FCV_CFG, :orbit, :penetrating_horizon, false)
    @test !bad.orbit.penetrating_horizon
    @test_throws ArgumentError mode_self_force(FCV_ORB, bad, FCV_M, FCV_HIST.levels[end],
                                               FCV_HIST.x_p)
end
