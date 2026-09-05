# reproduce/arXiv-2606.04998/fig_modes.jl — ESSF::fig-modes.
#
# arXiv:2606.04998 Fig. 5 (`fig:modes`, main.tex:691-696): "m-mode solutions for a scalar charge
# orbiting a black hole at the ISCO for three different black hole spins. Each plane shows an
# m-mode of the scalar field, Psi_m, up to m = 20, with the first few modes highlighted in colour.
# Plotted is their complex magnitude. Within the regularized region, the regular field
# Psi^R_m = Psi_m - Psi^P_m is plotted. The computational domain extends over z^2 in [0, 1]
# (equatorial symmetry) and to R = 10^4 M towards the right, but the plot is clipped at r = 2 r0."
# and main.tex:718: "All m-modes up to m = 20 are shown on a log scale, demonstrating the decay in
# power over mode number. Note that regular zero crossings in the solution of the regularized
# field within the worldtube display as sharp drops on the chosen log scale, particularly for the
# m = 0 mode. Also visible are features in the solution within the worldtube that are resolved by
# our adaptive grid."
#
# WHAT IS PLOTTED IS NOT A MODE-SUM CURVE. Each of the 63 surfaces is log10 of the complex
# magnitude of ONE m-mode field over the whole (r, z^2) domain — the same `sol.u` that Fig. 1
# draws (row 42's contract: `u` holds Psi^R inside the worldtube and Psi outside it, which is
# verbatim what the caption asks for) — NOT |F_r^m| against m. The per-mode self force is
# registered all the same (§ GOLDEN below), because it is one call away from the same data and
# it is what rows 46/53 consume, but it is not the figure's ink.
#
# THIN driver: zero physics. `kerr_isco` -> `CircularOrbit` -> `crp_paper_config` ->
# `block_domain` -> `mode_pool`/`solve_modes_parallel` (rows 43/47) and then only the SEALED
# read-outs `puncture_chart_position`, `worldtube_blocks`, `inside_worldtube`,
# `element_matrices`, `tensor_lin`, `mode_self_force`. The driver arranges, gates and draws.
#
#   julia --project=reproduce reproduce/arXiv-2606.04998/fig_modes.jl --check
#   julia --project=reproduce reproduce/arXiv-2606.04998/fig_modes.jl --plot-only   # redraw only
#
# THE THREE CONFIGURATIONS. The caption says "at the ISCO for three different black hole
# spins" and the published panel titles say `r0 = 9, 6, 2`; each published abscissa ends exactly
# at twice its title's radius (18, 12, 4 — the caption's "clipped at r = 2 r0"), so the panels
# were run at 9, 6 and 2, and only the middle one IS its spin's ISCO:
#
#     a = -0.998   r_ISCO = 8.9943744548 M   (title r0 = 9 — the ISCO to 0.06 %)   PLOTTED: ISCO
#     a =  0       r_ISCO = 6            M   (title r0 = 6 — identical)            PLOTTED: ISCO
#     a = +0.998   r_ISCO = 1.2369706552 M   (title r0 = 2 — 1.6 x the ISCO)       PLOTTED: 2 M
#
# So panels 1-2 follow the caption and panel 3 follows the PUBLISHED TITLE, which is the
# configuration the paper actually drew. `kerr_isco` supplies the first two radii; the third is
# the paper's own printed number. Every panel title carries the radius actually solved.
#
# r0 = 2 M WAS NOT CONSTRUCTIBLE UNTIL ESSF::row-57 FIXED IT (47f0ec0 + 7082f9a + 12b43bd).
# `CircularOrbit(a, 2.0).Lz` was `NaN`: `_angular_momentum_circular`
# (capabilities/Spacetime/src/Geodesics/constants.jl:13) divides by
# `((p - 2) p - a^2 (x^2 - 1))`, which on the EQUATOR is exactly `p (p - 2)` and vanishes at
# `p = 2` with its numerator, and `_carter_circular` (constants.jl:17) inherits it as `0 * NaN`
# — so `Lz` AND `Q` were 0/0 while `E`, which carries no such factor, was already right;
# `zc2 = 4 (Lz^2 + ...)` (main.tex:453) took the NaN and the puncture died in `legendre_P3`. It
# only bit where `r = 2 M` carries a circular geodesic at all, `a > 1/sqrt(2) = 0.70710678`;
# below that the negative radicand makes the `DomainError` those spins raise physically correct
# (at a = 0, `p = 2` is the horizon, where the pre-fix code returned a silent `E = 0.0`). The
# fix is bitwise-neutral away from the endpoint; a conditioning band survives around it (about
# `log10|p - 2|^-1` digits) and is that row's envelope, not this figure's business.
#
# THE HORIZON-AXIS LAYER, AND THE RULE THAT DECIDES THE HOLD (ESSF::fix-horizon-axis-layer).
# WHAT REGULARITY MEANS HERE, stated correctly: `eq:mode_sum` factors `sin^|m| theta` OUT of the
# field, so `Psi_m` is the variable that is ANALYTIC at the axis and approaches a finite,
# generically NON-ZERO limit there (main.tex:227: each mode "approaches a non-zero constant
# amplitude ... towards the poles"). A healthy plane therefore has an axis:equator ratio of
# ORDER ONE and a fitted `log|Psi_m|` vs `log(1 - z^2)` slope near ZERO — it must NOT fall like
# `(1 - z^2)^(m/2)`, which is `Phi`'s m-mode, never `Psi_m`.
# WHAT GOES WRONG: at the horizon the operator degenerates and its second radial branch is
# `Psi_m ~ (r - r_+)^(i nu)` with `nu = 2 m |Omega - Omega_H| (r_+^2 + a^2)/(r_+ - r_-)` — unit
# modulus, infinitely oscillatory, and excluded ONLY by the polynomial trial space. When hp-AMR
# refines the horizon element away, that exclusion weakens and the branch appears as a layer on
# the innermost radial node line, growing toward the axis. It is present at EVERY spin and is a
# question of SCALE, so the driver does not key on the spin: it MEASURES.
# THE RULE IS THE PACKAGE'S, not this driver's: sweep the panel with the default
# (`horizon_hold = false`, reference parity), read `horizon_axis_profile(dom, orbit, mesh,
# u).layer_ratio` for EVERY mode, and hand both sweeps to the sealed `select_hold_sweep`
# (`Solve/horizon_axis.jl`, shared with the error-budget figure). It publishes the held sweep
# exactly when the default's MAXIMUM OVER m crosses LAYER_GATE = 1.05 and the hold strictly
# reduces it, refuses a hold that does not, and reports a hole rather than substituting
# silently. The maximum, not the top mode: `layer_ratio` is not monotone in m (MEASURED at
# a = +0.998, r0 = 2: it peaks at 3.889 at m = 12 and falls to 3.6e-3 by m = 20). The plotted run's own
# `layer_ratio` is then GATED (GATE 7) and both numbers, with the setting, are stated on the
# panel and registered (`pK_layer_ratio`, `pK_layer_ratio_default_max`, `pK_horizon_hold`).
#
# [HOLE] GMRES STAGNATES AT m = 20 FOR a = +0.998, r0 = 2 M. Measured in this node's evidence in
# BOTH hold settings: `converged = false` with a final residual of 2.3e-5 (default) / 7.6e-5
# (held) against the yaml's 1e-12 tolerance, stopping at 71 / 137 iterations, far short of the
# 500 cap; every other grid of both sweeps converged, and the per-grid staircase GROWS with
# depth (1, 7, 7, 7, 11, 21, 61, 235 at m = 19). The mode STAYS in the figure — drawn dashed and
# fainter, named in the panel note with the plotted run's own residual — and `pK_converged`
# carries the verdict per mode. Handed to ESSF::fix-horizon-axis-layer's owner.
#
# THE WORLDTUBE RADIUS IS A DOCUMENTED CHOICE, NOT A PAPER NUMBER. The yaml is templated in
# `worldtube_radius` (anc/ScalarSelfForce.yaml:12) and the paper publishes a value only for its
# header example, r0 = 10 M with r_in = 2 M: `worldtube_radius = 4.8`. That number cannot be
# reused here — the lattice needs `r0 - 1.5 w_r > r_+` for a non-empty v block (yaml:55-57), and
# 4.8 puts the v->t transition at -1.2 M for r0 = 6. The driver therefore scales the worldtube
# with the radial gap the orbit actually has,
#
#     w_r = 0.6 (r0 - r_+),
#
# which (i) reproduces the paper's own published value EXACTLY at its own point
# (0.6 (10 - 2) = 4.8), (ii) keeps the v block at a fixed 10 % of [r_+, r0] at every spin, and
# (iii) is the only one-parameter family in reach that stays legal at all three ISCOs. It is an
# ASSUMPTION where the paper prints nothing, and every registered number depends on it.
#
# EIGHT AMR GRIDS, the paper's production resolution: `Amr.Iterations: 8` (yaml:153) counts
# spectre SOLVES and its trigger skips the initial guess, so the sweep is 7 refinements — the
# driver DERIVES `cfg.amr.iterations - cfg.amr.trigger_offset` exactly as fig_linsolv.jl does.
#
# GOLDEN. Per panel: the configuration (spin, r0, r_+, w_r, the six block bounds, x_p, the
# final element/dof counts) plus the rule's verdict (`horizon_hold`, `layer_ratio_default_max`)
# and, per m — the plotted amplitude `max_abs_psi`, the particle diagnostics `abs_psi_xp` /
# `abs_dpsi_r_xp` at the two highest AMR grids (the raw ingredients of `resolution_error`, which
# is deliberately NOT re-derived here: its ONE definition lives in `Solve/amr.jl` and rows
# 44/50/53 share it), the extracted `Fr`, the final grid's GMRES iteration count, the sealed
# row-58 instrument on the plotted run (`layer_ratio`, `axis_equator_ratio`, `axis_exponent`,
# `branch_index`), and the solver's own verdict (`converged`, `final_residual`). `Fr` comes from the SEALED raw-ingredient method
# `mode_self_force(orbit, m, psi, dpsi_r)`, whose third argument `d_z Psi^R` is the method's own
# default zero; that is bitwise what the `AmrHistory` method passes at an equatorial puncture in
# the z^2 chart (`_sf_dpsi_z(:z2, 0, d_w) = 2 sqrt(0) d_w`), and the driver GATES the two
# preconditions (`dom.chart === :z2`, `x_p[2] == 0`) instead of assuming them.
#
# Style: `_style/manuals/arXiv-2606.04998.md` sections 0 + 5. Rendering limitations under
# "RENDERING" below.

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()

# `--plot-only` redraws from the surface CACHE this driver writes beside the golden on `--regen`
# (63 nodal element lattices — far too much ink for a 200 KB golden), with no solve and no
# comparison. Same-Julia-version scratch artifact, untracked (`expected/.gitignore`); its absence
# is a HOLE, never a silent fallback.
const SURFACE_CACHE = joinpath(@__DIR__, "expected", "fig_modes_surface.jls")
const PLOT_ONLY = "--plot-only" in opts.extra

PLOT_ONLY && !isfile(SURFACE_CACHE) && ReproHarness.hole(
    "--plot-only needs the surface cache " *
    "$(relpath(SURFACE_CACHE, ReproHarness.REPRO_ROOT)); run --regen first")

PLOT_ONLY || ReproHarness.require_tier(opts, :heavy)   # measured; see the README

using Printf, LinearAlgebra, Statistics, Serialization
using GrTools
using GrTools.ScalarSelfForce
using GrTools.EllipticCore
using GrTools.Spacetime: kerr_isco

# ───────────────────────────── parameters ────────────────────────────────────────────────────
# (spin, radius): `:isco` is Spacetime's Bardeen closed form, a number is the paper's own
# printed panel title. Panels 1-2 are their spins' ISCOs (the caption); panel 3 is the PUBLISHED
# configuration r0 = 2 M, constructible since the p = 2 removable-singularity fix.
const PANELS  = ((-0.998, :isco), (0.0, :isco), (0.998, 2.0))
const LAYER_GATE = 1.05                # the measured `layer_ratio` a plotted panel may not pass
const WR_FRAC = 0.6                    # w_r = WR_FRAC (r0 - r_+); see the header
const WDZ2    = 0.25                   # yaml:12
const PORD    = 4                      # yaml:12
const LREF    = 0                      # yaml:12
const ROUT    = 1.0e4                  # yaml:11
const MMS     = CRP_PAPER_MMODES        # 0:20 (yaml has one MModeNumber per run)
const NWORK   = 5                      # <= 24 CPUs with 1 thread + 1 BLAS thread per worker

for f in (:solve_modes_parallel, :mode_pool, :mode_self_force, :puncture_chart_position)
    isdefined(ScalarSelfForce, f) || ReproHarness.hole(
        "ScalarSelfForce.$(f) is not available: Fig. 5 needs the m-mode fan-out " *
        "(ESSF::mmode-parallel-driver) over the hp-AMR driver (ESSF::amr-driver).")
end

panels = map(PANELS) do (a, spec)
    r0  = spec === :isco ? kerr_isco(a) : float(spec)
    orb = CircularOrbit(a, r0)
    wr  = WR_FRAC * (r0 - orb.rp)
    cfg = crp_paper_config(; spin = a, orbital_radius = r0, m_mode = first(MMS),
                           inner_radius = orb.rp, outer_radius = ROUT, worldtube_radius = wr,
                           worldtube_dz2 = WDZ2, P = PORD, L = LREF)
    dom = block_domain(orb, cfg)
    (a = a, r0 = r0, spec = spec, orb = orb, wr = wr, cfg = cfg, dom = dom,
     xp = puncture_chart_position(dom, orb), rmax = 2 * r0)
end
const ITER = panels[1].cfg.amr.iterations - panels[1].cfg.amr.trigger_offset   # 8 - 1 = 7

println("\n[fig_modes] three configurations, m = ", first(MMS), ":", last(MMS),
        ", ", ITER + 1, " AMR grids, P = ", PORD, ", L = ", LREF)
for p in panels
    @printf("[fig_modes] a = %+.4f  r₀ = %.10f (%s)  r₊ = %.10f  w_r = %.10f  x_p = %s\n",
            p.a, p.r0, p.spec === :isco ? "the ISCO" : "the published title",
            p.orb.rp, p.wr, string(p.xp))
    println("[fig_modes]    block bounds ", p.dom.bounds_r, "  clip at 2r₀ = ", p.rmax)
end
flush(stdout)

gatefails = String[]

# ───────────────────────────── the sweeps ────────────────────────────────────────────────────
# ONE pool for the whole figure (the row-47 contract: a fresh worker pays its entire JIT before
# its first solve, so a pool per panel would pay it three times).
struct Patch                      # one element's nodal lattice of ONE mode, as ParaView draws it
    r::Vector{Float64}            # the element's own collocation radii
    w::Vector{Float64}            # ... and z² nodes
    v::Matrix{Float64}            # |Ψ_m| (|Ψ^R_m| inside the worldtube) at those nodes
    tube::Bool                    # is this element inside the regularized region?
end

"the plotted quantity of row 58's rule: `layer_ratio` for every mode of a fanned-out sweep"
layer_ratios(p, s) = [haskey(s, m) ?
                      horizon_axis_profile(p.dom, p.orb, s[m].levels[end].mesh, s[m].u).layer_ratio :
                      NaN for m in MMS]

surfaces = Vector{Vector{Vector{Patch}}}()      # [panel][mode][element]
tables   = Vector{NamedTuple}()

if !PLOT_ONLY
    # a STABLE path (mktempdir() deletes itself at exit, taking the monitor log with it): the
    # compute-usage record has to outlive the run to be archivable evidence.
    monitor_log = joinpath(tempdir(), "fig_modes_load_monitor_$(getpid()).log")
    # The workers run `using ScalarSelfForce`, which the `reproduce` environment cannot serve:
    # it depends on the GrTools umbrella, not on the capability package. So the pool is launched
    # in the PACKAGE's own project, located from the loaded module rather than typed in. Both
    # environments resolve ScalarSelfForce to this same working tree (`[sources]` path deps), so
    # master and workers run identical code; only the master needs Plots.
    ssf_project = normpath(joinpath(dirname(pathof(ScalarSelfForce)), "..", "Project.toml"))
    isfile(ssf_project) || ReproHarness.hole(
        "the ScalarSelfForce project file is not where `pathof` says it is: $(ssf_project)")
    println("[fig_modes] worker project: ", ssf_project); flush(stdout)
    t_all = time()
    sweeps = mode_pool(NWORK; threads_per_worker = 1, blas_threads = 1,
                       project = ssf_project) do pool
        map(panels) do p
            # THE SWITCH RULE, MEASURED — never assumed, never keyed on the spin. Sweep with the
            # package DEFAULT (`horizon_hold = false`, reference parity), read the sealed row-58
            # instrument on every fanned-out mode, and re-sweep with the hold only if the
            # near-horizon strip carries more than LAYER_GATE of the rest of the plane anywhere.
            # THE MAXIMUM OVER m IS THE RULE'S INPUT, not the top mode: `layer_ratio` is NOT
            # monotone in m — MEASURED at a = +0.998, r0 = 2 it rises from m = 6, peaks at 3.889
            # at m = 12 and falls back to 3.6e-3 by m = 20, so a top-mode probe reads the panel
            # as healthy when it is not.
            t = time()
            sdef = solve_modes_parallel(p.orb, p.cfg, MMS; pool = pool, iterations = ITER,
                                        log = monitor_log)
            laydef = layer_ratios(p, sdef)
            @printf("[fig_modes] a = %+.4f: %d/%d modes in %.1f s (%d workers, default); layer_ratio max %.4g at m = %d\n",
                    p.a, length(sdef.results), length(MMS), time() - t, sdef.run.nworkers,
                    maximum(laydef), argmax(laydef) - 1)
            flush(stdout)
            # The DECISION is the sealed `select_hold_sweep` (row 58 / Solve/horizon_axis.jl),
            # not arithmetic this driver owns: called once against the default alone — which
            # publishes it when the ratio is inside the gate and otherwise reports the hole —
            # and, when that says a held sweep is needed, once more with both, so the helper
            # also refuses a hold that fails to reduce the ratio.
            sel = select_hold_sweep(sdef => laydef; gate = LAYER_GATE)
            if sel.hole
                t = time()
                sheld = solve_modes_parallel(p.orb, p.cfg, MMS; pool = pool, iterations = ITER,
                                             solve_kwargs = (; horizon_hold = true),
                                             log = monitor_log)
                @printf("[fig_modes] a = %+.4f: re-swept %d/%d modes in %.1f s with horizon_hold = true\n",
                        p.a, length(sheld.results), length(MMS), time() - t)
                sel = select_hold_sweep(sdef => laydef, sheld => layer_ratios(p, sheld);
                                        gate = LAYER_GATE)
            end
            println("[fig_modes] a = ", p.a, ": ", sel.reason,
                    " ⇒ horizon_hold = ", sel.hold)
            sel.hole && push!(gatefails, "GATE 7: a = $(p.a) — $(sel.reason)")
            flush(stdout)
            (set = sel.sweep, hold = sel.hold, probe_ratio = sel.default_ratio)
        end
    end
    sets = [q.set for q in sweeps]
    @printf("[fig_modes] all three panels: %.1f s; load monitor %s\n", time() - t_all, monitor_log)
    for (p, s) in zip(panels, sets)      # the compute-usage record, per sweep (types are the
        r = s.run                        # row-47 `ModeRunRecord`'s own — printed, never coerced)
        println("[fig_modes] a = ", p.a, ": ", r.nworkers, " workers x ", r.cpus_per_worker,
                " CPU, budget ", r.budget, ", guard ", r.load_guard,
                ", planned ", round(r.planned_share; digits = 4),
                "\n[fig_modes]    load before ", r.guard_before,
                "\n[fig_modes]    load after  ", r.guard_after,
                "\n[fig_modes]    load max    ", r.max_load,
                "\n[fig_modes]    startup ", round(r.startup_time; digits = 1), " s, wall ",
                round(r.wall_time; digits = 1), " s, samples ", r.samples)
    end
    flush(stdout)

    for (p, s) in zip(panels, sets)
        isempty(s.failures) || push!(gatefails,
            "GATE 1: a = $(p.a) lost modes $( [f.m for f in s.failures] ) — $(s.failures[1].message)")
        length(s.results) == length(MMS) || push!(gatefails,
            "GATE 1: a = $(p.a) returned $(length(s.results)) of $(length(MMS)) modes")
        # preconditions of the raw-ingredient extraction (never assumed)
        p.dom.chart === :z2 || push!(gatefails,
            "GATE 5: a = $(p.a) has chart $(p.dom.chart), not :z2 — ∂_zΨ^ℛ is then not the " *
            "raw-ingredient method's default zero")
        iszero(p.xp[2]) || push!(gatefails,
            "GATE 5: a = $(p.a) has x_p[2] = $(p.xp[2]) ≠ 0 — the puncture is off the equator")
    end

    for (p, s) in zip(panels, sets)
        psurf = Vector{Vector{Patch}}()
        maxpsi = Float64[]; fr = Float64[]
        psi_xp = Float64[]; dpsi_xp = Float64[]
        psi_xp_prev = Float64[]; dpsi_xp_prev = Float64[]
        ndofs = Int[]; ndofs0 = Int[]; nelem = Int[]; iters = Int[]
        ratio = Float64[]; expo = Float64[]; layer = Float64[]; nu = Float64[]
        conv = Bool[]; resid = Float64[]
        tube_lo = Float64[]; tube_hi = Float64[]; tube_rmin = Float64[]
        out_lo = Float64[]; out_hi = Float64[]
        for m in MMS
            res  = s[m]
            lvl  = res.levels[end]
            prev = res.levels[end-1]
            mesh = lvl.mesh
            pat  = Patch[]
            for e in eachindex(mesh.elements)
                el  = mesh.elements[e]
                blk = mesh.blocks[mesh.blockindex[el.block]]
                em  = element_matrices(Float64, blk, el, mesh.orders[e])
                nx, ny = mesh.orders[e]
                off = mesh.offsets[e]
                # the caption's clip, applied to the NODES: the wavezone element of a
                # compact ISCO domain starts inside 2 r₀ and runs to R = 10⁴ M, so a per-element
                # test would drag four decades of radius into the plane's statistics.
                keep = findall(<=(p.rmax), em.x[1])
                isempty(keep) && continue                            # clipped at r = 2 r₀
                rs = collect(em.x[1])[keep]
                ws = collect(em.x[2])
                v  = [abs(res.u[off + tensor_lin(i, j, nx)]) for i in keep, j in 1:ny]
                push!(pat, Patch(rs, ws, v, inside_worldtube(p.dom, el)))
            end
            push!(psurf, pat)
            # the plotted amplitude and the in-/out-of-worldtube dynamic ranges of THIS plane
            allv = reduce(vcat, [vec(q.v) for q in pat])
            tv = [x for q in pat if q.tube for x in q.v]
            ov = [x for q in pat if !q.tube for x in q.v]
            push!(maxpsi, maximum(allv))
            push!(tube_lo, minimum(tv)); push!(tube_hi, maximum(tv))
            push!(out_lo,  minimum(ov)); push!(out_hi,  maximum(ov))
            # radius of the deepest in-worldtube node of this plane (the "sharp drop")
            rmin = 0.0; vmin = Inf
            for q in pat, i in eachindex(q.r), j in eachindex(q.w)
                q.tube && q.v[i, j] < vmin && (vmin = q.v[i, j]; rmin = q.r[i])
            end
            push!(tube_rmin, rmin)
            f = mode_self_force(p.orb, m, lvl.psi, lvl.dpsi_r)   # ∂_zΨ^ℛ = 0: gated above
            push!(fr, f.Fr)
            push!(psi_xp, abs(lvl.psi));   push!(dpsi_xp, abs(lvl.dpsi_r))
            push!(psi_xp_prev, abs(prev.psi)); push!(dpsi_xp_prev, abs(prev.dpsi_r))
            push!(ndofs, mesh.ndofs); push!(ndofs0, res.levels[1].mesh.ndofs)
            push!(nelem, length(mesh.elements))
            push!(iters, lvl.record.iterations)
            # the sealed row-58 instrument on the plotted run, and the solver's own verdict
            hp = horizon_axis_profile(p.dom, p.orb, mesh, res.u)
            push!(ratio, hp.ratio); push!(expo, hp.exponent); push!(layer, hp.layer_ratio)
            push!(nu, horizon_branch_index(p.orb, m))
            push!(conv, res.converged); push!(resid, lvl.record.final_residual)
        end
        @printf("[fig_modes] a = %+.4f: sampled %d planes\n", p.a, length(psurf)); flush(stdout)
        push!(surfaces, psurf)
        push!(tables, (maxpsi = maxpsi, fr = fr, psi_xp = psi_xp, dpsi_xp = dpsi_xp,
                       psi_xp_prev = psi_xp_prev, dpsi_xp_prev = dpsi_xp_prev,
                       ndofs = ndofs, ndofs0 = ndofs0, nelem = nelem, iters = iters,
                       ratio = ratio, expo = expo, layer = layer, nu = nu,
                       conv = conv, resid = resid,
                       tube_lo = tube_lo, tube_hi = tube_hi, tube_rmin = tube_rmin,
                       out_lo = out_lo, out_hi = out_hi))
    end

    metas = [(a = p.a, r0 = p.r0, spec = p.spec, rp = p.orb.rp, wr = p.wr, rmax = p.rmax,
              bounds = collect(p.dom.bounds_r), xp = p.xp, omega = p.orb.Omega,
              hold = q.hold, layer_probe = q.probe_ratio) for (p, q) in zip(panels, sweeps)]
    if opts.mode === :regen
        open(SURFACE_CACHE, "w") do io
            serialize(io, (surfaces = surfaces, tables = tables, meta = metas))
        end
        println("[fig_modes] wrote the surface cache ",
                relpath(SURFACE_CACHE, ReproHarness.REPRO_ROOT), " (",
                filesize(SURFACE_CACHE), " bytes, untracked) for --plot-only")
    end
else
    CACHE = open(deserialize, SURFACE_CACHE)
    surfaces = CACHE.surfaces
    tables = CACHE.tables
    metas = CACHE.meta
    @printf("[fig_modes] --plot-only: %d panels x %d modes from %s\n", length(surfaces),
            length(first(surfaces)), relpath(SURFACE_CACHE, ReproHarness.REPRO_ROOT))
end

# ───────────────────────────── the node's gates and the golden ───────────────────────────────
const ZHI = 0.5      # the reference's vertical window: ticks 0 … −5 (manual § 5)
const ZLO = -5.5

# The paper's two claims about this figure are MEASURED per panel and REGISTERED, not asserted:
# a claim that fails is a golden row that says so (the fig_linsolv precedent), so `--check` keeps
# comparing numbers instead of aborting, and a later change that breaks a claim which currently
# HOLDS fails the comparison. `gatefails` is reserved for structural faults (a lost mode, a
# violated precondition of the extraction, an oversized artifact).
claims = NamedTuple[]
for (k, (mt, tb)) in enumerate(zip(metas, tables))
    tag = @sprintf("p%d", k)
    # ---- CLAIM A: "the decay in power over mode number" (main.tex:718) ----------------------
    # measured on the PLOTTED quantity: the largest |Ψ_m| of the drawn plane.
    A = tb.maxpsi
    mono2 = all(A[i+1] < A[i] for i in 3:length(A)-1)      # m ≥ 2 (A is 1-based in m)
    # the m values (not indices) where the decay breaks, over the CLAIMED range m >= 2; -1 = none.
    # `A` is 1-based in m, so index i is mode m = i - 1 and a break at i means A[i+1] >= A[i].
    breaks = [i - 1 for i in 3:length(A)-1 if !(A[i+1] < A[i])]
    firstbreak = isempty(breaks) ? -1 : first(breaks)
    lastbreak  = isempty(breaks) ? -1 : last(breaks)
    # ---- CLAIM B: "regular zero crossings … display as sharp drops … particularly for the
    #      m = 0 mode" (main.tex:718) — the in-worldtube dynamic range of each plane, against
    #      the same range OUTSIDE the tube where the field is smooth.  A hard zero at a node
    #      would make the count infinite; the floor is stated, not silent, and lies far below
    #      anything this figure resolves.
    fl(x) = max(x, 1.0e-300)
    tube_dec = log10.(tb.tube_hi ./ fl.(tb.tube_lo))
    out_dec  = log10.(tb.out_hi ./ fl.(tb.out_lo))
    # second instrument for the same claim, local to the particle instead of global over the
    # plane: |Ψ^ℛ(x_p)|, the ingredient row 44 turns into F_r^m.
    P = tb.psi_xp
    mono_xp = all(P[i+1] < P[i] for i in 3:length(P)-1)
    m0_drop = tube_dec[1] >= 2
    m0_deeper_than_neighbours = tube_dec[1] > tube_dec[2] && tube_dec[1] > tube_dec[3]
    top = log10(maximum(A))
    if !PLOT_ONLY
        register!(tag * "_spin", mt.a)
        register!(tag * "_r_isco", mt.r0)
        register!(tag * "_r_plus", mt.rp)
        register!(tag * "_worldtube_radius", mt.wr)
        register!(tag * "_block_bounds_r", collect(mt.bounds))
        register!(tag * "_x_p", collect(mt.xp))
        register!(tag * "_r_clip", mt.rmax)
        register!(tag * "_omega", mt.omega)
        register!(tag * "_final_elements", tb.nelem)
        register!(tag * "_final_ndofs", tb.ndofs)
        register!(tag * "_level0_ndofs", tb.ndofs0)
        register!(tag * "_final_gmres_iterations", tb.iters)
        register!(tag * "_max_abs_psi", A)
        register!(tag * "_log10_max_abs_psi", log10.(A))
        register!(tag * "_abs_psi_xp", tb.psi_xp)
        register!(tag * "_abs_dpsi_r_xp", tb.dpsi_xp)
        register!(tag * "_abs_psi_xp_prev_level", tb.psi_xp_prev)
        register!(tag * "_abs_dpsi_r_xp_prev_level", tb.dpsi_xp_prev)
        register!(tag * "_Fr", tb.fr)
        register!(tag * "_worldtube_decades", tube_dec)
        register!(tag * "_outside_decades", out_dec)
        register!(tag * "_worldtube_min_radius", tb.tube_rmin)
        register!(tag * "_log10_top_plane", top)
        register!(tag * "_decay_monotone_from_m2", mono2)
        register!(tag * "_psi_xp_monotone_from_m2", mono_xp)
        register!(tag * "_decay_last_break_m", lastbreak)
        register!(tag * "_decay_first_break_m", firstbreak)
        register!(tag * "_m0_drops_two_decades", m0_drop)
        register!(tag * "_m0_deeper_than_m1_m2", m0_deeper_than_neighbours)
        register!(tag * "_deepest_worldtube_mode", argmax(tube_dec) - 1)
        register!(tag * "_horizon_hold", mt.hold)
        register!(tag * "_layer_ratio_default_max", mt.layer_probe)
        register!(tag * "_layer_ratio", tb.layer)
        register!(tag * "_axis_equator_ratio", tb.ratio)
        register!(tag * "_axis_exponent", tb.expo)
        register!(tag * "_branch_index", tb.nu)
        register!(tag * "_converged", tb.conv)
        register!(tag * "_final_residual", tb.resid)
    end
    # GATE 7 (row 58's rule, on the PLOTTED run): the near-horizon strip must never carry more
    # than LAYER_GATE of the rest of the plane. The driver CHOOSES `horizon_hold` to satisfy it,
    # so a violation here means the rule did not work and the panel is not drawable as physics.
    maximum(tb.layer) <= LAYER_GATE || push!(gatefails,
        "GATE 7: a = $(mt.a) — plotted layer_ratio max " *
        @sprintf("%.4g", maximum(tb.layer)) * " > " * @sprintf("%.2f", LAYER_GATE) *
        " at m = $(argmax(tb.layer) - 1) (horizon_hold = $(mt.hold), default sweep " *
        @sprintf("%.4g", mt.layer_probe) * ")")
    @printf("\n[fig_modes] panel %d  a = %+.4f  r₀ = %.6f  (%d elements, %d dofs)\n",
            k, mt.a, mt.r0, tb.nelem[1], tb.ndofs[1])
    push!(claims, (mono2 = mono2, lastbreak = lastbreak, firstbreak = firstbreak,
                   mono_xp = mono_xp))
    println("[fig_modes]   CLAIM A max|Ψ_m| decays monotonically for m ≥ 2 : ",
            mono2 ? "YES" : "NO  (breaks first after m = $(firstbreak), last after m = $(lastbreak))",
            "   |Ψ^ℛ(x_p)|: ", mono_xp ? "YES" : "NO")
    @printf("[fig_modes]   CLAIM B m=0 drops %.2f decades in the worldtube (m=1 %.2f, m=2 %.2f); deepest plane m = %d\n",
            tube_dec[1], tube_dec[2], tube_dec[3], argmax(tube_dec) - 1)
    @printf("[fig_modes]   top plane at log₁₀|Ψ| = %+.3f\n", top)
    @printf("[fig_modes]   horizon_hold = %s (default sweep layer_ratio max %.4g, gate %.2f); plotted layer_ratio max %.4g at m = %d\n",
            string(mt.hold), mt.layer_probe, LAYER_GATE, maximum(tb.layer), argmax(tb.layer) - 1)
    let bad = [i - 1 for i in eachindex(tb.conv) if !tb.conv[i]]
        isempty(bad) ? println("[fig_modes]   every mode converged") :
            @printf("[fig_modes]   NOT CONVERGED: m = %s (worst final residual %.3e against the 1e-12 tolerance)\n",
                    join(bad, ", "), maximum(tb.resid[i] for i in eachindex(tb.conv) if !tb.conv[i]))
    end
    println("   m   log₁₀max|Ψ_m|   |Ψ^ℛ(x_p)|    F_r^m           tube[dec]  axis:eq    layer   its  conv")
    for i in eachindex(A)
        @printf("  %2d   %+8.4f      %.4e   %+.6e   %6.2f   %9.3e %8.3e %4d  %s\n", i - 1,
                log10(A[i]), tb.psi_xp[i], tb.fr[i], tube_dec[i], tb.ratio[i], tb.layer[i],
                tb.iters[i], tb.conv[i] ? "yes" : "NO")
    end
    flush(stdout)
end

if !PLOT_ONLY
    @repro_result "panel_count" length(metas)
    @repro_result "modes_per_panel" length(MMS)
    @repro_result "amr_grids" ITER + 1
    @repro_result "worldtube_radius_fraction" WR_FRAC
    @repro_result "polynomial_order" PORD
    @repro_result "initial_level" LREF
    # the lattice every number above belongs to: the production `ImposeEquatorialSymmetry: True`
    # (yaml:36) z² half-domain, NOT the full-z lattice of spectre's own regression input
    @repro_result "impose_equatorial_symmetry" panels[1].cfg.orbit.impose_equatorial_symmetry
end

# ───────────────────────────── the figure (presentation only) ────────────────────────────────
# RENDERING. The reference is a matplotlib 3D wireframe stack assembled in Keynote. Plots.jl/GR
# has no `plot_wireframe`, so this driver projects the (r, z², log₁₀|Ψ|) lattice itself with a
# fixed oblique axonometry and draws each plane as NaN-separated polylines through one series per
# mode — the same ink the reference lays down, in the same order (high m first, so the coloured
# low modes sit on top). Honest limitations, none structural:
#   * no depth sorting between planes: like matplotlib's own wireframe, an overlap is resolved by
#     draw order and transparency, not by occlusion;
#   * one plane is drawn per ELEMENT, from that element's own Legendre-Gauss nodes and nothing
#     else, so neighbouring patches do not touch — which is exactly how the reference's
#     ParaView/DG cells read (its curtain is visibly made of separate cell patches);
#   * the collocation lattice is subsampled to 3 lines per element per dimension for the five
#     coloured modes and 2 for the faint m ≥ 5 curtain: the 8th AMR grid carries orders up to
#     P + 7, and drawing every node line costs the 200 KB golden budget without adding a
#     visible feature.
# Every presentation constant below (box proportions, oblique offsets, colours, alphas, tick
# choices, label offsets) is a STYLE choice recorded here; none touches a registered value.
using Plots
using Plots: RGB, Shape, text
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

const PW  = 1.00     # box width  (r)
const PH  = 1.62     # box height (log₁₀|Ψ|)
const DXW = 0.17     # oblique offset of z² = 1, horizontal
const DYW = 0.30     # ... and vertical
const TIL = 0.05     # the r axis tilts down to the right, as the reference's does
const GAP = 0.58     # gap between panels
const PXU = 127      # px per canvas unit — the golden's size knob
const SC  = PXU / 127
# Plots' `font()` parser reads an INTEGER argument as the point size and ANY OTHER Real as a
# ROTATION IN DEGREES, so a Float64 size silently tilts the label (and leaves the size at its
# default). Every text size therefore goes through `FS`.
FS(pt) = max(4, round(Int, pt * SC))
# manual § 5's categorical palette: m = 0 brick red, 1 golden, 2 teal, 3 blue, 4 salmon, ≥5 gray.
# The grey ramp is the PUBLISHED one, `m/23` — measured off the paper's own fig:convergence lines
# by ESSF::fig-convergence and corrected in `_style/manuals/arXiv-2606.04998.md` in this commit
# (the manual's earlier `linspace(0, 0.8, 21)` was an eyeball reading).
const CMODE = (RGB(0.690, 0.255, 0.243), RGB(0.933, 0.765, 0.247), RGB(0.224, 0.706, 0.541),
               RGB(0.239, 0.490, 0.741), RGB(0.886, 0.412, 0.361))
graylevel(m) = m / 23        # the published ramp, MEASURED from fig:convergence's own lines
mcolor(m) = m <= 4 ? CMODE[m+1] : RGB(graylevel(m), graylevel(m), graylevel(m))
malpha(m) = m <= 4 ? 0.80 : 0.50 - 0.20 * (m - 5) / 15

"the panel's own screen map: r → [0,1], log₁₀|Ψ| → [0,1], z² as the oblique third axis"
struct Cam; x0::Float64; rlo::Float64; rhi::Float64; end
sx(c::Cam, r, w) = c.x0 + PW * (r - c.rlo) / (c.rhi - c.rlo) + DXW * w
sy(c::Cam, r, w, h) = PH * clamp((h - ZLO) / (ZHI - ZLO), 0.0, 1.0) + DYW * w -
                      TIL * (r - c.rlo) / (c.rhi - c.rlo)

sub(n, k) = n <= k ? collect(1:n) : unique(round.(Int, range(1, n; length = k)))

"nice round ticks inside [lo, hi]: the smallest 1/2/5 step leaving at most six (style only)"
function nticks(lo, hi; ints_hint::Bool = true)
    span = hi - lo
    base = 10.0^floor(log10(span))
    for s in (0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
        step = base * s
        ts = filter(t -> lo + 0.03span <= t <= hi, collect(step * ceil(lo / step):step:hi))
        length(ts) <= 6 && return (ts, step)
    end
    return ([lo + span / 2], span)
end

# label with exactly as many decimals as the chosen step needs (1.25 must not print as "1.2")
function fmt_tick(t, step)
    nd = 0
    while nd < 3 && abs(step * 10.0^nd - round(step * 10.0^nd)) > 1e-9
        nd += 1
    end
    nd == 0 && return string(Int(round(t)))
    nd == 1 && return @sprintf("%.1f", t)
    return @sprintf("%.2f", t)
end

plt = plot([NaN], [NaN]; label = "", legend = false, axis = false, grid = false,
           ticks = false, framestyle = :none, background_color = :white,
           aspect_ratio = :equal, dpi = 100)

cams = Cam[]
for (k, surf) in enumerate(surfaces)
    meta = metas[k]; tb = tables[k]
    cam = Cam((k - 1) * (PW + DXW + GAP), meta.rp, meta.rmax)
    push!(cams, cam)
    # the planes, high m first; a plane whose own solve did not converge is drawn DASHED and
    # fainter, and the panel note below names it — it stays in the figure, flagged, never hidden
    for m in reverse(0:(length(surf)-1))
        ok = tb.conv[m+1]
        nl = m <= 4 ? 3 : 2      # wireframe stride: the faint high-m curtain needs fewer lines
        xs = Float64[]; ys = Float64[]
        for q in surf[m+1]
            nx, ny = length(q.r), length(q.w)
            for j in sub(ny, nl)
                for i in 1:nx
                    v = q.v[i, j]
                    if q.r[i] <= meta.rmax && v > 0 && log10(v) > ZLO
                        push!(xs, sx(cam, q.r[i], q.w[j])); push!(ys, sy(cam, q.r[i], q.w[j], log10(v)))
                    else
                        push!(xs, NaN); push!(ys, NaN)
                    end
                end
                push!(xs, NaN); push!(ys, NaN)
            end
            for i in sub(nx, nl)
                q.r[i] > meta.rmax && continue
                for j in 1:ny
                    v = q.v[i, j]
                    if v > 0 && log10(v) > ZLO
                        push!(xs, sx(cam, q.r[i], q.w[j])); push!(ys, sy(cam, q.r[i], q.w[j], log10(v)))
                    else
                        push!(xs, NaN); push!(ys, NaN)
                    end
                end
                push!(xs, NaN); push!(ys, NaN)
            end
        end
        plot!(plt, xs, ys; linecolor = mcolor(m), linealpha = ok ? malpha(m) : 0.55 * malpha(m),
              linestyle = ok ? :solid : :dash,
              linewidth = (m <= 4 ? 0.70 : 0.45) * SC, label = "")
    end
    # the box
    corners = [(r, w, h) for h in (ZLO, ZHI) for w in (0.0, 1.0) for r in (cam.rlo, cam.rhi)]
    edges = [(1,2),(3,4),(1,3),(2,4),(5,6),(7,8),(5,7),(6,8),(1,5),(2,6),(3,7),(4,8)]
    for (i, j) in edges
        a1 = corners[i]; a2 = corners[j]
        plot!(plt, [sx(cam, a1[1], a1[2]), sx(cam, a2[1], a2[2])],
              [sy(cam, a1[1], a1[2], a1[3]), sy(cam, a2[1], a2[2], a2[3])];
              linecolor = RGB(0.25, 0.25, 0.25), linewidth = 0.9SC, label = "")
    end
    # vertical ticks 0 … −5 on both front verticals, labels outside
    for t in 0:-1:-5
        for (rr, ww, dxlab, al) in ((cam.rlo, 0.0, -0.045, :right),
                                    (cam.rhi, 1.0, 0.045, :left))
            x = sx(cam, rr, ww); y = sy(cam, rr, ww, float(t))
            plot!(plt, [x, x + sign(dxlab) * 0.022], [y, y]; linecolor = :black,
                  linewidth = 0.9SC, label = "")
            annotate!(plt, x + dxlab, y, text(string(t), FS(9), :black, al))
        end
    end
    # r ticks along the bottom front edge
    # up to six ticks whatever the step: the 1/2/5 family never needs more than one decimal
    # here, and six is what the published panels carry (5/10/15, 2…12, 1.5…4.0).
    rt, rstep = nticks(cam.rlo, cam.rhi)
    for t in rt
        x = sx(cam, t, 0.0); y = sy(cam, t, 0.0, ZLO)
        plot!(plt, [x, x], [y, y - 0.028]; linecolor = :black, linewidth = 0.9SC, label = "")
        annotate!(plt, x, y - 0.085, text(fmt_tick(t, rstep), FS(length(rt) >= 6 ? 6 : 8), :black, :center))
    end
    annotate!(plt, sx(cam, (cam.rlo + cam.rhi) / 2, 0.0), sy(cam, cam.rhi, 0.0, ZLO) - 0.21,
              text("\$r\$", FS(14), :black, :center))
    annotate!(plt, sx(cam, (cam.rlo + cam.rhi) / 2, 0.5), sy(cam, cam.rlo, 1.0, ZHI) + 0.315,
              text(@sprintf("\$a = %.3f, r_0 = %.3f\$", meta.a, meta.r0), FS(12),
                   :black, :center))
    # THE PANEL NOTE, every line of it read off this panel's own measurement: the row-58 rule's
    # verdict and the setting it chose, then any mode whose solve did not converge, then the
    # decay claim if it fails. Nothing here is typed.
    notes = String[]
    push!(notes, @sprintf("layer %.3g (def %.3g)  hold=%s", maximum(tb.layer), meta.layer_probe,
                          meta.hold ? "true" : "false"))
    let bad = [i - 1 for i in eachindex(tb.conv) if !tb.conv[i]]
        isempty(bad) || push!(notes, @sprintf("m=%s GMRES stalled at %.1e", join(bad, ","),
                                              maximum(tb.resid[i] for i in eachindex(tb.conv) if !tb.conv[i])))
    end
    claims[k].mono2 || push!(notes, "PRELIM: decay stops after m=$(claims[k].firstbreak)")
    for (i, nt) in enumerate(notes)
        annotate!(plt, sx(cam, (cam.rlo + cam.rhi) / 2, 0.5),
                  sy(cam, cam.rlo, 1.0, ZHI) + 0.225 - 0.077 * (i - 1),
                  text(nt, FS(6), RGB(0.35, 0.35, 0.35), :center))
    end
end

# the shared vertical label and the colour legend
annotate!(plt, cams[1].x0 - 0.30, PH / 2,
          text("\$\\log_{10}|\\Psi|\$", FS(14), :black, :center, rotation = 90))
lx = cams[end].x0 + PW + DXW + 0.20
for (i, lab) in enumerate(("\$m = 0\$", "\$m = 1\$", "\$m = 2\$", "\$m = 3\$", "\$m = 4\$",
                           "\$m \\geq 5\$"))
    m = i - 1
    y = PH * 0.78 - 0.185 * (i - 1)
    plot!(plt, Shape([lx, lx + 0.13, lx + 0.13, lx], [y, y, y + 0.058, y + 0.058]);
          fillcolor = (m <= 4 ? mcolor(m) : RGB(0.20, 0.20, 0.20)), linewidth = 0,
          fillalpha = 1.0, label = "")
    annotate!(plt, lx + 0.17, y + 0.033, text(lab, FS(11), :black, :left))
end

fx0 = cams[1].x0 - 0.44; fx1 = lx + 0.62
fy0 = -0.34; fy1 = PH + DYW + 0.45
plot!(plt; xlims = (fx0, fx1), ylims = (fy0, fy1),
      size = (round(Int, PXU * (fx1 - fx0)), round(Int, PXU * (fy1 - fy0))))

const GOLD_PNG = joinpath(@__DIR__, "expected", "fig_modes.png")
publish = opts.mode === :regen || PLOT_ONLY
out_png = publish ? GOLD_PNG : joinpath(mktempdir(), "fig_modes.png")
saved = save_figure!(plt, out_png; max_bytes = 200 * 1024, pdf = publish, id = "fig_modes")
saved.ok || push!(gatefails, "GATE 6: " * saved.message)
publish || println("[fig_modes] --check: rendered to a scratch path, ",
                   relpath(GOLD_PNG, ReproHarness.REPRO_ROOT), " untouched")

if PLOT_ONLY
    println("[fig_modes] --plot-only: redrawn from the surface cache, no physics recomputed")
    isempty(gatefails) || (println("\n[fig_modes] FINDING:"); foreach(f -> println("    ", f), gatefails))
    exit(isempty(gatefails) ? 0 : 1)
end

@repro_result "png_bytes_within_budget" saved.ok
@repro_result "gates_pass" isempty(gatefails)
if !isempty(gatefails)
    println("\n[fig_modes] FINDING — node gates violated:")
    foreach(f -> println("    ", f), gatefails)
end

code = ReproHarness.conclude!(opts; paper = "2606.04998", id = "fig_modes",
                              script = @__FILE__, rtol = 1e-8)
exit(isempty(gatefails) ? code : max(code, 1))
