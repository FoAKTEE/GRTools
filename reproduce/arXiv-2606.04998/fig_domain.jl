# reproduce/arXiv-2606.04998/fig_domain.jl — ESSF::fig-domain.
#
# arXiv:2606.04998 Fig. 1 (`fig:domain`, main.tex:526-530): "Computational domain for
# a = 0.5 M, r0 = 10 M, m = 1 composed of discontinuous Galerkin elements. The domain extends
# from the Kerr horizon r+ to the outer boundary at R = 10^4 M, and from z^2 = 0 (equator) to
# z^2 = 1 (poles) with equatorial symmetry (optional). Plotted is the complex magnitude of the
# field Psi_m. In the regularized region (red rectangle), the complex magnitude of the
# regularized field Psi^R_m = Psi_m - Psi^P_m is plotted instead. Null slicing is applied in
# the v and u regions with boost function H = +-1. The inset magnifies the region around the
# position of the puncture, x_p."
#
# THIN driver: zero physics. Every number comes from `ScalarSelfForce` / `EllipticCore` —
# `crp_paper_config` -> `block_domain` -> `solve_mode_amr` (the row-43 hp-AMR sweep) and then
# only the SEALED read-outs `field_at`, `reconstruct_full_field`, `inside_worldtube`,
# `puncture_chart_position`, `worldtube_blocks`, `element_matrices`, `map_forward`,
# `logical_interval`, `two_to_one_balanced`. The driver arranges, gates and draws.
#
# PARAMETERS ARE THE CAPTION'S, not the yaml header's: a = 0.5 M, r0 = 10 M, m = 1
# (main.tex:529). Everything the caption does NOT name is the production example of the
# paper's own input file (anc/ScalarSelfForce.yaml:8-12): R = 10^4 M, worldtube_radius
# = 4.8 M, worldtube_dz2 = 0.25, P = 4, L = 0, ImposeEquatorialSymmetry: True. The inner
# radius is r+(a = 0.5) = 1.8660254037844386 (block_domain REQUIRES bounds_r[1] == r+ for the
# `:horizon` no-BC contract — the row-39 landed gate).
#
# SIX h-REFINEMENT STEPS, from the paper's own sentence about THIS figure (main.tex:558):
# "The result of this h refinement strategy after six refinement steps is shown in
# \cref{fig:domain}, with the region around the puncture magnified in the inset." So the
# plotted mesh is `solve_mode_amr(...; iterations = 6)` — level 6 of the sweep, NOT the yaml's
# production `Amr.Iterations: 8`. The inset of the reference figure agrees independently: its
# element cascade shows six 2:1 steps down to the puncture corner, and the coarse elements of
# that inset carry a 9-cell (10-point) collocation grid, i.e. N = P + 6 = 10.
#
# THE ABSCISSA IS THE PAPER'S OWN VISUALISATION COORDINATE. A domain that runs to R = 10^4 M
# cannot be drawn against linear r. The paper's run writes
# `RadiallyCompressedCoordinates: {InnerRadius: *t_to_u_transition, OuterRadius: *outer_radius,
# Compression: *outer_shell_distribution}` (yaml:222-225, transcribed in `cfg.observers`), and
# spectre's implementation (ref-code/spectre `src/Domain/RadiallyCompressedCoordinates.cpp:27-52`
# with `src/Domain/CoordinateMaps/Interval.cpp:120-129`) defines, for r > r_in,
#
#     r_c(r) = 1/2 (A + B + (B - A)/(b - a) (a + b - 2 a b / r)),
#     A = a = r_in,  b = R,  B = r_in * log10(R)          (the "compressed outer radius")
#
# i.e. at the production point r_c(r) = 97.77904 - 1790.4487/r on [24.4, 10^4] -> [24.4, 97.6],
# and r_c(r) = r inside r_in. That map and its return trip are `radially_compressed` /
# `radially_compressed_inv` in `capabilities/ScalarSelfForce/src/Domain/config.jl` (gated
# against the C++ expressions in test_runtime_config.jl); this driver only CALLS them, so the
# 66 golden samples laid out uniformly in the drawn coordinate are pinned against the package,
# never against arithmetic the driver owns. FRAMING: the reference render itself stops at the
# u-block element boundary r = 48.68 M (measured: the published worldtube patch sits at
# 4.5-21.5 % of the drawn ribbon, against 5.6-21.9 % predicted by r_c cut at that element),
# so this figure is framed the same way — the outermost wavezone element runs off the right
# edge exactly as it does in the paper. Nothing is removed from the solve or the golden.
#
# WHAT IS PLOTTED IS `sol.u` EVERYWHERE. The row-42 contract is that a `ModeSolution`'s `u`
# holds Psi^R inside the worldtube and Psi outside it, which is verbatim what the caption
# asks for. `reconstruct_full_field` (Psi everywhere) is computed only for the GATES — it is
# log-divergent at x_p and is not drawn.
#
# THE ORDER MAP IS ANCHORED TO THE FIXED `hp_refine`. EllipticCore's `hp_refine` closes the
# caller's flags under the 2:1 rule (`balanced_flags`) before the p half, so an element the
# balancer splits KEEPS its polynomial order in the split dimension — spectre
# `Domain/Amr/UpdateAmrDecision.cpp:82-87` + `ParallelAlgorithms/Amr/Projectors/Mesh.cpp:30-34`.
# MEASURED here: the two balance-split children of block 1 at level 1 carry `(5, 4)`, not
# `(5, 5)`. `element_order_x/y`, `level_ndofs`, `level_map_sha256`, `level_gmres_iterations`
# and the sampled field values therefore all depend on that closure; `element_id_sha256`, the
# drawn element boundaries and the inset cascade do NOT (the h-refined element SET is the same
# either way). If this golden is ever checked against a tree WITHOUT the closure, those are
# exactly the rows that move.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2606.04998/fig_domain.jl --check
#
# Style: `_style/manuals/arXiv-2606.04998.md` sections 0 + 1. Honest Plots.jl-vs-ParaView
# limitations are listed under "RENDERING" below.

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()

# `--plot-only` redraws the figure from the surface CACHE this driver writes beside the golden
# on `--regen`, with no solve and no comparison — the fast loop for a purely cosmetic change.
# It cannot redraw from the golden itself the way fig_linsolv does: this figure's ink is 57
# element sample lattices (~8k values), which do not belong in a 200 KB golden, so the cache
# carries them instead. It is a same-Julia-version scratch artifact, untracked
# (`expected/.gitignore`), and its absence is a HOLE, never a silent fallback.
const SURFACE_CACHE = joinpath(@__DIR__, "expected", "fig_domain_surface.jls")
const PLOT_ONLY = "--plot-only" in opts.extra

PLOT_ONLY && !isfile(SURFACE_CACHE) && ReproHarness.hole(
    "--plot-only needs the surface cache " *
    "$(relpath(SURFACE_CACHE, ReproHarness.REPRO_ROOT)); run --regen first")

PLOT_ONLY || ReproHarness.require_tier(opts, :heavy)   # measured; see the README

using Printf, LinearAlgebra, SHA, Statistics, Serialization
using GrTools
using GrTools.ScalarSelfForce
using GrTools.EllipticCore

# ───────────────────────────── parameters ────────────────────────────────────────────────────
const SPIN = 0.5            # main.tex:529 (the caption)
const R0   = 10.0           # main.tex:529
const MM   = 1              # main.tex:529
const ITER = 6              # main.tex:558 — "after six refinement steps"
const ROUT = 1.0e4
const WR   = 4.8
const WDZ2 = 0.25
const PORD = 4
const LREF = 0

for f in (:solve_mode_amr, :puncture_chart_position, :reconstruct_full_field, :field_at)
    isdefined(ScalarSelfForce, f) || ReproHarness.hole(
        "ScalarSelfForce.$(f) is not available: Fig. 1 needs the hp-AMR driver " *
        "(ESSF::amr-driver, capabilities/ScalarSelfForce/src/Solve/amr.jl).")
end

orb = CircularOrbit(SPIN, R0)
const RIN = orb.rp                      # block_domain requires bounds_r[1] == r+
cfg = crp_paper_config(; spin = SPIN, orbital_radius = R0, m_mode = MM, inner_radius = RIN,
                       outer_radius = ROUT, worldtube_radius = WR, worldtube_dz2 = WDZ2,
                       P = PORD, L = LREF)
dom = block_domain(orb, cfg)
xp = puncture_chart_position(dom, orb)

println("\n[fig_domain] a = ", SPIN, "M  r₀ = ", R0, "M  m = ", MM, "  P = ", PORD,
        "  L = ", LREF, "  AMR steps = ", ITER)
println("[fig_domain] lattice ", dom.nr, "×", dom.nz, " blocks (", dom.chart, "), r₊ = ", RIN)
println("[fig_domain] radial block bounds ", dom.bounds_r)
println("[fig_domain] angular block bounds ", dom.bounds_z)
println("[fig_domain] worldtube block(s) ", worldtube_blocks(dom), "   x_p = ", xp)

# ───────────────────────── the drawn radial coordinate and the frame ─────────────────────────
# The paper's own visualisation abscissa, taken from `cfg.observers` through the SEALED
# `radially_compressed` / `radially_compressed_inv` (ScalarSelfForce `Domain/config.jl`,
# gated against spectre's `Interval` `Inverse` branch in test_runtime_config.jl). The driver
# owns no arithmetic here — it only names the coordinate it draws in.
rcc(r)     = radially_compressed(cfg.observers, r)
rcc_inv(x) = radially_compressed_inv(cfg.observers, x)

const RC_LO = rcc(dom.bounds_r[1]); const RC_HI = rcc(dom.bounds_r[end])
const WT_R = (dom.bounds_r[3], dom.bounds_r[4])       # the worldtube block's radial bounds
const WT_W = (dom.bounds_z[1], dom.bounds_z[2])       # …and its angular bounds
const DELTA = 1.0e-6

# The published render stops at the FIRST element boundary of the u block — the wavezone
# column's own 1/r map evaluated at its logical midpoint — so the frame is derived from the
# lattice rather than typed in: `r = 2 r_u R / (r_u + R)` at the production point, 48.681 M.
const U_BLOCK = only(filter(b -> block_region(dom, b) == :u && b < dom.nr, 0:(dom.nr*dom.nz-1)))
const U_SPLIT = map_forward(blockmap(dom.mesh.blocks[dom.mesh.blockindex[U_BLOCK]], 1), 0.0)
const R_FRAME = 1.068 * U_SPLIT   # +6.8 % so that boundary itself stays inside the frame


"one element's sample lattice: its own collocation nodes plus its two faces, per dimension."
struct ElemGrid; r::Vector{Float64}; w::Vector{Float64}; v::Matrix{Float64}; end

gatefails = String[]
grids = ElemGrid[]
on_horizon = Bool[]
inside_wt = Bool[]
UMAX = 0.0

if !PLOT_ONLY
t_solve = time()
hist = solve_mode_amr(orb, cfg, MM; iterations = ITER, dom = dom)
t_solve = time() - t_solve
@printf("[fig_domain] solve_mode_amr(iterations=%d): %.1f s\n", ITER, t_solve)
display(hist); println()

sol = hist.solution
mesh = sol.mesh
full = reconstruct_full_field(sol)          # Ψ_m everywhere — GATES ONLY, never drawn

# ───────────────────────────── GATE 1: the refinement-level / order map ──────────────────────
# The metadata rule of the sealed AMR node (row 43, `amr_flags` + `hp_refine` ≡
# `h_refine(p_refine(mesh; hold = closed), closed)` with `closed = balanced_flags(mesh, flags)`),
# asserted level by level:
#   (a) the criterion flags exactly the elements CONTAINING the puncture — 2 of them in the z²
#       chart — and they sit at h-level (L+1+k, L+1+k) with the order HELD at (P, P);
#   (b) every element of level k descends from exactly one element of level k−1 and carries
#       that ancestor's order incremented in each dimension the CLOSED flag set did NOT split:
#       an element the 2:1 balance splits keeps its order in the split dimension (spectre
#       `Domain/Amr/UpdateAmrDecision.cpp:82-87` + `ParallelAlgorithms/Amr/Projectors/Mesh.cpp:30-34`,
#       the EllipticCore `hp_refine` contract). MEASURED consequence at this point: the two
#       balance-split children of block 1 at level 1 carry (5, 4), not (5, 5);
#   (c) the element set is 2:1 balanced at every level;
#   (d) every order lies in [P, P+k] and the maximum is exactly P+k.
# (b) reads the closure back through the SEALED `balanced_flags`, so it is a consistency
# assertion between consecutive levels and the recorded flag set, not an independent
# re-derivation; (a), (c) and (d) are independent of it.
function ancestor_of(prev::Mesh2D, id::ElementId)
    haskey(prev.index, id) && return id
    for d in (1, 2)
        id.segments[d].level > 0 || continue
        p = parent_id(id, d)
        haskey(prev.index, p) && return p
    end
    if id.segments[1].level > 0 && id.segments[2].level > 0
        p = parent_id(parent_id(id, 1), 2)
        haskey(prev.index, p) && return p
    end
    return nothing
end

puncture_hlevels = Int[]
puncture_orders = Int[]
order_rule_ok = true
balance_ok = true
order_window_ok = true
for (k, lv) in enumerate(hist.levels)
    global order_rule_ok, balance_ok, order_window_ok
    step = lv.step
    two_to_one_balanced(lv.mesh.elements) || (balance_ok = false;
        push!(gatefails, "GATE 1c: level $step is not 2:1 balanced"))
    if !all(o -> PORD <= o[1] <= PORD + step && PORD <= o[2] <= PORD + step, lv.mesh.orders) ||
       maximum(maximum, lv.mesh.orders) != PORD + step
        order_window_ok = false
        push!(gatefails, "GATE 1d: level $step orders $(sort(unique(lv.mesh.orders))) leave " *
              "the window [$PORD, $(PORD + step)] or miss the maximum $(PORD + step)")
    end
    length(lv.punctures) == 2 || push!(gatefails,
        "GATE 1a: level $step flags $(length(lv.punctures)) puncture elements, expected 2 " *
        "(the z² chart puts x_p on the domain boundary)")
    for el in lv.punctures
        e = lv.mesh.index[el]
        for d in (1, 2)
            push!(puncture_hlevels, el.segments[d].level)
            el.segments[d].level == LREF + 1 + step || push!(gatefails,
                "GATE 1a: puncture element $el at level $step has h-level " *
                "$(el.segments[d].level) in dim $d, expected $(LREF + 1 + step)")
        end
        push!(puncture_orders, lv.mesh.orders[e][1], lv.mesh.orders[e][2])
        lv.mesh.orders[e] == (PORD, PORD) || push!(gatefails,
            "GATE 1a: puncture element $el at level $step has order " *
            "$(lv.mesh.orders[e]), expected $((PORD, PORD)) (the order is HELD)")
    end
    step == 0 && continue
    prev = hist.levels[k-1]
    closed = balanced_flags(prev.mesh, Dict(el => (true, true) for el in prev.punctures))
    for (e, id) in enumerate(lv.mesh.elements)
        anc = ancestor_of(prev.mesh, id)
        if anc === nothing
            order_rule_ok = false
            push!(gatefails, "GATE 1b: element $id of level $step has no ancestor at level $(step-1)")
            continue
        end
        o = prev.mesh.orders[prev.mesh.index[anc]]
        h = get(closed, anc, (false, false))
        want = (h[1] ? o[1] : o[1] + 1, h[2] ? o[2] : o[2] + 1)
        if lv.mesh.orders[e] != want
            order_rule_ok = false
            push!(gatefails, "GATE 1b: element $id of level $step has order " *
                  "$(lv.mesh.orders[e]), the closed hp rule gives $want " *
                  "(ancestor $anc, closed flags $h)")
        end
    end
end


samp_r = Float64[]; samp_w = Float64[]; samp_tag = String[]
addpt!(r, w, tag) = (push!(samp_r, float(r)); push!(samp_w, float(w)); push!(samp_tag, tag))
# (a) the equatorial ray through x_p, uniform in the drawn coordinate (40 points)
for k in 1:40
    addpt!(rcc_inv(RC_LO + (RC_HI - RC_LO) * (k - 0.5) / 40), 0.0, "ray_equator")
end
# (b) the polar ray at r = r₀ (20 points)
for k in 1:20
    addpt!(R0, (k - 0.5) / 20, "ray_polar")
end
# (c) the worldtube boundary, inside and outside (18 points)
# the tube stations are FRACTIONS of the worldtube block, so they follow `worldtube_dz2`
# and `worldtube_radius` instead of repeating them (1/8, 3/8, 3/4 of the angular extent;
# 1/8, 5/16, 23/32 across the radial one — the values the golden was first laid out on)
wt_w(f) = WT_W[1] + f * (WT_W[2] - WT_W[1])
wt_r(f) = WT_R[1] + f * (WT_R[2] - WT_R[1])
# the angular fractions are dyadic, so they are EXACT and are used as derived; the radial
# ones are not (5.2 + 0.125·9.6 is 6.4 + 4e-16), so the probe keeps the round figure and the
# tie to `worldtube_radius` is an assertion instead — the probe set must not drift with the
# last bit of a fraction, but it must not silently stop tracking the config either
@assert (wt_w(0.125), wt_w(0.375), wt_w(0.75)) == (0.03125, 0.09375, 0.1875)
@assert all(isapprox.((wt_r(0.125), wt_r(0.3125), wt_r(0.71875)), (6.4, 8.2, 12.1); rtol = 1e-12))
for w in (wt_w(0.125), wt_w(0.375), wt_w(0.75))
    addpt!(WT_R[1] + DELTA, w, "tube_inner_r_in"); addpt!(WT_R[1] - DELTA, w, "tube_outer_r_in")
    addpt!(WT_R[2] - DELTA, w, "tube_inner_r_out"); addpt!(WT_R[2] + DELTA, w, "tube_outer_r_out")
end
for r in (6.4, 8.2, 12.1)
    addpt!(r, WT_W[2] - DELTA, "tube_inner_w"); addpt!(r, WT_W[2] + DELTA, "tube_outer_w")
end
# (d) the horizon column and the poles (12 points)
for w in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0); addpt!(dom.bounds_r[1], w, "horizon"); end
for k in 1:6; addpt!(rcc_inv(RC_LO + (RC_HI - RC_LO) * (k - 0.5) / 6), 1.0, "pole"); end
# (e) a coarse lattice (40 points)
for i in 1:8, j in 1:5
    addpt!(rcc_inv(RC_LO + (RC_HI - RC_LO) * (i - 0.5) / 8), (j - 0.5) / 5, "lattice")
end

samp_psiR = [abs(field_at(sol, sol.u, samp_r[k], samp_w[k])) for k in eachindex(samp_r)]
samp_psi  = [abs(field_at(sol, full,  samp_r[k], samp_w[k])) for k in eachindex(samp_r)]

# ───────────────────────────── GATE 2: the field itself ──────────────────────────────────────
all(isfinite, samp_psiR) || push!(gatefails,
    "GATE 2: |Ψ^ℛ_m| is not finite at every sample (horizon column / poles included)")
all(isfinite, samp_psi) || push!(gatefails, "GATE 2: |Ψ_m| is not finite at every sample")
all(>(0), samp_psiR) || push!(gatefails, "GATE 2: |Ψ^ℛ_m| vanishes at a sample")
inner = findall(t -> startswith(t, "tube_inner"), samp_tag)
ratio_inner = samp_psiR[inner] ./ samp_psi[inner]
all(<(1), ratio_inner) || push!(gatefails,
    "GATE 2: |Ψ^ℛ_m| is NOT smaller than |Ψ_m| at every inner worldtube-boundary sample " *
    "(worst ratio $(maximum(ratio_inner)))")
xp == (R0, 0.0) || push!(gatefails, "GATE 2: puncture_chart_position gave $xp, expected $((R0, 0.0))")
# SMALLER, and SMOOTHER, measured on the elements the criterion flags — the two that touch
# x_p. Over the WHOLE tube a max/min span does not discriminate (MEASURED 11.88 for Ψ^ℛ against
# 11.52 for Ψ: no Gauss node sits on the puncture, so the retarded field never blows up at a
# node and the span is set by the tube's outer half in both cases — the planned "flatter over
# the tube" proxy is simply not an instrument, and is replaced here by the one that is). On the
# puncture elements the log divergence is resolved and the relative spread separates them.
tube_els = [e for e in eachindex(mesh.elements) if inside_worldtube(dom, mesh.elements[e])]
psiR_tube = Float64[]; psi_tube = Float64[]
for e in tube_els
    el = mesh.elements[e]
    blk = mesh.blocks[mesh.blockindex[el.block]]
    em = element_matrices(Float64, blk, el, mesh.orders[e])
    for x2 in em.x[2], x1 in em.x[1]
        push!(psiR_tube, abs(field_at(sol, sol.u, x1, x2; element = e)))
        push!(psi_tube, abs(field_at(sol, full, x1, x2; element = e)))
    end
end
tube_amplitude_ratio = maximum(psiR_tube) / maximum(psi_tube)
tube_amplitude_ratio < 1 || push!(gatefails,
    "GATE 2: max|Ψ^ℛ_m| over the tube is not smaller than max|Ψ_m| ($(tube_amplitude_ratio))")
relspread(v) = (maximum(v) - minimum(v)) / maximum(v)
punc_R = Float64[]; punc_P = Float64[]
for el in hist.levels[end].punctures
    e = mesh.index[el]
    blk = mesh.blocks[mesh.blockindex[el.block]]
    em = element_matrices(Float64, blk, el, mesh.orders[e])
    a = Float64[]; b = Float64[]
    for x2 in em.x[2], x1 in em.x[1]
        push!(a, abs(field_at(sol, sol.u, x1, x2; element = e)))
        push!(b, abs(field_at(sol, full, x1, x2; element = e)))
    end
    push!(punc_R, relspread(a)); push!(punc_P, relspread(b))
end
puncture_spread_R = maximum(punc_R); puncture_spread_P = minimum(punc_P)
puncture_spread_R < puncture_spread_P || push!(gatefails,
    "GATE 2: |Ψ^ℛ_m| is not smoother than |Ψ_m| on the puncture elements " *
    "(relative spread $(puncture_spread_R) vs $(puncture_spread_P))")

# ───────────────────────────── registered values (the golden) ────────────────────────────────
lvl = hist.levels[end]
lm_block = [el.block for el in mesh.elements]
lm_rlev  = [el.segments[1].level for el in mesh.elements]
lm_ridx  = [el.segments[1].index for el in mesh.elements]
lm_zlev  = [el.segments[2].level for el in mesh.elements]
lm_zidx  = [el.segments[2].index for el in mesh.elements]
lm_ox    = [o[1] for o in mesh.orders]
lm_oy    = [o[2] for o in mesh.orders]
lm_wt    = [inside_worldtube(dom, el) ? 1 : 0 for el in mesh.elements]
idbytes  = string(collect(zip(lm_block, lm_rlev, lm_ridx, lm_zlev, lm_zidx)))
mapbytes = string(collect(zip(lm_block, lm_rlev, lm_ridx, lm_zlev, lm_zidx, lm_ox, lm_oy)))

@repro_result "spin" SPIN
@repro_result "orbital_radius" R0
@repro_result "m_mode" MM
@repro_result "amr_steps" ITER
@repro_result "grid_points_per_dim" PORD
@repro_result "initial_refinement_L" LREF
@repro_result "inner_radius_rplus" RIN
@repro_result "outer_radius" ROUT
@repro_result "block_bounds_r" collect(dom.bounds_r)
@repro_result "block_bounds_z2" collect(dom.bounds_z)
@repro_result "worldtube_blocks" collect(worldtube_blocks(dom))
@repro_result "worldtube_bounds_r" collect(WT_R)
@repro_result "worldtube_bounds_z2" collect(WT_W)
@repro_result "puncture_position" collect(xp)
@repro_result "rcc_inner_radius" float(cfg.observers.rcc_inner_radius)
@repro_result "rcc_outer_radius" float(cfg.observers.rcc_outer_radius)
@repro_result "rcc_compressed_outer_radius" rcc(cfg.observers.rcc_outer_radius)
@repro_result "level_n_elements" [length(l.mesh.elements) for l in hist.levels]
@repro_result "level_ndofs" [l.mesh.ndofs for l in hist.levels]
@repro_result "level_min_order" [minimum(first, l.mesh.orders) for l in hist.levels]
@repro_result "level_max_order" [maximum(first, l.mesh.orders) for l in hist.levels]
@repro_result "level_n_puncture_elements" [length(l.punctures) for l in hist.levels]
@repro_result "level_gmres_iterations" [l.record.iterations for l in hist.levels]
@repro_result "puncture_h_levels" puncture_hlevels
@repro_result "puncture_orders" puncture_orders
@repro_result "element_block" lm_block
@repro_result "element_r_level" lm_rlev
@repro_result "element_r_index" lm_ridx
@repro_result "element_z_level" lm_zlev
@repro_result "element_z_index" lm_zidx
@repro_result "element_order_x" lm_ox
@repro_result "element_order_y" lm_oy
@repro_result "element_inside_worldtube" lm_wt
@repro_result "element_id_sha256" bytes2hex(sha256(idbytes))
@repro_result "level_map_sha256" bytes2hex(sha256(mapbytes))
@repro_result "sample_r" samp_r
@repro_result "sample_z2" samp_w
@repro_result "sample_abs_psi_regular" samp_psiR
@repro_result "sample_abs_psi_retarded" samp_psi
@repro_result "abs_psi_regular_max" maximum(abs, sol.u)
@repro_result "abs_psi_regular_min" minimum(abs, sol.u)
@repro_result "abs_psi_retarded_max" maximum(abs, full)
@repro_result "tube_boundary_ratio_max" maximum(ratio_inner)
@repro_result "tube_amplitude_ratio" tube_amplitude_ratio
@repro_result "puncture_element_spread_regular" punc_R
@repro_result "puncture_element_spread_retarded" punc_P
@repro_result "psi_regular_at_xp" abs(hist.levels[end].psi)
@repro_result "dpsi_r_regular_at_xp" abs(hist.levels[end].dpsi_r)
@repro_result "two_to_one_balanced_all_levels" balance_ok
@repro_result "order_window_all_levels" order_window_ok
@repro_result "hp_order_rule_all_levels" order_rule_ok

# per-element sample lattice: the element's OWN collocation nodes plus its two faces, so the
# drawn surface passes through the nodal values and the faint grid lines ARE the Legendre-Gauss
# grid the paper describes (main.tex:566-570)
for e in eachindex(mesh.elements)
    el = mesh.elements[e]
    blk = mesh.blocks[mesh.blockindex[el.block]]
    em = element_matrices(Float64, blk, el, mesh.orders[e])
    lo1, hi1 = logical_interval(Float64, el, 1); lo2, hi2 = logical_interval(Float64, el, 2)
    rs = vcat(map_forward(blockmap(blk, 1), lo1), collect(em.x[1]),
              map_forward(blockmap(blk, 1), hi1))
    ws = vcat(map_forward(blockmap(blk, 2), lo2), collect(em.x[2]),
              map_forward(blockmap(blk, 2), hi2))
    v = [abs(field_at(sol, sol.u, rs[i], ws[j]; element = e)) for i in eachindex(rs), j in eachindex(ws)]
    push!(grids, ElemGrid(rs, ws, v))
end

    UMAX = maximum(abs, sol.u)
    on_horizon = [mesh.elements[e].segments[1].level == 0 && grids[e].r[1] ≈ dom.bounds_r[1]
                  for e in eachindex(mesh.elements)]
    inside_wt = [inside_worldtube(dom, el) for el in mesh.elements]
    if opts.mode === :regen
        open(SURFACE_CACHE, "w") do io
            serialize(io, (grids = grids, umax = UMAX, on_horizon = on_horizon,
                           inside_wt = inside_wt))
        end
        println("[fig_domain] wrote the surface cache ",
                relpath(SURFACE_CACHE, ReproHarness.REPRO_ROOT), " (",
                filesize(SURFACE_CACHE), " bytes, untracked) for --plot-only")
    end

else
    # --plot-only: the scene comes back from the cache; no solve, no comparison
    _c = open(deserialize, SURFACE_CACHE)
    grids = _c.grids; on_horizon = _c.on_horizon; inside_wt = _c.inside_wt; UMAX = _c.umax
    @printf("[fig_domain] --plot-only: %d element lattices from %s\n",
            length(grids), relpath(SURFACE_CACHE, ReproHarness.REPRO_ROOT))
end

# ───────────────────────────── the figure (presentation only) ────────────────────────────────
# RENDERING. The reference is a ParaView render annotated in Keynote; Plots.jl/GR has no 3D
# surface primitive that can carry (i) a per-quad colour, (ii) a hard field discontinuity at
# the worldtube boundary and (iii) white element edges on the surface. So this driver does its
# own axonometric PERSPECTIVE projection and paints the surface as depth-sorted quadrilaterals
# (painter's algorithm) in a 2-D `:shape` series, exactly the way the reference's own renderer
# would. Honest limitations, none of them structural:
#   * the white element edges and the faint collocation-grid lines are drawn AFTER all quads,
#     so a line on a far element is not hidden by a nearer fold of the surface (the CRP
#     surface folds only at the wavezone dip, where the effect is sub-pixel);
#   * lighting is a single Lambert term with an ambient floor — the reference's ParaView
#     default shading is what makes its high-|Ψ| slopes read dark brown rather than tan;
#   * the raster PNG is the tracked golden; the vector PDF is written ONLY on `--regen` /
#     `--plot-only` and is kept out of the index by `expected/.gitignore` beside it (the root
#     `.gitignore` is deliberately untouched).
# Every presentation constant below (camera, box aspect, height scale, light, colour anchors,
# annotation offsets) is a STYLE choice recorded here; none of them touches a registered value.
using Plots
using Plots: RGB, cgrad, Shape, text, bbox
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

const LX = 2.05                 # scene length in the drawn radial coordinate
const LY = 1.00                 # scene width in z²  (the reference's ribbon is ~2:1)
const LZ = 0.78                 # height of |Ψ|_max above the base plane
const CAM_AZ = deg2rad(225.0)   # camera azimuth: +r to the upper right, +z² to the upper left
const CAM_EL = deg2rad(30.0)    # low elevation, as the manual records
const CAM_D  = 2.9              # perspective distance (the reference tapers ~0.65 end to end)
const R_FRAME = 52.0            # right edge of the frame: the reference render itself stops at
                                # the u-block element boundary r = 48.68 M (see the header)
const PX_PER_UNIT = 272         # raster resolution of the frame; the golden-size knob
# Point sizes are absolute, data-coordinate offsets are not, so the ONE thing that has to
# follow PX_PER_UNIT is type and stroke: `SC` keeps the layout invariant under the resolution
# knob (429 px/unit is the size the annotation offsets were laid out at). Without it, lowering
# the resolution to meet the golden budget silently enlarges every label until they collide.
const SC = PX_PER_UNIT / 429
const AMBIENT = 0.62
const LIGHT = (-0.42, -0.28, 0.86)
const RED = RGB(0.90, 0.09, 0.09)          # the paper's reserved worldtube accent
# colour anchors sampled from the reference render (`ssf_domain.pdf` at 150 dpi): deep navy →
# blue → teal → sage → khaki → tan → orange, the manual's "terrain-like blue→teal→sand→brown"
# ramp. The last anchor carries it on into orange (the reference's warmest pixels measure
# (198,168,108) and (186,144,78)); the maroon of the published render is that orange shaded.
const ANCHORS = [(0.031,0.031,0.345), (0.125,0.157,0.408), (0.188,0.282,0.502),
                 (0.188,0.376,0.533), (0.282,0.565,0.659), (0.408,0.627,0.627),
                 (0.502,0.659,0.627), (0.565,0.659,0.565), (0.627,0.627,0.471),
                 (0.706,0.612,0.400), (0.776,0.612,0.353), (0.800,0.565,0.259),
                 (0.831,0.494,0.180)]

sx3(r) = LX * (rcc(r) - RC_LO) / (RC_HI - RC_LO)
sz3(v) = LZ * v / UMAX
scene3(r, w, v) = (sx3(r), LY * w, sz3(v))
const CF = (-cos(CAM_AZ)*cos(CAM_EL), -sin(CAM_AZ)*cos(CAM_EL), -sin(CAM_EL))
const CRv = (-sin(CAM_AZ), cos(CAM_AZ), 0.0)
const CUv = (-cos(CAM_AZ)*sin(CAM_EL), -sin(CAM_AZ)*sin(CAM_EL), cos(CAM_EL))
const CCTR = (LX/2, LY/2, LZ/3)
const CPOS = (CCTR[1]-CAM_D*CF[1], CCTR[2]-CAM_D*CF[2], CCTR[3]-CAM_D*CF[3])
dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
function projs(p)
    v = (p[1]-CPOS[1], p[2]-CPOS[2], p[3]-CPOS[3]); t = dot3(v, CF); f = CAM_D / t
    return (dot3(v, CRv)*f, dot3(v, CUv)*f, t)
end
P2(r, w, v) = (p = projs(scene3(r, w, v)); (p[1], p[2]))

function cramp(t)
    t = clamp(t, 0.0, 1.0); s = t*(length(ANCHORS)-1)
    i = clamp(floor(Int, s)+1, 1, length(ANCHORS)-1); f = s-(i-1)
    a = ANCHORS[i]; b = ANCHORS[i+1]
    return (a[1]+(b[1]-a[1])*f, a[2]+(b[2]-a[2])*f, a[3]+(b[3]-a[3])*f)
end
# one categorical palette carrying (value, shade) so per-quad shaded colour goes through the
# single `fill_z` channel of ONE depth-ordered series
const NV = 64; const NS = 12
const PAL = [begin
                 c = cramp((iv-1)/(NV-1)); f = AMBIENT + (1-AMBIENT)*(is-1)/(NS-1)
                 RGB(clamp(c[1]*f,0,1), clamp(c[2]*f,0,1), clamp(c[3]*f,0,1))
             end for iv in 1:NV for is in 1:NS]
const PGRAD = cgrad(PAL; categorical = true)
palidx(v, s) = clamp(round(Int, v*(NV-1)), 0, NV-1)*NS + clamp(round(Int, s*(NS-1)), 0, NS-1) + 1


shp = Shape[]; dep = Float64[]; zz = Float64[]
for gd in grids
    nx, ny = length(gd.r), length(gd.w)
    for i in 1:nx-1, j in 1:ny-1
        c = ((gd.r[i],gd.w[j],gd.v[i,j]), (gd.r[i+1],gd.w[j],gd.v[i+1,j]),
             (gd.r[i+1],gd.w[j+1],gd.v[i+1,j+1]), (gd.r[i],gd.w[j+1],gd.v[i,j+1]))
        sc = map(q -> scene3(q...), c); pr = map(projs, sc)
        a = sc[2] .- sc[1]; b = sc[4] .- sc[1]
        n = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]); n[3] < 0 && (n = .-n)
        sh = max(0.0, dot3(n, LIGHT) / sqrt(dot3(n,n)*dot3(LIGHT,LIGHT)))
        push!(shp, Shape([p[1] for p in pr], [p[2] for p in pr]))
        push!(dep, sum(p -> p[3], pr)/4)
        push!(zz, palidx(sum(q -> q[3], c)/4 / UMAX, sh))
    end
end
ord = sortperm(dep; rev = true)
@printf("[fig_domain] surface: %d elements → %d depth-sorted quads\n", length(grids), length(shp))

function gridlines(gs)
    gx = Float64[]; gy = Float64[]
    for gd in gs
        nx, ny = length(gd.r), length(gd.w)
        for j in 2:ny-1
            for i in 1:nx; p = P2(gd.r[i], gd.w[j], gd.v[i,j]); push!(gx,p[1]); push!(gy,p[2]); end
            push!(gx, NaN); push!(gy, NaN)
        end
        for i in 2:nx-1
            for j in 1:ny; p = P2(gd.r[i], gd.w[j], gd.v[i,j]); push!(gx,p[1]); push!(gy,p[2]); end
            push!(gx, NaN); push!(gy, NaN)
        end
    end
    return gx, gy
end
function edgelines(gs)
    bx = Float64[]; by = Float64[]
    for gd in gs
        nx, ny = length(gd.r), length(gd.w)
        for (ii, jj) in ((1:nx, 1:1), (1:nx, ny:ny), (1:1, 1:ny), (nx:nx, 1:ny))
            for i in ii, j in jj; p = P2(gd.r[i], gd.w[j], gd.v[i,j]); push!(bx,p[1]); push!(by,p[2]); end
            push!(bx, NaN); push!(by, NaN)
        end
    end
    return bx, by
end

plt = plot(shp[ord]; fill_z = permutedims(zz[ord]), c = PGRAD, clims = (1, NV*NS),
           linewidth = 0, linealpha = 0, legend = false, colorbar = false, axis = false,
           grid = false, ticks = false, framestyle = :none, background_color = :white,
           aspect_ratio = :equal, dpi = 100)
gx, gy = gridlines(grids)
plot!(plt, gx, gy; linecolor = RGB(1,1,1), linealpha = 0.34, linewidth = 0.5SC, label = "")
bx, by = edgelines(grids)
plot!(plt, bx, by; linecolor = :white, linewidth = 1.6SC, label = "")

# the Kerr horizon edge r = r₊, drawn pole → equator so the rotated label reads left to right
hedge = sort([(e, grids[e]) for e in eachindex(grids) if on_horizon[e]]; by = t -> -t[2].w[1])
hx = Float64[]; hy = Float64[]
for (_, gd) in hedge, j in length(gd.w):-1:1
    p = P2(gd.r[1], gd.w[j], gd.v[1,j]); push!(hx, p[1]); push!(hy, p[2])
end
plot!(plt, hx, hy; linecolor = :black, linewidth = 5.0SC, label = "")
hang = atand(hy[end]-hy[1], hx[end]-hx[1])
annotate!(plt, (hx[1]+hx[end])/2-0.115, (hy[1]+hy[end])/2-0.105,
          text("Kerr horizon \$r_+\$", 14SC, :black, :center, rotation = hang))

# the v | t | u slicing arc along the pole edge (a quadratic Bézier above the surface)
pv = P2(dom.bounds_r[1], 1.0, 0.0); pu = P2(38.0, 1.0, 0.0)
A0 = (pv[1]-0.05, pv[2]+0.06); A2 = (pu[1], pu[2]+0.30)
A1 = (A0[1]-0.09, (A0[2]+A2[2])/2+0.10)
bez(t) = ((1-t)^2*A0[1]+2(1-t)*t*A1[1]+t^2*A2[1], (1-t)^2*A0[2]+2(1-t)*t*A1[2]+t^2*A2[2])
BP = bez.(range(0, 1; length = 80))
plot!(plt, [p[1] for p in BP], [p[2] for p in BP]; linecolor = :black, linewidth = 3.2SC,
      label = "", arrow = Plots.arrow(:closed, :head, 1.0, 1.0))
for (t, lab, hl) in ((0.02, "v", "\$H=-1\$"), (0.5, "t", ""), (0.97, "u", "\$H=+1\$"))
    p = bez(t)
    annotate!(plt, p[1]-0.055, p[2]+0.045, text("\$\\mathbf{$(lab)}\$", 18SC, :black, :center))
    isempty(hl) || annotate!(plt, p[1]-0.075, p[2]+0.175, text(hl, 11SC, :black, :right))
end

# the (z², r) direction arrows under the near corner
c0 = P2(dom.bounds_r[1], 0.0, 0.0); cz = P2(dom.bounds_r[1], 0.34, 0.0); cr = P2(7.4, 0.0, 0.0)
aoy = -0.26
plot!(plt, [c0[1], cz[1]], [c0[2]+aoy, cz[2]+aoy]; linecolor = :black, linewidth = 2.4SC,
      label = "", arrow = Plots.arrow(:closed, :both, 0.8, 0.8))
plot!(plt, [c0[1], cr[1]], [c0[2]+aoy, cr[2]+aoy]; linecolor = :black, linewidth = 2.4SC,
      label = "", arrow = Plots.arrow(:closed, :head, 0.8, 0.8))
annotate!(plt, cz[1]-0.03, cz[2]+aoy-0.10, text("\$z^2 \\equiv \\cos(\\theta)^2\$", 13SC, :right))
annotate!(plt, cr[1]+0.05, cr[2]+aoy-0.02, text("\$r\$", 14SC, :left))

# the red worldtube rectangle, its label, and x_p
zrect = 0.10 * UMAX
rcorn = [P2(WT_R[1],WT_W[1],zrect), P2(WT_R[2],WT_W[1],zrect),
         P2(WT_R[2],WT_W[2],zrect), P2(WT_R[1],WT_W[2],zrect)]
plot!(plt, [p[1] for p in rcorn[[1,2,3,4,1]]], [p[2] for p in rcorn[[1,2,3,4,1]]];
      linecolor = RED, linewidth = 4.5SC, label = "")
# the label is placed RELATIVE to the rectangle's top edge (offsets measured on the render),
# so it stays in the white band between the surface and the rectangle whatever the frame does
tmid = ((rcorn[3][1] + rcorn[4][1]) / 2, (rcorn[3][2] + rcorn[4][2]) / 2)
annotate!(plt, tmid[1] + 0.52, tmid[2] + 0.17,
          text("Regularized field \$\\Psi^{R}_{m} = \\Psi_{m} - \\Psi^{P}_{m}\$", 14SC, RED, :center))
epx = findfirst(e -> inside_wt[e] && grids[e].r[end] ≈ R0 && grids[e].w[1] ≈ 0.0,
                eachindex(grids))
pxp = P2(R0, 0.0, grids[epx].v[end, 1])
scatter!(plt, [pxp[1]], [pxp[2]]; markercolor = :black, markersize = 6.5SC,
         markerstrokewidth = 0, label = "")
annotate!(plt, pxp[1]+0.045, pxp[2]-0.045, text("\$\\mathbf{x}_p\$", 15SC, :black, :left))

# the inset: the puncture cascade, magnified (same camera, a window around x_p)
const IN_R = (9.0, 11.4); const IN_W = (0.0, 0.055)
insel = [grids[e] for e in eachindex(grids)
         if grids[e].r[1] < IN_R[2] && grids[e].r[end] > IN_R[1] && grids[e].w[1] < IN_W[2]]
icor = [P2(r, w, 0.006) for r in IN_R for w in IN_W]
ixl = (minimum(p->p[1], icor), maximum(p->p[1], icor))
iyl = (minimum(p->p[2], icor), maximum(p->p[2], icor))
imx = 0.02*(ixl[2]-ixl[1])
plot!(plt; inset = (1, bbox(0.70, 0.05, 0.29, 0.34, :left, :bottom)), subplot = 2)
ishp = Shape[]; idep = Float64[]; izz = Float64[]
for gd in insel
    nx, ny = length(gd.r), length(gd.w)
    for i in 1:nx-1, j in 1:ny-1
        c = ((gd.r[i],gd.w[j],gd.v[i,j]), (gd.r[i+1],gd.w[j],gd.v[i+1,j]),
             (gd.r[i+1],gd.w[j+1],gd.v[i+1,j+1]), (gd.r[i],gd.w[j+1],gd.v[i,j+1]))
        pr = map(q -> projs(scene3(q...)), c)
        push!(ishp, Shape([p[1] for p in pr], [p[2] for p in pr]))
        push!(idep, sum(p->p[3], pr)/4); push!(izz, palidx(sum(q->q[3],c)/4/UMAX, 1.0))
    end
end
iord = sortperm(idep; rev = true)
plot!(plt, ishp[iord]; subplot = 2, fill_z = permutedims(izz[iord]), c = PGRAD,
      clims = (1, NV*NS), linewidth = 0, linealpha = 0, legend = false, colorbar = false,
      ticks = false, grid = false, framestyle = :box, foreground_color_border = :black,
      background_color_inside = :white, aspect_ratio = :equal, label = "",
      xlims = (ixl[1]-imx, ixl[2]+imx), ylims = (iyl[1]-imx, iyl[2]+imx))
igx, igy = gridlines(insel)
plot!(plt, igx, igy; subplot = 2, linecolor = RGB(1,1,1), linealpha = 0.5, linewidth = 0.7SC, label = "")
ibx, iby = edgelines(insel)
plot!(plt, ibx, iby; subplot = 2, linecolor = :white, linewidth = 2.6SC, label = "")

allx = Float64[]; ally = Float64[]
for s in shp; append!(allx, s.x); append!(ally, s.y); end
append!(allx, filter(isfinite, bx)); append!(ally, filter(isfinite, by))
append!(allx, [p[1] for p in BP]); append!(ally, [p[2] for p in BP])
append!(allx, [c0[1], cz[1], cr[1]]); append!(ally, [c0[2]+aoy, cz[2]+aoy, cr[2]+aoy])
xcut = P2(R_FRAME, 0.0, 0.0)[1]
fx0 = minimum(allx) - 0.34; fx1 = xcut + 0.62   # left margin: the H = −1 label; right: the inset
fy0 = minimum(ally) - 0.09; fy1 = maximum(ally) + 0.14
# the canvas is DERIVED from the frame, so `aspect_ratio = :equal` letterboxes nothing and the
# whole raster carries picture. PX_PER_UNIT is then the one knob for the golden's size budget.
plot!(plt; subplot = 1, xlims = (fx0, fx1), ylims = (fy0, fy1),
      size = (round(Int, PX_PER_UNIT * (fx1 - fx0)), round(Int, PX_PER_UNIT * (fy1 - fy0))))

# `--check` must never touch the tracked artifact: it renders to a SCRATCH path so the
# plotting path is still exercised, measures it against the golden budget, and leaves the
# committed PNG alone. Only `--regen`/`--plot-only` publish. The size verdict is folded into
# the gate list BEFORE `gates_pass` and `conclude!`, so a figure that outgrows the budget is a
# FAIL row with a full comparison table, never a stacktrace.
const GOLD_PNG = joinpath(@__DIR__, "expected", "fig_domain.png")
publish = opts.mode === :regen || PLOT_ONLY
out_png = publish ? GOLD_PNG : joinpath(mktempdir(), "fig_domain.png")
saved = save_figure!(plt, out_png; max_bytes = 200 * 1024, pdf = publish, id = "fig_domain")
saved.ok || push!(gatefails, "GATE 6: " * saved.message)
publish || println("[fig_domain] --check: rendered to a scratch path, ",
                   relpath(GOLD_PNG, ReproHarness.REPRO_ROOT), " untouched")

if PLOT_ONLY
    println("[fig_domain] --plot-only: redrawn from the surface cache, no physics recomputed")
    isempty(gatefails) || (println("\n[fig_domain] FINDING:"); foreach(f -> println("    ", f), gatefails))
    exit(isempty(gatefails) ? 0 : 1)
end

# ───────────────────────────── the node's own gates ──────────────────────────────────────────
@repro_result "png_bytes_within_budget" saved.ok
@repro_result "gates_pass" isempty(gatefails)
if !isempty(gatefails)
    println("\n[fig_domain] FINDING — node gates violated:")
    foreach(f -> println("    ", f), gatefails)
end

code = ReproHarness.conclude!(opts; paper = "2606.04998", id = "fig_domain",
                              script = @__FILE__, rtol = 1e-10)
exit(isempty(gatefails) ? code : max(code, 1))
