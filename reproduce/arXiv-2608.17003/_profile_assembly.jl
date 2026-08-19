# reproduce/arXiv-2608.17003/_profile_assembly.jl — PRESENTATION-SIDE assembly of the radial
# profiles that Figs. 3 and 4 of arXiv:2608.17003 plot.  Included by `fig03.jl` and `fig04.jl`
# (one file, so the assembly is never written twice).
#
# ZERO PHYSICS.  Every physical quantity is produced by `AccretionDisc.SpiralWaves`:
#
#   SW.build_grid          the App.-C adaptive grid of one azimuthal mode
#   SW.solve_mode          the implicit-midpoint solve of the master system (47)/(64)
#   SW.amf                 F_{J,m}(r), Eq. (72a), Fourier weight (73) already inside
#   SW.torque_density      √−g⟨𝒯⟩_m(r), the integrand of (72b), same weight convention
#   SW.cumulative_torque   𝔗_m(r) = ∫ dr √−g⟨𝒯⟩_m, Eq. (72b)
#   NumericsCore.CubicSplineInterpolant   the interpolant `SpiralWaves` itself uses to put modes
#                                         solved on DIFFERENT adaptive grids onto a common one
#                                         (see SpiralWaves/image_synthesis.jl)
#
# What this file adds is bookkeeping only:
#
#  1. the m-loop (threaded; each task owns its grid/solution and writes disjoint slots, and the
#     reduction over m is SERIAL in a fixed order, so the result is independent of the thread count);
#  2. the reduction of each mode to CELL AVERAGES on a common plotting grid of edges x_j:
#
#         ⟨d𝔗/dr⟩_j = [𝔗_m(x_{j+1}) − 𝔗_m(x_j)] / (x_{j+1} − x_j),
#         ⟨dF_J/dr⟩_j = [F_{J,m}(x_{j+1}) − F_{J,m}(x_j)] / (x_{j+1} − x_j),
#
#     i.e. the same conservative discrete derivative that `SpiralWaves.balance_residual` forms on the
#     solver's own cells, evaluated instead on the FIGURE's cells so that the two curves of Fig. 3 are
#     compared at the resolution at which they are drawn.  The corotation deposition therefore appears
#     in ⟨dF_J/dr⟩ as a single-cell spike whose HEIGHT is set by the width of the cell that contains
#     r_c (it is a Landau-regularised near-delta: its plotted height is a property of the plotting
#     grid, in the published figure exactly as here) — only its position and its absence from
#     ⟨d𝔗/dr⟩ are structural.
#
#  3. the mode sum of (73): the per-mode weights (2π for m = 0, 4π for m ≥ 1) are already inside
#     `amf`/`torque_density`, so the sum over the solved m ≥ 1 modes is a plain sum.
#
# TRUNCATION CONVENTION (explicit, and the only modelling statement in this file): outside its own
# adaptive domain [r_in, r_out] a mode contributes NOTHING to either cell average.  This is realised
# by clamping the evaluation points into [r_in, r_out], which makes 𝔗_m and F_{J,m} constant there —
# the exact continuation of the paper's own model, whose domain criterion (C2)-(C3) places both edges
# in the free-wave zone where the torque density has died away and F_{J,m} is constant.  The residual
# torque density at the edges is measured and returned as `edge_density_fraction`; a large value
# would mean the domain was cut too early and the assembly is reported as unconverged.
#
# [FUTURE] the natural home of `profile_sum` is `AccretionDisc.SpiralWaves.diagnostics` (it is the
# mode-summed form of `balance_residual`); it is kept here while the figure lane is separate from the
# package lane, and is a candidate for promotion with its own package tests.

module ProfileAssembly

using GrTools
using GrTools.NumericsCore: CubicSplineInterpolant

const SW = GrTools.AccretionDisc.SpiralWaves

export profile_sum, graded_edges, centered_edges, cell_centers

"""
    graded_edges(blocks) -> Vector{Float64}

Concatenate the half-open blocks `((lo, hi, step), ...)` into one strictly increasing vector of cell
edges (the last block contributes its upper end).  Used to refine the plotting grid around
corotation without inflating the golden: away from r_c the profiles are smooth and a coarse cell is
indistinguishable from a fine one, at r_c the cell width sets the visible height of the deposition
spike.
"""
function graded_edges(blocks)
    xs = Float64[]
    for (k, (lo, hi, step)) in enumerate(blocks)
        n = round(Int, (hi - lo) / step)
        n ≥ 1 || error("graded_edges: block $(k) has no cells")
        for i in 0:(n - 1)
            push!(xs, lo + i * step)
        end
        k == length(blocks) && push!(xs, hi)
    end
    all(diff(xs) .> 0) || error("graded_edges: the assembled edges are not strictly increasing")
    xs
end

"""
    centered_edges(x0, lo, hi, step) -> Vector{Float64}

Uniform cell edges of width `step` covering `[lo, hi]`, ANCHORED so that `x0` is a cell CENTRE
(edges at `x0 ± step/2`).  Corotation is a near-delta in the AMF density, so where it falls relative
to the cells decides whether the published single-point spike is rendered as one cell or split
between two; anchoring the plotting grid on `x0 = r_c` samples the profile AT corotation, which is
what the published figures show.  The realised range is `[lo, hi]` rounded outward to whole cells.
"""
function centered_edges(x0::Real, lo::Real, hi::Real, step::Real)
    step > 0 || throw(ArgumentError("centered_edges: need step > 0"))
    lo < x0 < hi || throw(ArgumentError("centered_edges: need lo < x0 < hi"))
    kmin = floor(Int, (lo - x0) / step - 0.5)
    kmax = ceil(Int, (hi - x0) / step - 0.5)
    xs = [x0 + (k + 0.5) * step for k in kmin:kmax]
    all(diff(xs) .> 0) || error("centered_edges: the assembled edges are not strictly increasing")
    xs
end

"Cell centres of an edge vector."
cell_centers(edges::AbstractVector) = [(edges[i] + edges[i+1]) / 2 for i in 1:(length(edges) - 1)]

# one mode: solve and reduce to the two cell-average vectors (pure — same inputs, same bits)
function _one_mode(bg, m::Integer, redges::Vector{Float64}, grid_kwargs, solve_kwargs)
    t0 = time()
    g = SW.build_grid(bg, m; grid_kwargs...)
    sol = SW.solve_mode(bg, m; grid = g, threads = false, solve_kwargs...)
    F = SW.amf(bg, sol)
    D = SW.torque_density(bg, sol)
    𝔗 = SW.cumulative_torque(bg, sol; density = D)
    splF = CubicSplineInterpolant(sol.r, F)
    splT = CubicSplineInterpolant(sol.r, 𝔗)
    lo, hi = first(sol.r), last(sol.r)
    n = length(redges) - 1
    Fv = [splF(clamp(x, lo, hi), 1) for x in redges]
    Tv = [splT(clamp(x, lo, hi), 1) for x in redges]
    dF = Vector{Float64}(undef, n); dT = Vector{Float64}(undef, n)
    @inbounds for j in 1:n
        w = redges[j+1] - redges[j]
        dF[j] = (Fv[j+1] - Fv[j]) / w
        dT[j] = (Tv[j+1] - Tv[j]) / w
    end
    dmax = maximum(abs, D)
    # only the edges the WKB criterion (C2)-(C3) MOVED are truncations; the fiducial edges (C1)
    # r_ISCO and 5a• are physical boundaries of the model and are excluded from the check
    fid_in, fid_out = SW.r_isco(bg), 5 * bg.aEMRI
    cut_in = g.r_in > fid_in * (1 + 1e-10)
    cut_out = g.r_out < fid_out * (1 - 1e-10)
    edge_frac = iszero(dmax) ? 0.0 :
        max(cut_in ? abs(D[1]) / dmax : 0.0, cut_out ? abs(D[end]) / dmax : 0.0)
    (dT = dT, dF = dF, T = SW.total_torque(bg, sol; density = D),
     T_L = SW.lindblad_torque(bg, sol; flux = F), nodes = length(sol.r),
     residual = sol.residual, edge_frac = edge_frac, peak = maximum(abs, dT),
     seconds = time() - t0)
end

"""
    profile_sum(bg, ms, redges; grid_kwargs = (;), solve_kwargs = (;), threads = true,
                verbose = true) -> NamedTuple

Solve every azimuthal mode in `ms` on its own App.-C grid and return the mode sum (73) of the two
cell-averaged radial derivatives on the cells of `redges` (RADII, not r/a•):

  * `dT_dr`   Σ_m ⟨d𝔗/dr⟩   — the torque density of (72b), the blue curve of Fig. 3 and Fig. 4;
  * `dF_dr`   Σ_m ⟨dF_J/dr⟩ — the AMF density of (72a), the red curve of Fig. 3;
  * `T_total`, `T_L_total` — the mode sums of (74) and (75) on the same modes (fixed-order);
  * `per_mode_peak` — max_j |⟨d𝔗/dr⟩| of each mode, the mode-sum convergence diagnostic;
  * `edge_density_fraction` — the largest |√−g⟨𝒯⟩| at a WKB-CUT domain edge relative to that mode's
    own maximum (the truncation-convention check of the file header; the fiducial edges r_ISCO and
    5a• of (C1) are physical model boundaries and do not enter);
  * `nodes`, `residual_max`, `seconds`, `budget`.

Deterministic: the modes are reduced serially in the order of `ms`, so the result does not depend on
the thread count.
"""
function profile_sum(bg, ms::AbstractVector{<:Integer}, redges::AbstractVector{<:Real};
                     grid_kwargs = NamedTuple(), solve_kwargs = NamedTuple(),
                     threads::Bool = true, verbose::Bool = true)
    isempty(ms) && throw(ArgumentError("profile_sum: `ms` must be non-empty"))
    re = collect(Float64, redges)
    n = length(re) - 1
    budget = SW.compute_budget()
    verbose && println("[profile_sum] ", budget.note)
    t0 = time()
    out = Vector{Any}(undef, length(ms))
    if threads
        Threads.@threads :static for i in eachindex(ms)
            @inbounds out[i] = _one_mode(bg, ms[i], re, grid_kwargs, solve_kwargs)
        end
    else
        for i in eachindex(ms)
            @inbounds out[i] = _one_mode(bg, ms[i], re, grid_kwargs, solve_kwargs)
        end
    end
    dT = zeros(Float64, n); dF = zeros(Float64, n)
    Ttot = 0.0; TLtot = 0.0; nodes = 0; resmax = 0.0; edgemax = 0.0
    peaks = Float64[]
    for o in out                                   # FIXED-ORDER serial reduction
        dT .+= o.dT
        dF .+= o.dF
        Ttot += o.T; TLtot += o.T_L
        nodes += o.nodes
        resmax = max(resmax, o.residual)
        edgemax = max(edgemax, o.edge_frac)
        push!(peaks, o.peak)
    end
    (dT_dr = dT, dF_dr = dF, T_total = Ttot, T_L_total = TLtot, per_mode_peak = peaks,
     edge_density_fraction = edgemax, nodes = nodes, residual_max = resmax,
     seconds = time() - t0, budget = budget)
end

end # module ProfileAssembly
