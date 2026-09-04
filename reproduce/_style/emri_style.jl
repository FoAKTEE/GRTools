# Shared plot style for EXACT figure reproduction (phys-agentic-loop/_common/figure_reproduction_criterion.md).
# The PHYSICS DATA is always computed by the migrated Julia packages/* — this module only governs PRESENTATION.
#
# Per-paper colormaps/fonts are EXTRACTED FROM THE VLM SPEC of each published figure (Rule 1) and added here as
# named constants when first needed — do NOT guess them up front (AHA: a wrong colormap is worse than none).
# The emri-disk papers render with LaTeX (Hirata 2011 = eps; Duque&Sberna 2025 / Dyson&D'Orazio 2026 = matplotlib),
# so the default font is serif/Computer Modern. matplotlib colormap constants (viridis/inferno/RdBu/...) are added
# per-figure from the spec.
using Plots
using Plots: RGB, cgrad

# LaTeX-paper defaults: serif/CM font, boxed frame, no grid unless the spec shows one.
function emri_defaults!()
    default(fontfamily = "Computer Modern", grid = false, framestyle = :box,
            titlefontsize = 14, guidefontsize = 13, tickfontsize = 11,
            colorbar_titlefontsize = 12, foreground_color_legend = nothing)
    return nothing
end

# ── Per-paper colormap registry (populated from VLM specs in PHASE 1) ───────────────────────────────────────
# Example shape (fill from the published figure's spec, e.g. matplotlib "viridis" sampled at 17 points):
#   const VIRIDIS = cgrad([RGB(r,g,b) for (r,g,b) in _VIRIDIS_RGB])
# (none yet — added as figures are speced)

# ── figure output: one place for savefig + the golden-size budget ────────────────────────────
"""
    save_figure!(plt, png; max_bytes = 200*1024, pdf = false, id = "figure", io = stdout)
        -> (; path, bytes, pdf_path, pdf_bytes, ok, message)

Write `plt` to `png` (creating its directory), optionally also to the sibling `.pdf`, and
REPORT the raster size against the reproduction budget instead of throwing.

The budget is the reproduce/ rule for a tracked golden artifact: a figure that quietly grows
is a repository regression, but it is a COSMETIC one, so it must reach the driver's own gate
list and become a FAIL row — not a stacktrace that aborts before `conclude!` and leaves the
run with no comparison table at all. Hence `ok`/`message` rather than an `error`: the caller
pushes `message` onto its gate list when `ok` is false.

`pdf` is opt-in because a vector copy of a dense figure is large and is never the golden;
drivers write it on `--regen`/`--plot-only` only, and keep it out of the index with a local
`expected/.gitignore`.
"""
function save_figure!(plt, png::AbstractString; max_bytes::Integer = 200 * 1024,
                      pdf::Bool = false, id::AbstractString = "figure", io::IO = stdout)
    mkpath(dirname(png))
    savefig(plt, png)
    nb = filesize(png)
    pdf_path = replace(png, r"\.png$" => ".pdf")
    pdf_bytes = 0
    if pdf
        savefig(plt, pdf_path)
        pdf_bytes = filesize(pdf_path)
    end
    ok = nb <= max_bytes
    msg = ok ? "" :
        "$(id): $(basename(png)) is $(nb) bytes > $(max_bytes) — reduce the canvas, the dpi " *
        "or the drawn detail (the tracked golden must stay inside the reproduce/ budget)"
    println(io, "[", id, "] wrote ", basename(png), " (", nb, " bytes",
            ok ? "" : ", OVER BUDGET", ")",
            pdf ? " and $(basename(pdf_path)) ($(pdf_bytes) bytes, untracked)" : "")
    return (path = png, bytes = nb, pdf_path = pdf_path, pdf_bytes = pdf_bytes,
            ok = ok, message = msg)
end

# ── log-axis decade ticks: one definition for every driver ───────────────────────────────────
"""
    decade_ticks(lo, hi; step = 1, anchor = 0, inside = true) -> Vector{Float64}

Powers of ten to label a log axis spanning `[lo, hi]`: one tick every `step` decades, the
sequence anchored so that exponent `anchor` is always on it (`step = 3, anchor = -1` gives
1e-1, 1e-4, 1e-7, … whatever the data range).

`inside = true` (the default) keeps only ticks that lie INSIDE `[lo, hi]`. That clamp is
load-bearing: GR draws an out-of-range tick pinned to the axis end, which on a shared-x panel
stack shows up as a stray label row between the panels. `inside = false` instead widens to the
enclosing decades and then one decade further on each side — the right choice when the caller
sets its own `ylims` wider than the data and wants the ticks to reach them.

Hoisted here from two drivers that had grown the same function twice
(`arXiv-2606.04998/fig_linsolv.jl`, `arXiv-2606.04998/fig_convergence.jl`); both call sites are
bitwise unchanged.
"""
function decade_ticks(lo, hi; step::Integer = 1, anchor::Integer = 0, inside::Bool = true)
    e0, e1 = inside ? (ceil(Int, log10(lo)), floor(Int, log10(hi))) :
                      (floor(Int, log10(lo)) - 1, ceil(Int, log10(hi)) + 1)
    return [10.0^e for e in e0:e1 if mod(e - anchor, step) == 0]
end
