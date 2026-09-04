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
