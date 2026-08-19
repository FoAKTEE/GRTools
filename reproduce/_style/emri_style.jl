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
