# reproduce/_style — presentation layer for reproduction figures

Physics data always comes from the GrTools capability packages; the files here govern
**presentation only** (fonts, colormaps, frame style), so that a reproduced figure can be
VLM-compared against the published one.

| file | origin | use |
|---|---|---|
| `emri_style.jl` | emri-disk `results/emri-disk/figures/emri_style.jl` (header kept) | `emri_defaults!()` — serif/Computer Modern, boxed frame, no grid; matplotlib-paper defaults (Plots.jl/GR) |
| `mma_style.jl` | emri-disk `results/gw-flux/figures/mma_style.jl` | exact Mathematica `ColorData["Rainbow"]` / `["TemperatureMap"]` gradients for the Mathematica-rendered gw-flux papers |
| `manuals/arXiv-2606.04998.md` | emri-disk `progress/elliptic-ssf/FIGURE_STYLE_MANUAL.md` | VLM-inspected per-figure style manual (Vu et al., scalar self-force with NR methods) |
| `manuals/arXiv-2608.17003.md` | emri-disk `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md` | VLM-inspected per-figure style manual (Dittmann, Hegade K. R., Fairbairn, EMRIs in gaseous disks) |

The manuals are **copies**: the canonical, evolving versions live in the emri-disk mission
repository under `progress/<mission>/FIGURE_STYLE_MANUAL.md` (elliptic-ssf → 2606.04998,
emri-gasdisk → 2608.17003). Re-copy them here whenever the upstream manual changes
(compare with `diff`); a driver that produces a figure must read the relevant manual
section first and then follow the recurring VLM procedure: (1) read the manual, (2) VLM-check
the new figure, (3) VLM-compare new vs reference and decide whether the plotting script needs
another iteration.

Usage in a driver:

```julia
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))   # defines emri_defaults!()
emri_defaults!()
```

Per-paper colormap constants are added to these files only when a VLM spec of the
published figure demands them (never guessed up front).
