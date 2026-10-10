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
| `manuals/arXiv-2207.14637.md` | authored here (a559b4e); indexed by emri-disk `progress/emri-eccentric/FIGURE_STYLE_MANUAL.md` | VLM-inspected per-figure style manual (Fairbairn & Rafikov 2022, FR22, density waves excited by eccentric planets, linear theory) |
| `manuals/arXiv-2407.20398.md` | authored here (a559b4e); indexed by emri-disk `progress/emri-eccentric/FIGURE_STYLE_MANUAL.md` | VLM-inspected per-figure style manual (Fairbairn & Rafikov 2024, FR24 arXiv v1, eccentric planet-disc interactions: migration and eccentricity evolution) |
| `manuals/arXiv-2506.19917.md` | authored here (a559b4e); indexed by emri-disk `progress/emri-eccentric/FIGURE_STYLE_MANUAL.md` | VLM-inspected per-figure style manual (Fairbairn & Dittmann 2025, FD25, pushing the limits of eccentricity in planet-disc interactions) |

The 2606.04998 and 2608.17003 manuals are **copies**: the canonical, evolving versions live in
the emri-disk mission repository under `progress/<mission>/FIGURE_STYLE_MANUAL.md` (elliptic-ssf →
2606.04998, emri-gasdisk → 2608.17003). Re-copy them here whenever the upstream manual changes
(compare with `diff`). The three eccentric planet-disc manuals (2207.14637, 2407.20398,
2506.19917) are canonical here; the emri-disk file `progress/emri-eccentric/FIGURE_STYLE_MANUAL.md`
is only their index plus the cross-paper conventions. A driver that produces a figure must read
the relevant manual section first and then follow the recurring VLM procedure: (1) read the
manual, (2) VLM-check the new figure, (3) VLM-compare new vs reference and decide whether the
plotting script needs another iteration.

Usage in a driver:

```julia
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))   # defines emri_defaults!()
emri_defaults!()
```

Per-paper colormap constants are added to these files only when a VLM spec of the
published figure demands them (never guessed up front).
