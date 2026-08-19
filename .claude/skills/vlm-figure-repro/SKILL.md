---
name: vlm-figure-repro
description: 'VLM-verified figure reproduction workflow: read the per-paper FIGURE_STYLE_MANUAL, generate the figure, VLM-check the new visualization, then VLM-compare new vs reference and decide whether the plotting script must be updated; iterate. Use whenever producing or checking any reproduction figure (arXiv:2606.04998 elliptic-ssf, arXiv:2608.17003 emri-gasdisk, or any later ref-paper).'
---

# VLM-verified figure reproduction (recurring skill)

Recurring 3-step procedure for EVERY reproduction figure in this mission,
mandated by the mission addon prompt (`progress/prompt/emri-disk-addon.md`,
"Visualization"): preserve the presentation style of the reference paper;
(1) read the manual, (2) VLM-check the new visualization, (3) VLM-compare
the new figure against the reference and decide whether the plotting
script needs updating — make a new figure if needed. Never generate or sign
off a figure without completing all three steps.

## Figure-style manuals (one per paper — pick the one for your figure)

| paper | manual | figures covered |
|---|---|---|
| arXiv:2606.04998 (Vu et al., self-force with NR methods; ledger ELLIPTIC-SSF) | `progress/elliptic-ssf/FIGURE_STYLE_MANUAL.md` | Figs 1–7 (ssf_domain, matrix, precond, ssf_convergence, ssf_modes, linsolv, ssf_error_budget) |
| arXiv:2608.17003 (Dittmann, Hegade K. R., Fairbairn, EMRIs in gaseous disks; ledger EMRI-GASDISK, DAG prefix EGD::) | `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md` | Figs 1–14 (spiral_grid, spirals_two_rows, dT_dr_vs_dFdr, dT_dr_versus_ap, T_GR_versus_T_Newt, chi_rp, rp_h0, fig_pq_m95h0_violin, num_vs_analytical, torque_reversal_in_extreme_disks, convergence, pole_errs, soft_differences, m95_vs_b) |

When a new reference paper is added, first write its manual by VLM
inspection of every reference figure (same structure: §0 shared
conventions — typography, palette, line widths, markers, grid, legends,
figure size, colorbars — then one section per figure with type, layout,
verbatim axis labels, scales/ranges, series encoding, annotations/insets,
colormap, plotted data + parameters, and a VLM MATCH CHECKLIST), add a row
to this table, and commit both.

## Step 1 — Read the manual FIRST

Before writing a single line of plotting code, Read the relevant section of
the paper's manual: §0 (shared rcParams / palette block) AND the per-figure
section. It records plot type, panel layout + shared axes, axis
quantities/scales/ranges, colormaps and series encodings (e.g. 2608:
`#0f4c81/#e2583e/...` categorical cycle, plasma(0.1–0.9) for h0, viridis
for a•, Set2 for the violin, RdBu ±0.77 for δU maps), legend/annotation
placement, typography (usetex CM serif, `\mathfrak` p,q), and a matplotlib
reproduction recipe. Follow the recipe; do not improvise a style from
memory.

## Step 2 — Generate with the recipe, then VLM-CHECK the new figure alone

Write the plotting script using the recipe. Render to a file (PNG or PDF).
Use REAL computed data from the GRTools pipeline — never synthetic
stand-ins presented as results. Then Read the produced image with the Read
tool and check it on its own against the manual's per-figure attributes:
fonts render (no missing-LaTeX fallbacks), ticks in/all sides, right
colours and markers, axis labels verbatim, ranges, legend style/placement,
panel stamps, inset/colorbar present. Fix script-level defects here before
comparing with the reference.

## Step 3 — VLM-COMPARE new vs reference, decide, iterate

1. Read BOTH images visually with the Read tool: the reference figure (its
   PDF/PNG in `ref-paper/<paper>/src/figures/`, or the pre-rendered PNGs)
   and the newly generated figure.
2. Compare STRUCTURE, not vibes, walking through the manual's "VLM MATCH
   CHECKLIST" for that figure: same panels and arrangement, same axes and
   scales/ranges, same trends/slopes/crossings, same characteristic
   features (spikes, plateaus, ordering of series, violin spreads, zero
   crossings, ringing, ISCO truncation…), plus series encoding and
   legend/annotation placement.
3. DECIDE: (a) style/script mismatch → update the plotting script,
   regenerate, go back to 3.1; (b) physics/data mismatch → do NOT touch the
   plot; record it as an error-ledger trial with expected/observed and
   escalate to the owning DAG node; (c) all checklist items match →
   record the verdict. Repeat until the structure matches or the mismatch
   is escalated.

## HARD RULES

- Require the SAME visual structure between reference and reproduction
  before claiming a match. Similar-ish is not a match.
- NEVER rationalize a mismatch as a "paper limit", rendering artifact, or
  "close enough". A mismatch is a finding; record it and fix or escalate.
- NEVER force agreement using the known answer (no fudging data, rescaling
  to hit the reference, or plotting the paper's values as your own).
- Record the VLM verdict in the commit body as a `verify:` object, e.g.:

  ```
  verify: {method: vlm-compare, manual: progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md#8,
           reference: ref-paper/arXiv-2608.17003/src/figures/fig_pq_m95h0_violin.pdf,
           produced: progress/emri-gasdisk/figs/m95_violin_repro.png,
           verdict: match|mismatch, notes: "<checklist items passed/failed>"}
  ```

  A `verdict: match` is only permitted after steps 2 and 3 were actually
  executed in this session (produced image Read alone, then both images
  Read and compared); paste-free claims are inadmissible under
  alignment.md §0.
