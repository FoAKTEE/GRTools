---
name: vlm-figure-repro
description: 'VLM-verified figure reproduction workflow: read the mission FIGURE_STYLE_MANUAL, generate the figure, then VLM-compare against the reference figure and iterate. Use whenever producing or checking any reproduction figure.'
---

# VLM-verified figure reproduction

Recurring workflow for every reproduction figure in this mission
(arXiv:2606.04998 and any other ref-paper). Never generate or sign off a
figure without completing all three steps.

## Step 1 — Read the manual FIRST

Before writing a single line of plotting code, Read the relevant section of
the mission's figure-style manual:

- `progress/elliptic-ssf/FIGURE_STYLE_MANUAL.md` (arXiv:2606.04998 figures)

It records, per figure: plot type, panel layout + shared axes, axis
quantities/scales/ranges, colormaps and series encodings, legend/annotation
placement, typography, and a matplotlib reproduction recipe. Follow the
recipe; do not improvise a style from memory.

## Step 2 — Generate with the recipe

Write the plotting script using the recipe (layout, colormap, scales,
grayscale/line-style hierarchy, red-accent semantics, CM-serif fonts).
Render to a file (PNG or PDF). Use REAL computed data from the pipeline —
never synthetic stand-ins presented as results.

## Step 3 — VLM-compare and iterate

1. Read BOTH images visually with the Read tool: the reference figure
   (its PDF in `ref-paper/<paper>/src/`) and the newly generated figure.
2. Compare STRUCTURE, not vibes: same panels and arrangement, same axes and
   scales/ranges, same trends/slopes, same characteristic features
   (restart segments, tolerance floors, stacked planes, budget ordering,
   worldtube dips, etc.), plus series encoding and legend/annotation
   placement. Use the per-figure "Match check" list in the manual.
3. If ANY structural mismatch: state it explicitly, update the plotting
   script (or fix the underlying data/pipeline if the mismatch is physical),
   regenerate, and go back to 3.1. Repeat until the structure matches or
   the mismatch is escalated as a physics discrepancy in the error ledger.

## HARD RULES

- Require the SAME visual structure between reference and reproduction
  before claiming a match. Similar-ish is not a match.
- NEVER rationalize a mismatch as a "paper limit", rendering artifact, or
  "close enough". A mismatch is a finding; record it and fix or escalate.
- NEVER force agreement using the known answer (no fudging data, rescaling
  to hit the reference, or plotting the paper's values as your own).
- Record the VLM verdict in the commit body as a `verify:` object, e.g.:

  ```
  verify: {method: vlm-compare, reference: ref-paper/arxiv-2606.04998/src/linsolv.pdf,
           produced: progress/elliptic-ssf/figs/linsolv_repro.png,
           verdict: match|mismatch, notes: "<structural comparison summary>"}
  ```

  A `verdict: match` is only permitted after step 3 was actually executed
  in this session (both images Read and compared); paste-free claims are
  inadmissible under alignment.md §0.
