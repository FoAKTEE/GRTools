# Provenance of the vendored agent skills

| skill | origin | captured |
|---|---|---|
| `compute-usage` | `emri-disk/.claude/skills/compute-usage/SKILL.md` (verbatim, md5 `c62a2e990916415d63dc7363187902ac`, source mtime 2026-07-28T07:22:32Z) | 2026-08-19 |
| `vlm-figure-repro` | `emri-disk/.claude/skills/vlm-figure-repro/SKILL.md` (verbatim, md5 `08b7945396a8ab0025400a28e9b5c4f5`, source mtime 2026-08-19T10:05:41Z) | 2026-08-19 |
| `repro-verify` | written here for GRTools (no upstream copy) | 2026-08-19 |

Notes:

- The two copied skills are **verbatim** snapshots of the mission repo
  (`/data/haiyangw/claude/emri-disk`). They may be revised there concurrently;
  re-copy and update the md5s above rather than editing the copies by hand.
- `vlm-figure-repro` points at the per-paper `FIGURE_STYLE_MANUAL.md` files
  (`progress/elliptic-ssf/`, `progress/emri-gasdisk/`) — those paths are
  **relative to the mission repo**, not to GRTools. Read the manual there; the workflow (read manual → generate with
  the recipe → VLM-compare structure → iterate) is repo-independent, and the
  per-paper manuals for GRTools reproductions live under `reproduce/_style/manuals/`.
- `compute-usage` limits (≤ 24 CPUs for our jobs, total system CPU < 90 %, no
  GPU) are host-owner policy and apply to every job launched from this repo.

Re-sync log: the `vlm-figure-repro` source was revised in the mission repo on
2026-08-19 (per-paper manual table, explicit VLM-check-then-compare steps) while
this vendoring was in flight; the copy above is the post-revision snapshot.
