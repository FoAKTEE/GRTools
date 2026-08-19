---
name: repro-verify
description: 'Before any completion claim about a numerical result in GRTools, run the reproduction check: julia --project=reproduce reproduce/run_all.jl --check [--paper <id>] and paste the output; regenerate goldens only with --regen and a commit that explains why.'
---

# repro-verify — the reproduction gate for numerical claims

Binding for every GRTools node whose deliverable is a *number*, *table*, or
*figure* reproduced from a paper. The rule is alignment.md §0: the agent
proposes, the verifier admits. Here the verifier is the reproduction harness
under `reproduce/`, and its **pasted output** is the only admissible evidence.

## Step 1 — Run the check (never skip, never summarize from memory)

```bash
cd /data/haiyangw/claude/GRTools
julia --project=reproduce reproduce/run_all.jl --check                 # whole suite
julia --project=reproduce reproduce/run_all.jl --check --paper 2608.17003   # one paper
```

- `--check` recomputes every registered target and compares it against the
  committed **golden** values at the tolerance each target declares. Exit 0 =
  reproduced; nonzero = a real discrepancy.
- Run it **before** writing the completion claim, not after. If you did not run
  it in this session, you have no result to claim.
- The check is a computation: obey the `compute-usage` skill (below) when it is
  large.

## Step 2 — Paste the output

Paste the verbatim tail of the run (target ids, observed vs golden, tolerance,
verdict, exit code) into

- the commit body, as a `verify:` object, and
- the node's evidence directory
  (`results/<paper-slug>/paper_<PAPER>/evidence/<node>/` in the mission repo,
  emri-disk), as a log file that the knowledge-ledger row points at.

A completion claim with no pasted verifier output is inadmissible — do not
write one.

## Step 3 — Goldens change only deliberately

```bash
julia --project=reproduce reproduce/run_all.jl --regen [--paper <id>]
```

- `--regen` **overwrites** the goldens. It is never part of "making the check
  pass". Use it only when the *reference* changed: a new paper target, a fixed
  bug in the physics, a deliberate tolerance or convention change.
- A `--regen` commit must state, in the body, **why** the golden moved: the
  old value, the new value, and the root cause of the difference, tagged
  `[SOLID|PRELIMINARY|HOLE|FUTURE]`. A regen with no explanation is a silently
  rewritten result — the exact failure this gate exists to prevent.
- Never regen to hide a failing check, and never hand-edit a golden file.

## Compute discipline

Before launching a full `--check` (or any parallel sweep it triggers), follow
the **`compute-usage`** skill: `nproc; uptime`, keep our jobs at
**≤ 24 CPUs**, keep total system CPU **< 90 %**, never GPU. Prefer
`julia --project=reproduce -t4` for BigFloat-heavy targets (measured: `-t16` is
net slower than `-t4` under BigFloat GC contention); scale with processes.

## Commit requirement

Every commit that lands or changes a numerical result carries a `verify:` body
object naming the check and its verdict, e.g.

```
- verify: julia --project=reproduce reproduce/run_all.jl --check --paper 2608.17003
          -> 14/14 targets within tol, exit 0 (log: .../evidence/<node>/repro_check.log)
```

and a `result:` object with a claim tag. `verdict: pass` is only permitted if
the command was actually executed in this session and its output pasted.

## Status of the harness

`reproduce/run_all.jl` is landed by the reproduction-harness node of the
GRTOOLS-REPRO DAG. If it is not present in your checkout, the node you depend
on has not landed yet — say so and block; do not invent an alternative check
and do not claim a reproduction.
