---
name: compute-usage
description: 'Compute-usage self-check before and after launching any large computation job (Julia test suites, solver runs, parallel m-mode sweeps, BigFloat builds). Use whenever launching parallel workers, long-running solves, or multi-process jobs.'
---

# Compute-usage self-check (owner policy, binding)

Owner limits for this host (256 cores, shared multi-user):

- **CPU**: use no more than 50% of all CPUs in principle; in practice use *middle,
  non-occupied* cores and **at most 24 CPUs** for our jobs.
- **Whole-system guard**: before AND after launching jobs, total CPU usage of ALL
  users must stay **below 90%**.
- **GPU**: banned by admin. Never use CUDA/AMDGPU paths; test only on CPUs.

## Before launching

```bash
nproc; uptime
ps -eo pcpu --no-headers | awk '{s+=$1} END {printf "total CPU%%: %.0f of %d00\n", s, '"$(nproc)"'}'
```

- If `total/(nproc*100) > 0.90 - (planned_job_share)`, do NOT launch; shrink the job
  or wait.
- Cap explicitly: `julia -t N` / `-p N` with `N ≤ 24`; prefer `N ≤ 4` threads for
  BigFloat-heavy work (GC contention makes `-t16` net slower than `-t4` — measured);
  scale BigFloat via processes, not threads.
- For embarrassingly parallel sweeps (e.g. independent m-modes), size the worker pool
  as `min(24, modes)` and leave it there; never auto-scale up.

## After launching

Re-run the same check. If total system usage exceeds 90%, reduce the pool (kill by
exact PID — never `pkill -f` a pattern contained in your own command line).

## Record

Any commit whose `run:` object involved a parallel/large job states the core count
used and that the 90% guard was checked (before/after).
