# Eccentric EMRI regression drivers

`circular_regression.jl` executes the existing arXiv:2608.17003 figure drivers
with the adiabatic eccentric response. Run from the GRTools root, with one
process at a time and an external timeout:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 10800 prlimit --cpu=7200:7210 julia -t1 --startup-file=no \
  --project=reproduce reproduce/_emri_eccentric/circular_regression.jl fig10 /absolute/evidence/fig10
```

Repeat for `fig01` through `fig14`. No `GRTOOLS_FIG*` environment overrides,
`--regen`, or `--plot-only` are accepted. All arithmetic is Float64. The driver
records per-mode CPU and wall cost in `mode_costs.jsonl`, figure cost/status in
`figNN_status.json`, and every computed/golden key with its original tolerance
in `figNN.json`. A failed or interrupted figure is not a regression pass.
For the manifest invocation, set `ECC_REGRESSION_EVIDENCE` to an absolute output
directory and use `--check --tier heavy`; this runs all fourteen sequentially.

The numerical configuration, radial grid, mode ceilings, image synthesis,
profile reductions, printed checks and per-key tolerances remain owned by the
original drivers. `Base.include(mapexpr, ...)` adapts only calls to the response
and diagnostic entry points listed in `BACKEND_CALLS`, the final comparison,
and presentation/cache writes. Nested includes apply this transformation only
inside `reproduce/arXiv-2608.17003/`. No package method is replaced. Original
driver files and goldens are read only, and all original expected files are
hashed before and after each run. Figure 5's explicitly Newtonian denominator
continues to use the Newtonian solver; every GR numerator uses the eccentric
response. Configuration and reduction logic are not copied into another driver.
The adapter preserves the original `linsolve=:auto` selection (sparse for
Float64). The general eccentric solver's `:banded` default is unchanged.

Each response constructs an `Adiabatic` `DiscModel`, an `EquatorialOrbit(e=0)`,
a fixed softening equal to the original EGD softening, and `ModeContext(m,k=0)`.
The frozen EGD radial nodes become an `EccentricGrid`. Surveys and wakes use
`solve_modes(context_builder, perturber, modes; bank, grid_builder, transform)`.
One `source_radial_bank` serves all modes, even when their corotation grids and
domains differ. The profile drivers' individual calls share one bank covering
m=1:256 on the fiducial domain. Only one such bank is cached at a time, owned by
this driver and reset between figures. No circular source or response is used.

The bank uses the landed Chebyshev interpolation and direct phase projection
(`method=:direct`, `Nq=2`). The circular carrier is constant in radial phase;
fast phase sums are unnecessary here. Every bank is independently checked at
seven radii plus the orbital radius, for every requested mode and all fourteen
source components/derivatives, against `1e-8 |ref| + 1e-12 U`. Failure aborts the
figure. `bank_costs.jsonl` records this check and the bank's CPU/wall cost. This
source budget does not relax any golden tolerance or bound the BVP error.

`solve_modes` calls the existing eccentric context solver with C-66 radiation
rows. Its `transform` callback reduces each solution immediately; torque
surveys retain only reductions, and wake synthesis retains only radial states.
Torque and flux use the source retained by each eccentric solution. Explicit-bank contexts and grids are built one at a time. Per-mode
JSONL costs include solves and, after the first mode, grid construction. The
first grid is prepared before the bank; profile callers also build their grids
before entering the adapter. Their grid costs, shared-bank work and diagnostic
reductions are included in figure totals; figure CPU cost includes both. Profile-call costs include the first
bank build. `Nq=16` resolves the retained k=-2:2 sideband/limit checks. This is a
regression on the frozen discretization, not certification of finite-e source
quadrature, adaptive-grid selection or exterior truncation.

Known HOLE (lander review of lane w15-g38). The C-66 rows are eigenvectors of the
UNREGULARIZED edge operator; the EGD (66)-(68) rows are eigenvectors of the
Landau-shifted edge matrix (`SpiralWaves/solver.jl`, `solve_mode(bg, m)`). Where a
domain edge lies inside the Landau profile (small a near the ISCO, narrow high-m
domains) they define different e=0 problems, at order eps(r_edge)/|omega_d(r_edge)|:
up to 1e-5 per mode and 1.3e-7 in the Fig. 6 total at a = 11 M. The regression
therefore does not yet meet its contract; the fix is tracked as obligation
gr-o-e0-boundary-convention. The test file pins the convention near corotation.

The original drivers are loaded only from this checkout's reproduction tree.
No mirrored source artifact is executed. Their arguments are temporarily set
to `--check --tier heavy` and then restored; tolerances, reductions and numerical
configurations remain unchanged. Evidence must use a fresh, empty directory outside the reproduction tree;
its parent must already exist. Symlinks into the reproduction tree and nonempty
attempt directories are rejected before any write. Run each figure in a separate process with
the external CPU limit above. A signal/timeout leaves measured partial costs
in the command receipt and completed-mode log and is a cost HOLE.

Shared-bank tests verify bitwise identity on identical samples, interpolation
against direct source jets and responses, input order, and bank isolation.
The package test is run by direct include (no change to `test/runtests.jl`):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 600 julia -t1 --startup-file=no \
  --project=capabilities/AccretionDisc -e 'include("capabilities/AccretionDisc/test/test_eccentric_regression.jl")'
```

Figure 10 keeps its landed qualification: the published isolated-marker torque
deficit is 5–18%, corresponding to a 2–9% driven source-amplitude difference.
The landed EGD node is preliminary. Its earlier claim of convergence to the
published values by m=40 was withdrawn because those markers were clipped by
the frame. This regression compares with the unchanged EGD numerical golden;
it does not resolve that paper-to-calculation residual or tune the source.
