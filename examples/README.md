# GrTool examples

Self-contained, runnable scripts demonstrating the suite. Each prints labeled
results with built-in cross-checks against exact limits or identities.

From the repository root (first run resolves and precompiles the suite):

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=. examples/01_kerr_geodesics.jl
```

| script | what it shows |
|---|---|
| `01_kerr_geodesics.jl` | Kerr bound-orbit constants of motion (E, Lz, Q), fundamental frequencies, ISCO/separatrix, BigFloat inheritance — `Spacetime` |
| `02_spheroidal_harmonics.jl` | Spin-weighted spheroidal eigenvalues (two backends), S(θ,φ) evaluation, spherical limit, complex spheroidicity — `Harmonics` |
| `03_wigner_coupling.jl` | 3j/6j/Clebsch–Gordan with orthogonality checks, 3-harmonic coupling, exact BigInt core — `Wigner` |
| `04_teukolsky_fluxes.jl` | Circular and eccentric point-particle Teukolsky modes: GW energy/angular-momentum fluxes, the J̇ = (m/ω)Ė identity, the Newtonian 32/5 limit — `Perturbation` |
| `05_quasinormal_modes.jl` | Kerr QNM frequencies across spin (Leaver continued fraction vs :fast backend), overtones — `Perturbation` |
| `06_elliptic_core.jl` | Spectral reference elements (exact quadrature/differentiation) and complex-generic matrix-free GMRES — `EllipticCore` |
| `07_bigfloat_precision.jl` | End-to-end BigFloat precision inheritance: 256- vs 512-bit self-convergence across three packages |

Conventions throughout: geometric units `G = c = 1`, black-hole mass `M = 1`;
spin `a` (or `χ`) dimensionless; fields `∝ e^{i(mφ − ωt)}`. See
`docs/src/topics/conventions.md` for the full statement.
