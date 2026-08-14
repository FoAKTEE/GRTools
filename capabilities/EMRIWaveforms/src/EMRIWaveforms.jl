"""
    EMRIWaveforms

GrTool tier-5 capability: the EMRI trajectory + observed-waveform pipeline (design §4.8, §6.8),
fusing FastEMRIWaveforms and KerrBackwards. HDF5 HOME #2.

# Submodules (pipeline order)
- `Trajectory`  — inspiral ODE evolution (vendored DOPR853 integrator + Spacetime
  geodesics/separatrix + Roots/NLsolve parameter-maps). Single `run_inspiral(...; backend=:few|:kb)`.
- `FluxModels`  — flux ODE RHS over `Data` flux grids (Schwarz-ecc / Kerr-ecc-eq).
- `AmplitudeModels` — RomanAmplitude NN (via `Kernels` amp GEMM) + `Data` amplitude-grid bicubic.
  Single `RomanAmplitude` + `get_amplitudes(...; backend=:few|:kb)`.
- `ModeSelect`  — power-thresholded harmonic mode selection.
- `ModeSum`     — waveform assembly `h = Σ A·sYlm·e^{-iΦ}` (TD / FD-SPA / AAK via `Kernels`).
- `AAK`         — augmented-analytic-kludge basis + single `aak_summation(...; backend=:few|:kb)` /
  `aak_waveform`.
- `Mappings`    — orbital-parameter coordinate maps (reuse `NumericsCore` conventions).
- `Generate`    — top-level observed-waveform generation.

# Single-API surface (§4.8, §6.8 [MILESTONE] CLOSED)
The three duplicated FastEMRIWaveforms / KerrBackwards forks — the `Trajectory` inspiral driver,
`AmplitudeModels.RomanAmplitude`, and `AAK.aak_summation`/`aak_waveform` — are each ONE canonical API:
a `backend::Symbol` selector (`:few` default, `:kb`) dispatches to the shared kernel and reproduces both
legacy forks bit-for-bit (the redundant fork scaffolding is deleted). `FluxModels`, `ModeSelect`, and
`ModeSum` retain `KerrBackwards`/`FEW` submodules for genuinely DISTINCT backends (Schwarz-ecc vs
Kerr-ecc-eq flux grids; KB (l,m,n) vs FEW selection; time-domain vs FD-SPA mode-sum) — these are separate
physics, not duplicated scaffolding.

The 4 dead FEW import edges (Teukolsky / TeukolskyPN / TeukolskySpinFluxes / SpinningOrbit) are
DELETED → this capability is Symbolics-clean. It is Float64/GPU-only: it CALLS the BigFloat-generic
lower tiers at Float64 and imposes no precision constraint upward.
"""
module EMRIWaveforms

# HDF5 HOME #2 ext trigger (design §5.4 / S4-b §7.2/§10): ALL grid reader BODIES live in Data's HDF5
# weakdep-extension; this direct import (HDF5 is in our [deps]) is what activates DataHDF5Ext, so
# Data.read_kerrampgrid / Data.read_roman_nn / Data.read_ampgrid / Data.read_fluxgrid have bodies.
# No EMRIWaveforms source performs direct HDF5 reads anymore.
import HDF5

# --- pipeline submodules (stubs during scaffold; filled tests-first) ---
include("Mappings.jl")
include("Trajectory.jl")
include("FluxModels.jl")
include("AmplitudeModels.jl")
include("ModeSelect.jl")
include("ModeSum.jl")
include("AAK.jl")
include("Generate.jl")

using .Mappings
using .Trajectory
using .FluxModels
using .AmplitudeModels
using .ModeSelect
using .ModeSum
using .AAK
using .Generate

end # module EMRIWaveforms
