# Hybrid

Tier-5 GrTool capability (design §4.9): the **WaSABI self-force / post-Newtonian hybrid**
EMRI waveform engine.

## Physics owned

The SF-PN-hybrid inspiral-and-waveform pipeline. A self-force-accurate forcing model is fused
with a post-Newtonian completion to drive a spinning EMRI inspiral, and the resulting trajectory
is mapped to a binary gravitational waveform. This is the WaSABI-C scheme: a 9-variable
SF-PN-hybrid forcing (Fω/Fν/Fm/Fχ) with δx-coupled trajectory evolution and closed-form
shift channels.

## Legacy package absorbed

- `WaSABI` — entire SF-PN-hybrid waveform package (forcing, evolution, binding-energy
  interpolants, model registry, hybrid binary waveform, model citation).

The legacy `WaSABI` package stays green as the characterization oracle and is **never** added as
a live test dependency (it pulls a heavy SF-PN + ODE closure). Tests instead ground against
wolframscript runs of WaSABI's WL forcing/interpolation/NDSolve sources, captured fixtures /
baked literals from the legacy `.m`/`.wl` data, and analytic identities.

## Submodules / public API

- `Models` — `BinaryInspiralModel` container, waveform/trajectory accessors.
- `BindingEnergy` — interpolated binding-energy / flux channels.
- `HybridC` — WaSABI-C 9-var SF-PN-hybrid forcing and δx-coupled trajectory codegen.
- `Drivers` — N-var `SpinInspiralSolution` inspiral evolution (OrdinaryDiffEqVerner / SciMLBase).
- `Registry` — model registry.
- `BinaryWaveform` — hybrid binary waveform assembly.
- `Introspection` — `cite_model`.

The upstream-undefined `inspiral_evaluate` path is carried as an honest stub-`[HOLE]`.

## Dependencies

GrTool: `NumericsCore` (interpolation), `BigSpecial` (`digamma`), `Spacetime`, `Data`
(`.wl`/`.jld2` artifacts). External: `OrdinaryDiffEqVerner`, `SciMLBase`, `JLD2`.

This capability is a top-of-DAG leaf, so the heavy ODE/SciML precompile floor is confined here
and never reaches the core numeric/GPU waveform paths.
