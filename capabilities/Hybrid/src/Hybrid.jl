"""
    Hybrid

WaSABI self-force / post-Newtonian HYBRID waveform capability (design §4.9).

Owns the SF-PN-hybrid EMRI waveform pipeline absorbed from the legacy `WaSABI` package:

  - `Models`        — `BinaryInspiralModel` container + waveform/trajectory accessors.
  - `BindingEnergy` — interpolated binding-energy / flux channels.
  - `HybridC`       — the WaSABI-C 9-var SF-PN-hybrid forcing (Fω/Fν/Fm/Fχ) and the
                      δx-coupled trajectory evolution / closed-form shift channels.
  - `Drivers`       — N-var `SpinInspiralSolution` inspiral evolution over
                      OrdinaryDiffEqVerner / SciMLBase.
  - `Registry`      — model registry.
  - `BinaryWaveform`— the hybrid binary waveform assembly.
  - `Introspection` — `cite_model`.

Routes `digamma` to `BigSpecial`, interpolation to `NumericsCore`, and `.wl`/`.jld2`
artifacts to `Data`. Top-of-DAG leaf: the heavy ODE/SciML precompile is confined here.
"""
module Hybrid

import JLD2  # Data JLD2 weakdep-ext trigger (§5.4): DataJLD2Ext activates only when the consumer loads the trigger

# Submodule includes (populated TDD per round). Layering: the interpolation foundation feeds
# BindingEnergy, both feed the Models forcing/amplitude leaves; Drivers/Registry/BinaryWaveform
# sit on top.
include("Interpolation.jl")
include("BindingEnergy.jl")
include("Models.jl")
include("Drivers.jl")
include("HybridC.jl")
include("Registry.jl")
include("BinaryWaveform.jl")
include("Introspection.jl")

end # module Hybrid
