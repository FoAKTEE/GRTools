"""
    Hybrid.BinaryWaveform

High-level hybrid binary-waveform assembly (design §4.9): build a [`BinaryInspiralModel`](@ref)
from initial conditions, integrate via the registry-dispatched forcing model, and emit the strain
modes `h_{l,m}(t) = A_{l,m}(Ω(t)) exp(-i m φ(t))` along the trajectory, plus the `(l, m)` mode
list and the orbital-trajectory accessors.
"""
module BinaryWaveform

using ..Models: Amp0PASchwarz, Amp1PAT1, Amp0PAKerr, Amp1PAT1ea,
    Forcing0PASchwarz, Forcing1PAT1, Forcing0PAKerr,
    h0PA_ampΩ, h1PA_ampΩ, real_amp22, mode_keys_0PAKerr
using ..Drivers: InspiralSolution, int_inspiral, omega_at, phi_at, r0_at
import ..Drivers: duration
using ..Registry: get_inspiral, get_amplitudes, list_inspiral_models

include("BinaryWaveform/waveform.jl")
include("BinaryWaveform/binary_inspiral.jl")

export mode_list, waveform_mode, BinaryInspiralModel, binary_inspiral, waveform, trajectory,
    duration, initial_conditions

end # module BinaryWaveform
