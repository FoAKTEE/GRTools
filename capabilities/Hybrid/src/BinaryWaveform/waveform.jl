# Strain amplitudes and spin-weighted spherical-harmonic mode list.
# Port of ref-code/WaSABI/Kernel/Waveform.m (WaSABI`Waveform` context).

"""
    mode_list(model)

Return the sorted list of available `(l, m)` modes (positive `m`). Port of `ModeList`, which
returns `Keys[GetAmplitudes[model]]`: here the `(l, m)` keys of the amplitude model.
"""
function mode_list end

mode_list(amp::Amp0PASchwarz) = sort(collect(keys(amp.Z)))
mode_list(amp::Amp1PAT1) = sort(collect(keys(amp.base.Z)))
mode_list(::Amp0PAKerr) = mode_keys_0PAKerr()
mode_list(amp::Amp1PAT1ea) = sort(collect(keys(amp.base.base.Z)))

"""
    waveform_mode(model, l, m)

Return a spin-weighted spherical-harmonic mode amplitude of the gravitational-wave strain.
Port of `WaveformMode` (decl-only generic upstream).
"""
function waveform_mode end
