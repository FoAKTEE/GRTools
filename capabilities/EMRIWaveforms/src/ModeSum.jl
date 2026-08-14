"""
    EMRIWaveforms.ModeSum

Time-domain observed-strain reduction `h = Σ A·sYlm·e^{-iΦ}` (design §6.8: FEW
`summation/{directmodesum,interpolatedmodesum,fdinterpolatedmodesum}.jl`; KerrBackwards
`{summation,fdsummation}.jl`). Thin convention adapters delegating to the single consolidated,
already-characterized `Kernels.direct_mode_sum` / `Kernels.interpolated_mode_sum`. v0.1 keeps the two
forks behind the qualified deep namespaces `ModeSum.KerrBackwards` / `ModeSum.FEW` (§4.8).

- KerrBackwards fork: `(teuk_modes, ylms, Φ_φ, Φ_r, m_arr, n_arr)` is exactly the Kernels signature ⇒
  direct delegation (bit-identical by construction).
- FEW fork: `(…, Φ_θ, l_arr, m_arr, k_arr, n_arr)` carries the polar phase `Φ_θ`/index `k` and the
  partner sign `(−1)^l`. The KerrEccEq waveform pipeline is equatorial (`k = 0`), so the adapter folds
  `(−1)^l` into the `m < 0` partner harmonics and delegates to the consolidated Kernels kernel; a nonzero
  `k_arr` is rejected (would need the not-yet-consolidated generic-Kerr polar term).
- FD-SPA (`fd_interpolated_mode_sum` / `fd_flip_average`): the frequency-domain stationary-phase
  mode-sum, hoisted to `Kernels.fd_interpolated_mode_sum` (S4-a, §5.5 'Time-domain + FD-SPA + AAK
  mode-sum live in Kernels' / §6.8 single copy — bit-identical to the pre-hoist local kernel per the
  baked Kernels fixture). The FEW names here are thin facades; the SPA primitives
  (`_spa_func`/`_cube_roots`/`_fd_seg_freq_range`) stay addressable as const aliases (§7.7 facade guard).
"""
module ModeSum

import Kernels

# =====================================================================================
# KerrBackwards fork — Kernels signature verbatim ⇒ direct delegation.
# =====================================================================================
module KerrBackwards

import Kernels
export direct_mode_sum, interpolated_mode_sum

"""
    direct_mode_sum(teuk_modes, ylms, Phi_phi, Phi_r, m_arr, n_arr; threaded=false)

Direct ±m time-domain mode-sum (KB `summation.mode_sum`). Delegates to `Kernels.direct_mode_sum`.
"""
direct_mode_sum(teuk_modes, ylms, Phi_phi, Phi_r, m_arr, n_arr; threaded::Bool = false) =
    Kernels.direct_mode_sum(teuk_modes, ylms, Phi_phi, Phi_r, m_arr, n_arr; threaded = threaded)

"""
    interpolated_mode_sum(t, teuk_modes, ylms, Phi_phi, Phi_r, m_arr, n_arr; dt=10.0, out_len=nothing, threaded=false)

Interpolated time-domain mode-sum (KB `summation.interpolated_mode_sum`). Delegates to
`Kernels.interpolated_mode_sum`.
"""
interpolated_mode_sum(t, teuk_modes, ylms, Phi_phi, Phi_r, m_arr, n_arr;
                      dt::Real = 10.0, out_len = nothing, threaded::Bool = false) =
    Kernels.interpolated_mode_sum(t, teuk_modes, ylms, Phi_phi, Phi_r, m_arr, n_arr;
                                  dt = dt, out_len = out_len, threaded = threaded)

end # module KerrBackwards

# =====================================================================================
# FEW fork — fold (−1)^l partner sign, require equatorial k=0, delegate to Kernels.
# =====================================================================================
module FEW

import Kernels
export direct_mode_sum, interpolated_mode_sum, fd_interpolated_mode_sum, fd_flip_average

# Fold the FEW (−1)^l partner-harmonic sign into the m<0 ylm block so the consolidated Kernels kernel
# (whose partner term is ylms[K+j]·conj(A)·e^{+iψ}) reproduces FEW's (−1)^l·ylms[K+j]·conj(A)·e^{+iψ}.
function _fold_partner_sign(ylms::AbstractVector, l_arr, K)
    out = collect(ylms)
    @inbounds for j in 1:K
        out[K + j] = (-1)^l_arr[j] * out[K + j]
    end
    return out
end

_require_equatorial(k_arr) = all(iszero, k_arr) ||
    throw(ArgumentError("FEW.direct_mode_sum: nonzero k_arr (generic-Kerr polar term) is not supported \
                         by the consolidated equatorial Kernels mode-sum"))

"""
    direct_mode_sum(teuk_modes, ylms, Phi_phi, Phi_theta, Phi_r, l_arr, m_arr, k_arr, n_arr; threaded=false)

Direct mode-sum (FEW `DirectModeSum.sum`). Requires equatorial `k_arr = 0`; folds the `(−1)^l` partner
sign into the `m < 0` harmonics and delegates to `Kernels.direct_mode_sum`.
"""
function direct_mode_sum(teuk_modes::AbstractMatrix, ylms::AbstractVector,
                         Phi_phi, Phi_theta, Phi_r, l_arr, m_arr, k_arr, n_arr; threaded::Bool = false)
    _require_equatorial(k_arr)
    K = size(teuk_modes, 2)
    ylms2 = _fold_partner_sign(ylms, l_arr, K)
    return Kernels.direct_mode_sum(teuk_modes, ylms2, Phi_phi, Phi_r, m_arr, n_arr; threaded = threaded)
end

"""
    interpolated_mode_sum(t, teuk_modes, ylms, Phi_phi, Phi_theta, Phi_r, l_arr, m_arr, k_arr, n_arr;
                          dt=10.0, out_len=nothing, threaded=false)

Interpolated mode-sum (FEW `InterpolatedModeSum.sum`). Equatorial `k_arr = 0`; folds `(−1)^l` and
delegates to `Kernels.interpolated_mode_sum`.
"""
function interpolated_mode_sum(t, teuk_modes::AbstractMatrix, ylms::AbstractVector,
                               Phi_phi, Phi_theta, Phi_r, l_arr, m_arr, k_arr, n_arr;
                               dt::Real = 10.0, out_len = nothing, threaded::Bool = false)
    _require_equatorial(k_arr)
    K = size(teuk_modes, 2)
    ylms2 = _fold_partner_sign(ylms, l_arr, K)
    return Kernels.interpolated_mode_sum(t, teuk_modes, ylms2, Phi_phi, Phi_r, m_arr, n_arr;
                                         dt = dt, out_len = out_len, threaded = threaded)
end

# --- FD-SPA frequency-domain stationary-phase mode-sum ------------------------------
# Hoisted to Kernels (S4-a; single canonical copy per design §5.5/§6.8). The delegation is
# bit-identical: Kernels' gather-form kernel replicates the legacy scatter's per-bin accumulation
# order exactly (characterized in Kernels against the baked pre-hoist fixture fdspa_ref.jl captured
# from THIS module's local kernel @ f409fc7). The internal SPA primitives remain addressable under
# their FEW names as const facades (§7.7 facade/sweep guard — the test-referenced symbols survive).

const _spa_func = Kernels.spa_func
const _cube_roots = Kernels._cube_roots
const _fd_seg_freq_range = Kernels._fd_seg_freq_range

"""
    fd_flip_average(waveform, frequency) -> (hplus, hcross)

Convert the raw `±f` SPA spectrum into `h₊(f)`/`h₍×₎(f)` via the flip-average symmetry
(`fdinterp.py`). Facade for [`Kernels.fd_flip_average`].
"""
const fd_flip_average = Kernels.fd_flip_average

"""
    fd_interpolated_mode_sum(t, teuk_modes, ylms, Phi_phi, Phi_r, f_phi, f_r, m_arr, n_arr;
                             f_arr, include_minus_m=true, separate_modes=false) -> (hplus, hcross)

Frequency-domain stationary-phase mode-sum (FEW `FDInterpolatedModeSum.sum`). Thin adapter keeping the
FEW-fork signature; delegates to the single consolidated [`Kernels.fd_interpolated_mode_sum`]
(bit-identical to the pre-hoist local kernel). `separate_modes=true` returns the raw per-mode `±f`
spectrum (no flip-average; `include_minus_m` forced off); else `(h₊,h₍×₎)`.
"""
fd_interpolated_mode_sum(t::AbstractVector, teuk_modes::AbstractMatrix,
                         ylms::AbstractVector, Phi_phi::AbstractVector,
                         Phi_r::AbstractVector, f_phi::AbstractVector,
                         f_r::AbstractVector, m_arr::AbstractVector{<:Integer},
                         n_arr::AbstractVector{<:Integer};
                         f_arr::AbstractVector, include_minus_m::Bool = true,
                         separate_modes::Bool = false) =
    Kernels.fd_interpolated_mode_sum(t, teuk_modes, ylms, Phi_phi, Phi_r, f_phi, f_r,
                                     m_arr, n_arr; f_arr = f_arr,
                                     include_minus_m = include_minus_m,
                                     separate_modes = separate_modes)

end # module FEW

end # module ModeSum
