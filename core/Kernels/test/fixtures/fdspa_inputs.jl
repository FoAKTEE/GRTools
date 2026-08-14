# fdspa_inputs.jl — deterministic input generators for the FD-SPA mode-sum characterization.
#
# Shared VERBATIM between the one-shot pre-hoist capture script (scratchpad/s4a_capture_fdspa_fixture.jl,
# run against the committed EMRIWaveforms.ModeSum.FEW.fd_interpolated_mode_sum BEFORE the kernel moved to
# Kernels) and test_modesum.jl — so the baked fdspa_ref.jl literals are grounded on the PRE-CHANGE
# behavior with bit-identical inputs (no hand transcription).
#
# Case 1: the EMRIWaveforms test_modesum.jl configuration (nt=5, K=2, m=[2,2] — all mode frequencies
#         positive, sources strictly above the zero bin; has SPA-turnover NaN bins).
# Case 2: a zero-crossing configuration (mode 1: m=1,n=-1 with f_phi−f_r crossing 0 inside segment 3),
#         so single (mode,segment) bin-ranges straddle the zero bin and the ±m partner scatter lands
#         direct+partner contributions on the SAME target bins — exercising the direct-before-partner,
#         partner-before-direct AND the i==zidx interleaved accumulation orders of the legacy kernel.

function fdspa_case1_inputs()
    nt = 5; K = 2
    t = collect(0.0:100.0:400.0)
    teuk = ComplexF64[complex(0.5 * cos(0.01 * i * j), 0.3 * sin(0.013 * i * j)) for i in 1:nt, j in 1:K]
    ylms = ComplexF64[1.0 + 0.0im, 0.6 - 0.2im, 0.8 + 0.1im, 0.4 + 0.3im]
    Phi_phi = [0.02 * (i - 1) * 100 for i in 1:nt]
    Phi_r   = [0.013 * (i - 1) * 100 for i in 1:nt]
    f_phi = fill(2.0e-3, nt) .+ [1e-5 * (i - 1) for i in 1:nt]
    f_r   = fill(1.3e-3, nt) .+ [0.7e-5 * (i - 1) for i in 1:nt]
    m_arr = [2, 2]; n_arr = [0, 1]
    df = 2.0e-5; nf = 601
    f_arr = [df * (i - (nf + 1) ÷ 2) for i in 1:nf]
    return (t = t, teuk = teuk, ylms = ylms, Phi_phi = Phi_phi, Phi_r = Phi_r,
            f_phi = f_phi, f_r = f_r, m_arr = m_arr, n_arr = n_arr, f_arr = f_arr)
end

function fdspa_case2_inputs()
    nt = 6; K = 3
    t = collect(0.0:100.0:500.0)
    teuk = ComplexF64[complex(0.4 * cos(0.21 * i + 0.3 * j) + 0.05 * j,
                              0.3 * sin(0.17 * i - 0.2 * j)) for i in 1:nt, j in 1:K]
    ylms = ComplexF64[1.0 + 0.0im, 0.6 - 0.2im, -0.35 + 0.45im,
                      0.8 + 0.1im, 0.4 + 0.3im, 0.25 - 0.55im]
    Phi_phi = [0.66 * (i - 1) + 0.037 * (i - 1)^2 for i in 1:nt]
    Phi_r   = [0.44 * (i - 1) + 0.011 * (i - 1)^2 for i in 1:nt]
    # mode 1 frequency f_phi − f_r crosses 0 near t = 250 (inside segment 3); mild sinusoidal
    # curvature keeps the per-segment frequency cubic away from the degenerate-linear limit.
    f_phi = [1.05e-3 + 1.2e-6 * (100 * (i - 1) - 250) + 2.0e-5 * sin(0.5 * (i - 1)) for i in 1:nt]
    f_r   = fill(1.05e-3, nt)
    m_arr = [1, 2, 2]; n_arr = [-1, 0, 1]
    df = 2.5e-5; nf = 321
    f_arr = [df * (i - (nf + 1) ÷ 2) for i in 1:nf]
    return (t = t, teuk = teuk, ylms = ylms, Phi_phi = Phi_phi, Phi_r = Phi_r,
            f_phi = f_phi, f_r = f_r, m_arr = m_arr, n_arr = n_arr, f_arr = f_arr)
end
