# arXiv:1301.1591 — CONFIGURATION shared by the a = 0.8 mode-scattering drivers (T1.1 analogue).
#
# ZERO PHYSICS: run configuration (legacy `results/dg-teukolsky-spin/scripts/run_scatter.jl`:
# s = −2, a = 0.8, K = 24, N = 8, filter_alpha = 36, tf = 3000, ℓ = 2..6 coupled, ID0 for the m = 0
# rows and ID1 for the m = 2 rows of the paper's non-axisymmetric table) plus the PRINTED decay-rate
# tables the measurements are compared against. All estimators are package functions.

using GrTools
using GrTools.Perturbation: TimeDomain as TDS

const S_SPIN, A_SCAT, K_SCAT, TF_SCAT = -2, 0.8, 24, 3000.0
const ELLS_SCAT = collect(2:6)
const P4_WINDOW_SCAT = (700.0, 1300.0)

"Printed s = −2 decay rates (finite | scri) [arXiv:1301.1591, Tab. s2_LPI, ID0 m = 0 block]; NaN = paper marks ×."
const T_ID0_M0 = Dict(
    2 => Dict(2 => (7, 6), 3 => (8, 7), 4 => (9, 8), 5 => (10, 9)),
    3 => Dict(2 => (7, 6), 3 => (8, 7), 4 => (9, 8), 5 => (10, 9)),
    4 => Dict(2 => (8, 7), 3 => (9, 7), 4 => (10, 8), 5 => (11, 9)),
    5 => Dict(2 => (9, 8), 3 => (10, 8), 4 => (11, 8), 5 => (NaN, 9)))
"Printed s = −2, m = 2 rates [arXiv:1301.1591, Tab. m_neq_0_rates_sm2], ID1 seeds."
const T_ID1_M2 = Dict(
    2 => Dict(2 => (7, 6), 3 => (8, 7), 4 => (9, 8), 5 => (10, 9)),
    3 => Dict(2 => (7, 6), 3 => (8, 7), 4 => (9, 8), 5 => (10, 9)))

"""
    run_scatter(lp, m, id) -> (τs, Ψ, maxu, dt)

Seed multipole `lp` (via `init_mode`) inside the coupled window ℓ = 2..6 at a = 0.8 and record
every projected mode at the three probes.
"""
function run_scatter(lp::Int, m::Int, id::Symbol)
    pol = TDS.spin_excise_policy(1.0, A_SCAT)
    p = TDS.SpinModesProblem(S_SPIN, ELLS_SCAT, m; N = 8, K = K_SCAT, M = 1.0, a = A_SCAT,
                             excise = pol.excise, filter_alpha = 36.0)
    @assert abs(p.mesh.edges[1] - (pol.Rplus - pol.excise)) < 1e-12
    dt = TDS.spin_cfl_dt(p; C = 0.6)
    rec = max(1, round(Int, 1.0 / dt))
    τs, Ψ, maxu = TDS.solve_spin_modes(p; tfinal = TF_SCAT, dt = dt, record_every = rec, id = id,
                                       init_mode = findfirst(==(lp), ELLS_SCAT),
                                       ic = (amp = 1.0, R0 = 0.75, w = 1500.0))
    return (τs, Ψ, maxu, dt)
end

"""
    register_scatter!(lp, m, id, τs, Ψ, maxu, dt, table)

Register one row block of the scattering matrix: for each measured ℓ at the mid and scri probes,
the raw P4-window fit, the envelope fit, the sliding plateau/min/max and the ADAPTIVE-window fit
(the per-entry statistic of the legacy table: the least floor-biased uniform choice), plus the
printed law and the deviation. `NaN` laws (paper marks the entry ×) are carried as NaN, not
silently dropped.
"""
function register_scatter!(lp::Int, m::Int, id::Symbol, τs, Ψ, maxu, dt, table)
    pre = "lp$(lp)_m$(m)_$(String(id))"
    @repro_result "$(pre)_dt" dt
    @repro_result "$(pre)_maxu_end" maxu[end]
    @repro_result "$(pre)_maxu_finite" all(isfinite, maxu)
    for l in 2:5, (ip, pn, li) in ((2, "mid", 1), (3, "scri", 2))
        q = findfirst(==(l), ELLS_SCAT)
        ψ = abs.(Ψ[ip, q, :])
        law = table[lp][l][li]
        @repro_result "$(pre)_$(pn)_l$(l)_law" Float64(isnan(law) ? NaN : -law)
        @repro_result "$(pre)_$(pn)_l$(l)_raw" TDS.fit_tail_window(τs, ψ, P4_WINDOW_SCAT)
        nenv, npk = TDS.fit_tail_envelope(τs, ψ, P4_WINDOW_SCAT)
        @repro_result "$(pre)_$(pn)_l$(l)_env" nenv
        @repro_result "$(pre)_$(pn)_l$(l)_npeaks" npk
        @repro_result "$(pre)_$(pn)_l$(l)_adaptive" TDS.fit_tail_adaptive(τs, ψ; t0 = 800.0)
        pl, nmin, nmax, nwin = TDS.fit_tail_sliding(τs, ψ)
        @repro_result "$(pre)_$(pn)_l$(l)_plateau" pl
        @repro_result "$(pre)_$(pn)_l$(l)_nmin" nmin
        @repro_result "$(pre)_$(pn)_l$(l)_nmax" nmax
        @repro_result "$(pre)_$(pn)_l$(l)_nwin" nwin
    end
    return nothing
end
