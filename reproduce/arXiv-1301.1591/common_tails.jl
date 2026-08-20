# arXiv:1301.1591 — CONFIGURATION shared by the vacuum s = −2 Kerr tail drivers.
#
# ZERO PHYSICS: every number below is a run configuration (the P4-proven gate config of
# `Perturbation/test/test_timedomain_spin.jl` gate C+E) or a fit window; every computation is a
# call into `Perturbation.TimeDomain` (`SpinModesProblem`, `solve_spin_modes`, `spin_cfl_dt`,
# `spin_excise_policy`, `fit_tail_window`, `fit_tail_envelope`, `fit_tail_sliding`).
#
# Config provenance: legacy driver `results/dg-teukolsky-spin/scripts/run_spin_tails.jl`
# (N = 8, filter_alpha = 36 for |s| = 2, ID0 seed (ℓ=2, m=0), R0 = 0.75 for a ≠ 0, w = 1500,
# CFL C = 0.6, ~1M sample cadence, ℓ = 2..4 coupled).

using GrTools
using GrTools.Perturbation: TimeDomain as TDS

const S_SPIN   = -2
const M_AZIM   = 0
const ELLS     = [2, 3, 4]
const N_ORDER  = 8
const FILTER_A = 36.0
const CFL_C    = 0.6
const IC_W     = 1500.0
const TFINAL   = 3000.0
"P4 committed gate window — pre-declared, never tuned per case [arXiv:1301.1591 campaign A]."
const P4_WINDOW = (700.0, 1300.0)
"Empirical ID0 tail laws [arXiv:1301.1591, Eqs. (28),(29)] for a seed ℓ′ = 2: (finite, scri)."
const LAW_ID0_LP2 = Dict(2 => (-7, -6), 3 => (-8, -7), 4 => (-9, -8))
const PROBES = ("hor", "mid", "scri")

"Spin tag used in registered result names: 0.9 → `a0p9`."
spin_tag(a) = "a" * replace(string(a), "." => "p")

"""
    run_tail(a, K; tfinal = TFINAL) -> (τs, Ψ, maxu, dt, excise)

One vacuum s = −2 tail run at spin `a` with `K` elements. The excision width comes from the
package policy `TDS.spin_excise_policy`; the mesh is asserted to have RETAINED both the excised
inner edge and R₊ (build_mesh snaps a forced coordinate to the nearest uniform edge, so a too
small excision is silently deleted).
"""
function run_tail(a, K; tfinal = TFINAL)
    pol = TDS.spin_excise_policy(1.0, a)
    p = TDS.SpinModesProblem(S_SPIN, ELLS, M_AZIM; N = N_ORDER, K = K, M = 1.0, a = a,
                             excise = pol.excise, filter_alpha = FILTER_A)
    Rin = pol.Rplus - pol.excise
    @assert abs(p.mesh.edges[1] - Rin) < 1e-12 "excision lost (raise K): mesh starts at $(p.mesh.edges[1]) ≠ $Rin"
    @assert any(e -> abs(e - pol.Rplus) < 1e-12, p.mesh.edges) "R₊ is not a mesh edge"
    dt = TDS.spin_cfl_dt(p; C = CFL_C)
    rec = max(1, round(Int, 1.0 / dt))
    τs, Ψ, maxu = TDS.solve_spin_modes(p; tfinal = tfinal, dt = dt, record_every = rec,
                                       id = :ID0, init_mode = 1,
                                       ic = (amp = 1.0, R0 = (a == 0 ? 0.8 : 0.75), w = IC_W))
    return (τs, Ψ, maxu, dt, pol.excise)
end

"""
    register_tail_fits!(tag, τs, Ψ, maxu, dt; ells = ELLS, series_stride = 0, series_ells = (2,))

Register every measured number of one tail run: the raw P4-window exponent, the envelope
exponent, the sliding-window plateau/min/max per probe and ℓ, plus the run's `dt`, terminal
field norm and boundedness flag. With `series_stride > 0` a decimated |ψ| record is registered
too (the plotted array of the published-figure analogue), for `series_ells`.
"""
function register_tail_fits!(tag, τs, Ψ, maxu, dt; ells = ELLS, series_stride = 0,
                             series_ells = (2,))
    @repro_result "$(tag)_dt" dt
    @repro_result "$(tag)_maxu_end" maxu[end]
    @repro_result "$(tag)_maxu_finite" all(isfinite, maxu)
    for (ip, pn) in enumerate(PROBES), (q, l) in enumerate(ells)
        ψ = abs.(Ψ[ip, q, :])
        @repro_result "$(tag)_$(pn)_l$(l)_raw" TDS.fit_tail_window(τs, ψ, P4_WINDOW)
        nenv, npk = TDS.fit_tail_envelope(τs, ψ, P4_WINDOW)
        @repro_result "$(tag)_$(pn)_l$(l)_env" nenv
        @repro_result "$(tag)_$(pn)_l$(l)_npeaks" npk
        pl, nmin, nmax, nwin = TDS.fit_tail_sliding(τs, ψ)
        @repro_result "$(tag)_$(pn)_l$(l)_plateau" pl
        @repro_result "$(tag)_$(pn)_l$(l)_nmin" nmin
        @repro_result "$(tag)_$(pn)_l$(l)_nmax" nmax
        @repro_result "$(tag)_$(pn)_l$(l)_nwin" nwin
    end
    if series_stride > 0
        idx = 1:series_stride:length(τs)
        @repro_result "$(tag)_tau" collect(τs[idx])
        for (ip, pn) in enumerate(PROBES), l in series_ells
            q = findfirst(==(l), ells)
            @repro_result "$(tag)_$(pn)_l$(l)_abs" abs.(Ψ[ip, q, idx])
        end
    end
    return nothing
end
