# arXiv:1301.1591 — Campaign-A vacuum s = −2 Kerr tail at a = 0.9 M (the DG figure-set analogues
# `F-S1.4_Tails_Kerr_s-2_a0.9_finite` and `F-S1.4s_Tails_Kerr_s-2_a0.9_scri` of the legacy
# reproduction, emri-disk `results/dg-teukolsky-spin/`).
#
# WHAT IS REPRODUCED: the empirical ID0 late-time tail laws of the paper, Eqs. (28),(29): a seed
# ℓ′ = 2, m = 0, s = −2 perturbation decays as τ^-7 at finite radius and τ^-6 at scri (ℓ = 2), with
# the ℓ-ladder −8/−7 and −9/−8. The measured exponents are registered at all three probes
# (horizon R₊, mid-domain, scri) for ℓ = 2,3,4, with the raw P4-gate-window fit, the
# oscillation-robust envelope fit and the sliding-window plateau/scatter.
#
# HONEST BOUNDARY (measured, not tuned): the committed |s| = 2 defaults carry a κ-mediated numerical
# oscillation that rides every s = −2 tail, so a raw fit lands within ~0.7 of the law rather than on
# it; that is the P4 gate tolerance (±1.0) and it is why the envelope estimator exists. No knob is
# tuned per case — the window (700,1300) and the config are the committed gate's.
#
#   julia --project=reproduce -t4 reproduce/arXiv-1301.1591/tails_a0p9.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

include(joinpath(@__DIR__, "common_tails.jl"))
using Printf
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))

const A = 0.9
const K = 24

t0 = time()
τs, Ψ, maxu, dt, excise = run_tail(A, K)
println("[repro] a=$(A) K=$(K) excise=$(round(excise, digits=5)) dt=$(dt) wall=$(round(time()-t0, digits=1)) s")

@repro_result "excise" excise
@repro_result "n_samples" length(τs)
@repro_result "tau_end" τs[end]
register_tail_fits!(spin_tag(A), τs, Ψ, maxu, dt; series_stride = 10, series_ells = (2,))
# the two published law anchors of this figure pair (ℓ = 2): finite −7, scri −6
@repro_result "law_finite_l2" LAW_ID0_LP2[2][1]
@repro_result "law_scri_l2"   LAW_ID0_LP2[2][2]

# ---- figure: log10|ψ_ℓ| vs log10 τ at the finite and scri probes -------------------------------
# Presentation follows the legacy renders `F-S1.4_Tails_Kerr_s-2_a0.9_{finite,scri}.png`
# (VLM-compared): same quantity and frame (log10-log10), same ψ₂/ψ₃/ψ₄ colour-and-dash scheme, same
# "(expected law, measured envelope LPI)" legend. The legacy inset (running envelope LPI vs log10 τ)
# is NOT reproduced here — the running LPI is a diagnostic of the same envelope fit that the golden
# already stores per ℓ.
if opts.mode === :regen
    emri_defaults!()
    tag = spin_tag(A)
    styles = ((:red, :solid), (:black, :dash), (:blue, :dashdot))
    sel = findall(t -> t >= 60.0, τs)
    lt = log10.(τs[sel])
    plt = plot(layout = (1, 2), size = (980, 400), legend = :topright, legendfontsize = 7,
               left_margin = 6Plots.mm, bottom_margin = 6Plots.mm, top_margin = 4Plots.mm,
               grid = true, gridalpha = 0.25, gridstyle = :dot)
    for (k, (ip, pname, panel)) in enumerate(((2, "mid", "finite radius (R_mid)"), (3, "scri", "scri+")))
        for (q, l) in enumerate(ELLS)
            c, ls = styles[q]
            lp = log10.(max.(abs.(Ψ[ip, q, sel]), 1e-20))
            law = (l == 2 ? LAW_ID0_LP2[2] : l == 3 ? LAW_ID0_LP2[3] : LAW_ID0_LP2[4])[k == 1 ? 1 : 2]
            env, _ = TDS.fit_tail_envelope(τs, abs.(Ψ[ip, q, :]), P4_WINDOW)
            plot!(plt[k], lt, lp, c = c, ls = ls, lw = 0.8,
                  label = @sprintf("psi_%d  (%d, %.3f)", l, law, env))
        end
        plot!(plt[k], title = "a = 0.9, s = -2, ID0 (2,0) seed: " * panel,
              titlefontsize = 10, xlabel = "log10 tau",
              ylabel = k == 1 ? "log10 |psi_l(tau, R)|" : "", ylim = (-20, 2))
    end
    png = joinpath(@__DIR__, "tails_a0p9.png")
    savefig(plt, png)
    println("[repro] wrote ", png, " (", filesize(png), " bytes)")
end

exit(ReproHarness.conclude!(opts; paper = "1301.1591", id = "tails_a0p9", script = @__FILE__,
                            rtol = 1e-8))
