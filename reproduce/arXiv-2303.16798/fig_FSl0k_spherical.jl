# arXiv:2303.16798 fig:fluxes_spherical, panels FSl10k / FSl20k / FSl30k (paper.tex:829): the
# linear-in-spin SPHERICAL-orbit energy-flux corrections log10|F^E_{S,lm0k}| as an (l, k) heatmap
# for m = 1, 2, 3 at a = 0.9 M, p = 10, I = 30°, e = 0.
#
# PUBLISHED MESSAGE: the global maximum sits at k = l − m (lowest available l) and a warm ridge
# follows k ≈ l − m, with a deep low-flux lobe at large l and k < 0.
#
# THIN driver: every mode is `Perturbation.TeukolskySpinModeSphericalCorrectionAnalytical`
# (the s∥ convention consistent with the equatorial circular flux); the script bins and plots.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2303.16798/fig_FSl0k_spherical.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Perturbation: TeukolskySpinModeSphericalCorrectionAnalytical
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const A, P = 0.9, 10.0
const X = cos(30 * pi / 180)
const KS = -11:11
const LS = Dict(1 => 2:9, 2 => 2:9, 3 => 3:10)
const XT = Dict(1 => 2:2:8, 2 => 2:2:8, 3 => 4:2:10)

jobs = [(m, l, k) for m in 1:3 for l in LS[m] for k in KS]
vals = Vector{Float64}(undef, length(jobs))
Threads.@threads for i in eachindex(jobs)
    m, l, k = jobs[i]
    F = try
        md = TeukolskySpinModeSphericalCorrectionAnalytical(l, m, k, A, P, X)
        real(md.FluxesCorrection.Energy.I + md.FluxesCorrection.Energy.H)
    catch                                  # continued-fraction limit at large |k|: negligible mode
        0.0
    end
    vals[i] = log10(abs(F) + 1e-300)
end

ridge_frac = Float64[]
peak_l = Int[]; peak_k = Int[]
for m in 1:3
    ls = LS[m]
    Z = fill(NaN, length(KS), length(ls))       # Z[k_idx, l_idx] ⇒ heatmap(ls, ks, Z): x = l, y = k
    sub = [(i, jobs[i]) for i in eachindex(jobs) if jobs[i][1] == m]
    for (i, (_, l, k)) in sub
        Z[findfirst(==(k), KS), findfirst(==(l), ls)] = vals[i]
    end
    ok = 0
    for l in ls
        col = [(k, vals[i]) for (i, (_, ll, k)) in sub if ll == l]
        kbest = col[argmax([v for (_, v) in col])][1]
        ok += (kbest == l - m) ? 1 : 0
    end
    push!(ridge_frac, ok / length(ls))
    bi = argmax([vals[i] for (i, _) in sub])
    push!(peak_l, sub[bi][2][2]); push!(peak_k, sub[bi][2][3])
    println("  m=", m, "  GLOBAL PEAK (l,k) = (", peak_l[end], ",", peak_k[end],
            ")  expect (", first(ls), ",", first(ls) - m, ")   ridge k=l−m in ", ok, "/", length(ls), " columns")
    flush(stdout)

    hm = heatmap(collect(ls), collect(KS), Z, xlabel = "l", ylabel = "k",
                 title = raw"$\log_{10}\,|\mathcal{F}^E_{S,l" * "$(m)" * raw"0k}|$",
                 c = MMA_TEMPERATUREMAP, clims = (-30, -8), colorbar_ticks = -30:5:-10,
                 xticks = XT[m], yticks = -10:5:10, aspect_ratio = :auto)
    savefig(hm, joinpath(@__DIR__, "FSl$(m)0k.png"))
    @repro_result "log10absF_m$(m)" Z
end

@repro_result "k"               collect(KS)
@repro_result "peak_l"          peak_l
@repro_result "peak_k"          peak_k
@repro_result "ridge_fraction"  ridge_frac
# the published claim, mode by mode: the global peak is at (l = max(2,m), k = l − m)
@repro_result "peak_at_lowest_l_and_k_equals_l_minus_m" all(m -> peak_l[m] == first(LS[m]) &&
                                                                 peak_k[m] == first(LS[m]) - m, 1:3)

exit(ReproHarness.conclude!(opts; paper = "2303.16798", id = "fig_FSl0k_spherical",
                            script = @__FILE__, rtol = 1e-8))
