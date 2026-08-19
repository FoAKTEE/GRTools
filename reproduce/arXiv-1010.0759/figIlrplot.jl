# arXiv:1010.0759 fig:ilrplot — Schwarzschild m = 1, 2, 3 inner-Lindblad resonances:
#   LEFT  : locations r₁(r₀) (roots of D(r) = m(Ω−Ω₀) − κ)
#   RIGHT : strengths log₁₀|N|(r₀) from the FULL relativistic Teukolsky overlap
#           (eq:P-ILR + eq:SP + eq:Tnorm), with the tab:ilrS markers overlaid.
#
# THIN driver: `AccretionDisc.lindblad_locations` + `AccretionDisc.lindblad_torque` (the ℓ-sum runs
# through the GrTools Perturbation MST engine and the Harmonics spheroidal eigenvalue).
#
# GRID: 12 points/curve on 8 ≤ r₀ ≤ 250 (the published figure shows ~16). The strength evaluation
# costs ~5 s per (m, r₀), so a denser grid belongs to the heavy tier.
#
#   julia --project=reproduce -t4 reproduce/arXiv-1010.0759/figIlrplot.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, lindblad_locations, lindblad_torque
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
using Plots

const D0 = RelativisticDisc(a = 0.0)
const MS = (1, 2, 3)
r0s = exp10.(range(log10(8.0), log10(250.0); length = 12))
@repro_result "r0_grid" r0s

locs = Dict(m => [first(lindblad_locations(D0, r0, m)) for r0 in r0s] for m in MS)
for m in MS
    @repro_result "r1_m$(m)" locs[m]
end

Ns = Dict{Int,Vector{Float64}}()
for m in MS
    t0 = time()
    Ns[m] = [lindblad_torque(D0, r0, m; kind = :ILR) for r0 in r0s]
    @repro_result "N_m$(m)" Ns[m]
    @printf("m=%d strengths done in %.1f s\n", m, time() - t0)
end

# tab:ilrS golden markers
gold_r0   = [8.0, 10.0, 20.0, 50.0, 100.0, 250.0]
gold_r1m1 = [6.48, 7.13, 10.14, 16.91, 25.25, 43.25]
gold_r1m2 = [6.97, 8.19, 14.50, 33.44, 64.97, 159.48]
gold_Nm1  = [7.11e-1, 1.10e-1, 6.12e-3, 4.86e-4, 8.33e-5, 8.35e-6]
gold_Nm2  = [4.69e1, 1.74e1, 5.37, 3.22, 2.75, 2.50]

emri_defaults!()
lstyle = Dict(3 => :solid, 2 => :dash, 1 => :dot)
mshape = Dict(3 => :cross, 2 => :xcross, 1 => :star5)
pL = plot(title = "ILR locations in Schwarzschild", xlabel = "r_0 (secondary)", ylabel = "r_1 (resonance)",
          xscale = :log10, yscale = :log10, xlims = (7, 300), ylims = (6, 200), legend = :topleft, color = :black)
for m in (3, 2, 1)
    plot!(pL, r0s, locs[m]; line = (lstyle[m], :black), marker = (mshape[m], 3, :black), label = "m=$m")
end
scatter!(pL, gold_r0, gold_r1m2; marker = (:xcross, 5, :black), label = "")
scatter!(pL, gold_r0, gold_r1m1; marker = (:star5, 5, :black), label = "")
pR = plot(title = "ILR strengths in Schwarzschild", xlabel = "r_0 (secondary)", ylabel = "log_{10} |N|",
          xscale = :log10, xlims = (7, 300), ylims = (-5, 2), legend = :topright, color = :black)
for m in (3, 2, 1)
    plot!(pR, r0s, log10.(abs.(Ns[m])); line = (lstyle[m], :black), marker = (mshape[m], 3, :black), label = "m=$m")
end
scatter!(pR, gold_r0, log10.(gold_Nm2); marker = (:xcross, 5, :black), label = "")
scatter!(pR, gold_r0, log10.(gold_Nm1); marker = (:star5, 5, :black), label = "")
png = joinpath(@__DIR__, "expected", "figIlrplot.png")
savefig(plot(pL, pR; layout = (1, 2), size = (1100, 480), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm), png)
@printf("wrote %s (%d bytes)\n", relpath(png, @__DIR__), filesize(png))

# published asymptotics: N(m=2) → −2.36, N(m=3) → −7.50
@repro_result "N_m2_at_r0_250" lindblad_torque(D0, 250.0, 2; kind = :ILR)
@repro_result "N_m3_at_r0_250" lindblad_torque(D0, 250.0, 3; kind = :ILR)
@printf("asymptotics r₀=250: N(m=2)=%.3f (paper −2.36), N(m=3)=%.3f (paper −7.50)\n",
        lindblad_torque(D0, 250.0, 2; kind = :ILR), lindblad_torque(D0, 250.0, 3; kind = :ILR))

exit(ReproHarness.conclude!(opts; paper = "1010.0759", id = "figIlrplot",
                            script = @__FILE__, rtol = 1e-8))
