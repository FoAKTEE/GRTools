# arXiv:2303.16798 fig:comparison_generic, panels error_eccentric_l2m2 (top) and _l5m4 (bottom)
# (paper.tex:587): the per-mode relative difference ΔF^E_{S,lmnk} = |1 − F_{lmn0}(I)/F_{lmn0}(0)|
# vs inclination I ∈ {0.5°, 1°, 2°, 4°, 8°}, for n = 0..4, k = 0, at a = 0.9 M, p = 12, e = 0.3 —
# the O(I²) convergence of the inclined mode flux to its equatorial limit (log-log, slope-2 guides).
#
# THIN driver, both fluxes from Perturbation.SpinFluxes with the SAME ω-derivative step h:
#   F(I) = TeukolskySpinModeCorrectionAnalytical(l, m, n, 0, a, p, e, cos I)  (s∥ convention)
#   F(0) = equatorial_spin_flux(l, m, n, a, p, e)
# h = 1e-5 (not the 3e-4 default): the default floors the ΔF ≈ 1e-5…1e-4 relative differences with
# the linearisation error of the frequency derivative.
#
# KNOWN GAP (measured in the legacy reproduction, error_eccentric.spec.md Rule 3, and unchanged
# here): the published panels show the n-curves spread over ~22×, while the validated flux gives a
# nearly n-independent ΔF. The inclined per-mode flux is oracle-matched against the BHPToolkit
# Mathematica source at (2,2,0,0) and (2,2,1,0), so the spread comes from an undocumented per-n
# figure-construction step, not from the flux. The slope-2 message reproduces exactly.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2303.16798/fig_error_eccentric.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Perturbation: equatorial_spin_flux, TeukolskySpinModeCorrectionAnalytical
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const A, P, E = 0.9, 12.0, 0.3
const IDEG = (0.5, 1.0, 2.0, 4.0, 8.0)
const HFD = 1e-5
const NCOLS = (MMA_BLUE, MMA_GOLD, MMA_GREEN, MMA_RED, MMA_PURPLE)
const PANELS = ((2, 2, "l2m2"), (5, 4, "l5m4"))

Feq(l, m, n) = (f = equatorial_spin_flux(l, m, n, A, P, E; h = HFD); real(f.dE1I + f.dE1H))
FI(l, m, n, Ideg) = (md = TeukolskySpinModeCorrectionAnalytical(l, m, n, 0, A, P, E,
                                                               cos(Ideg * pi / 180); h = HFD);
                     real(md.FluxesCorrection.Energy.I + md.FluxesCorrection.Energy.H))

Irad = [Id * pi / 180 for Id in IDEG]
for (l, m, tag) in PANELS
    f0 = [Feq(l, m, n) for n in 0:4]
    jobs = [(n, Id) for n in 0:4 for Id in IDEG]
    dF = Vector{Float64}(undef, length(jobs))
    Threads.@threads for j in eachindex(jobs)
        n, Id = jobs[j]
        dF[j] = abs(1 - FI(l, m, n, Id) / f0[n + 1])
    end

    plt = plot(xlabel = "I", ylabel = raw"$\Delta\mathcal{F}^E_{S,lmnk}$",
               xscale = :log10, yscale = :log10, legend = :topleft,
               xticks = ([0.01, 0.02, 0.05, 0.10], ["0.01", "0.02", "0.05", "0.10"]))
    Iline = range(0.006, 0.16, length = 50)
    for C in (0.03, 0.3, 3.0, 30.0, 300.0)
        plot!(plt, collect(Iline), C .* Iline .^ 2, line = (:dash, :gray, 1), label = "")
    end
    slopes = Float64[]
    for n in 0:4
        y = [dF[j] for j in eachindex(jobs) if jobs[j][1] == n]
        scatter!(plt, Irad, y, marker = (:circle, 5), color = NCOLS[n + 1], label = "n = $(n)")
        s = (log(y[end]) - log(y[1])) / (log(Irad[end]) - log(Irad[1]))
        push!(slopes, s)
        println("  ", tag, " n=", n, "  ΔF=", y, "  slope=", s); flush(stdout)
        @repro_result "dF_$(tag)_n$(n)" y
    end
    ylims!(plt, 1e-6, 5e-2)
    savefig(plt, joinpath(@__DIR__, "error_eccentric_$(tag).png"))

    @repro_result "F_equatorial_$(tag)" f0
    @repro_result "slopes_$(tag)"       slopes
    @repro_result "slopes_near_two_$(tag)" all(s -> isapprox(s, 2.0; atol = 0.05), slopes)
    # the measured per-n spread at I = 8° — the published panels show ≈ 22×, see the gap note above
    y8 = [dF[j] for j in eachindex(jobs) if jobs[j][2] == last(IDEG)]
    @repro_result "per_n_spread_at_8deg_$(tag)" maximum(y8) / minimum(y8)
end

@repro_result "I_rad" Irad
@repro_result "h_fd"  HFD

exit(ReproHarness.conclude!(opts; paper = "2303.16798", id = "fig_error_eccentric",
                            script = @__FILE__, rtol = 1e-8))
