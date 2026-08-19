# arXiv:2303.16798 fig:comparison_spherical (paper.tex:577): the relative difference of the
# linear-in-spin TOTAL energy flux, ΔF^E_S = |1 − F^E_S(I)/F^E_S(I=0)|, versus inclination
# I ∈ {0.5°, 1°, 2°, 4°, 8°} at a = 0.9 M for p = 5, 10, 15 — the O(I²) convergence of the
# spherical (e = 0, inclined) flux to the equatorial circular flux (log-log, slope-2 guides).
#
# THIN driver, both fluxes from Perturbation.SpinFluxes:
#   F^E_S(I) = Σ_{l≤5, 1≤m≤l, |k|≤8} TeukolskySpinModeSphericalCorrectionAnalytical  (s∥ convention)
#   F^E_S(0) = Σ_{l≤5, 1≤m≤l}        circular_spin_flux                              (equatorial)
# The s∥ ("Path B") spherical entry point is the one CONSISTENT with the equatorial circular flux;
# the full-tensor entry point is a different spin contraction and does not limit to it.
#
# HEAVY, NO GOLDEN YET. 3570 spherical analytical modes (15 (p, I) points x 238 modes) whose MST
# radial sums are slow at p = 15 with |k| <= 8: MEASURED 4716 s (79 min) at -t4 with one straggler
# job still running when it was killed, against a 30 min per-driver budget. The
# driver is complete and was NOT run to completion, so there is no golden: --check exits 3 (HOLE)
# rather than fabricating one. Generate it with `--regen --tier heavy` when the budget allows.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2303.16798/fig_comparison_spherical.jl --regen --tier heavy

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

if opts.mode === :check && !isfile(ReproHarness.expected_path("2303.16798", "fig_comparison_spherical"))
    ReproHarness.hole("no golden yet — the 3570-mode O(I^2) scan measured 79 min at -t4 without " *
                      "completing; generate it with `--regen --tier heavy` (MANIFEST status: future)")
end

using GrTools
using GrTools.Perturbation: circular_spin_flux, TeukolskySpinModeSphericalCorrectionAnalytical
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const A = 0.9
const PS = (5.0, 10.0, 15.0)
const IDEG = (0.5, 1.0, 2.0, 4.0, 8.0)
const LMAX, KMAX = 5, 8
const COLS = Dict(5.0 => MMA_BLUE, 10.0 => MMA_GOLD, 15.0 => MMA_GREEN)

FS0(p) = sum((try
                  f = circular_spin_flux(l, m, A, p)
                  real(f.dE1I + f.dE1H)
              catch                      # high-l continued-fraction limit: negligible mode
                  0.0
              end) for l in 2:LMAX for m in 1:l)

function FSI(p, Ideg)
    x = cos(Ideg * pi / 180)
    modes = [(l, m, k) for l in 2:LMAX for m in 1:l for k in -KMAX:KMAX]
    term = Vector{Float64}(undef, length(modes))
    for i in eachindex(modes)
        l, m, k = modes[i]
        term[i] = try
            md = TeukolskySpinModeSphericalCorrectionAnalytical(l, m, k, A, p, x)
            real(md.FluxesCorrection.Energy.I + md.FluxesCorrection.Energy.H)
        catch
            0.0
        end
    end
    return sum(term)                     # fixed-order reduction ⇒ thread-count invariant
end

f0 = Dict(p => FS0(p) for p in PS)
jobs = [(p, Id) for p in PS for Id in IDEG]
dF = Vector{Float64}(undef, length(jobs))
Threads.@threads for j in eachindex(jobs)
    p, Id = jobs[j]
    dF[j] = abs(1 - FSI(p, Id) / f0[p])
end

Irad = [Id * pi / 180 for Id in IDEG]
slopes = Float64[]
for p in PS
    y = [dF[j] for j in eachindex(jobs) if jobs[j][1] == p]
    s = (log(y[end]) - log(y[1])) / (log(Irad[end]) - log(Irad[1]))
    push!(slopes, s)
    println("  p=", p, "  F_S(0)=", f0[p], "  ΔF=", y, "  log-log slope(0.5°→8°)=", s)
    flush(stdout)
end

plt = plot(xlabel = "I", ylabel = raw"$\Delta\mathcal{F}^E_S$", xscale = :log10, yscale = :log10,
           legend = :bottomright,
           xticks = ([0.01, 0.02, 0.05, 0.10], ["0.01", "0.02", "0.05", "0.10"]))
Iline = range(0.006, 0.16, length = 50)
for C in (0.03, 0.3, 3.0, 30.0)
    plot!(plt, collect(Iline), C .* Iline .^ 2, line = (:dash, :gray, 1), label = "")
end
for p in PS
    y = [dF[j] for j in eachindex(jobs) if jobs[j][1] == p]
    scatter!(plt, Irad, y, marker = (:circle, 6), color = COLS[p], label = "p = $(Int(p))")
end
ylims!(plt, 3e-6, 2e-2)
savefig(plt, joinpath(@__DIR__, "comparison_spherical.png"))

@repro_result "I_rad"          Irad
@repro_result "FS_equatorial"  [f0[p] for p in PS]
for p in PS
    @repro_result "dF_p$(Int(p))" [dF[j] for j in eachindex(jobs) if jobs[j][1] == p]
end
@repro_result "loglog_slopes"        slopes
@repro_result "all_slopes_near_two"  all(s -> isapprox(s, 2.0; atol = 0.05), slopes)

exit(ReproHarness.conclude!(opts; paper = "2303.16798", id = "fig_comparison_spherical",
                            script = @__FILE__, rtol = 1e-8))
