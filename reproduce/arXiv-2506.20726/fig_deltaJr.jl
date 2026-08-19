# arXiv:2506.20726 fig:deltaJr (paper.tex:663-676): the linear-in-spin radial-action correction at
# fixed orbital parameters, δJ_r|_{p,e,x}, vs eccentricity e (log-log) for four orbits — all O(e²).
#
# THIN driver: δJ comes from `Spacetime.SpinningSecondary.spin_actions_pex`, promoted in this
# migration out of the legacy figure-side helper `results/gw-flux/figures/spin_actions.jl`.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2506.20726/fig_deltaJr.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Spacetime.SpinningSecondary: spin_actions_pex
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

# (label, a, p, x, colour, linestyle) — the four published curves
const ORBITS = [("(0.99,8,0.8)", 0.99, 8.0, 0.8, MMA_BLUE,  :dashdot),
                ("(0.99,8,0.1)", 0.99, 8.0, 0.1, MMA_GREEN, :dash),
                ("(0.99,3,0.8)", 0.99, 3.0, 0.8, :orange,   :dot),
                ("(0.10,8,0.8)", 0.10, 8.0, 0.8, MMA_RED,   :solid)]

es = 10 .^ range(-5, log10(0.5), length = 40)

data = Dict{String,Vector{Float64}}()
for (lbl, a, p, x, _, _) in ORBITS
    y = fill(NaN, length(es))
    for (i, e) in enumerate(es)
        y[i] = try
            abs(spin_actions_pex(a, p, e, x)[1])
        catch
            NaN
        end
    end
    data[lbl] = y
    i1 = argmin(abs.(es .- 0.1)); i3 = argmin(abs.(es .- 1e-3))
    slope = log(y[i1] / y[i3]) / log(es[i1] / es[i3])
    println("  ", lbl, "  δJr(e=0.1)=", y[i1], "  log-log slope(1e-3→0.1)=", slope)
end

plt = plot(xlabel = "e", ylabel = raw"$\delta J_r$", xscale = :log10, yscale = :log10,
           legend = :bottomright, xlims = (1e-5, 0.5), ylims = (1e-11, 1e-1), legendfontsize = 7)
for (lbl, _, _, _, c, ls) in ORBITS
    plot!(plt, es, data[lbl], color = c, linestyle = ls, lw = 2, label = lbl)
end
savefig(plt, joinpath(@__DIR__, "deltaJr.png"))

slopes = Float64[]
@repro_result "e" collect(es)
for (lbl, _, _, _, _, _) in ORBITS
    @repro_result "deltaJr_$(lbl)" data[lbl]
    i1 = argmin(abs.(es .- 0.1)); i3 = argmin(abs.(es .- 1e-3))
    push!(slopes, log(data[lbl][i1] / data[lbl][i3]) / log(es[i1] / es[i3]))
end
@repro_result "loglog_slopes" slopes
@repro_result "all_slopes_are_two" all(s -> isapprox(s, 2.0; atol = 0.02), slopes)

exit(ReproHarness.conclude!(opts; paper = "2506.20726", id = "fig_deltaJr", script = @__FILE__,
                            rtol = 1e-9))
