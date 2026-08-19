# arXiv:2603.13482 fig:delta_Omega (main.tex:476): the linear-in-spin Boyer–Lindquist frequency
# corrections Ω¹_i in FOUR spin gauges (C̃ = fixed shifted constants, p = fixed turning points,
# C = fixed constants of motion, J = fixed actions) as a function of p − p_s near the separatrix,
# for a = 0.95 M, e = 0.5, x = 0.5. Three stacked log-log panels (Ω¹_r, Ω¹_z, Ω¹_φ).
#
# PUBLISHED MESSAGE: the C, J and p gauges diverge as the separatrix is approached (Υ_r → 0
# degenerates the constants/action Jacobians and dΩ/dp ∼ (p−p_s)^(−1/2)), while the C̃ gauge stays
# FINITE — the regularity statement of the paper.
#
# THIN driver: all four gauges come from Spacetime.SpinningSecondary
# (`spin_linear_frequencies`, `spin_frequencies_turning_point_gauge`,
# `spin_frequencies_constants_gauge`, `spin_frequencies_action_gauge`) and the separatrix from
# `Spacetime.separatrix` — the legacy figure script reconstructed the p and C gauges from
# finite-difference Jacobians inline and bisected its own separatrix.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2603.13482/fig_delta_Omega.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Spacetime: separatrix
using GrTools.Spacetime.SpinningSecondary: spin_linear_frequencies,
                                           spin_frequencies_turning_point_gauge,
                                           spin_frequencies_constants_gauge,
                                           spin_frequencies_action_gauge
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const A, E, X = 0.95, 0.5, 0.5

ps = separatrix(A, E, X)
println("  separatrix p_s = ", ps)

dps = 10 .^ range(-4, 1, length = 60)

const GAUGES = (("Ctilde", spin_linear_frequencies,                 MMA_BLUE,  :solid,   "~C"),
                ("p",      spin_frequencies_turning_point_gauge,    MMA_GREEN, :dot,     "p"),
                ("C",      spin_frequencies_constants_gauge,        MMA_GOLD,  :dash,    "C"),
                ("J",      spin_frequencies_action_gauge,           MMA_RED,   :dashdot, "J"))

data = Dict{String,Matrix{Float64}}()          # name => 3 × length(dps) of |Ω¹_i|
for (name, f, _, _, _) in GAUGES
    M = fill(NaN, 3, length(dps))
    for (i, dp) in enumerate(dps)
        v = try
            collect(f(A, ps + dp, E, X))
        catch
            [NaN, NaN, NaN]
        end
        M[:, i] = abs.(v)
    end
    data[name] = M
    println("  ", rpad(name, 7), " |Ω¹| @ p−p_s=1e-4 : ", M[:, 1], "   @ 10 : ", M[:, end])
end

ylabs = (raw"$\Omega_r^1$", raw"$\Omega_z^1$", raw"$\Omega_\phi^1$")
panels = Any[]
for c in 1:3
    pl = plot(xscale = :log10, yscale = :log10, xlabel = c == 3 ? raw"$p-p_s$" : "",
              ylabel = ylabs[c], legend = c == 1 ? :topright : false,
              xlims = (1e-4, 10), ylims = (1e-4, 1e6))
    for (name, _, col, ls, lab) in GAUGES
        y = data[name][c, :]
        msk = isfinite.(y) .& (y .> 0)
        plot!(pl, dps[msk], y[msk], color = col, linestyle = ls, lw = 2, label = lab)
    end
    push!(panels, pl)
end
savefig(plot(panels..., layout = (3, 1), size = (560, 760)),
        joinpath(@__DIR__, "delta_Omega.png"))

@repro_result "p_separatrix" ps
@repro_result "dp"           collect(dps)
for (name, _, _, _, _) in GAUGES
    @repro_result "absOmega1_$(name)" data[name]
end
# the published message, as MEASURED magnitudes + the growth ratios near/far over the scanned range.
# The C̃ gauge is BOUNDED: max |Ω¹| = 0.049 over the whole scan (its near/far ratio 10.3 is the mild
# bump at p − p_s ≈ 1e-4…1e-3 that the published panel also shows), while the C and J gauges reach
# 1e5…1e6 and the turning-point gauge 35.
growth = Dict(name => maximum(data[name][:, 1]) / maximum(data[name][:, end]) for (name, _, _, _, _) in GAUGES)
maxabs = Dict(name => maximum(data[name]) for (name, _, _, _, _) in GAUGES)
for (name, _, _, _, _) in GAUGES
    @repro_result "growth_near_over_far_$(name)" growth[name]
    @repro_result "max_absOmega1_$(name)"        maxabs[name]
end
@repro_result "Ctilde_bounded"           all(isfinite, data["Ctilde"]) && maxabs["Ctilde"] < 1
@repro_result "C_and_J_diverge_strongly" maxabs["C"] > 1e4 && maxabs["J"] > 1e4
@repro_result "p_gauge_milder_than_C"    maxabs["p"] < maxabs["C"]
@repro_result "Ctilde_bounded_while_others_diverge" maxabs["Ctilde"] < 1 && maxabs["C"] > 1e4 &&
                                                    maxabs["J"] > 1e4 && maxabs["p"] > 1

exit(ReproHarness.conclude!(opts; paper = "2603.13482", id = "fig_delta_Omega",
                            script = @__FILE__, rtol = 1e-8))
