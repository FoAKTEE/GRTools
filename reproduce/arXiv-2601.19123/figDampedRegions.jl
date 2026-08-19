# arXiv:2601.19123 fig:DampedRegions (F3) — "Excitation Regions for r_p = 15M": for each spin
# (a = 0.01, 0.60, 0.99 and the Newtonian limit) the m = 1…15 Lindblad intervals [r_ILR, r_OLR],
# drawn as overlapping semi-transparent bands that darken toward r_p.
#
# THIN driver: `AccretionDisc.lindblad_locations` (roots of m(Ω^k_r − Ω^k_{r_p}) = ±κ_r). The
# Newtonian row is the analytic limit of the same condition, ς_m^∓ = ((m∓1)/m)^{2/3} — the limit
# checked to 3e-6 by reproduce/arXiv-1010.0758/check_RT8_GT80.jl (sigma_ILR/OLR).
#
#   julia --project=reproduce -t4 reproduce/arXiv-2601.19123/figDampedRegions.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, lindblad_locations
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
using Plots

const RP = 15.0
const MS = 1:15
rows = [("\$a=0.01\$", 0.01, RGB(0.31, 0.51, 0.55)),
        ("\$a=0.60\$", 0.60, RGB(0.78, 0.56, 0.35)),
        ("\$a=0.99\$", 0.99, RGB(0.48, 0.36, 0.55))]

emri_defaults!()
plt = plot(; title = "Excitation Regions for \$r_p=15M\$", xlabel = "\$r/M\$",
           xlims = (5, 29), ylims = (0, 4), yticks = false, legend = :topright,
           xticks = 6:2:28, size = (720, 360))

for (i, (lbl, a, col)) in enumerate(rows)
    d = RelativisticDisc(a = a)
    rin = Float64[]; rout = Float64[]
    for m in MS
        rI, rO = lindblad_locations(d, RP, m)
        push!(rin, rI); push!(rout, rO)
    end
    @repro_result "rILR_a$(a)" rin
    @repro_result "rOLR_a$(a)" rout
    y0 = 4 - i
    for m in MS
        ri = isnan(rin[m]) ? 5.0 : rin[m]
        isnan(rout[m]) && continue
        plot!(plt, Shape([ri, rout[m], rout[m], ri], [y0+0.05, y0+0.05, y0+0.9, y0+0.9]);
              fillcolor = col, fillalpha = 0.13, linealpha = 0, label = "")
    end
    plot!(plt, [NaN], [NaN]; seriestype = :scatter, mc = col, ms = 5, label = lbl)
    @printf("%-12s band (m=1) = [%.3f, %.3f];  (m=15) = [%.3f, %.3f]\n",
            lbl, rin[1], rout[1], rin[15], rout[15])
end

# Newtonian row: the analytic ς_m^∓ limit of the same condition
newt_in  = [m == 1 ? NaN : ((m-1)/m)^(2/3)*RP for m in MS]
newt_out = [((m+1)/m)^(2/3)*RP for m in MS]
@repro_result "rILR_newtonian" newt_in
@repro_result "rOLR_newtonian" newt_out
let col = RGB(0.50, 0.55, 0.60), y0 = 0
    for m in MS
        ri = isnan(newt_in[m]) ? 5.0 : newt_in[m]
        plot!(plt, Shape([ri, newt_out[m], newt_out[m], ri], [y0+0.05, y0+0.05, y0+0.9, y0+0.9]);
              fillcolor = col, fillalpha = 0.13, linealpha = 0, label = "")
    end
    plot!(plt, [NaN], [NaN]; seriestype = :scatter, mc = col, ms = 5, label = "Newton")
end
vline!(plt, [RP]; lc = :black, ls = :dash, lw = 2, label = "")
annotate!(plt, 6.0, 0.5, text("\$r_{LR}^{(m=1)}=0 \\rightarrow\$", 9, :left))
png = joinpath(@__DIR__, "expected", "figDampedRegions.png")
savefig(plt, png)
@printf("wrote %s (%d bytes)\n", relpath(png, @__DIR__), filesize(png))
@repro_result "r_p" RP

exit(ReproHarness.conclude!(opts; paper = "2601.19123", id = "figDampedRegions",
                            script = @__FILE__, rtol = 1e-9))
