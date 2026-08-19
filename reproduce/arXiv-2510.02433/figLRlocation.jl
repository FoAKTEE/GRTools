# arXiv:2510.02433 fig:LRlocations — the resonance offset |r_LR − r_s| vs r_s for m = 2, 3, 4, 8,
# outer (solid) and inner (dashed) Lindblad resonances, at a/M = 0 and a/M = 0.9.
#
# THIN driver: `AccretionDisc.lindblad_locations` (roots of m(Ω_φ−Ω_s) = ±κ on the exact Kerr
# equatorial geodesic frequencies), `Spacetime.kerr_isco` for the inner edge.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2510.02433/figLRlocation.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, lindblad_locations
using GrTools.Spacetime: kerr_isco
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
using Plots

const MS = (2, 3, 4, 8)
const COL = Dict(2 => RGB(0.282, 0.471, 0.816), 3 => RGB(0.933, 0.522, 0.290),
                 4 => RGB(0.839, 0.373, 0.373), 8 => RGB(0.584, 0.424, 0.706))
emri_defaults!()

panels = Any[]
for a in (0.0, 0.9)
    d = RelativisticDisc(a = a)
    rsg = exp10.(range(log10(kerr_isco(a) + 0.2), log10(100.0); length = 60))
    @repro_result "rs_grid_a$(a)" rsg
    plt = plot(xscale = :log10, yscale = :log10, xlabel = "\$r_s/M\$",
               ylabel = "\$|r_\\mathrm{LR}-r_s|/M\$", title = "a/M = $a", legend = :bottomright)
    for m in MS
        dout = Float64[]; din = Float64[]
        for rs in rsg
            rI, rO = lindblad_locations(d, rs, m)
            push!(dout, rO - rs); push!(din, rs - rI)
        end
        @repro_result "dout_a$(a)_m$(m)" dout
        @repro_result "din_a$(a)_m$(m)" din
        plot!(plt, rsg, dout; color = COL[m], lw = 2, ls = :solid, label = "m=$m (OLR)")
        plot!(plt, rsg, din;  color = COL[m], lw = 2, ls = :dash, label = "")
    end
    push!(panels, plt)
    @printf("a=%.1f: ISCO=%.4f; at r_s=20 the m=2 offsets are OLR %.4f / ILR %.4f\n",
            a, kerr_isco(a), lindblad_locations(d, 20.0, 2)[2] - 20, 20 - lindblad_locations(d, 20.0, 2)[1])
end
png = joinpath(@__DIR__, "expected", "figLRlocation.png")
savefig(plot(panels...; layout = (1, 2), size = (1100, 450), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm), png)
@printf("wrote %s (%d bytes)\n", relpath(png, @__DIR__), filesize(png))
@repro_result "isco_a0" kerr_isco(0.0)
@repro_result "isco_a09" kerr_isco(0.9)

exit(ReproHarness.conclude!(opts; paper = "2510.02433", id = "figLRlocation",
                            script = @__FILE__, rtol = 1e-9))
