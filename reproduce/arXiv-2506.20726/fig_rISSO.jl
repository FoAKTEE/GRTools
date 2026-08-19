# arXiv:2506.20726 fig:rISSO (paper.tex:531), both panels:
#   LEFT  `rISSOgplot`  — the geodesic ISSO radius r_ISSOg/M vs inclination x, for a = 0.3/0.6/0.9 M.
#   RIGHT `drISSOplot`  — the secondary-spin ISSO shift δr_ISSO/(s∥ M) vs x, same three spins.
#
# THIN driver: the ISSO root-finder is `Spacetime.kerr_isso` (promoted out of the legacy figure
# scripts in this migration) and the O(s∥) radial shift is the CTilde-gauge δp of
# `Spacetime.SpinningSecondary.spin_linear_orbital_params` at the ISSO radius. Zero physics here.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2506.20726/fig_rISSO.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Spacetime: kerr_isso
using GrTools.Spacetime.SpinningSecondary: spin_linear_orbital_params
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const SPINS = (0.3, 0.6, 0.9)
xs = collect(range(-1, 1, length = 200))          # 200 points: never lands exactly on the polar x=0

risso  = Dict(a => [kerr_isso(a, x) for x in xs] for a in SPINS)
drisso = Dict(a => [spin_linear_orbital_params(a, kerr_isso(a, x), 0.0, x)[1] for x in xs]
              for a in SPINS)

for a in SPINS
    println("  a=", a, "M  r_ISSO(x=-1)=", risso[a][1], "  r_ISSO(x=+1)=", risso[a][end],
            "  δr(x=-1)=", drisso[a][1], "  δr(x=+1)=", drisso[a][end])
end

styles = Dict(0.3 => (:red, :solid), 0.6 => (:blue, :dash), 0.9 => (:green, :dashdot))

pl = plot(xlabel = "x", ylabel = raw"$r_{\rm ISSOg}/M$", legend = :bottomleft,
          xticks = -1:0.5:1, yticks = 2:1:9, ylims = (2, 9))
for a in SPINS
    c, ls = styles[a]
    plot!(pl, xs, risso[a], color = c, linestyle = ls, lw = 2.5, label = "a=$(a)M")
end
savefig(pl, joinpath(@__DIR__, "rISSOgplot.png"))

pr = plot(xlabel = "x", ylabel = raw"$\delta r_{\rm ISSO}/(s_\parallel M)$",
          legend = :bottomright, xticks = -1:0.5:1)
for a in SPINS
    c, ls = styles[a]
    plot!(pr, xs, drisso[a], color = c, linestyle = ls, lw = 2.5, label = "a=$(a)M")
end
savefig(pr, joinpath(@__DIR__, "drISSOplot.png"))

@repro_result "x"              xs
for a in SPINS
    @repro_result "rISSO_a$(a)"  risso[a]
    @repro_result "drISSO_a$(a)" drisso[a]
end
# the published endpoints (legacy figure-script values, quoted in the ISSO promotion commit)
@repro_result "rISSO_a0.9_prograde"   kerr_isso(0.9, 1.0)
@repro_result "rISSO_a0.9_retrograde" kerr_isso(0.9, -1.0)
@repro_result "rISSO_a0.3_prograde"   kerr_isso(0.3, 1.0)
# both curves are monotone decreasing in x (prograde orbits sit deeper)
@repro_result "rISSO_monotone_in_x"   all(a -> all(diff(risso[a]) .< 0), SPINS)

exit(ReproHarness.conclude!(opts; paper = "2506.20726", id = "fig_rISSO", script = @__FILE__,
                            rtol = 1e-9))
