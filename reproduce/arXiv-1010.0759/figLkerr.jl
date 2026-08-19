# arXiv:1010.0759 fig:lkerr (r₀ = 50M) and fig:lkerr20 (r₀ = 20M) — Lindblad resonance LOCATIONS
# vs spin a for m = 1, 2, 3 (inner) and m = 3, 2, 1 (outer).
#
# THIN driver: `AccretionDisc.lindblad_locations`. Anchors at a = 0 are the tab:ilrS rows
# (r₀=50: m1 ILR 16.91, m2 ILR 33.44, m2 OLR 64.83, m1 OLR 78.33;
#  r₀=20: 10.14 / 14.50 / 25.47 / 30.65).
#
#   julia --project=reproduce -t4 reproduce/arXiv-1010.0759/figLkerr.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, lindblad_locations
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
using Plots

spins = collect(range(-0.9, 0.9; length = 37))
@repro_result "spin_grid" spins

function panel(r0)
    plt = plot(yscale = :log10, xlabel = "spin a", ylabel = "r_1 (resonance)",
               title = "Lindblad locations, r_0=$(Int(r0))M", legend = :right)
    for (lab, m, kind) in (("m=1 ILR", 1, :I), ("m=2 ILR", 2, :I), ("m=3 ILR", 3, :I),
                           ("m=3 OLR", 3, :O), ("m=2 OLR", 2, :O), ("m=1 OLR", 1, :O))
        ys = [begin
                  rI, rO = lindblad_locations(RelativisticDisc(a = a), r0, m)
                  kind === :I ? rI : rO
              end for a in spins]
        @repro_result "r0$(Int(r0))_m$(m)_$(kind)" ys
        plot!(plt, spins, ys; label = lab)
    end
    plt
end

p50 = panel(50.0)
p20 = panel(20.0)
emri_defaults!()
png = joinpath(@__DIR__, "expected", "figLkerr.png")
savefig(plot(p50, p20; layout = (1, 2), size = (1000, 420)), png)
@printf("wrote %s (%d bytes)\n", relpath(png, @__DIR__), filesize(png))

# a = 0 anchors vs tab:ilrS
for (r0, tab) in ((50.0, (16.91, 33.44, 64.83, 78.33)), (20.0, (10.14, 14.50, 25.47, 30.65)))
    d = RelativisticDisc(a = 0.0)
    i1, o1 = lindblad_locations(d, r0, 1); i2, o2 = lindblad_locations(d, r0, 2)
    @printf("r₀=%.0f a=0: m1ILR=%.2f(%.2f) m2ILR=%.2f(%.2f) m2OLR=%.2f(%.2f) m1OLR=%.2f(%.2f)\n",
            r0, i1, tab[1], i2, tab[2], o2, tab[3], o1, tab[4])
    @repro_result "anchor_a0_r0$(Int(r0))" [i1, i2, o2, o1]
end

exit(ReproHarness.conclude!(opts; paper = "1010.0759", id = "figLkerr",
                            script = @__FILE__, rtol = 1e-9))
