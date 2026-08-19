# arXiv:2603.12189 fig:comparepvalues_leadingharmonics (fig. 7, main.tex:889): the error of three
# APPROXIMATE spinning trajectories against the exact one, |δF_exact − δF_approx|, for the single
# (l,m,n,k) = (2,2,0,0) mode versus p ∈ [6, 15] at a = 0.9 M, e = 0.2, x = √3/2, in four channels
# (E at infinity / at the horizon, Lz at infinity / at the horizon).
#
# The three approximations, all built from `Perturbation.KerrSpinOrbitCorrection`:
#   SG          — the shifted geodesic: the exact nmax = kmax = 8 orbit correction with the
#                 LIBRATION zeroed (rS = zS = UrS = UzS = 0), i.e. the secular shift only
#   n^s = 4     — the orbit correction truncated at nmax = kmax = 4
#   n^s = 6     — the orbit correction truncated at nmax = kmax = 6
# The published message is the ordering SG > n^s = 4 > n^s = 6 at every p, i.e. the spinning
# trajectory converges harmonic by harmonic and the shifted geodesic is the crudest member.
#
# THIN driver: the flux of each trajectory is `Perturbation._spin_mode_correction(oc, 2,2,0,0)`;
# the only thing built here is the zeroed-libration wrapper around the exact correction closure,
# which is a CHOICE OF INPUT to the package function, not physics.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2603.12189/fig07.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

using GrTools
using GrTools.Perturbation: KerrSpinOrbitCorrection, _spin_mode_correction
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const A, E, X = 0.9, 0.2, sqrt(3) / 2
const PS = collect(6.0:0.5:15.0)
const APPROX = (:SG, :h4, :h6)
const CHAN = (:EI, :EH, :LI, :LH)

zero_libration(oc) = merge(oc, (; OrbitCorrection = (wr, wz) -> begin
    c = oc.OrbitCorrection(wr, wz)
    (ΔtS = c.ΔtS, rS = 0.0, zS = 0.0, ΔφS = c.ΔφS, UrS = 0.0, UzS = 0.0)
end))

sflux(oc) = let md = _spin_mode_correction(oc, 2, 2, 0, 0)
    (EI = real(md.dE1I), EH = real(md.dE1H), LI = real(md.dL1I), LH = real(md.dL1H))
end

D = Dict(v => Dict(c => fill(NaN, length(PS)) for c in CHAN) for v in APPROX)
Threads.@threads for i in eachindex(PS)
    p = PS[i]
    try
        oc8 = KerrSpinOrbitCorrection(A, p, E, X, 8, 8)
        fe = sflux(oc8)
        for (v, o) in ((:SG, zero_libration(oc8)),
                       (:h4, KerrSpinOrbitCorrection(A, p, E, X, 4, 4)),
                       (:h6, KerrSpinOrbitCorrection(A, p, E, X, 6, 6)))
            fa = sflux(o)
            for c in CHAN
                D[v][c][i] = abs(getfield(fe, c) - getfield(fa, c))
            end
        end
    catch err
        println("  p=", p, " FAILED: ", typeof(err)); flush(stdout)
    end
end

i10 = argmin(abs.(PS .- 10))
for c in CHAN, v in APPROX
    println("  ", rpad(string(c), 3), " ", rpad(string(v), 3), "  @p=6 ", D[v][c][1],
            "  @p=10 ", D[v][c][i10], "  @p=15 ", D[v][c][end])
end
flush(stdout)

sty = Dict(:SG => (MMA_RED, :dash, :circle), :h4 => (:orange, :dashdot, :utriangle),
           :h6 => (MMA_BLUE, :dot, :cross))
lab = Dict(:SG => "Shifted geodesic", :h4 => raw"$n^s_{\max}=k^s_{\max}=4$",
           :h6 => raw"$n^s_{\max}=k^s_{\max}=6$")
ylab = Dict(:EI => raw"$|\delta\mathcal{F}^{E,\infty}_{2200}$ exact $-$ approx$|$",
            :EH => raw"$|\delta\mathcal{F}^{E,H}_{2200}$ exact $-$ approx$|$",
            :LI => raw"$|\delta\mathcal{F}^{L_z,\infty}_{2200}$ exact $-$ approx$|$",
            :LH => raw"$|\delta\mathcal{F}^{L_z,H}_{2200}$ exact $-$ approx$|$")
panels = Any[]
for c in CHAN
    pl = plot(xlabel = "p", ylabel = ylab[c], yscale = :log10,
              legend = (c === :EI ? :bottomleft : false), legendfontsize = 6, ytickfontsize = 6)
    for v in APPROX
        col, ls, mk = sty[v]
        plot!(pl, PS, D[v][c], color = col, linestyle = ls, marker = mk, ms = 3, lw = 1.5, label = lab[v])
    end
    push!(panels, pl)
end
savefig(plot(panels..., layout = (2, 2), size = (1050, 740)), joinpath(@__DIR__, "fig07.png"))

@repro_result "p" PS
for v in APPROX, c in CHAN
    @repro_result "err_$(v)_$(c)" D[v][c]
end
# the published message, per point: SG > n^s=4 > n^s=6
ok = [count(i -> D[:SG][c][i] > D[:h4][c][i] > D[:h6][c][i], eachindex(PS)) for c in CHAN]
@repro_result "n_points"                 length(PS)
@repro_result "ordering_SG_gt_h4_gt_h6"  ok
@repro_result "ordering_holds_everywhere" all(ok .== length(PS))

exit(ReproHarness.conclude!(opts; paper = "2603.12189", id = "fig07", script = @__FILE__,
                            rtol = 1e-8))
