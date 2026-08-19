# arXiv:2506.20726 fig:trajectory_r_z_plane (paper.tex:501): a near-spherical trajectory projected
# into the r–z plane for s∥ = −M, −0.5M, 0, 0.5M, M at fixed turning points (a = 0.95M, p = 8M,
# x = 0.5). The spherical geodesic is r = p, z = z_max cos w_z with z_max = √(1−x²); the secondary
# spin bends it by s∥·(rS, zS) from the pex-gauge spherical orbit correction.
#
# THIN driver: (rS, zS) come from `Perturbation.KerrSpinOrbitCorrectionSpherical(...).OrbitCorrection`.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2506.20726/fig_trajectory_rz.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Perturbation: KerrSpinOrbitCorrectionSpherical
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const A, P, X = 0.95, 8.0, 0.5
const ZMAX = sqrt(1 - X^2)

oc  = KerrSpinOrbitCorrectionSpherical(A, P, 0.0, X, 8)
wzs = collect(range(0, 2pi, length = 400))
RS  = [oc.OrbitCorrection(wz).rS for wz in wzs]
ZS  = [oc.OrbitCorrection(wz).zS for wz in wzs]
ZG  = ZMAX .* cos.(wzs)

const SPINS = [(-1.0, :blue, :solid, "s∥ = -M"), (-0.5, :orange, :dash, "s∥ = -0.5M"),
               (0.0, :green, :dot, "s∥ = 0"), (0.5, :red, :dashdot, "s∥ = 0.5M"),
               (1.0, :purple, :dash, "s∥ = M")]

plt = plot(xlabel = "r/M", ylabel = "z", legend = :outerright, xticks = [7.99, 8.0, 8.01])
for (s, c, ls, lab) in SPINS
    r = P .+ s .* RS
    z = ZG .+ s .* ZS
    plot!(plt, r, z, color = c, linestyle = ls, lw = 2, label = lab)
    println("  ", lab, "  r ∈ [", minimum(r), ", ", maximum(r), "]  max|r-8| = ",
            maximum(abs.(r .- P)))
end
savefig(plt, joinpath(@__DIR__, "trajectory_r_z_plane.png"))

@repro_result "w_z"       wzs
@repro_result "rS"        RS
@repro_result "zS"        ZS
@repro_result "z_geodesic" ZG
@repro_result "max_abs_rS" maximum(abs.(RS))
@repro_result "max_abs_zS" maximum(abs.(ZS))
# fixed turning points: the polar correction vanishes at the turning points w_z = 0, π, 2π
@repro_result "zS_at_turning_points" [ZS[1], ZS[200], ZS[end]]

exit(ReproHarness.conclude!(opts; paper = "2506.20726", id = "fig_trajectory_rz",
                            script = @__FILE__, rtol = 1e-9, atol = 1e-14))
