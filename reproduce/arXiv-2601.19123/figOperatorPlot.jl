# arXiv:2601.19123 fig:OperatorPlot (F2) — "Radial Profile of Potential Function": the effective
# potential B_m(r) of the reduced enthalpy-wave equation (eq:reduceeq) at the fiducial parameters
# r_p = 15M, h = 0.1, a = 0.6, for m = 1…15.
#
# THIN driver: `AccretionDisc.B_m_geometric` — the LEADING GEOMETRIC potential B^{(−2)}_m/h², whose
# zeros are exactly the Lindblad turning points (`lindblad_locations`) and whose 1/(r−r_p)²
# corotation pole coefficient is `B_m_corotation`.
#
# SCOPE / documented gap: the published figure also shows the warm +∞ pressure walls bracketing r_p;
# those need the unpublished height-averaged B^{(0)} term ([HOLE: H-DISC-FULLOP] in the package).
# The blue excitation family, the −∞ corotation pole and the Lindblad turning points reproduce here.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2601.19123/figOperatorPlot.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, B_m_geometric, B_m_corotation, lindblad_locations, risco_disc
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
using Plots

const A = 0.6; const RP = 15.0; const H = 0.1; const MMAX = 15
d = RelativisticDisc(a = A, h0 = H)
rs = range(risco_disc(d) + 0.05, 40.0; length = 800)
@repro_result "r_grid" collect(rs)

emri_defaults!()
plt = plot(; xlabel = "\$r/M\$", ylabel = "\$B_m(r)\$", title = "Radial profile of the potential (a=0.6, r_p=15M, h=0.1)",
           legend = false, ylims = (-2e4, 2e4), size = (760, 480))
rILRs = Float64[]; rOLRs = Float64[]
for m in 1:MMAX
    B = [abs(r - RP) < 1e-6 ? NaN : B_m_geometric(d, r, RP, m) for r in rs]
    m in (1, 3, 8, 15) && @repro_result "B_m$(m)" B
    rI, rO = lindblad_locations(d, RP, m)
    push!(rILRs, rI); push!(rOLRs, rO)
    plot!(plt, rs, B; lw = 1.2, color = RGB(0.2, 0.35 + 0.04*m, 0.9 - 0.03*m))
end
@repro_result "rILR_m1to15" rILRs
@repro_result "rOLR_m1to15" rOLRs
vline!(plt, [RP]; lc = :black, ls = :dash, lw = 1.5)
png = joinpath(@__DIR__, "expected", "figOperatorPlot.png")
savefig(plt, png)
@printf("wrote %s (%d bytes)\n", relpath(png, @__DIR__), filesize(png))

# structural assertions registered as numbers
@repro_result "B_at_ILR_m3" B_m_geometric(d, rILRs[3], RP, 3)
@repro_result "B_at_OLR_m3" B_m_geometric(d, rOLRs[3], RP, 3)
pole = [1e-6^2 * B_m_geometric(d, RP + 1e-6, RP, m) * H^2 for m in 1:MMAX]
poleref = [1e-6^2 * B_m_corotation(d, RP + 1e-6, RP, m) for m in 1:MMAX]
@repro_result "corotation_pole_coefficient" pole
@repro_result "corotation_pole_max_reldiff" maximum(abs.(pole .- poleref) ./ abs.(poleref))
@printf("Lindblad zeros: |B(r_ILR,m=3)|=%.3e, |B(r_OLR,m=3)|=%.3e ; corotation pole matches B_m_corotation to %.2e\n",
        abs(B_m_geometric(d, rILRs[3], RP, 3)), abs(B_m_geometric(d, rOLRs[3], RP, 3)),
        maximum(abs.(pole .- poleref) ./ abs.(poleref)))

exit(ReproHarness.conclude!(opts; paper = "2601.19123", id = "figOperatorPlot",
                            script = @__FILE__, rtol = 1e-8,
                            tolerances = Dict("corotation_pole_max_reldiff" => (rtol = 0.5, atol = 1e-12),
                                              "B_at_ILR_m3" => (rtol = 1e-6, atol = 1e-6),
                                              "B_at_OLR_m3" => (rtol = 1e-6, atol = 1e-6))))
