# arXiv:2303.16798 fig:fourier_coefficients (paper.tex:532): the radial spin-correction Fourier
# coefficients δχ^S_{r,n} vs n for a = 0.9 M, p = 15, I = 15°, e ∈ {0.2, 0.4, 0.6} (nmax = 16,
# kmax = 8). The published message is the least-squares truncation FLOOR: past a certain |n| the
# coefficients stop decreasing and the floor RISES with eccentricity.
#
# THIN driver: the coefficients are `Spacetime.SpinningSecondary._numeric_correction(...).δχrSIm`
# (the 2-torus DFT least-squares engine); the script only takes |·|, plots and registers.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2303.16798/fig_fourier_coefficients.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Spacetime.SpinningSecondary: _numeric_correction
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const A, P, NMAX, KMAX = 0.9, 15.0, 16, 8
const X = cos(15 * pi / 180)
const ES = (0.2, 0.4, 0.6)
const COLS = (MMA_BLUE, MMA_GOLD, MMA_GREEN)

ns = vcat(-NMAX:-1, 1:NMAX)          # δχrSIm index order (n = 0 excluded)

plt = plot(xlabel = "n", ylabel = raw"$\delta\chi^S_{r,n}$", yscale = :log10,
           legend = :right, xticks = -15:5:15)

coeffs = Dict{Float64,Vector{Float64}}()
for (e, c) in zip(ES, COLS)
    y = abs.(_numeric_correction(A, P, e, X, NMAX, KMAX).δχrSIm)
    coeffs[e] = y
    scatter!(plt, ns, y, marker = (:circle, 5), color = c, label = "e = $(e)")
    ipk = findfirst(==(1), ns)
    println("  e=", e, "  peak |δχ_{r,+1}| = ", y[ipk],
            "   floor |δχ_{r,+16}| = ", y[end], "   |δχ_{r,-16}| = ", y[1])
    flush(stdout)
end
ylims!(plt, 1e-22, 3e-2)
savefig(plt, joinpath(@__DIR__, "fourier_coefficients.png"))

ipk = findfirst(==(1), ns)
peaks  = [coeffs[e][ipk] for e in ES]
floors = [max(coeffs[e][1], coeffs[e][end]) for e in ES]

@repro_result "n"        ns
for e in ES
    @repro_result "abs_dchi_rS_e$(e)" coeffs[e]
end
@repro_result "peak_at_n_plus1"  peaks
@repro_result "floor_at_n_16"    floors
# The published message is "the truncation floor RISES with e". MEASURED HERE (and identically in
# the legacy figure, see fourier_coefficients.spec.md Rule 3): only the e = 0.6 floor is resolvable
# in Float64 — the e = 0.2 and e = 0.4 floors sit at the LEAST-SQUARES CONDITIONING floor of the
# Float64 solve, MEASURED at 1.7e4 … 3.3e6 times ε·peak (registered below), i.e. ~1e-11 relative,
# while the published figure reaches ~1e-20 (it was solved in EXTENDED precision). The ordering
# between the two low-e floors is therefore solver noise, not physics: registered as measured,
# never as a monotone claim.
@repro_result "floor_over_eps_times_peak" floors ./ (eps(Float64) .* peaks)
@repro_result "floor_e06_is_highest"      floors[3] == maximum(floors)
@repro_result "floor_below_peak_by_1e6"   all(floors ./ peaks .< 1e-6)
@repro_result "peaks_grow_linearly_in_e"  maximum(abs.(peaks ./ collect(ES) ./ (peaks[1] / ES[1]) .- 1))

exit(ReproHarness.conclude!(opts; paper = "2303.16798", id = "fig_fourier_coefficients",
                            script = @__FILE__, rtol = 1e-8))
