# arXiv:1301.1591 — Campaign-C self-convergence panels of the s = −2 RT-chart sector (the DG
# figure-set analogue `F-S1.1_Convergence_s-2_RT`, itself the analogue of the DG repo's
# F1.1 super-convergence figure).
#
# WHAT IS MEASURED: the rel-L2-in-time error of the ℓ = 2 probe series against the finest-K run at
# the same polynomial order N, on a SHARED fixed dt = 1e-3 and sampling cadence (so the τ grids are
# identical and the comparison is grid-to-grid), at the horizon / mid-domain / scri probes. The
# fitted order is the negated log-log slope of E(K) (`fit_convergence_order`). Expected classes:
# N+1 in the interior, 2N+1 at the compactified scri boundary.
#
# CONFIG (legacy `run_convergence.jl`): a = 0.9, s = −2, m = 0, ℓ = 2..4 coupled, ID0 seed (ℓ=2)
# with R0 = 0.75, w = 600, tfinal = 50 (the pre-junk epoch: at tf = 100 the κ-oscillation junk
# floors the ladders), constraint_kappa = 0.3, and the modal filter OFF — filter_alpha = 0 is
# required here because the |s| = 2 filter is a long-run stabiliser that would contaminate the
# spectral order being measured; over tf = 50 the documented noise growth (≤ 0.034/M) amplifies
# round-off by ≤ e^1.7, far below the smallest measured error.
#
# tier HEAVY: three ladders + three reference runs, ≈ 8 min at -t4.
#
#   julia --project=reproduce -t4 reproduce/arXiv-1301.1591/fig_convergence.jl --check --tier heavy

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

using GrTools
using GrTools.Perturbation: TimeDomain as TDS
using Printf
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))

const S_SPIN, A, M_AZIM, ELLS = -2, 0.9, 0, [2, 3, 4]
const TF, DT, REC = 50.0, 1.0e-3, 500
const KAPPA = 0.3
const IC = (amp = 1.0, R0 = 0.75, w = 600.0)
const LADDERS = ((2, [16, 24, 32, 48, 64, 96, 128], 256),
                 (4, [8, 12, 16, 24, 32, 48, 64], 128),
                 (6, [6, 8, 12, 16, 24, 32], 96))
const PROBES = ("hor", "mid", "scri")

"One ℓ = 2 probe series at order N, K elements — shared dt/cadence ⇒ identical τ grids."
function probe_series(N, K)
    pol = TDS.spin_excise_policy(1.0, A)
    p = TDS.SpinModesProblem(S_SPIN, ELLS, M_AZIM; N = N, K = K, M = 1.0, a = A,
                             excise = pol.excise, filter_alpha = 0.0, constraint_kappa = KAPPA)
    @assert abs(p.mesh.edges[1] - (pol.Rplus - pol.excise)) < 1e-12
    cfl = TDS.spin_cfl_dt(p; C = 0.6)
    @assert DT < cfl "shared dt=$(DT) above CFL $(cfl) at K=$(K) N=$(N)"
    τs, Ψ, mx = TDS.solve_spin_modes(p; tfinal = TF, dt = DT, record_every = REC, id = :ID0,
                                     init_mode = 1, ic = IC)
    @assert all(isfinite, mx)
    return τs, Ψ[:, 1, :]
end

relL2(u, v) = sqrt(sum(abs2, u .- v) / sum(abs2, v))

errs = Dict{Int,Matrix{Float64}}()      # N => (K × 3) errors
for (N, Ks, Kref) in LADDERS
    t0 = time()
    τr, Ψr = probe_series(N, Kref)
    E = zeros(length(Ks), 3)
    for (i, K) in enumerate(Ks)
        τs, Ψ = probe_series(N, K)
        @assert length(τs) == length(τr)
        for ip in 1:3
            E[i, ip] = relL2(Ψ[ip, :], Ψr[ip, :])
        end
        println("[repro] N=$(N) K=$(K): ", join((@sprintf("%s=%.3e", PROBES[ip], E[i, ip]) for ip in 1:3), " "))
        flush(stdout)
    end
    errs[N] = E
    @repro_result "N$(N)_K" Float64.(Ks)
    for ip in 1:3
        @repro_result "N$(N)_E$(PROBES[ip])" E[:, ip]
        @repro_result "N$(N)_rate_$(PROBES[ip])" TDS.fit_convergence_order(Ks, E[:, ip])
        @repro_result "N$(N)_class_$(PROBES[ip])" (PROBES[ip] == "scri" ? 2N + 1 : N + 1)
    end
    println("[repro] N=$(N) ladder wall=$(round(time()-t0, digits=1)) s"); flush(stdout)
end

# ---- figure: 1×3 log-log self-convergence panels with per-N slope legends ------------------------
# Presentation follows the legacy render `results/dg-teukolsky-spin/figures/F-S1.1_Convergence_s-2_RT.png`
# (VLM-compared: same panel structure, ladders, rates and ordering; colours/markers/legend text matched).
if opts.mode === :regen
    emri_defaults!()
    styles = Dict(2 => (:red, :star5), 4 => (:blue, :circle), 6 => (:magenta, :xcross))
    titles = ("R+ (horizon)", "R_mid (interior)", "R = 1 = scri+")
    plt = plot(layout = (1, 3), size = (1000, 340), legend = :bottomleft, legendfontsize = 6,
               left_margin = 6Plots.mm, bottom_margin = 6Plots.mm, top_margin = 3Plots.mm,
               grid = true, gridalpha = 0.25, gridstyle = :dot)
    for ip in 1:3, (N, Ks, _) in LADDERS
        E = errs[N][:, ip]
        r = TDS.fit_convergence_order(Ks, E)
        cls = ip == 3 ? 2N + 1 : N + 1
        c, mk = styles[N]
        plot!(plt[ip], Ks, E, xscale = :log10, yscale = :log10, c = c, m = mk, ms = 4, lw = 1.2,
              label = @sprintf("N=%d (p=%.2f, class %d)", N, r, cls))
        plot!(plt[ip], title = titles[ip], xlabel = "number of subdomains K",
              ylabel = ip == 1 ? "rel. L2-in-time error, l = 2 probe series" : "",
              xticks = ([10, 100], ["10^1", "10^2"]))
    end
    plot!(plt, plot_title = "Convergence, s = -2 RT-chart vacuum sector (a = 0.9, coupled l = 2..4)",
          plot_titlefontsize = 11)
    png = joinpath(@__DIR__, "convergence.png")
    savefig(plt, png)
    println("[repro] wrote ", png, " (", filesize(png), " bytes)")
end

exit(ReproHarness.conclude!(opts; paper = "1301.1591", id = "convergence", script = @__FILE__,
                            rtol = 1e-8))
