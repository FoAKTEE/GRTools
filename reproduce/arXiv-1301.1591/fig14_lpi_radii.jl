# arXiv:1301.1591 Fig. 14 analogue — local power-law indices (LPIs) of the ℓ = 2 projection of an
# s = −2, a = 0.9, ID0, ℓ′ = 4 perturbation at the paper's extraction radii (R_H = R₊, R1..R9,
# R_S = 1); the DG figure-set artefact `F-S1.14_LPI_radii_s-2_a0.9_lp4`.
#
# PAPER EXPECTATION [Tab. s2_LPI, row ℓ′ = 4]: μ₂ = 8|7 — LPI → −8 at finite radii (horizon
# included) and −7 at scri.
#
# MEASURED (legacy campaign, reproduced here): at Float64 on the committed stack the ℓ′ = 4 → ℓ = 2
# DOWN-MODE is buried in the junk floor before it reaches its asymptotic index; the envelope LPIs
# stall near −5.5. The paper's Fig. 14 required QUADRUPLE precision. This driver therefore
# reproduces a documented NO-GO, and the golden stores the measured (non-matching) LPIs — it is a
# faithful record of what the committed Float64 sector produces, never a fabricated match.
#
# tier HEAVY: one K = 32, tf = 2200 run with 12 probes, ≈ 7 min at -t4 (legacy wall 358 s).
#
#   julia --project=reproduce -t4 reproduce/arXiv-1301.1591/fig14_lpi_radii.jl --check --tier heavy

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

using GrTools
using GrTools.Perturbation: TimeDomain as TDS

const S_SPIN, A, M_AZIM, K, TF = -2, 0.9, 0, 32, 2200.0
const ELLS = collect(2:6)
const LP = 4
"Paper extraction radii R1..R9 [arXiv:1301.1591, Fig. 14]."
const REXT = [0.8229, 0.9114, 0.9544, 0.9740, 0.9883, 0.9919, 0.9971, 0.9987, 0.9997]
const FIT_WINDOW = (1200.0, 2000.0)

pol = TDS.spin_excise_policy(1.0, A)
p = TDS.SpinModesProblem(S_SPIN, ELLS, M_AZIM; N = 8, K = K, M = 1.0, a = A,
                         excise = pol.excise, filter_alpha = 36.0, extract_R = REXT)
@assert abs(p.mesh.edges[1] - (pol.Rplus - pol.excise)) < 1e-12
dt = TDS.spin_cfl_dt(p; C = 0.6)
rec = max(1, round(Int, 1.0 / dt))
println("[repro] fig14: dt=$(dt) probes=$(length(p.probe_R))"); flush(stdout)

t0 = time()
τs, Ψ, maxu = TDS.solve_spin_modes(p; tfinal = TF, dt = dt, record_every = rec, id = :ID0,
                                   init_mode = findfirst(==(LP), ELLS),
                                   ic = (amp = 1.0, R0 = 0.75, w = 1500.0))
println("[repro] wall=$(round(time()-t0, digits=1)) s  maxu_end=$(maxu[end])"); flush(stdout)

q2 = findfirst(==(2), ELLS)
names = vcat(["RH", "Rmid", "RS"], ["R$(i)" for i in 1:length(REXT)])

@repro_result "dt" dt
@repro_result "maxu_end" maxu[end]
@repro_result "maxu_finite" all(isfinite, maxu)
@repro_result "probe_R" collect(p.probe_R)
@repro_result "law_finite" -8      # paper Tab. s2_LPI, ℓ′=4 → ℓ=2, finite radii
@repro_result "law_scri" -7        # ... and at scri
for (ip, nm) in enumerate(names)
    ψ = abs.(Ψ[ip, q2, :])
    @repro_result "lpi_$(nm)_raw" TDS.fit_tail_window(τs, ψ, FIT_WINDOW)
    nenv, npk = TDS.fit_tail_envelope(τs, ψ, FIT_WINDOW)
    @repro_result "lpi_$(nm)_env" nenv
    @repro_result "lpi_$(nm)_npeaks" npk
end
idx = 1:10:length(τs)
@repro_result "tau" collect(τs[idx])
@repro_result "psi_l2_RH" abs.(Ψ[1, q2, idx])
@repro_result "psi_l2_RS" abs.(Ψ[3, q2, idx])

exit(ReproHarness.conclude!(opts; paper = "1301.1591", id = "fig14_lpi_radii",
                            script = @__FILE__, rtol = 1e-8))
