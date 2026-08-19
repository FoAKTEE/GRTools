# arXiv:1010.0759 — the THREE convention reconciliations between the GrTools Teukolsky engine and
# Hirata's relativistic Lindblad amplitudes, and the full relativistic reconstruction of c₁₃:
#   (1) eq:cplex1 closure  |c₁₃|² + α = (2ω)⁶ |ℵ|²/|C|²  with the Teukolsky–Press ε = Γ/2 horizon
#       factor.  THIS IS THE GATE that fixes α: Hirata's PRINTED bracket (ϖ²+Γ²)(ϖ²+4Γ²) misses
#       the closure by up to 9.1e-5 (ERRATA #8); the package uses (ϖ²+(Γ/2)²)(ϖ²+Γ²).
#   (2) the Newtonian limit c₁₃ ∝ (2ω)^{1−ℓ} and eq:cplex1 with α → 0 (Schwarzschild);
#   (3) the Teukolsky–Starobinsky |C|² of the paper (l.655) ≡ the package flux |C|²;
#   (4) the SWSH sphere normalization ∫|S_pkg|²dΩ = 1 ⇒ Z^out = Z^∞_pkg/√(2π);
#   (5) end-to-end: the package energy fluxes ≡ Hirata's P_out / P_down with those conversions.
#
# THIN driver: `AccretionDisc.Modes.c13_alpha` returns (c13, c14, aleph, alpha, lambda, Csq);
# `Perturbation` supplies the point-particle amplitudes and fluxes, `Harmonics` the SWSH.
#
#   julia --project=reproduce -t4 reproduce/arXiv-1010.0759/check_c13_alpha.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: c13_alpha
using GrTools.Perturbation: TeukolskyPointParticleModeEccentric
using GrTools.Harmonics: SpinWeightedSpheroidalHarmonicS
using GrTools.NumericsCore: gl_nodes_weights

# (1) eq:cplex1 closure at six (ℓ, m, a, ω)
println("(1) eq:cplex1  |c13|² + α = (2ω)⁶|ℵ|²/|C|²")
worst = 0.0
for (l, m, a, w) in ((2, 2, 0.5, 0.12), (2, 1, 0.5, 0.08), (3, 2, 0.7, 0.15),
                     (2, 2, 0.0, 0.10), (4, 3, 0.9, 0.20), (3, 3, 0.3, 0.18))
    cc = c13_alpha(l, m, a, w; Nm = 300)
    rhs = (2w)^6/cc.Csq*abs2(cc.aleph)
    rel = abs(abs2(cc.c13) + cc.alpha - rhs)/abs(rhs)
    global worst = max(worst, rel)
    @printf("  (l,m,a,ω)=(%d,%d,%.1f,%.2f)  |c13|=%.7g  α=%.6g  relerr=%.3e\n",
            l, m, a, w, abs(cc.c13), cc.alpha, rel)
end
@repro_result "cplex1_worst_relerr" worst
@repro_result "cplex1_closes" worst < 1e-11

# a representative amplitude bundle (regression on the numbers themselves)
cc = c13_alpha(2, 2, 0.5, 0.12; Nm = 300)
@repro_result "c13_l2m2_a05_w012" cc.c13
@repro_result "alpha_l2m2_a05_w012" cc.alpha
@repro_result "Csq_l2m2_a05_w012" cc.Csq
@repro_result "lambda_l2m2_a05_w012" cc.lambda

# (2) Newtonian limit: |c13(ω₁)|/|c13(ω₂)| → (ω₁/ω₂)^{1−ℓ}, and α → 0 closure at a = 0
println("\n(2) Newtonian limit  c13 ∝ (2ω)^{1−ℓ}  and  |c13|² = (2ω)⁶|ℵ|²/|C|²  (α→0, a=0)")
pl = Float64[]; nl = Float64[]
for l in 2:6
    w1, w2 = 0.004, 0.002
    r = abs(c13_alpha(l, 2, 0.0, w1; Nm = 300).c13)/abs(c13_alpha(l, 2, 0.0, w2; Nm = 300).c13)
    push!(pl, abs(r - (w1/w2)^(1-l))/(w1/w2)^(1-l))
    @printf("  ℓ=%d  ratio=%.8g  (ω1/ω2)^{1−ℓ}=%.8g  reldiff=%.3e\n", l, r, (w1/w2)^(1-l), pl[end])
end
for l in 2:5
    w = 0.003; c = c13_alpha(l, 2, 0.0, w; Nm = 300)
    Csq_schw = c.lambda^2*(c.lambda + 2)^2 + 144*w^2
    push!(nl, abs(abs2(c.c13) - (2w)^6/Csq_schw*abs2(c.aleph))/((2w)^6/Csq_schw*abs2(c.aleph)))
    @printf("  ℓ=%d  Newtonian eq:cplex1 relerr = %.3e   (Schwarzschild |C|² = %.10g)\n", l, nl[end], Csq_schw)
end
@repro_result "newtonian_c13_powerlaw_max_reldiff" maximum(pl)
@repro_result "newtonian_cplex1_max_relerr" maximum(nl)

# (3) the paper's |C|² (l.655) vs the package flux |C|²
println("\n(3) Teukolsky–Starobinsky |C|²: paper form ≡ package flux form")
tsmax = 0.0
for (l, m, a, w) in ((2, 2, 0.5, 0.12), (3, 2, 0.7, 0.15), (4, 3, 0.9, 0.20), (2, 2, 0.0, 0.10))
    c = c13_alpha(l, m, a, w; Nm = 300); lam = c.lambda
    Ck = ((lam + 2)^2 + 4a*m*w - 4a^2*w^2)*(lam^2 + 36*m*a*w - 36*a^2*w^2) +
         (2lam + 3)*(96*a^2*w^2 - 48*m*a*w) + 144*w^2*(1 - a^2)
    global tsmax = max(tsmax, abs(c.Csq - Ck)/abs(Ck))
    @printf("  (l,m,a,ω)=(%d,%d,%.1f,%.2f)  |C|²=%.10g  reldiff=%.3e\n", l, m, a, w, c.Csq, abs(c.Csq - Ck)/abs(Ck))
end
@repro_result "TS_Csq_max_reldiff" tsmax

# (4) SWSH normalization ∫|S|² dΩ = 1 (Gauss–Legendre in cos θ)
S = SpinWeightedSpheroidalHarmonicS(-2, 2, 2, 0.5)
x, wq = gl_nodes_weights(Float64, 200)
norm2 = 2π*sum(wq[i]*abs2(S(acos(x[i]), 0.0)) for i in eachindex(x))
@repro_result "swsh_sphere_norm" norm2
@printf("\n(4) ∫|S_pkg|² dΩ = %.12f  ⇒  Z^out = Z^∞_pkg/√(2π)\n", norm2)

# (5) end-to-end: package fluxes ≡ Hirata P_out / P_down under the pinned conversions
println("\n(5) P_out = |Z^out|²/(2ω²) ≡ EnergyFlux.I ;  P_down = α|Z^down|²/(2ω²) ≡ EnergyFlux.H")
eI = Float64[]; eH = Float64[]
for (l, m, n, a, p, e) in ((2, 2, 1, 0.5, 10.0, 0.1), (3, 3, 1, 0.6, 12.0, 0.15), (2, 2, -1, 0.7, 9.0, 0.2))
    md = TeukolskyPointParticleModeEccentric(l, m, n, a, p, e; s = -2)
    w = md.omega
    al = c13_alpha(l, m, a, w; Nm = 300).alpha
    Pout = abs2(md.Zinf/sqrt(2π))/(2w^2); Pdown = al*abs2(md.ZH/sqrt(2π))/(2w^2)
    push!(eI, abs(Pout - md.EnergyFlux.I)/abs(md.EnergyFlux.I))
    push!(eH, abs(Pdown - md.EnergyFlux.H)/abs(md.EnergyFlux.H))
    @printf("  (l,m,n,a,p,e)=(%d,%d,%+d,%.1f,%.1f,%.2f) ω=%.5g  P_out rel=%.3e  P_down rel=%.3e\n",
            l, m, n, a, p, e, w, eI[end], eH[end])
end
@repro_result "flux_identity_max_relerr_I" maximum(eI)
@repro_result "flux_identity_max_relerr_H" maximum(eH)

exit(ReproHarness.conclude!(opts; paper = "1010.0759", id = "check_c13_alpha",
                            script = @__FILE__, rtol = 1e-8,
                            tolerances = Dict("cplex1_worst_relerr" => (rtol = 0.5,),
                                              "newtonian_cplex1_max_relerr" => (rtol = 0.5,),
                                              "TS_Csq_max_reldiff" => (rtol = 1.0, atol = 1e-13),
                                              "flux_identity_max_relerr_I" => (rtol = 0.5,),
                                              "flux_identity_max_relerr_H" => (rtol = 0.5,))))
