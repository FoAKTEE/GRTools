# arXiv:1010.0759 (Hirata II) — RT-N1: the PUBLISHED Newtonian-limit normalized torque strengths
#   N = −2.36 (m=2 ILR), −7.50 (m=3 ILR), +0.19 (m=1 OLR), +8.62 (m=2 OLR)   (l.1008 / eq:smolr)
# together with the RT-A2 Laplace-coefficient identities the closed forms rest on.
#
# THIN driver: AccretionDisc.Modes supplies every number — `lindblad_strength_newtonian` (the
# r₀-independent closed form), `lindblad_strength` (the same kernel at finite r₀, exact Kerr
# kinematics), `laplace_b` / `laplace_b_deriv`.
#
# TWO PUBLISHED TYPOS are located by this check and are carried in the package docstrings:
#   (i)  l.1008 prints N = ∓(π/6) ς² |·|² — the factor m/(m∓1) is missing (it would give
#        −1.18 / −5.00 instead of the tabulated −2.36 / −7.50);
#   (ii) eq:smolr prints b′(ς⁻) in the OLR bracket; the correct argument is ς⁺ (m=2 OLR is 8.62,
#        not 5e-4).
#
#   julia --project=reproduce -t4 reproduce/arXiv-1010.0759/check_RT_N1.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, lindblad_strength, lindblad_strength_newtonian,
                             laplace_b, laplace_b_deriv

const PUBLISHED = ((2, :ILR, -2.36), (3, :ILR, -7.50), (1, :OLR, 0.19), (2, :OLR, 8.62))

println("kind m   ς        bracket-free N (closed form)   published")
for (m, kind, pub) in PUBLISHED
    N = lindblad_strength_newtonian(m; kind = kind)
    @repro_result "N_$(kind)_m$(m)" N
    @printf("%-4s %d   %+10.6f   %+8.4f   (rel %.2e)\n", kind, m, N, pub, abs(N - pub)/abs(pub))
end

# r₀-independence: the closed form is the r₀ → ∞ limit of the finite-r₀ kernel (a = 0)
d0 = RelativisticDisc(a = 0.0)
for (m, kind, _) in PUBLISHED
    Nlim = lindblad_strength_newtonian(m; kind = kind)
    for r0 in (1.0e4, 1.0e6)
        Nr = lindblad_strength(d0, r0, m; kind = kind)
        @printf("  m=%d %s  r0=%-8g N=%+.8f   |N/N_lim − 1| = %.2e\n", m, kind, r0, Nr, abs(Nr/Nlim - 1))
    end
    @repro_result "N_$(kind)_m$(m)_r0_1e6" lindblad_strength(d0, 1.0e6, m; kind = kind)
end

# RT-A2 identities: the reciprocal relation b^{(m)}(1/ς) = ς b^{(m)}(ς) and the m=1 indirect
# coefficient b^{(1)′}(0) = 1 (the gauge equivalence −ς⁺ ≡ −4(ς⁺)^{-2} at ς⁺₁ = 2^{2/3}).
recip = [abs(laplace_b(0.5, m, 1/s) - s*laplace_b(0.5, m, s)) for (m, s) in ((1, 0.63), (2, 0.6), (3, 0.5))]
@repro_result "laplace_reciprocal_maxabs" maximum(recip)
@repro_result "b1_prime_at_0" laplace_b_deriv(0.5, 1, 1.0e-6)
s1 = 2.0^(2/3)
@repro_result "m1_OLR_indirect_gauge_gap" abs(-s1 - (-4*s1^(-2)))
@printf("\nLaplace reciprocal max|b(1/ς)−ς b(ς)| = %.2e ; b^{(1)′}(0) = %.10f ; |−ς⁺ + 4(ς⁺)^{-2}| = %.2e\n",
        maximum(recip), laplace_b_deriv(0.5, 1, 1.0e-6), abs(-s1 + 4*s1^(-2)))

exit(ReproHarness.conclude!(opts; paper = "1010.0759", id = "check_RT_N1",
                            script = @__FILE__, rtol = 1e-10))
