# arXiv:2510.02433 — the three closed-form scaling statements of the torque discussion:
#   (C1) Newtonian inner-disc torque scaling  T ∝ r^{n_r},  n_r = 3/2 − 2Σ_p  for Σ ∝ r^{−Σ_p},
#        h ∝ r^{(2Σ_p−1)/4}  ⇒  radiation-pressure inner disc (Σ_p = −3/2) gives n_r = 9/2 (l.605).
#        Verified BOTH as exact rational arithmetic and against the package NT/SS profiles.
#   (C2) the pressure-modified Lindblad resonance SATURATES: |r_LR − r_s| → 2H/3 as m → ∞ (l.514,
#        Artymowicz 1993) — `AccretionDisc.lindblad_locations_pressure`.
#   (C3) the suppression onset ξ = m h ≳ 1 ⟺ m ≳ 1/h, and the cutoff `pressure_cutoff`.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2510.02433/check_satscale.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, NovikovThorneDisc, lindblad_locations_pressure,
                             pressure_cutoff, ss_surface_density, ss_aspect_ratio, powerlaw_fit

# (C1) exponent arithmetic  n_r = 3/2 − 2Σ_p
println("(C1) T ∝ r^{n_r}, n_r = 3/2 − 2Σ_p")
nr(Sp) = 3//2 - 2*Sp
for Sp in (-3//2, -1//1, 0//1, 1//2)
    @printf("  Σ_p=%-6s Σ∝r^{%-6s} h∝r^{%-6s} ⇒ n_r = %s\n", string(Sp), string(-Sp), string((2Sp-1)//4), string(nr(Sp)))
end
@repro_result "n_r_radiation_pressure" Float64(nr(-3//2))
@repro_result "n_r_equals_9_over_2" nr(-3//2) == 9//2
# the SS profile in the package IS the Σ_p = −3/2 case: log-slopes 3/2 and −1
nt0 = NovikovThorneDisc(a = 0.0)
sl(f, r) = (log(f(r*1.001)) - log(f(r/1.001)))/(log(r*1.001) - log(r/1.001))
@repro_result "SS_Sigma_logslope" sl(r -> ss_surface_density(nt0, r), 100.0)
@repro_result "SS_h_logslope" sl(r -> ss_aspect_ratio(nt0, r), 100.0)
@printf("  package SS profile: dlnΣ/dlnr = %.12f (3/2), dlnh/dlnr = %.12f (−1)  ⇒ Σ_p = %.12f\n",
        sl(r -> ss_surface_density(nt0, r), 100.0), sl(r -> ss_aspect_ratio(nt0, r), 100.0),
        -sl(r -> ss_surface_density(nt0, r), 100.0))

# (C2) saturation of the pressure-modified resonance at 2H/3
println("\n(C2) pressure-modified resonance saturation (h = 0.01, r_s = 50M): |r_LR − r_s| → 2H/3")
h0 = 0.01; rp = 50.0; H = h0*rp
dN = RelativisticDisc(a = 0.0, h0 = h0)
csN(r) = h0*r*r^(-3/2)                      # Newtonian c_s = H Ω (the paper's Keplerian limit)
mout = Float64[]; min_ = Float64[]
for m in (5, 10, 50, 200, 1000, 5000)
    rI, rO = lindblad_locations_pressure(dN, rp, m; csound = csN)
    push!(mout, rO - rp); push!(min_, rp - rI)
    @printf("  m=%5d : r_LR,out−r_s=%+.5f   r_s−r_LR,in=%+.5f   (2H/3 = %.5f)\n", m, rO - rp, rp - rI, 2H/3)
end
@repro_result "sat_outer_offsets" mout
@repro_result "sat_inner_offsets" min_
@repro_result "sat_target_2H_over_3" 2H/3
@repro_result "sat_outer_relerr_m5000" abs(mout[end] - 2H/3)/(2H/3)

# (C3) suppression onset ξ = m h
println("\n(C3) suppression onset ξ = m h ≳ 1 ⟺ m ≳ 1/h")
for h in (0.1, 0.03, 0.01)
    @printf("  h=%.3f : ξ=1 at m=1/h=%.0f ; pressure_cutoff(1)=%.6f ; cutoff(3)=%.6f\n",
            h, 1/h, pressure_cutoff(1.0), pressure_cutoff(3.0))
end
@repro_result "cutoff_at_xi_1" pressure_cutoff(1.0)
@repro_result "cutoff_at_xi_3" pressure_cutoff(3.0)
@repro_result "cutoff_at_xi_0" pressure_cutoff(0.0)
cutoff_half_width() = begin
    lo, hi = 0.0, 5.0
    for _ in 1:200
        mid = (lo + hi)/2
        pressure_cutoff(mid) > 0.5 ? (lo = mid) : (hi = mid)
    end
    (lo + hi)/2
end
@repro_result "cutoff_half_width" cutoff_half_width()
@printf("  cutoff drops below 1/2 at ξ = %.6f\n", cutoff_half_width())

exit(ReproHarness.conclude!(opts; paper = "2510.02433", id = "check_satscale",
                            script = @__FILE__, rtol = 1e-9))
