# arXiv:2510.02433 fig:discProfile (RT4) — Novikov–Thorne relativistic thin-disc surface density
# Σ(r) [g/cm²] and aspect ratio h(r) = H/r for a/M = 0 and 0.9, with the Shakura–Sunyaev power-law
# overlay (Σ ∝ r^{3/2}, h ∝ 1/r) normalized to coincide with the a=0 curve in the Newtonian limit.
#
# THIN driver: `AccretionDisc.nt_surface_density` / `nt_aspect_ratio` / `ss_surface_density` /
# `ss_aspect_ratio` (Novikov–Thorne inner region, Abramowicz–Fragile Eqs. 99 with the Page–Thorne
# correction factors), `Spacetime.kerr_isco` for the inner edge.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2510.02433/figDiscProfile.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: NovikovThorneDisc, nt_factors, nt_surface_density, nt_aspect_ratio,
                             ss_surface_density, ss_aspect_ratio, nt_isco
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
using Plots

nt0 = NovikovThorneDisc(a = 0.0)
nt9 = NovikovThorneDisc(a = 0.9)
grid(nt) = exp10.(range(log10(nt_isco(nt) + 0.05), log10(1000.0); length = 200))
g0 = grid(nt0); g9 = grid(nt9)
@repro_result "r_grid_a0" g0
@repro_result "r_grid_a09" g9
S0 = nt_surface_density.(Ref(nt0), g0); S9 = nt_surface_density.(Ref(nt9), g9)
h0 = nt_aspect_ratio.(Ref(nt0), g0);    h9 = nt_aspect_ratio.(Ref(nt9), g9)
Sss = ss_surface_density.(Ref(nt0), g0); hss = ss_aspect_ratio.(Ref(nt0), g0)
@repro_result "Sigma_NT_a0" S0
@repro_result "Sigma_NT_a09" S9
@repro_result "h_NT_a0" h0
@repro_result "h_NT_a09" h9
@repro_result "Sigma_SS" Sss
@repro_result "h_SS" hss
@repro_result "isco_a0" nt_isco(nt0)
@repro_result "isco_a09" nt_isco(nt9)

# published anchor + structural landmarks
@repro_result "h_a0_at_1000M" nt_aspect_ratio(nt0, 1000.0)
for (lab, nt, g, S, h) in (("a=0", nt0, g0, S0, h0), ("a=0.9", nt9, g9, S9, h9))
    iS = argmin(S); ih = argmax(h)
    @repro_result "Sigma_min_$(lab)" S[iS]
    @repro_result "Sigma_min_radius_$(lab)" g[iS]
    @repro_result "h_max_$(lab)" h[ih]
    @repro_result "h_max_radius_$(lab)" g[ih]
    @printf("%-5s : Σ_min=%.4e g/cm² at r=%.3f M ; h_max=%.5f at r=%.3f M\n", lab, S[iS], g[iS], h[ih], g[ih])
end
for r in (50.0, 250.0, 1000.0)
    @printf("  r=%6.0f  Σ_NT(a0)/Σ_SS=%.4f  Σ_NT(a.9)/Σ_SS=%.4f  h_NT(a0)/h_SS=%.4f  h_NT(a.9)/h_SS=%.4f\n",
            r, nt_surface_density(nt0, r)/ss_surface_density(nt0, r),
               nt_surface_density(nt9, r)/ss_surface_density(nt0, r),
               nt_aspect_ratio(nt0, r)/ss_aspect_ratio(nt0, r),
               nt_aspect_ratio(nt9, r)/ss_aspect_ratio(nt0, r))
end
f = nt_factors(0.9, 10.0)
@repro_result "factors_a09_r10" [f.A, f.B, f.C, f.D, f.E, f.Q]

emri_defaults!()
const CblueSS = RGB(0.30, 0.45, 0.90); const Corange = RGB(0.95, 0.55, 0.10); const Cred = RGB(0.86, 0.18, 0.18)
pS = plot(; xlabel = "\$r/M\$", ylabel = "\$\\Sigma\\ [{\\rm g/cm^2}]\$", xscale = :log10, yscale = :log10,
          legend = :topleft, xlims = (3, 1000))
plot!(pS, g0, Sss; color = CblueSS, ls = :dash, lw = 2, label = "Shakura - Sunyaev (\$\\Sigma\\propto r^{3/2}\$)")
plot!(pS, g0, S0; color = Corange, lw = 2.2, label = "Novikov - Thorne, \$a/M=0\$")
plot!(pS, g9, S9; color = Cred, lw = 2.2, label = "Novikov - Thorne, \$a/M=0.9\$")
pH = plot(; xlabel = "\$r/M\$", ylabel = "\$h\$", xscale = :log10, yscale = :log10,
          legend = :topright, xlims = (3, 1000))
plot!(pH, g0, hss; color = CblueSS, ls = :dash, lw = 2, label = "Shakura - Sunyaev (\$h\\propto 1/r\$)")
plot!(pH, g0, h0; color = Corange, lw = 2.2, label = "Novikov - Thorne, \$a/M=0\$")
plot!(pH, g9, h9; color = Cred, lw = 2.2, label = "Novikov - Thorne, \$a/M=0.9\$")
png = joinpath(@__DIR__, "expected", "figDiscProfile.png")
savefig(plot(pS, pH; layout = (1, 2), size = (1100, 430), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm), png)
@printf("wrote %s (%d bytes)\n", relpath(png, @__DIR__), filesize(png))

exit(ReproHarness.conclude!(opts; paper = "2510.02433", id = "figDiscProfile",
                            script = @__FILE__, rtol = 1e-10))
