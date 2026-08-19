# Teukolsky point-particle fluxes: gravitational-wave energy flux from circular orbits.
#
# Run from the repository root:
#     julia --project=. examples/04_teukolsky_fluxes.jl
#
# TeukolskyPointParticleModeCircular(s, l, m, a, p) solves the radial Teukolsky
# equation with a circular-orbit point-particle source and returns amplitudes and
# fluxes (EnergyFlux/AngularMomentumFlux, each split into infinity .I and horizon .H).

using Perturbation
using Printf

include(joinpath(@__DIR__, "_check.jl"))   # check / check_close: every printed cross-check is asserted

a, p = 0.5, 8.0

Ωφ = 1 / (p^1.5 + a)          # circular equatorial orbital frequency (mode omega = m * Omega_phi)

println("Dominant modes, circular orbit a = $a, p = $p (mass ratio scaled out):")
println("  (l, m)    omega        Edot_inf        Edot_hor        Ldot_inf")
total = 0.0
for l in 2:4, m in 1:l
    mode = TeukolskyPointParticleModeCircular(-2, l, m, a, p)
    ef, jf = mode.EnergyFlux, mode.AngularMomentumFlux
    # +/- m symmetry: the m<0 partner carries the same energy flux.
    global total += 2 * (ef.I + ef.H)
    @printf "  (%d,%2d)  %9.6f  %14.6e  %14.6e  %14.6e\n" l m mode.omega ef.I ef.H jf.I
    check_close("mode ($l,$m): omega == m * Omega_phi", mode.omega, m * Ωφ; rtol = 1e-10)
    check("mode ($l,$m): infinity-side energy flux is positive", ef.I > 0)
    # a = 0.5 prograde, m > 0: omega < m*Omega_H, so every mode here is superradiant and the
    # horizon flux is NEGATIVE (the hole loses energy to the wave).
    check("mode ($l,$m): superradiant horizon flux is negative (Edot_hor = $(ef.H))", ef.H < 0)
    check_close("mode ($l,$m): Ldot_inf == (m/omega) Edot_inf", jf.I, (m / mode.omega) * ef.I; rtol = 1e-8)
end
@printf "\nTotal Edot (l<=4, both signs of m) = %.8e\n" total
check("total energy flux is positive and finite", isfinite(total) && total > 0)

# Consistency: the angular-momentum flux obeys Ldot = (m/omega) Edot mode by mode.
mode22 = TeukolskyPointParticleModeCircular(-2, 2, 2, a, p)
ratio22 = mode22.AngularMomentumFlux.I / (2 / mode22.omega * mode22.EnergyFlux.I)
@printf "J/E identity (2,2): Ldot / (m/omega * Edot) = %.12f (exact 1)\n" ratio22
check_close("J/E identity Ldot = (m/omega) Edot for the (2,2) mode", ratio22, 1.0; atol = 1e-10)

# Newtonian check: at large p the (2,+-2) pair approaches the quadrupole 32/5 p^-5.
pfar = 80.0
e22 = TeukolskyPointParticleModeCircular(-2, 2, 2, 0.0, pfar).EnergyFlux.I
e2m2 = TeukolskyPointParticleModeCircular(-2, 2, -2, 0.0, pfar).EnergyFlux.I
newt = (32 / 5) * pfar^-5
ratio_newt = (e22 + e2m2) / newt
@printf "\nSchwarzschild p=%.0f quadrupole: Edot(2,2)+Edot(2,-2) = %.6e vs 32/5 p^-5 = %.6e (ratio %.4f)\n" pfar (e22 + e2m2) newt ratio_newt
# At p = 80 the leading 1-PN correction is O(1/p) ~ few %, so the ratio must sit just under 1.
check("(2,+-2) flux approaches the Newtonian 32/5 p^-5 limit (ratio = $ratio_newt)",
      0.90 < ratio_newt < 1.02)

# Eccentric orbits expose the (l, m, n) radial-harmonic structure.
md = TeukolskyPointParticleModeEccentric(2, 2, 1, a, p, 0.2)
@printf "\nEccentric (l,m,n)=(2,2,1), e=0.2:  omega = %.8f, Edot_inf = %.6e\n" md.omega md.EnergyFlux.I
# Eccentric radial harmonic n = 1 sits one Omega_r above the n = 0 (circular-limit) frequency.
check("eccentric (2,2,1) mode has a positive, finite energy flux",
      isfinite(md.EnergyFlux.I) && md.EnergyFlux.I > 0)
check("eccentric (2,2,1) frequency exceeds the circular (2,2,0) frequency",
      md.omega > mode22.omega)
