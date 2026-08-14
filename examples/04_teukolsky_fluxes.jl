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

a, p = 0.5, 8.0

println("Dominant modes, circular orbit a = $a, p = $p (mass ratio scaled out):")
println("  (l, m)    omega        Edot_inf        Edot_hor        Ldot_inf")
total = 0.0
for l in 2:4, m in 1:l
    mode = TeukolskyPointParticleModeCircular(-2, l, m, a, p)
    ef, jf = mode.EnergyFlux, mode.AngularMomentumFlux
    # +/- m symmetry: the m<0 partner carries the same energy flux.
    global total += 2 * (ef.I + ef.H)
    @printf "  (%d,%2d)  %9.6f  %14.6e  %14.6e  %14.6e\n" l m mode.omega ef.I ef.H jf.I
end
@printf "\nTotal Edot (l<=4, both signs of m) = %.8e\n" total

# Consistency: the angular-momentum flux obeys Ldot = (m/omega) Edot mode by mode.
mode22 = TeukolskyPointParticleModeCircular(-2, 2, 2, a, p)
@printf "J/E identity (2,2): Ldot / (m/omega * Edot) = %.12f (exact 1)\n" mode22.AngularMomentumFlux.I / (2 / mode22.omega * mode22.EnergyFlux.I)

# Newtonian check: at large p the (2,+-2) pair approaches the quadrupole 32/5 p^-5.
pfar = 80.0
e22 = TeukolskyPointParticleModeCircular(-2, 2, 2, 0.0, pfar).EnergyFlux.I
e2m2 = TeukolskyPointParticleModeCircular(-2, 2, -2, 0.0, pfar).EnergyFlux.I
newt = (32 / 5) * pfar^-5
@printf "\nSchwarzschild p=%.0f quadrupole: Edot(2,2)+Edot(2,-2) = %.6e vs 32/5 p^-5 = %.6e (ratio %.4f)\n" pfar (e22 + e2m2) newt ((e22 + e2m2) / newt)

# Eccentric orbits expose the (l, m, n) radial-harmonic structure.
md = TeukolskyPointParticleModeEccentric(2, 2, 1, a, p, 0.2)
@printf "\nEccentric (l,m,n)=(2,2,1), e=0.2:  omega = %.8f, Edot_inf = %.6e\n" md.omega md.EnergyFlux.I
