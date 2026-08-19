# Kerr quasinormal modes: the (s=-2, l=2, m=2, n=0) frequency across black-hole spin.
#
# Run from the repository root:
#     julia --project=. examples/05_quasinormal_modes.jl
#
# qnm_frequency(s, l, m, n, chi) returns the complex QNM frequency M*omega
# (Re = ringing frequency, Im < 0 = damping). Two backends: :leaver (Leaver
# continued fraction) and :fast; both are continued-fraction-residual roots.

using Perturbation
using Printf

include(joinpath(@__DIR__, "_check.jl"))   # check / check_close: every printed cross-check is asserted

s, l, m, n = -2, 2, 2, 0

# Schwarzschild fundamental mode — the classic 0.373672 - 0.088962i.
ω0 = qnm_frequency(s, l, m, n, 0.0)
@printf "Schwarzschild (2,2,0):  M*omega = %.6f %+.6fi\n\n" real(ω0) imag(ω0)
check_close("Schwarzschild (2,2,0) == the tabulated 0.373672 - 0.088962i",
            ω0, 0.373672 - 0.088962im; atol = 2e-6)

println("Spin sequence (seeded continuation up the chi ladder):")
println("   chi      Re(M omega)   Im(M omega)   |leaver - fast|")
ωL, ωF = ω0, qnm_frequency(s, l, m, n, 0.0; backend=:fast)
for χ in 0.1:0.1:0.9
    global ωL = qnm_frequency(s, l, m, n, χ; seed=ωL, backend=:leaver)
    global ωF = qnm_frequency(s, l, m, n, χ; seed=ωF, backend=:fast)
    @printf "  %.2f    %10.6f    %10.6f     %.2e\n" χ real(ωL) imag(ωL) abs(ωL - ωF)
    check_close("chi = $χ: :leaver and :fast backends agree", ωL, ωF; atol = 1e-12)
    check("chi = $χ: the mode is damped (Im M omega < 0)", imag(ωL) < 0)
end
# Prograde spin-up: the (2,2,0) ringing frequency rises and the damping weakens.
check("Re(M omega) rises with spin (chi = 0 -> 0.9)", real(ωL) > real(ω0))
check("|Im(M omega)| falls with spin (chi = 0 -> 0.9)", abs(imag(ωL)) < abs(imag(ω0)))

# Overtone: n = 1 damps faster than the fundamental.
ω1 = qnm_frequency(s, l, m, 1, 0.0)
@printf "\nSchwarzschild (2,2,1) overtone: M*omega = %.6f %+.6fi (|Im| grows with n)\n" real(ω1) imag(ω1)
check("overtone n=1 damps faster than n=0", abs(imag(ω1)) > abs(imag(ω0)))
