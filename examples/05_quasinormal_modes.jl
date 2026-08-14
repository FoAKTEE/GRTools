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

s, l, m, n = -2, 2, 2, 0

# Schwarzschild fundamental mode — the classic 0.373672 - 0.088962i.
ω0 = qnm_frequency(s, l, m, n, 0.0)
@printf "Schwarzschild (2,2,0):  M*omega = %.6f %+.6fi\n\n" real(ω0) imag(ω0)

println("Spin sequence (seeded continuation up the chi ladder):")
println("   chi      Re(M omega)   Im(M omega)   |leaver - fast|")
ωL, ωF = ω0, qnm_frequency(s, l, m, n, 0.0; backend=:fast)
for χ in 0.1:0.1:0.9
    global ωL = qnm_frequency(s, l, m, n, χ; seed=ωL, backend=:leaver)
    global ωF = qnm_frequency(s, l, m, n, χ; seed=ωF, backend=:fast)
    @printf "  %.2f    %10.6f    %10.6f     %.2e\n" χ real(ωL) imag(ωL) abs(ωL - ωF)
end

# Overtone: n = 1 damps faster than the fundamental.
ω1 = qnm_frequency(s, l, m, 1, 0.0)
@printf "\nSchwarzschild (2,2,1) overtone: M*omega = %.6f %+.6fi (|Im| grows with n)\n" real(ω1) imag(ω1)
