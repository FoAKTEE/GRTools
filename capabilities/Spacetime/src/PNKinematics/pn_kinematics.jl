# PNKinematics — Schwarzschild–Kepler circular-orbit kinematic identities (M = 1).
# Oracle: WaSABI 0PASchwarz.m — Ω = √(1/r0³) (InspiralModels:13), r0 = Ω^(-2/3) (AmplitudeModels:11).

"""
    kepler_omega(r0)

Schwarzschild circular-orbit orbital frequency Ω(r0) = √(1/r0³), in units M = 1.
"""
kepler_omega(r0) = sqrt(1 / r0^3)

"""
    r0_of_omega(Ω)

Areal radius r0(Ω) = Ω^(-2/3) of a Schwarzschild circular orbit (M = 1), the inverse of
[`kepler_omega`](@ref).
"""
r0_of_omega(Ω) = Ω^(-2 // 3)

"""
    pn_x(Ω)

Post-Newtonian frequency variable x = Ω^(2/3) (so r0 = 1/x on a Kepler orbit), M = 1.
"""
pn_x(Ω) = Ω^(2 // 3)

# Unicode aliases (DRY: keep both spellings for callers on either convention).
const kepler_Ω = kepler_omega
const r0_of_Ω = r0_of_omega
