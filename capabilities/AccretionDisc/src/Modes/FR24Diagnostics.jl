# FR24 v1: ref-paper/arXiv-2407.20398/src/Sections/2_linear_framework.tex:268-317.
# NewtonianDisc units: GM = a = Omega_p = Sigma_p = 1 for the FR24 reproduction.

"""
    fr24_corotation_radius(d, mode)

FR24 (35), the O(h_p²) radius in the disc's length units. `mode.omega` is the
TOTAL frequency; Omega(r_C) = mode.omega / mode.m. This is an approximation,
not the grid's corotation index. Use `bisect(r -> Omega(d,r)-mode.omega/m,...)`
for an independent exact-root diagnostic. No background profile is duplicated.
"""
function fr24_corotation_radius(d::NewtonianDisc, mode::HarmonicMode)
    mode.m > 0 && mode.omega > 0 || throw(ArgumentError("positive m and frequency required"))
    frequency, _ = promote(mode.omega, d.h0)
    r0 = (frequency / mode.m)^(-2//3)
    r0 * (1 - d.h0^2 * (d.p + d.q) * r0^(1-d.q) / 3)
end
