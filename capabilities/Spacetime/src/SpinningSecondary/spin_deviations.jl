# Linear-in-spin worldline deviations (spin-curvature coupling) δt3/δr3/δz3/δφ3 — consolidates
# SpinningOrbit/src/spin_deviations.jl (faithful 1:1 of SpinningTrajectory.wl L272-275). These are
# the O(s) position shifts of the spinning secondary's worldline from the geodesic, evaluated on
# the (geodesic) trajectory (r,z) and four-velocity (Ur,Uz) with constants (En,Jz,K) and spin a.
# Δ=r²−2r+a², Σ=r²+a²z². Pure arithmetic ⇒ type-generic {T}, BigFloat inherited.

"""`_δt3(r,z,Ur,Uz,K,a)` — linear-in-spin t-deviation (SpinningTrajectory.wl L272). R2."""
_δt3(r, z, Ur, Uz, K, a) =
    (r * (a^2 + r^2) * Ur - a^2 * (a^2 + r * (-2 + r)) * z * Uz) /
    (sqrt(K) * (a^2 + r * (-2 + r)) * (r^2 + a^2 * z^2))

"""`_δr3(r,z,En,Jz,K,a)` — linear-in-spin r-deviation (SpinningTrajectory.wl L273). R2."""
_δr3(r, z, En, Jz, K, a) =
    (r * (a^2 * En - a * Jz + En * r^2)) / (sqrt(K) * (r^2 + a^2 * z^2))

"""`_δz3(r,z,En,Jz,K,a)` — linear-in-spin z-deviation (SpinningTrajectory.wl L274). R2."""
_δz3(r, z, En, Jz, K, a) =
    (a * z * (Jz - a * En * (1 - z^2))) / (sqrt(K) * (r^2 + a^2 * z^2))

"""`_δφ3(r,z,Ur,Uz,K,a)` — linear-in-spin φ-deviation (SpinningTrajectory.wl L275). R2."""
_δφ3(r, z, Ur, Uz, K, a) =
    (a * (r * (1 - z^2) * Ur - (a^2 + r * (-2 + r)) * z * Uz)) /
    (sqrt(K) * (a^2 + r * (-2 + r)) * (1 - z^2) * (r^2 + a^2 * z^2))
