# Geometry — Schwarzschild (a = 0) tortoise coordinate, effective-potential, and
# circular-orbit constants of motion (M = 1). Consolidates the H1Lorenz/DiscModes
# lorenz_metric orbit fork (crosswalk row 358). Pure arithmetic ⇒ BigFloat inherited;
# every value is built at T = element type of the input, never down-cast to Float64.

"Schwarzschild lapse-squared f(r) = 1 − 2M/r (M = 1)."
@inline schwarzschild_lapse(r) = 1 - 2 / r

"Radial derivative f'(r) = 2M/r² (M = 1)."
@inline schwarzschild_dlapse(r) = 2 / r^2

"""
    bls_potential(r, l)

Barack–Lousto–Sago effective potential Vₗ(r) = ¼ f(r) (f'(r)/r + l(l+1)/r²) of the
Schwarzschild Regge–Wheeler/Lorenz-gauge radial equation (M = 1).
"""
@inline bls_potential(r, l) =
    (one(r) / 4) * schwarzschild_lapse(r) * (schwarzschild_dlapse(r) / r + l * (l + 1) / r^2)

"""
    tortoise_schw(r)

Schwarzschild tortoise coordinate r*(r) = r + 2 ln(r/2 − 1) (M = 1), r > 2.
"""
@inline tortoise_schw(r) = r + 2 * log(r / 2 - 1)

"""
    tortoise_schw_inverse(rs; tol, itmax) -> r

Inverse tortoise: solve r* = r + 2 ln(r/2 − 1) for r > 2 (M = 1), r* ∈ (−∞, ∞).
Newton iteration is performed in the cancellation-free variable `q = log(r/2-1)`, where
`q + exp(q) = r*/2 - 1`. If the exact exterior radius is closer to 2 than the output type can
represent, the nearest representable exterior point `nextfloat(2)` is returned. Type-generic — the
tolerance inherits precision through `eps`.
"""
function tortoise_schw_inverse(rs::Real; tol=nothing, itmax::Int = 200)
    rss = float(rs)
    T = typeof(rss)
    isfinite(rss) || throw(ArgumentError("tortoise_schw_inverse: rs must be finite (got $rs)"))
    tolv = tol === nothing ? 8 * eps(T) : T(tol)
    tolv >= zero(T) || throw(ArgumentError("tortoise_schw_inverse: tol must be nonnegative"))
    c = rss / 2 - one(T)
    q = c < zero(T) ? c : log1p(c)
    for _ in 1:itmax
        y = exp(q)
        step = (q + y - c) / (one(T) + y)
        q -= step
        if abs(step) <= tolv * (one(T) + abs(q))
            r = 2 * (one(T) + exp(q))
            return r > T(2) ? r : nextfloat(T(2))
        end
    end
    error("tortoise_schw_inverse did not converge in $itmax iterations")
end

"Circular-orbit specific energy E = (1 − 2/p)/√(1 − 3/p) (Schwarzschild, M = 1)."
@inline schw_orbit_energy(p) = (1 - 2 / p) / sqrt(1 - 3 / p)

"Circular-orbit specific angular momentum L = √p/√(1 − 3/p) (Schwarzschild, M = 1)."
@inline schw_orbit_ang_mom(p) = sqrt(p) / sqrt(1 - 3 / p)

"Circular-orbit azimuthal frequency Ω_φ = p^(−3/2) (Schwarzschild Kepler law, M = 1)."
@inline schw_orbit_phi_freq(p) = p^(-3 // 2)
