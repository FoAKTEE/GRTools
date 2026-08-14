# Kerr (and a=0) background quantities for the rescaled scalar Teukolsky equation: the mode
# potential V_{ℓm}(r), the layer-coupling f, and the imaginary friction μ. M=1 by default;
# the `_`-shims preserve the explicit-M surface via r*(r;M,a) = M·r*(r/M;1,a/M),
# Δ(r;M,a) = M²·Δ(r/M;1,a/M).
#
# CONSTRAINT — the inverse tortoise stays LOCAL (not Spacetime.tortoise_kerr_inverse): at
# deep-horizon nodes r−r₊ underflows and each solver lands on its own numeric floor
# (Spacetime: r₊+eps; here: r₊(1+1e-12)); V_{ℓm}, f ∝ Δ inherit the ~650× floor difference,
# and the baked characterization fixtures resolve it — changing the floor breaks them.

"Kerr Δ = r² - 2Mr + a² (M-explicit shim over the M=1 `Spacetime.delta_kerr` convention)."
@inline _Δ(r, M, a) = r^2 - 2M * r + a^2

"Kerr horizons (r₊, r₋) = M ± √(M²−a²) (M-explicit shim over the M=1 `Spacetime.horizons`)."
@inline _horizons(M, a) = (M + sqrt(M^2 - a^2), M - sqrt(M^2 - a^2))

"tortoise r*(r; M, a) via `Spacetime.tortoise_kerr` (M = 1) with the exact M-scaling."
@inline _tortoise(r; M=1.0, a=0.0) = M * tortoise_kerr(r / M, a / M)

"""
inverse tortoise: solve `_tortoise(r) = rs` for `r ∈ (r₊, ∞)`. Newton from a good guess,
bisection fallback — the repo `Backgrounds.inverse_tortoise` VERBATIM (its near-horizon floor
`r₊(1+1e-12)` and `tol=1e-13` are load-bearing for the baked characterization; see header),
with the r* closed form delegated to `Spacetime.tortoise_kerr`.
"""
function _inverse_tortoise(rs; M=1.0, a=0.0, tol=1e-13, maxit=200)
    rp, _ = _horizons(M, a)
    # initial guess: large r ⇒ r≈rs; near horizon ⇒ exponential approach
    r = rs > 4M ? rs : rp + 2M * exp((rs - rp) / (2M) - 1)
    r = max(r, rp * (1 + 1e-12))
    drdrs(r) = _Δ(r, M, a) / (r^2 + a^2)
    for _ in 1:maxit
        f = _tortoise(r; M=M, a=a) - rs
        abs(f) < tol * (1 + abs(rs)) && return r
        step = f * drdrs(r)                 # Newton: r -= f / (drs/dr)
        rnew = r - step
        if !(rnew > rp) || !isfinite(rnew)  # bisection-ish safeguard
            rnew = (r + rp) / 2 + (r - rp) / 4 * sign(-f)
            rnew = max(rnew, rp * (1 + 1e-13))
        end
        abs(rnew - r) < tol * (1 + r) && return rnew
        r = rnew
    end
    return r
end

"""
    potential_V(r, ℓ, m; M=1, a=0) -> V_{ℓm}

Rescaled scalar-Teukolsky mode potential (Paper 1 eq. potential):
`V = 3r²Δ²/(r²+a²)⁴ − 2rΔ(r−M)/(r²+a²)³ − Δ²/(r²+a²)³ + a²m²/(r²+a²)² − ℓ(ℓ+1)Δ/(r²+a²)²`.
At `a=0` this equals `−(1−2M/r)[ℓ(ℓ+1)/r² + 2M/r³]` (standard Schwarzschild scalar potential,
negated to the paper's sign convention `−ψ̈+ψ''+Vψ=0`).
"""
function potential_V(r, ℓ, m; M=1.0, a=0.0)
    d = _Δ(r, M, a)
    s = r^2 + a^2
    return 3r^2 * d^2 / s^4 - 2r * d * (r - M) / s^3 - d^2 / s^3 +
           a^2 * m^2 / s^2 - ℓ * (ℓ + 1) * d / s^2
end

"""
    coupling_f(r; M=1, a=0)

`f = Δ a²/(r²+a²)²` (Kerr inter-mode coupling coefficient, Paper 1 §II.E; 0 for a=0).
"""
@inline coupling_f(r; M=1.0, a=0.0) = _Δ(r, M, a) * a^2 / (r^2 + a^2)^2

"""
    friction_mu(r, m; M=1, a=0)

`μ = 4 i m M a r/(r²+a²)²` (imaginary friction of the m≠0 Kerr sector, Paper 1 §II;
0 for a=0 or m=0).
"""
@inline friction_mu(r, m; M=1.0, a=0.0) = 4im * m * M * a * r / (r^2 + a^2)^2
