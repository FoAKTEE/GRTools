# Hypergeometric layer for the MST radial sums: exact ₂F₁ (In solution) and Tricomi-U (Up
# solution) per-n terms + their up/down recurrence arrays. ₂F₁/U via HypergeometricFunctions.jl
# (research-grade complex/branch handling).
#
# AHA (design §4.3 dedup discipline): the MST parameter wrappers (_aF/_bF/_cF/_aU) and the
# _safe_2f1 |z|>6 / near-integer branch-SELECTION heuristic are a PER-CALLER decision and stay
# LOCAL here. The generic ₂F₁ pieces — the no-Γ/ψ Maclaurin series and the z→∞ connection — are
# NOT duplicated: the |z|≤6 Pfaff dispatcher delegates DOWN to BigSpecial.safe_2f1 and the large-|z|
# route to BigSpecial.f21_zinf (the canonical primitives, already characterization-gated to this same
# Teukolsky 2F1 in BigSpecial's core suite). Γ of a Complex{BigFloat} routes through BigSpecial.gamma_big.

import HypergeometricFunctions

# Teukolsky ₂F₁ / U parameter sets (Sasaki–Tagoshi §2.1) — LOCAL MST wrappers (AHA kept local).
_aF(s, nu, tau, eps) = nu + 1 - im * tau
_bF(s, nu, tau, eps) = -nu - im * tau
_cF(s, nu, tau, eps) = 1 - s - im * (eps + tau)
_aU(s, nu, tau, eps) = nu + s + 1 - im * eps

"""
    _safe_2f1(a, b, c, z) -> Complex

Robust Gauss ₂F₁: byte-identical to `HypergeometricFunctions._₂F₁` for Float64/ComplexF64 (it just
calls it) but does NOT fail for Complex{BigFloat} at |z|≳½ — there upstream's analytic-continuation
branches reach a `loggamma`/`digamma` of a Complex{BigFloat} (unimplemented) and throw. On that
throw, for the BigFloat path the LOCAL branch heuristic selects:
  • |z|>6 and a−b away from an integer ⇒ the z→∞ connection (FAST: the only viable route for the
    κ-tail r→1.5e4 where the Pfaff w→1 needs ~1e5 terms) — delegated to `BigSpecial.f21_zinf`;
  • otherwise the whole |z|≤6 Pfaff-transform dispatcher (z→w=z/(z−1), |w|<1 for Re z<½, ∀ MST
    In-coordinate x≤0, + the no-Γ/ψ Maclaurin series) is delegated VERBATIM to the canonical
    `BigSpecial.safe_2f1` — no longer re-copied here (DRY).
The |z|≤6 region (incl. every REF_HYP oracle row) stays on the Pfaff route.
"""
function _safe_2f1(a, b, c, z)
    try
        return HypergeometricFunctions._₂F₁(a, b, c, z)
    catch err
        (err isa DomainError || err isa MethodError) || rethrow(err)
        T = promote_type(typeof(a), typeof(b), typeof(c), typeof(z))   # working (BigFloat) precision
        aT = T(a); bT = T(b); cT = T(c); zT = T(z)
        ab = aT - bT; near_int = abs(real(ab) - round(real(ab))) < 1e-3 && abs(imag(ab)) < 1e-3
        if abs(zT) > 6 && !near_int
            return f21_zinf(aT, bT, cT, zT)                            # FAST large-|z| path (BigSpecial)
        end
        return safe_2f1(aT, bT, cT, zT)   # |z|≤6: canonical BigSpecial Pfaff+Maclaurin dispatcher
    end
end

"Exact ₂F₁ term `₂F₁[n+aF, bF-n, cF, x]` (Sasaki–Tagoshi §2.6 H2F1Exact)."
function _h2f1_exact(n::Integer, s, nu, tau, eps, x)
    a = n + _aF(s, nu, tau, eps); b = _bF(s, nu, tau, eps) - n; c = _cF(s, nu, tau, eps)
    return _safe_2f1(a, b, c, complex(x))
end

"Exact U term `cⁿ U[n+aU, 2n+b, c]`, b=2ν+2, c=−2i·zhat (Sasaki–Tagoshi §2.8 HUExact)."
function _hu_exact(n::Integer, s, nu, tau, eps, zhat)
    b = 2 * nu + 2
    c = -2 * im * zhat
    return c^n * HypergeometricFunctions.U(n + _aU(s, nu, tau, eps), 2n + b, c)
end

"Exact dHU/dzhat term (§2.9 dHUExact) = −2i(cⁿ⁻¹ n U[n+aU,2n+b,c] − cⁿ(n+aU) U[1+n+aU,1+2n+b,c])."
function _dhu_exact(n::Integer, s, nu, tau, eps, zhat)
    a = _aU(s, nu, tau, eps); b = 2 * nu + 2; c = -2 * im * zhat
    return -2 * im * (c^(n - 1) * n * HypergeometricFunctions.U(n + a, 2n + b, c) -
                      c^n * (n + a) * HypergeometricFunctions.U(1 + n + a, 1 + 2n + b, c))
end

"""
    _hu_array(s, nu, tau, eps, zhat, Nmax) -> Dict n=>HU[n]

HU[n] for n=−Nmax..Nmax via the up/down RECURRENCE (§2.8), seeded by exact HU[0],HU[1]. The
recurrence is REQUIRED for accuracy: exact `cⁿU` per-n degrades at the contributing n (the
cancellation that gave ~1e-4 R^up error; the recurrence is what BHPToolkit uses).
"""
function _hu_array(s, nu, tau, eps, zhat, Nmax::Integer)
    a = _aU(s, nu, tau, eps); b = 2 * nu + 2; c = -2 * im * zhat
    HU = Dict{Int,typeof(complex(_hu_exact(0, s, nu, tau, eps, zhat)))}()  # generic (F64/BigFloat)
    HU[0] = _hu_exact(0, s, nu, tau, eps, zhat)
    HU[1] = _hu_exact(1, s, nu, tau, eps, zhat)
    for n in 2:Nmax            # HUUp (§2.8): uses HU[n-1], HU[n-2]
        t1 = (-2 - a + b + n) * (-2 + b + 2n) * HU[n-2]
        t2 = (-3 + b + 2n) * (8 + (b + 2n)^2 + 2 * (a + n) * c - (b + 2n) * (6 + c)) * HU[n-1] / c
        HU[n] = (t1 + t2) / ((-1 + a + n) * (-4 + b + 2n))
    end
    for n in -1:-1:-Nmax       # HUDown (§2.8): uses HU[n+1], HU[n+2]
        t1 = -((1 + b + 2n) * (b^2 + 4n * (1 + n) + b * (2 + 4n - c) + 2 * a * c) * HU[n+1]) / c
        t2 = (1 + a + n) * (b + 2n) * HU[n+2]
        HU[n] = (t1 + t2) / ((-a + b + n) * (2 + b + 2n))
    end
    return HU
end

"""
    _dhu_array(s, nu, tau, eps, zhat, HU, Nmax) -> Dict n=>dHU[n]

dHU[n]/dzhat for n=−Nmax..Nmax via the §2.9 up/down recurrence (seeds dHU[0],dHU[1] exact;
consumes the HU array). Required for accuracy at large n, like `_hu_array`.
"""
function _dhu_array(s, nu, tau, eps, zhat, HU, Nmax::Integer)
    a = _aU(s, nu, tau, eps); b = 2 * nu + 2; c = -2 * im * zhat
    dHU = Dict{Int,typeof(complex(_dhu_exact(0, s, nu, tau, eps, zhat)))}()  # generic (F64/BigFloat)
    dHU[0] = _dhu_exact(0, s, nu, tau, eps, zhat)
    dHU[1] = _dhu_exact(1, s, nu, tau, eps, zhat)
    for n in 2:Nmax            # dHUUp (§2.9)
        t1 = (-2 - a + b + n) * (-2 + b + 2n) * dHU[n-2] / (-4 + b + 2n)
        t2 = (-3 + b + 2n) * (8 + b^2 + 4 * (-3 + n) * n + b * (-6 + 4n - c) + 2 * a * c) * dHU[n-1] / ((-4 + b + 2n) * c)
        t3 = 2 * im * (-3 + b + 2n) * (-2 + b + 2n) * HU[n-1] / c^2
        dHU[n] = (t1 + t2 + t3) / (-1 + a + n)
    end
    for n in -1:-1:-Nmax       # dHUDown (§2.9)
        t1 = (1 + b + 2n) * c * (b^2 + 4n * (1 + n) + b * (2 + 4n - c) + 2 * a * c) * dHU[n+1]
        t2 = (b + 2n) * (-(1 + a + n) * c^2 * dHU[n+2])
        t3 = (b + 2n) * (2 * im * (1 + b + 2n) * (2 + b + 2n) * HU[n+1])
        dHU[n] = (t1 + t2 + t3) / ((a - b - n) * (2 + b + 2n) * c^2)
    end
    return dHU
end
