# Separatrix — ONE canonical last-stable-orbit separatrix p_sep(a,e,x) + the Kerr ISCO.
# Stein–Warburton 2020 (arXiv:1912.07609) separatrix polynomials (transcribed from
# KerrBackwards/Utility.cc); the bracketed real root is found with NumericsCore's vendored A42
# (a42) — NO external Roots dependency. The FEW `get_separatrix` / KerrBackwards `separatrix`
# become thin aliases of this; the e=0 equatorial branch and the 4 legacy ISCO forks collapse
# to the single Bardeen `kerr_isco` closed form (AHA: same decision).

"""
    kerr_isco(a) -> r_isco

Innermost stable circular (equatorial) orbit radius (M = 1) for spin `a ∈ [-1,1]`; `a>0`
prograde, `a<0` retrograde. Bardeen–Press–Teukolsky 1972 (ApJ 178, 347) closed form.
Limits: `kerr_isco(0)=6`, `kerr_isco(1)=1` (prograde extremal), `kerr_isco(-1)=9` (retrograde).
Canonical replacement for the WaSABI/DiscModes/DiscEMRIForce/FEW ISCO forks.
"""
function kerr_isco(a::Real)
    Z1 = 1 + cbrt(1 - a^2) * (cbrt(1 - a) + cbrt(1 + a))
    Z2 = sqrt(3*a^2 + Z1^2)
    3 + Z2 - sign(a) * sqrt((3 - Z1) * (3 + Z1 + 2*Z2))
end

# Stein–Warburton 2020 separatrix polynomials; p_sep is the bracketed real root in p.
function _sep_poly_equat(p, a, e)
    a^4 * (-3 - 2e + e^2)^2 + p^2 * (-6 - 2e + p)^2 -
        2*a^2 * (1 + e) * p * (14 + 2*e^2 + 3*p - e*p)
end

function _sep_poly_polar(p, a, e)
    a^6 * (-1 + e)^2 * (1 + e)^4 + p^5 * (-6 - 2e + p) +
        a^2 * p^3 * (-4*(-1 + e)*(1 + e)^2 + (3 + e*(2 + 3e))*p) -
        a^4 * (1 + e)^2 * p * (6 + 2*e^3 + 2*e*(-1 + p) - 3*p - 3*e^2*(2 + p))
end

function _sep_poly_full(p, a, e, x)
    (-4*(3 + e)*p^11 + p^12 +
        a^12 * (-1 + e)^4 * (1 + e)^8 * (-1 + x)^4 * (1 + x)^4 -
        4*a^10 * (-3 + e)*(-1 + e)^3*(1 + e)^7 * p * (-1 + x^2)^4 -
        4*a^8 * (-1 + e)*(1 + e)^5 * p^3 * (-1 + x)^3 * (1 + x)^3 *
            (7 - 7*x^2 - e^2*(-13 + x^2) + e^3*(-5 + x^2) + 7*e*(-1 + x^2)) +
        8*a^6 * (-1 + e)*(1 + e)^3 * p^5 * (-1 + x^2)^2 *
            (3 + e + 12*x^2 + 4*e*x^2 + e^3*(-5 + 2*x^2) + e^2*(1 + 2*x^2)) -
        8*a^4 * (1 + e)^2 * p^7 * (-1 + x) * (1 + x) *
            (-3 + e + 15*x^2 - 5*e*x^2 + e^3*(-5 + 3*x^2) + e^2*(-1 + 3*x^2)) +
        4*a^2 * p^9 * (-7 - 7*e + e^3*(-5 + 4*x^2) + e^2*(-13 + 12*x^2)) +
        2*a^8 * (-1 + e)^2*(1 + e)^6 * p^2 * (-1 + x^2)^3 *
            (2*(-3 + e)^2*(-1 + x^2) + a^2*(e^2*(-3 + x^2) - 3*(1 + x^2) + 2*e*(1 + x^2))) -
        2*p^10 * (-2*(3 + e)^2 + a^2*(-3 + 6*x^2 + e^2*(-3 + 2*x^2) + e*(-2 + 4*x^2))) +
        a^6 * (1 + e)^4 * p^4 * (-1 + x^2)^2 *
            (-16*(-1 + e)^2*(-3 - 2*e + e^2)*(-1 + x^2) +
             a^2*(15 + 6*x^2 + 9*x^4 + e^2*(26 + 20*x^2 - 2*x^4) +
                  e^4*(15 - 10*x^2 + x^4) + 4*e^3*(-5 - 2*x^2 + x^4) -
                  4*e*(5 + 2*x^2 + 3*x^4))) -
        4*a^4 * (1 + e)^2 * p^6 * (-1 + x) * (1 + x) *
            (-2*(11 - 14*e^2 + 3*e^4)*(-1 + x^2) +
             a^2*(5 - 5*x^2 - 9*x^4 + 4*e^3*x^2*(-2 + x^2) +
                  e^4*(5 - 5*x^2 + x^4) + e^2*(6 - 6*x^2 + 4*x^4))) +
        a^2 * p^8 *
            (-16*(1 + e)^2*(-3 + 2*e + e^2)*(-1 + x^2) +
             a^2*(15 - 36*x^2 + 30*x^4 + e^4*(15 - 20*x^2 + 6*x^4) +
                  4*e^3*(5 - 12*x^2 + 6*x^4) + 4*e*(5 - 12*x^2 + 10*x^4) +
                  e^2*(26 - 72*x^2 + 44*x^4))))
end

"""
    separatrix(a, e, x) -> p_sep

Last-stable-orbit separatrix `p_sep(a, e, x)` (Stein–Warburton 2020, M = 1, `x = cos I`).
Branches: a=0 → `6+2e`; e=0,|x|=1 → `kerr_isco`; |x|=1 → root of the equatorial polynomial;
x=0 → root of the polar polynomial; otherwise the generic-polynomial root bracketed between
the equatorial and polar separatrices (sign of x selects the bracket).
"""
function separatrix(a::Real, e::Real, x::Real)
    T = float(promote_type(typeof(a), typeof(e), typeof(x)))
    a, e, x = T(a), T(e), T(x)
    if a == 0
        return 6 + 2e
    elseif e == 0 && abs(x) == 1
        return kerr_isco(x > 0 ? a : -a)
    elseif x == 1
        return a42(p -> _sep_poly_equat(p, a, e), 1 + e, 6 + 2e)
    elseif x == -1
        return a42(p -> _sep_poly_equat(p, a, e), 6 + 2e, 5 + e + 4*sqrt(1 + e))
    else
        polar_lo = 1 + sqrt(T(3)) + sqrt(3 + 2*sqrt(T(3)))
        polar_p = a42(p -> _sep_poly_polar(p, a, e), polar_lo, T(8))
        x == 0 && return polar_p
        if x > 0
            equat_p = a42(p -> _sep_poly_equat(p, a, e), 1 + e, 6 + 2e)
            return a42(p -> _sep_poly_full(p, a, e, x), equat_p, polar_p)
        else
            return a42(p -> _sep_poly_full(p, a, e, x), polar_p, T(12))
        end
    end
end

"""
    get_separatrix(a, e, x) -> p_sep

FEW-convention alias of [`separatrix`](@ref) (`few.utils.utility.get_separatrix`).
"""
@inline get_separatrix(a, e, x) = separatrix(a, e, x)

# ── ISSO: the innermost STABLE SPHERICAL orbit (the e=0, inclined generalisation of the ISCO) ─────
# Promoted from the arXiv:2506.20726 figure scripts (`repro_rISSOg.jl` / `repro_drISSO.jl`), where
# the root-finder was figure-local. A spherical (r = const) geodesic of radius `r` and inclination
# `x = cos I` carries the constants (E, Lz, Q) = constants_of_motion(a, r, 0, x); it is radially
# STABLE while the radial potential obeys R''(r) < 0 and marginally stable (the ISSO) at R''(r) = 0.
# With A = (Lz − aE)² + Q the potential expands to
#   R(r) = (E² − 1) r⁴ + 2 r³ + (2a²E² − 2aELz − A − a²) r² + 2A r + a²(a²E² − 2aELz + Lz² − A),
# so R''(r) = 12(E² − 1) r² + 12 r + 2(2a²E² − 2aELz − A − a²) is CLOSED FORM (the figure scripts
# used a h=1e-4 finite difference; the analytic form removes that O(h²) bias). |x| = 1 reduces to
# the Bardeen ISCO and a = 0 gives 6 for every inclination — both asserted in the test suite.

"""
    isso_radial_curvature(a, r, x) -> R''(r)

Second derivative of the Kerr radial potential `R` at the spherical-orbit radius `r`, evaluated
with the constants of motion `(E, Lz, Q)` of the spherical orbit `(a, r, e=0, x)` held fixed
(M = 1, `x = cos I`). Negative ⇒ the spherical orbit is radially stable; the zero is the ISSO.
"""
function isso_radial_curvature(a, r, x)
    En, Lz, Q = constants_of_motion(a, r, zero(r), x)
    A = (Lz - a * En)^2 + Q
    return 12 * (En^2 - 1) * r^2 + 12 * r + 2 * (2 * a^2 * En^2 - 2 * a * En * Lz - A - a^2)
end

"""
    kerr_isso(a, x; rmax = 15) -> r_isso

Innermost stable SPHERICAL orbit radius (M = 1) for spin `a` and inclination `x = cos I`: the
radius where [`isso_radial_curvature`](@ref) changes sign (`R'' > 0` inside, `R'' < 0` outside).
The equatorial limits reduce to the Bardeen ISCO, `kerr_isso(a, ±1) == kerr_isco(±a)`, and
`kerr_isso(0, x) == 6` for every inclination. Bracketed bisection from just outside the outer
horizon to `rmax`, guarded against the inner region where no bound spherical orbit exists
(a failure of the constants there is treated as "unstable", i.e. inside the ISSO).

Canonical home of the ISSO root-finder that was figure-local in the arXiv:2506.20726 scripts.
"""
function kerr_isso(a::Real, x::Real; rmax::Real = 15)
    T = float(promote_type(typeof(a), typeof(x), typeof(rmax)))
    a, x = T(a), T(x)
    lo = 1 + sqrt(1 - a^2) + T(3) / 10        # just outside the outer horizon (R'' > 0 there)
    hi = T(rmax)
    curv(r) = try
        T(isso_radial_curvature(a, r, x))
    catch
        one(T)                                 # no bound spherical orbit ⇒ treat as unstable
    end
    curv(hi) < 0 || throw(ArgumentError("kerr_isso: rmax=$(rmax) is not outside the ISSO (a=$a, x=$x)"))
    for _ in 1:200
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        curv(mid) < 0 ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end
