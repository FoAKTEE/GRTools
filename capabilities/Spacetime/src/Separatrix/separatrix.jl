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
