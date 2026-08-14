"""
    EMRIWaveforms.Mappings

Orbital-parameter coordinate maps for the EMRI pipeline (design §6.8: FEW
`mappings/{kerrecceq,pn,jacobian,schwarzecc}.jl`, `geodesic_elq.jl`; KerrBackwards
`inclination.jl`). Reuses the `NumericsCore` conventions (`p_to_y`) and `Spacetime` geodesics
(`get_separatrix`, `constants_of_motion`, `geodesic_energy/angular_momentum`); the 4 dead FEW imports
(Teukolsky / TeukolskyPN / TeukolskySpinFluxes / SpinningOrbit) and the dead FEW→TeukolskyPN edge are
dropped → Symbolics-clean.

Surface:
- `pn_energy` / `pn_angular_momentum` / `pn_carter` — PN constants of motion (FEW `_PN_E/_PN_L/_PN_C`).
- `xI_to_Y` — exact `x = cos I → Y = cos ι = Lz/√(Lz²+Q)` over `Spacetime.constants_of_motion`.
- `Y_to_xI` — numeric inverse of `xI_to_Y` (KB `inclination.jl`, Brent); `pn_Y_to_xI` is the FEW PN fork.
- `kerrecceq_forward_map` / `kerrecceq_backward_map` — apex `(a,p,e,xI)` ↔ grid `(u,w,y,z)`.
- `ELdot_to_PEdot_Jacobian` — equatorial `(Ė,L̇z) → (ṗ,ė)` flux Jacobian.
- `p_to_y` — re-exported from `NumericsCore.conventions` (no new body).
"""
module Mappings

import NumericsCore: p_to_y
using Spacetime: get_separatrix, constants_of_motion, geodesic_energy, geodesic_angular_momentum
using Roots: find_zero, Brent

export pn_energy, pn_angular_momentum, pn_carter,
       xI_to_Y, Y_to_xI, pn_Y_to_xI,
       kerrecceq_forward_map, kerrecceq_backward_map,
       ELdot_to_PEdot_Jacobian, p_to_y

# =====================================================================================
# PN constants-of-motion kernels (FEW few/utils/mappings/pn.jl) — pure arithmetic, no deps.
# By construction pn_carter = L²(1/Y² − 1), so Y = L/√(L²+Q) holds exactly for the PN constants.
# =====================================================================================

@inline function _check_pn_params(q, p, e, Y)
    (abs(q) >= 1.0 || abs(e) >= 1.0 || p <= 1.0 || abs(Y) > 1.0) &&
        throw(ArgumentError("PN mapping parameters out of range"))
    return nothing
end

function _d1(q, p, e, Y)
    _check_pn_params(q, p, e, Y)
    p2 = p * p; q2 = q * q
    inv_ep4 = 1.0 / ((1.0 + e)^4)
    e == 0.0 ? (p2 + q2 - 2.0 * p) * p2 :
        p2 * ((q2 * e * e) + (2.0 * q2 - 2.0 * p) * e + p2 + q2 - 2.0 * p) * inv_ep4
end

function _f1(q, p, e, Y)
    _check_pn_params(q, p, e, Y)
    p2 = p * p; q2 = q * q
    inv_ep4 = 1.0 / ((1.0 + e)^4)
    e == 0.0 ? p2 * p2 + p2 * q2 + 2.0 * p * q2 :
        2.0 * p * ((1.0 + e)^2 * (e + p / 2.0 + 1.0) * q2 + (p2 * p) / 2.0) * inv_ep4
end

function _g1(q, p, e, Y)
    _check_pn_params(q, p, e, Y)
    e == 0.0 ? 2.0 * q * p : 2.0 * q * p / (1.0 + e)
end

function _h1(q, p, e, Y)
    _check_pn_params(q, p, e, Y)
    q2 = q * q; p2 = p * p; Y2 = Y * Y; y2 = 1.0 - Y2
    inv_ep2 = 1.0 / ((1.0 + e)^2)
    if e == 0.0 && Y == 1.0
        return 2.0 * q * p
    elseif e == 0.0
        return (p2 - 2.0 * p + y2 * q2) / Y2
    elseif Y == 1.0
        return (-2.0 * p * e + p2 - 2.0 * p) * inv_ep2
    else
        return (y2 * q2 * e * e + (2.0 * y2 * q2 - 2.0 * p) * e + y2 * q2 + p2 - 2.0 * p) * inv_ep2 / Y2
    end
end

function _d2(q, p, e, Y)
    _check_pn_params(q, p, e, Y)
    q2 = q * q; p2 = p * p
    inv_em4 = 1.0 / ((-1.0 + e)^4)
    e == 0.0 ? 4.0 * p2 * p + 2.0 * p * q2 - 6.0 * p2 :
        p2 * ((q2 * e * e) + (-2.0 * q2 + 2.0 * p) * e + p2 + q2 - 2.0 * p) * inv_em4
end

function _f2(q, p, e, Y)
    _check_pn_params(q, p, e, Y)
    q2 = q * q; p2 = p * p
    inv_em4 = 1.0 / ((-1.0 + e)^4)
    e == 0.0 ? (2.0 * p + 2.0) * q2 + 4.0 * p2 * p :
        -2.0 * ((e - p / 2.0 - 1.0) * (-1.0 + e)^2 * q2 - (p2 * p) / 2.0) * p * inv_em4
end

function _g2(q, p, e, Y)
    _check_pn_params(q, p, e, Y)
    e == 0.0 ? 2.0 * q : -2.0 * q * p / (-1.0 + e)
end

function _h2(q, p, e, Y)
    _check_pn_params(q, p, e, Y)
    q2 = q * q; p2 = p * p; Y2 = Y * Y; y2 = 1.0 - Y2
    inv_em2 = 1.0 / ((-1.0 + e)^2)
    if e == 0.0 && Y == 1.0
        return 2.0 * p - 2.0
    elseif e == 0.0
        return (2.0 * p - 2.0) / Y2
    elseif Y == 1.0
        return (2.0 * p * e + p2 - 2.0 * p) * inv_em2
    else
        return (y2 * q2 * e * e + (-2.0 * y2 * q2 + 2.0 * p) * e + y2 * q2 + p2 - 2.0 * p) * inv_em2 / Y2
    end
end

"`pn_energy(q,p,e,Y)` — PN orbital energy `E` (FEW `_PN_E`)."
function pn_energy(q, p, e, Y)
    d1, f1, g1, h1 = _d1(q, p, e, Y), _f1(q, p, e, Y), _g1(q, p, e, Y), _h1(q, p, e, Y)
    d2, f2, g2, h2 = _d2(q, p, e, Y), _f2(q, p, e, Y), _g2(q, p, e, Y), _h2(q, p, e, Y)
    kappa = d1 * h2 - d2 * h1
    epsilon = d1 * g2 - d2 * g1
    rho = f1 * h2 - f2 * h1
    eta = f1 * g2 - f2 * g1
    sigma = g1 * h2 - g2 * h1
    rhs_sqrt = sigma * (sigma * epsilon * epsilon + rho * epsilon * kappa - eta * kappa * kappa)
    rhs_numer = kappa * rho + 2.0 * epsilon * sigma
    rhs_denom = rho * rho + 4.0 * eta * sigma
    grade = xor(Y < 0, q < 0) ? -1 : 1   # retrograde if (Y<0) xor (q<0)
    E_square = (rhs_numer - grade * 2.0 * sqrt(rhs_sqrt)) / rhs_denom
    return sqrt(E_square)
end

"`pn_angular_momentum(q,p,e,Y)` — PN axial angular momentum `Lz` (FEW `_PN_L`)."
function pn_angular_momentum(q, p, e, Y)
    E = pn_energy(q, p, e, Y)
    h2 = _h2(q, p, e, Y)
    rhs_first_term = -_g2(q, p, e, Y) * E / h2
    rhs_sqrt = rhs_first_term^2 + (_f2(q, p, e, Y) * E * E - _d2(q, p, e, Y)) / h2
    sgnY = Y > 0 ? 1 : -1
    return rhs_first_term + sgnY * sqrt(rhs_sqrt)
end

"`pn_carter(q,p,e,Y)` — PN Carter constant `Q = L²(1/Y² − 1)` (FEW `_PN_C`); `0` at `Y = 1`."
function pn_carter(q, p, e, Y)
    L = pn_angular_momentum(q, p, e, Y)
    Y == 1.0 ? 0.0 : L * L * ((1.0 / (Y * Y)) - 1.0)
end

# FEW PN-based Y→xI inversion (few/utils/mappings/pn._Y_to_xI_kernel_inner).
"`pn_Y_to_xI(a,p,e,Y)` — FEW PN closed-form inclination inverse (fork of [`Y_to_xI`])."
function pn_Y_to_xI(a, p, e, Y)
    (abs(Y) == 1 || Y == 0.0) && return Y
    E = pn_energy(a, p, e, Y); Lz = pn_angular_momentum(a, p, e, Y); Q = pn_carter(a, p, e, Y)
    a2 = a * a; Lz2 = Lz * Lz
    E2m1 = (E - 1) * (E + 1.0)
    QpLz2ma2omE2 = Q + Lz2 + a2 * E2m1
    denomsqr = QpLz2ma2omE2 + sqrt(QpLz2ma2omE2 * QpLz2ma2omE2 - 4.0 * Lz2 * a2 * E2m1)
    return sqrt(2.0) * Lz / sqrt(denomsqr)
end

# =====================================================================================
# Inclination maps x_I = cos I ↔ Y = cos ι (KerrBackwards inclination.jl).
# =====================================================================================

"""
    xI_to_Y(a, p, e, x) -> Y

Exact `x = cos I → Y = cos ι = Lz/√(Lz²+Q)` over `Spacetime.constants_of_motion`
(machine precision; `few.utils.utility.xI_to_Y`).
"""
function xI_to_Y(a, p, e, x)
    c = constants_of_motion(a, p, e, x)
    return c.Lz / sqrt(c.Lz^2 + c.Q)
end

const _Y_TO_XI_LIM = 0.998

"""
    Y_to_xI(a, p, e, Y) -> x

Numeric inverse of [`xI_to_Y`]: `Y ↦ x` is monotonic with `x = ±1 ⇒ Y = ±1`. For
`|Y| > 0.998` the orbit is treated as equatorial (`x = Y`); otherwise a bracketed Brent solve of
`xI_to_Y(a,p,e,x) − Y` over `[Y−0.15, Y+0.15]` clamped to `±0.998` (KB `inclination.jl`, `YLIM=0.998`).
"""
function Y_to_xI(a, p, e, Y)
    -1 ≤ Y ≤ 1 || throw(DomainError(Y, "Y_to_xI: Y = cos ι must lie in [-1, 1]"))
    abs(Y) > _Y_TO_XI_LIM && return float(Y)
    x_lo = max(Y - 0.15, -_Y_TO_XI_LIM)
    x_hi = min(Y + 0.15, _Y_TO_XI_LIM)
    return find_zero(x -> xI_to_Y(a, p, e, x) - Y, (x_lo, x_hi), Brent())
end

# =====================================================================================
# KerrEccEq apex ↔ interpolation-coordinate maps (FEW few/utils/mappings/kerrecceq.jl).
# Two charts: region A (p ≤ pLSO+DELTAPMAX) and region B. Exact inverses on each chart interior.
# =====================================================================================

const XMIN = 0.05
const AMAX = 0.999
const AMIN = -AMAX
const DELTAPMIN = 0.001
const DELTAPMAX = 9 + DELTAPMIN
const EMAX = 0.9
const ESEP = 0.25
const ALPHA_FLUX = 1.0 / 2.0
const BETA_FLUX = 2.0
const ALPHA_AMP = 1.0 / 3.0
const BETA_AMP = 3.0
const DPC_REGIONB = DELTAPMAX - 0.001
const PMAX_REGIONB = 200
const AMAX_REGIONB = 0.999
const AMIN_REGIONB = -AMAX_REGIONB
const EMAX_REGIONB = 0.9
const DELTAPMIN_REGIONB = 9

_asign(a) = a == 0 ? one(float(a)) : sign(a)

# Region A scalar kernels
u_of_p(p, pLSO, alpha) = begin
    check_term = (log(p - pLSO + DELTAPMAX - 2 * DELTAPMIN) - log(DELTAPMAX - DELTAPMIN)) / log(2)
    sgn = sign(check_term)
    sgn * (sgn * check_term)^alpha
end
p_of_u(u, pLSO, alpha) =
    (pLSO + DELTAPMIN) + (DELTAPMAX - DELTAPMIN) * (exp(u^(1 / alpha) * log(2)) - 1)
y_of_x(x) = (x - XMIN) / (1 - XMIN)
x_of_y(y) = y * (1 - XMIN) + XMIN
chi_of_a(a) = (1 - a)^(1 / 3)
z_of_a(a) = (chi_of_a(a) - chi_of_a(AMAX)) / (chi_of_a(AMIN) - chi_of_a(AMAX))
a_of_chi(chi) = 1 - chi^3
a_of_z(z) = a_of_chi(chi_of_a(AMAX) + z * (chi_of_a(AMIN) - chi_of_a(AMAX)))
Secc_of_uz(u, z, beta) = begin
    check_part = z + u^beta * (1 - z)
    sgn = sign(check_part)
    ESEP + (EMAX - ESEP) * sgn * sqrt(sgn * check_part)
end
w_of_euz(e, u, z, beta) = e / Secc_of_uz(u, z, beta)
e_of_uwz(u, w, z, beta) = Secc_of_uz(u, z, beta) * w

function _uwyz_of_apex_kernel(a, p, e, x, pLSO, alpha, beta)
    u = u_of_p(p, pLSO, alpha); y = y_of_x(x); z = z_of_a(a); w = w_of_euz(e, u, z, beta)
    return u, w, y, z
end

function apex_of_uwyz(u, w, y, z; alpha = ALPHA_FLUX, beta = BETA_FLUX)
    a = a_of_z(z); x = x_of_y(y); e = e_of_uwz(u, w, z, beta)
    pLSO = get_separatrix(abs(a), e, _asign(a))
    p = p_of_u(u, pLSO, alpha)
    return a, p, e, x
end

# Region B scalar kernels
U_of_p_flux(p, pLSO) =
    ((DELTAPMIN_REGIONB)^(-0.5) - (p - pLSO)^(-0.5)) /
    ((DELTAPMIN_REGIONB)^(-0.5) - (PMAX_REGIONB - pLSO)^(-0.5))
p_of_U_flux(U, pLSO) =
    ((DELTAPMIN_REGIONB)^(-0.5) - U * ((DELTAPMIN_REGIONB)^(-0.5) - (PMAX_REGIONB - pLSO)^(-0.5)))^(-2) + pLSO
U_of_p_amplitude(p, pLSO) = begin
    pc = pLSO + DPC_REGIONB
    (pc^(-0.5) - p^(-0.5)) / (pc^(-0.5) - (PMAX_REGIONB + pc)^(-0.5))
end
p_of_U_amplitude(U, pLSO) = begin
    pc = pLSO + DPC_REGIONB
    (pc^(-0.5) - U * (pc^(-0.5) - (PMAX_REGIONB + pc)^(-0.5)))^(-2)
end
W_of_e(e) = e / EMAX_REGIONB
e_of_W(w) = w * EMAX_REGIONB
Z_of_a(a) = (chi_of_a(a) - chi_of_a(AMAX_REGIONB)) / (chi_of_a(AMIN_REGIONB) - chi_of_a(AMAX_REGIONB))
a_of_Z(z) = a_of_chi(chi_of_a(AMAX_REGIONB) + z * (chi_of_a(AMIN_REGIONB) - chi_of_a(AMAX_REGIONB)))

_UWYZ_of_apex_kernel(a, p, e, x, pLSO, is_flux) =
    (is_flux ? U_of_p_flux(p, pLSO) : U_of_p_amplitude(p, pLSO), W_of_e(e), y_of_x(x), Z_of_a(a))

function apex_of_UWYZ(u, w, y, z, is_flux)
    a = a_of_Z(z); x = x_of_y(y); e = e_of_W(w)
    pLSO = get_separatrix(abs(a), e, _asign(a))
    p = is_flux ? p_of_U_flux(u, pLSO) : p_of_U_amplitude(u, pLSO)
    return a, p, e, x
end

function _alpha_beta_isflux(kind)
    kind == "flux" && return (ALPHA_FLUX, BETA_FLUX, true)
    kind == "amplitude" && return (ALPHA_AMP, BETA_AMP, false)
    throw(ArgumentError("kind must be \"flux\" or \"amplitude\""))
end

"""
    kerrecceq_forward_map(a, p, e, xI; pLSO=nothing, return_mask=false, kind="flux")

Apex `(a, p, e, xI)` → interpolation coordinates `(u, w, y, z)` (FEW scalar path). Only `xI = 1`.
`return_mask=true` appends `near` (region A flag).
"""
function kerrecceq_forward_map(a, p, e, xI; pLSO = nothing, return_mask = false, kind = "flux")
    xI == 1 || throw(ArgumentError("Only xI = 1 is supported."))
    alpha, beta, is_flux = _alpha_beta_isflux(kind)
    pl = pLSO === nothing ? get_separatrix(abs(a), e, _asign(a) * xI) : pLSO
    near = p <= pl + DELTAPMAX
    u, w, y, z = near ? _uwyz_of_apex_kernel(a, p, e, xI, pl, alpha, beta) :
                        _UWYZ_of_apex_kernel(a, p, e, xI, pl, is_flux)
    return return_mask ? (u, w, y, z, near) : (u, w, y, z)
end

"""
    kerrecceq_backward_map(u, w, y, z; regionA=true, kind="flux")

Interpolation coordinates `(u, w, y, z)` → apex `(a, p, e, xI)` (FEW scalar path). Only `y = 1`.
"""
function kerrecceq_backward_map(u, w, y, z; regionA = true, kind = "flux")
    y == 1 || throw(ArgumentError("Only xI = 1 is supported."))
    alpha, beta, is_flux = _alpha_beta_isflux(kind)
    return regionA ? apex_of_uwyz(u, w, y, z; alpha = alpha, beta = beta) :
                     apex_of_UWYZ(u, w, y, z, is_flux)
end

# =====================================================================================
# Equatorial flux Jacobian (E,Lz) → (p,ecc) (FEW few/utils/mappings/jacobian.jl).
# =====================================================================================

"`_schwarz_jac_kernel(p,e,Edot,Ldot)` — closed-form Schwarzschild `(ṗ,ė)`."
function _schwarz_jac_kernel(p, e, Edot, Ldot)
    pdot = (-2 * (Edot * sqrt((4 * e^2 - (-2 + p)^2) / (3 + e^2 - p)) * (3 + e^2 - p) * p^1.5 +
            Ldot * (-4 + p)^2 * sqrt(-3 - e^2 + p))) / (4 * e^2 - (-6 + p)^2)
    if e > 0
        edot = -((Edot * sqrt((4 * e^2 - (-2 + p)^2) / (3 + e^2 - p)) * p^1.5 *
                  (18 + 2 * e^4 - 3 * e^2 * (-4 + p) - 9 * p + p^2) +
                  (-1 + e^2) * Ldot * sqrt(-3 - e^2 + p) * (12 + 4 * e^2 - 8 * p + p^2)) /
                 (e * (4 * e^2 - (-6 + p)^2) * p))
    else
        edot = 0.0
    end
    return pdot, edot
end

"""
    ELdot_to_PEdot_Jacobian(a, p, ecc, xI, Edot, Lzdot) -> (pdot, eccdot)

Transform equatorial fluxes `(Ė, L̇z) → (ṗ, ė)`; `|xI| = 1`. `a == 0` uses the closed-form
Schwarzschild kernel, otherwise the Kerr constants of motion (FEW `ELdot_to_PEdot_Jacobian`).
"""
function ELdot_to_PEdot_Jacobian(a, p, ecc, xI, Edot, Lzdot)
    a == 0 && return _schwarz_jac_kernel(p, ecc, Edot, Lzdot)
    E = geodesic_energy(a, p, ecc, xI)
    Lz = geodesic_angular_momentum(a, p, ecc, xI)
    Lz2pQ = Lz * Lz
    LzmaE = Lz - a * E
    omE2 = (1.0 - E) * (1.0 + E)
    if ecc == 0.0
        u2 = 1.0 / p
        denom = u2 * (-6 * omE2 + u2 * (6.0 + u2 * (-a * a * omE2 - Lz2pQ)))
        numerEdot = -4.0 * E + u2 * u2 * (-2.0 * a * a * E + u2 * (2.0 * a * LzmaE))
        numerLzdot = u2 * u2 * (2.0 * Lz + u2 * (-2 * LzmaE))
        pdot = (numerEdot * Edot + numerLzdot * Lzdot) / denom
        eccdot = 0.0
    else
        ra = p / (1.0 - ecc); rp = p / (1.0 + ecc)
        ua2 = 1.0 / ra; up2 = 1.0 / rp
        denoma = ua2 * (-4.0 * omE2 + ua2 * (6 + ua2 * (-2 * Lz2pQ - 2 * a * a * omE2 +
                 ua2 * (2.0 * Lz2pQ - 4.0 * a * E * Lz + 2.0 * a * a * E * E))))
        denomp = up2 * (-4.0 * omE2 + up2 * (6 + up2 * (-2 * Lz2pQ - 2 * a * a * omE2 +
                 up2 * (2.0 * Lz2pQ - 4.0 * a * E * Lz + 2.0 * a * a * E * E))))
        numeraEdot = -2.0 * (E + ua2 * ua2 * (a * a * E + ua2 * (-2.0 * a * LzmaE)))
        numeraLzdot = 2.0 * ua2 * ua2 * (Lz - 2.0 * ua2 * LzmaE)
        numerpEdot = -2.0 * (E + up2 * up2 * (a * a * E + up2 * (-2.0 * a * LzmaE)))
        numerpLzdot = 2.0 * up2 * up2 * (Lz - 2.0 * up2 * LzmaE)
        radot = (numeraEdot * Edot + numeraLzdot * Lzdot) / denoma
        rpdot = (numerpEdot * Edot + numerpLzdot * Lzdot) / denomp
        pdot = 0.5 * (1 + ecc * (-2.0 + ecc)) * radot + 0.5 * (1.0 + ecc * (2.0 + ecc)) * rpdot
        eccdot = (1.0 + ecc) * ((1.0 + ecc * (-2.0 + ecc)) * radot + (ecc * ecc - 1.0) * rpdot) / (2.0 * p)
        abs(eccdot) < 3.0e-14 && (eccdot = 0.0)
    end
    return pdot, eccdot
end

end # module Mappings
