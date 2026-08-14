# MST radial assembly: the raw (norm=1) homogeneous solutions R^in / R^up and their r-derivatives
# = prefactor × dual-test-truncated Σ_n f_n · {2F1,U}_n. Consumes mst.jl (f_n) + hypergeom.jl
# (the 2F1/U terms + HU/dHU recurrence arrays). norm=1 gives the raw MST solution (the public
# unit-transmission normalization is applied later by TeukolskyRadial). AHA: the f_Up Pochhammer
# rescale delegates to BigSpecial.pochhammer (signature-identical to the legacy local copy).

"In-solution prefactor (Sasaki–Tagoshi §2.5). x=(rp−r)/(2κ), x<0 for r>rp so −x>0, 1−x>0."
_prefac_in(s, eps, tau, kappa, x) =
    (-x)^(-s - im * (eps + tau) / 2) * (1 - x)^(im * (eps - tau) / 2) * exp(im * eps * kappa * x)

"Up-solution prefactor (§2.5). zhat>0 and zhat−εκ>0 for r>rp (real positive bases)."
_prefac_up(s, eps, kappa, tau, nu, zhat) =
    2.0^nu * exp(-pi * eps) * exp(-im * pi * (nu + 1)) * exp(im * zhat) *
    zhat^(nu + im * (eps + tau) / 2) * (zhat - eps * kappa)^(-im * (eps + tau) / 2) *
    exp(-im * pi * s) * (zhat - eps * kappa)^(-s)

"Dual-test truncated MST sum: term(n), upward from 0 and downward from −1 (§2.10)."
function _mst_sum(term, Nmax::Integer, atol::Real, rtol::Real)
    resUp = zero(complex(term(0))); n = 0     # generic (ComplexF64 or Complex{BigFloat})
    while n <= Nmax
        t = term(n); prev = resUp; resUp += t
        (resUp == prev || abs(t) <= atol + abs(resUp) * rtol) && break
        n += 1
    end
    resDown = zero(resUp); n = -1
    while n >= -Nmax
        t = term(n); prev = resDown; resDown += t
        (resDown == prev || abs(t) <= atol + abs(resDown) * rtol) && break
        n -= 1
    end
    return resUp + resDown
end

# f_Up rescale (§2.4): fUp[n] = (−1)ⁿ Poch[ν+1+s−iε,n]/Poch[ν+1−s+iε,n] · f_n.
_fup_factor(nu, s, eps, n) =
    (iseven(n) ? 1 : -1) * pochhammer(nu + 1 + s - im * eps, n) / pochhammer(nu + 1 - s + im * eps, n)

"""
    _mst_radial_in(s, m, a, omega, nu, lambda, r; norm=1.0, Nmax=80) -> Complex

MST "In" radial solution R^in(r) (§2.10): `prefacIn/norm · Σ_n f_n · ₂F₁[n](x)`, x=(rp−r)/(2κ).
"""
function _mst_radial_in(s::Integer, m::Integer, a::Real, omega::Real, nu, lambda, r::Real;
                        norm=1.0, Nmax::Integer=80, atol::Real=1e-16, rtol::Real=1e-16)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa
    rp = 1 + kappa; x = (rp - r) / (2 * kappa)
    farr, _, _, _ = _mst_fn_array(s, m, a, omega, nu, lambda, Nmax)
    prefac = _prefac_in(s, eps, tau, kappa, x) / norm
    term = n -> prefac * farr[n] * _h2f1_exact(n, s, nu, tau, eps, x)
    return _mst_sum(term, Nmax, atol, rtol)
end

"""
    _mst_radial_in_deriv(s, m, a, omega, nu, lambda, r; norm=1.0, Nmax=80) -> Complex

First r-derivative dR^in/dr (§2.11): `dxdr/norm · Σ_n f_n(dprefacIn·₂F₁[n] + prefacIn·d₂F₁[n])`,
dxdr=−1/(2κ). dprefacIn via the log-derivative; d₂F₁ via the exact Gauss shift.
"""
function _mst_radial_in_deriv(s::Integer, m::Integer, a::Real, omega::Real, nu, lambda, r::Real;
                              norm=1.0, Nmax::Integer=80, atol::Real=1e-16, rtol::Real=1e-16)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa
    rp = 1 + kappa; x = (rp - r) / (2 * kappa); dxdr = -1 / (2 * kappa)
    farr, _, _, _ = _mst_fn_array(s, m, a, omega, nu, lambda, Nmax)
    aF = _aF(s, nu, tau, eps); bF = _bF(s, nu, tau, eps); cF = _cF(s, nu, tau, eps)
    p = -s - im * (eps + tau) / 2; qx = im * (eps - tau) / 2; rr = im * eps * kappa
    prefacIn = (-x)^p * (1 - x)^qx * exp(rr * x)
    dprefacIn = prefacIn * (p / x - qx / (1 - x) + rr)        # ∂prefacIn/∂x (log-derivative)
    dh2f1 = n -> ((n + aF) * (bF - n) / cF) *
                 _safe_2f1(n + aF + 1, bF - n + 1, cF + 1, complex(x))   # BigFloat-robust
    term = n -> farr[n] * (dprefacIn * _h2f1_exact(n, s, nu, tau, eps, x) + prefacIn * dh2f1(n))
    return (dxdr / norm) * _mst_sum(term, Nmax, atol, rtol)
end

"""
    _mst_radial_up(s, m, a, omega, nu, lambda, r; norm=1.0, Nmax=80) -> Complex

MST "Up" radial solution R^up(r) (§2.12): `prefacUp/norm · Σ_n fUp[n] · U[n](zhat)`,
zhat = ε(r−rm)/2; HU[n] via the recurrence array (required for accuracy).
"""
function _mst_radial_up(s::Integer, m::Integer, a::Real, omega::Real, nu, lambda, r::Real;
                        norm=1.0, Nmax::Integer=80, atol::Real=1e-16, rtol::Real=1e-16)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa
    rm = 1 - kappa; zhat = eps * r / 2 - eps * rm / 2
    farr, _, _, _ = _mst_fn_array(s, m, a, omega, nu, lambda, Nmax)
    HU = _hu_array(s, nu, tau, eps, zhat, Nmax)          # recurrence — required for accuracy
    prefac = _prefac_up(s, eps, kappa, tau, nu, zhat) / norm
    term = n -> prefac * _fup_factor(nu, s, eps, n) * farr[n] * HU[n]
    return _mst_sum(term, Nmax, atol, rtol)
end

"""
    _mst_radial_up_deriv(s, m, a, omega, nu, lambda, r; norm=1.0, Nmax=80) -> Complex

First r-derivative dR^up/dr (§2.13): `dzhatdr/norm · Σ_n fUp[n](dprefacUp·HU[n] + prefacUp·dHU[n])`,
dzhatdr=ε/2. dprefacUp via the log-derivative; dHU via the §2.9 recurrence array.
"""
function _mst_radial_up_deriv(s::Integer, m::Integer, a::Real, omega::Real, nu, lambda, r::Real;
                              norm=1.0, Nmax::Integer=80, atol::Real=1e-16, rtol::Real=1e-16)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa
    rm = 1 - kappa; zhat = eps * r / 2 - eps * rm / 2; dzhatdr = eps / 2
    farr, _, _, _ = _mst_fn_array(s, m, a, omega, nu, lambda, Nmax)
    HU = _hu_array(s, nu, tau, eps, zhat, Nmax)
    dHU = _dhu_array(s, nu, tau, eps, zhat, HU, Nmax)
    prefacUp = _prefac_up(s, eps, kappa, tau, nu, zhat)
    dprefacUp = prefacUp * (im + (nu + im * (eps + tau) / 2) / zhat +
                            (-im * (eps + tau) / 2 - s) / (zhat - eps * kappa))     # ∂prefacUp/∂zhat
    term = n -> _fup_factor(nu, s, eps, n) * farr[n] * (dprefacUp * HU[n] + prefacUp * dHU[n])
    return (dzhatdr / norm) * _mst_sum(term, Nmax, atol, rtol)
end
