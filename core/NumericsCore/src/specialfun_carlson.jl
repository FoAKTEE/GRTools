# Carlson symmetric elliptic integrals + complete/incomplete Legendre K/E/Pi + Jacobi am/sn/cn/dn.
#
# Canonical, type-generic (real|complex) consolidation of the two legacy forks: KerrGeodesics'
# complex Mathematica-PARAMETER (m=k^2) Carlson layer and FastEMRIWaveforms' real modulus wrappers
# are the SAME decision (the FEW layer is literally `real(KerrGeodesics._R*(complex(args)))`), so
# they collapse to one primitive set whose return type FOLLOWS the input type (real in -> real out,
# complex in -> complex out). Carlson 1995 / Numerical Recipes 6.12 duplication-theorem algorithm.
#
# Complete integrals (DLMF 19.25.1-2, alpha^2=n, k^2=m):
#   ellipK(m)    = R_F(0, 1-m, 1)
#   ellipE(m)    = R_F(0, 1-m, 1) - (m/3) R_D(0, 1-m, 1)
#   ellipPi(n,m) = R_F(0, 1-m, 1) + (n/3) R_J(0, 1-m, 1, 1-n)
# returned on Mathematica's principal branch for real m>1 (m approached as m-i0+ via _below_cut).

# 5th-order series => truncation error ~ errtol^6. Preserve the legacy Float64 threshold, but
# tighten it when the active type has enough precision that 1e-6 would become the accuracy floor.
const _ELL_ERRTOL = 1.0e-6
@inline _ell_errtol(::Type{T}) where {T<:AbstractFloat} =
    min(T(_ELL_ERRTOL), eps(T)^(one(T) / 6))

# Hard iteration cap for the duplication loops. Legitimate inputs reach _ELL_ERRTOL in <~20 iterations
# (each step shrinks the argument spread by ~1/4). The cap fires ONLY for the mathematically-divergent
# degenerate cases — R_F/R_D with >=2 vanishing arguments (e.g. R_F(0,0,z)=+inf, the m->1 limit of
# ellipK/ellipE), whose duplication residual `(ave-xt)/ave` is pinned at 1 and never drops below the
# tolerance. Capping makes such a stray call return (Inf) instead of looping forever; the genuine
# special values (e.g. E(1)=1) are handled analytically by the callers below.
const _ELL_MAXITER = 200

"""`_below_cut(v)` — signed-zero promotion of a real `1-m`/`1-n` argument onto Mathematica's
principal branch: real m>1 is the limit m-i0+, so `1-m` approaches from +i0+ and `sqrt(neg)`
takes the +i side."""
@inline _below_cut(v::Real) = complex(float(v), zero(float(v)))

function _rc(x::Complex{T}, y::Complex{T}) where {T<:AbstractFloat}
    xt, yt = x, y
    ave = zero(x); s = zero(x)
    errtol = _ell_errtol(T)
    for _ in 1:_ELL_MAXITER
        lam = 2*sqrt(xt)*sqrt(yt) + yt
        xt = (xt + lam)/4; yt = (yt + lam)/4
        ave = (xt + yt + yt)/3
        s = (yt - ave)/ave
        abs(s) < errtol && break
    end
    return (1 + s*s*(3//10 + s*(1//7 + s*(3//8 + s*(9//22)))))/sqrt(ave)
end

function _rf(x::Complex{T}, y::Complex{T}, z::Complex{T}) where {T<:AbstractFloat}
    xt, yt, zt = x, y, z
    ave = zero(x); dx = dy = dz = zero(x)
    errtol = _ell_errtol(T)
    for _ in 1:_ELL_MAXITER
        sqx, sqy, sqz = sqrt(xt), sqrt(yt), sqrt(zt)
        lam = sqx*(sqy + sqz) + sqy*sqz
        xt = (xt + lam)/4; yt = (yt + lam)/4; zt = (zt + lam)/4
        ave = (xt + yt + zt)/3
        dx = (ave - xt)/ave; dy = (ave - yt)/ave; dz = (ave - zt)/ave
        max(abs(dx), abs(dy), abs(dz)) < errtol && break
    end
    e2 = dx*dy - dz*dz
    e3 = dx*dy*dz
    return (1 + (e2/24 - 1//10 - 3*e3/44)*e2 + e3/14)/sqrt(ave)
end

function _rd(x::Complex{T}, y::Complex{T}, z::Complex{T}) where {T<:AbstractFloat}
    xt, yt, zt = x, y, z
    summ = zero(x); fac = one(T); ave = zero(x); dx = dy = dz = zero(x)
    errtol = _ell_errtol(T)
    for _ in 1:_ELL_MAXITER
        sqx, sqy, sqz = sqrt(xt), sqrt(yt), sqrt(zt)
        lam = sqx*(sqy + sqz) + sqy*sqz
        summ += fac/(sqz*(zt + lam))
        fac /= 4
        xt = (xt + lam)/4; yt = (yt + lam)/4; zt = (zt + lam)/4
        ave = (xt + yt + 3*zt)/5
        dx = (ave - xt)/ave; dy = (ave - yt)/ave; dz = (ave - zt)/ave
        max(abs(dx), abs(dy), abs(dz)) < errtol && break
    end
    ea = dx*dy; eb = dz*dz; ec = ea - eb; ed = ea - 6*eb; ee = ed + ec + ec
    C1, C2, C3, C4 = 3//14, 1//6, 9//22, 3//26
    C5 = C3/4; C6 = 3*C4/2
    return 3*summ + fac*(1 + ed*(-C1 + C5*ed - C6*dz*ee) +
                         dz*(C2*ee + dz*(-C3*ec + dz*C4*ea)))/(ave*sqrt(ave))
end

function _rj(x::Complex{T}, y::Complex{T}, z::Complex{T}, p::Complex{T}) where {T<:AbstractFloat}
    xt, yt, zt, pt = x, y, z, p
    summ = zero(x); fac = one(T); ave = zero(x); dx = dy = dz = dp = zero(x)
    errtol = _ell_errtol(T)
    for _ in 1:_ELL_MAXITER
        sqx, sqy, sqz = sqrt(xt), sqrt(yt), sqrt(zt)
        lam = sqx*(sqy + sqz) + sqy*sqz
        alpha = (pt*(sqx + sqy + sqz) + sqx*sqy*sqz)^2
        beta  = pt*(pt + lam)^2
        summ += fac*_rc(alpha, beta)
        fac /= 4
        xt = (xt + lam)/4; yt = (yt + lam)/4; zt = (zt + lam)/4; pt = (pt + lam)/4
        ave = (xt + yt + zt + 2*pt)/5
        dx = (ave - xt)/ave; dy = (ave - yt)/ave; dz = (ave - zt)/ave; dp = (ave - pt)/ave
        max(abs(dx), abs(dy), abs(dz), abs(dp)) < errtol && break
    end
    ea = dx*(dy + dz) + dy*dz
    eb = dx*dy*dz
    ec = dp*dp
    ed = ea - 3*ec
    ee = eb + 2*dp*(ea - ec)
    C1, C2, C3, C4 = 3//14, 1//3, 3//22, 3//26
    C5 = 3*C3/4; C6 = 3*C4/2; C7 = C2/2; C8 = C3 + C3
    return 3*summ + fac*(1 + ed*(-C1 + C5*ed - C6*ee) + eb*(C7 + dp*(-C8 + dp*C4)) +
                         dp*ea*(C2 - dp*C3) - C2*dp*ec)/(ave*sqrt(ave))
end

# --- public type-following Carlson symmetric forms ---

@inline _anycomplex(args...) = any(a -> a isa Complex, args)
@inline _ctype(args...) = float(real(promote_type(map(typeof, args)...)))

"""`RC(x, y)` — degenerate Carlson `R_C(x,y)=R_F(x,y,y)`. Real in -> real out; complex preserved."""
function RC(x::Number, y::Number)
    T = _ctype(x, y); v = _rc(Complex{T}(x), Complex{T}(y))
    return _anycomplex(x, y) ? v : real(v)
end

"""`RF(x, y, z)` — Carlson symmetric integral of the 1st kind `R_F`."""
function RF(x::Number, y::Number, z::Number)
    T = _ctype(x, y, z); v = _rf(Complex{T}(x), Complex{T}(y), Complex{T}(z))
    return _anycomplex(x, y, z) ? v : real(v)
end

"""`RD(x, y, z)` — degenerate Carlson `R_D(x,y,z)=R_J(x,y,z,z)` (2nd kind)."""
function RD(x::Number, y::Number, z::Number)
    T = _ctype(x, y, z); v = _rd(Complex{T}(x), Complex{T}(y), Complex{T}(z))
    return _anycomplex(x, y, z) ? v : real(v)
end

"""`RJ(x, y, z, p)` — Carlson symmetric integral of the 3rd kind `R_J`."""
function RJ(x::Number, y::Number, z::Number, p::Number)
    T = _ctype(x, y, z, p); v = _rj(Complex{T}(x), Complex{T}(y), Complex{T}(z), Complex{T}(p))
    return _anycomplex(x, y, z, p) ? v : real(v)
end

# --- complete Legendre integrals (Mathematica parameter m, principal branch, complex-valued) ---

"""`ellipK(m)` — complete elliptic integral of the 1st kind `EllipticK[m]` (m=k^2).
The endpoint `m=1` is the logarithmic divergence `+Inf`."""
function ellipK(m::Real)
    mf = float(m)
    isone(mf) && return complex(oftype(mf, Inf), zero(mf))
    return _rf(complex(zero(mf)), _below_cut(1 - mf), complex(one(mf)))
end

"""`ellipE(m)` — complete elliptic integral of the 2nd kind `EllipticE[m]`. At `m=1` the Carlson form
`R_F(0,0,1) - (1/3) R_D(0,0,1)` is an inf - inf limit, so the analytic value `E(1)=1` is returned."""
function ellipE(m::Real)
    mf = float(m)
    isone(mf) && return complex(one(mf), zero(mf))   # E(1)=1 (R_F(0,0,1)/R_D(0,0,1) both diverge)
    z0 = complex(zero(mf)); z1 = complex(one(mf)); mc = _below_cut(1 - mf)
    return _rf(z0, mc, z1) - (mf/3)*_rd(z0, mc, z1)
end

"""`ellipPi(n, m)` — complete elliptic integral of the 3rd kind `EllipticPi[n,m]`."""
function ellipPi(n::Real, m::Real)
    z0 = complex(zero(float(m))); z1 = complex(one(float(m))); mc = _below_cut(1 - m)
    return _rf(z0, mc, z1) + (n/3)*_rj(z0, mc, z1, _below_cut(1 - n))
end

"""`ellipPi(n::Complex, m)` — `EllipticPi[n,m]` for complex characteristic n; `1-n` is off the real
axis so no branch-cut nudge is needed."""
function ellipPi(n::Complex, m::Real)
    z0 = complex(zero(float(m))); z1 = complex(one(float(m))); mc = _below_cut(1 - m)
    return _rf(z0, mc, z1) + (n/3)*_rj(z0, mc, z1, oneunit(z0) - n)
end

# --- incomplete Legendre integrals (amplitude phi, parameter m), real regime (DLMF 19.25.5-7) ---

"""`_phi_reduce(phi)` — split `phi` into `(j, phi - j*pi)` for the period reduction of the
incomplete integrals (continuous, odd in phi)."""
@inline _phi_reduce(phi::Real) = (j = round(phi / pi); (j, phi - j*pi))

"""`ellipFinc(phi, m)` — incomplete elliptic integral of the 1st kind `EllipticF[phi,m]`."""
function ellipFinc(phi::Real, m::Real)
    j, ph = _phi_reduce(phi)
    s = sin(ph); c = cos(ph)^2
    val = s*_rf(complex(c), complex(1 - m*s^2), complex(one(c)))
    return real(val) + 2*j*real(ellipK(m))
end

"""`ellipEinc(phi, m)` — incomplete elliptic integral of the 2nd kind `EllipticE[phi,m]`."""
function ellipEinc(phi::Real, m::Real)
    j, ph = _phi_reduce(phi)
    s = sin(ph); c = cos(ph)^2; y = complex(1 - m*s^2); o = complex(one(c))
    val = s*_rf(complex(c), y, o) - (m/3)*s^3*_rd(complex(c), y, o)
    return real(val) + 2*j*real(ellipE(m))
end

"""`ellipPiinc(n, phi, m)` — incomplete elliptic integral of the 3rd kind `EllipticPi[n,phi,m]`."""
function ellipPiinc(n::Real, phi::Real, m::Real)
    j, ph = _phi_reduce(phi)
    s = sin(ph); c = cos(ph)^2; y = complex(1 - m*s^2); o = complex(one(c))
    val = s*_rf(complex(c), y, o) + (n/3)*s^3*_rj(complex(c), y, o, complex(1 - n*s^2))
    return real(val) + 2*j*real(ellipPi(n, m))
end

# --- Jacobi amplitude / elliptic functions (real u, 0<=m<=1) via descending-Landen / AGM ---
#
# Abramowitz & Stegun 16.4, 17.6. Type-generic: the AGM scale sequence carries the element type T
# (= promote_type(u,m)), so BigFloat precision is inherited (legacy KerrGeodesics hard-Float64'd it).

"""`jacobiAM(u, m)` — Jacobi amplitude `am(u|m)` (real u, 0<=m<=1)."""
function jacobiAM(u::Real, m::Real)
    T = float(promote_type(typeof(u), typeof(m)))
    uu = T(u); mm = T(m)
    mm == 0 && return uu
    as = T[one(T)]; cs = T[sqrt(mm)]
    a = one(T); b = sqrt(1 - mm); c = sqrt(mm); n = 0
    tol = 4*eps(T)
    while abs(c) > tol && n < 80
        a, b, c = (a + b)/2, sqrt(a*b), (a - b)/2
        push!(as, a); push!(cs, c); n += 1
    end
    phi = (T(2)^n) * as[n+1] * uu
    for k in n:-1:1
        phi = (phi + asin(clamp((cs[k+1]/as[k+1])*sin(phi), -one(T), one(T))))/2
    end
    return phi
end

"""`jacobiSN(u, m)` — `sn(u|m)`."""
jacobiSN(u::Real, m::Real) = sin(jacobiAM(u, m))

"""`jacobiCN(u, m)` — `cn(u|m)`."""
jacobiCN(u::Real, m::Real) = cos(jacobiAM(u, m))

"""`jacobiDN(u, m)` — `dn(u|m)`."""
function jacobiDN(u::Real, m::Real)
    s = jacobiSN(u, m)
    return sqrt(1 - oftype(s, m)*s^2)
end
