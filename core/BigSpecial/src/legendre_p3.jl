# Associated Legendre function of TYPE 3 — P_n^m on the cut z > 1 (DLMF §14.3(ii) /
# Mathematica LegendreP[n,m,3,z]) at general (half-integer) degree. This is the m-mode
# puncture-integral kernel of the CRP elliptic scheme (arXiv:2606.04998 main.tex:467-486):
# ∫ Rⁿ e^{-imΔφ} dΔφ ∝ ₃P_{n/2}^m(z), n odd, so HALF-INTEGER degree is the primary use case.
# Built on the existing BigSpecial primitives (f21_maclaurin, gamma_big, digamma_big) —
# type-generic throughout, so BigFloat precision is inherited.

"""
    legendre_P3(n, m, z) -> real

Associated Legendre function of type 3, `₃P_n^m(z)`, defined for real `z > 1` (branch cut
`(-∞, 1]`); degree `n` may be any real (half-integers are the CRP puncture's case), order `m`
any real (integer in the physics). Evaluated per arXiv:2606.04998 main.tex:480 as

    ₃P_n^m(z) = Γ(n+m+1)/Γ(n−m+1) · x^{m/2} · (2/(z+1))^{n+1} · ₂𝐅₁(n+m+1, n+1; m+1; x),

with `x = (z−1)/(z+1)` and `₂𝐅₁` the REGULARIZED Gauss function `₂F₁/Γ(c)` (valid for
`n ≥ 0, m ≥ 0`; here used for the mapped `n ≥ −1/2, m ≥ 0`, where all of `a,b,c > 0`).
Negative degree/order are routed through the exact connections (ibid.)

    ₃P_n^m = ₃P_{−n−1}^m  (applied for n < −1/2),
    ₃P_n^{−m} = Γ(n−m+1)/Γ(n+m+1) · ₃P_n^m  (m > 0),

applied in that order. Near the worldline (`x ≥ 0.995`, i.e. `z ≥ 399`) the regularized ₂F₁ is
evaluated by direct summation of DLMF 15.8.10 / 15.8.12 instead of the Maclaurin series — the
CRP's own switch; it requires `2n+1 ∈ ℤ` (any integer or half-odd-integer degree). For the
out-of-contract case `x ≥ 0.995` with `2n+1 ∉ ℤ` the Maclaurin route is kept (slow convergence,
still correct). Type-generic: BigFloat arguments give BigFloat-accurate results in BOTH branches
(the DLMF ψ/log sums use closed-form digamma at (half-)integer points, no Float64 cap).
"""
function legendre_P3(n::Real, m::Real, z::Real)
    T = float(promote_type(typeof(n), typeof(m), typeof(z)))
    nT = T(n); mT = T(m); zT = T(z)
    zT > 1 || throw(DomainError(z, "legendre_P3 is defined for z > 1"))
    nT < -T(1) / 2 && (nT = -nT - 1)               # ₃P_n^m = ₃P_{−n−1}^m
    if mT < 0                                      # ₃P_n^{−μ} = Γ(n−μ+1)/Γ(n+μ+1) ₃P_n^μ
        return gamma_big(nT + mT + 1) / gamma_big(nT - mT + 1) * _p3_positive(nT, -mT, zT)
    end
    return _p3_positive(nT, mT, zT)
end

"""
    legendre_P3_dz(n, m, z) -> real

Derivative `d ₃P_n^m/dz` via the exact recurrence (arXiv:2606.04998 main.tex:486)

    d₃P_n^m/dz = [(n−m+1) ₃P_{n+1}^m(z) − (n+1) z ₃P_n^m(z)] / (z²−1),

valid for all real `n, m` (the degree/order connections are applied inside `legendre_P3`).
"""
function legendre_P3_dz(n::Real, m::Real, z::Real)
    T = float(promote_type(typeof(n), typeof(m), typeof(z)))
    nT = T(n); mT = T(m); zT = T(z)
    return ((nT - mT + 1) * legendre_P3(nT + 1, mT, zT) - (nT + 1) * zT * legendre_P3(nT, mT, zT)) /
           (zT^2 - 1)
end

# ₃P_n^m for the mapped domain n ≥ −1/2, m ≥ 0. w = 2/(z+1) = 1−x is computed DIRECTLY (exact
# formula): forming 1−x from the rounded x would carry an absolute eps error into w, i.e. a
# relative error eps/w that the near-branch w^{−M} amplifies M-fold (≈1e-11 at z=1e4, M=10).
function _p3_positive(n::T, m::T, z::T) where {T<:Real}
    x = (z - 1) / (z + 1)
    w = 2 / (z + 1)
    pre = gamma_big(n + m + 1) / gamma_big(n - m + 1) * x^(m / 2) * w^(n + 1)
    a = n + m + 1; b = n + 1; c = m + 1
    Mr = 2n + 1                                    # c = a + b − M with M = 2n+1 ≥ 0 here
    if x >= T(995) / 1000 && isinteger(Mr)
        return pre * _reg2f1_near1(a, b, Int(Mr), w)
    end
    return pre * _reg2f1_far(a, b, c, x)
end

# Regularized ₂F₁(a,b;c;x)/Γ(c), far-from-worldline branch: the no-Γ Maclaurin primitive
# (all terms positive for a,b,c,x > 0 — no cancellation) divided by Γ(c).
_reg2f1_far(a::T, b::T, c::T, x::T) where {T<:Real} = f21_maclaurin(a, b, c, x) / gamma_big(c)

# Regularized ₂F₁(a,b;a+b−M;x)/Γ(a+b−M) for M = 0,1,2,… near x = 1, by direct summation of
# DLMF 15.8.10 (M=0) / 15.8.12 (M≥1) in w = 1−x (passed exactly by the caller):
#   𝐅 = Γ(M)/(Γ(a)Γ(b)) w^{−M} Σ_{k=0}^{M−1} (a−M)_k (b−M)_k / (k! (1−M)_k) w^k              [M≥1]
#     − (−1)^M/(Γ(a−M)Γ(b−M)) Σ_{k≥0} (a)_k (b)_k / (k! (k+M)!) w^k
#           · [ln w − ψ(k+1) − ψ(k+M+1) + ψ(a+k) + ψ(b+k)].
# Truncation follows the CRP: stop once a term is below eps relative to the running total
# (two consecutive small terms, guarding an isolated zero of the ψ/log bracket). The ψ values
# are seeded once via `_psi0` and advanced by the exact recurrence ψ(x+1) = ψ(x) + 1/x.
function _reg2f1_near1(a::T, b::T, M::Int, w::T) where {T<:Real}
    lw = log(w)
    tol = eps(T)
    total = zero(T)
    if M > 0
        t = one(T); s1 = t
        for k in 0:(M - 2)
            t *= (a - M + k) * (b - M + k) / ((1 - M + k) * (k + 1)) * w
            s1 += t
        end
        total = gamma_big(T(M)) / (gamma_big(a) * gamma_big(b)) * w^(-M) * s1
    end
    pre2 = -(iseven(M) ? one(T) : -one(T)) / (gamma_big(a - M) * gamma_big(b - M))
    f = one(T) / factorial_big(T, M)               # (a)_k (b)_k / (k! (k+M)!) w^k at k = 0
    psi = lw - _psi0(one(T)) - _psi0(T(M + 1)) + _psi0(a) + _psi0(b)
    small = 0
    for k in 0:100_000
        contrib = pre2 * f * psi
        total += contrib
        small = abs(contrib) <= tol * abs(total) ? small + 1 : 0
        small >= 2 && break
        f *= (a + k) * (b + k) / ((k + 1) * (k + M + 1)) * w
        psi += -one(T) / (k + 1) - one(T) / (k + M + 1) + one(T) / (a + k) + one(T) / (b + k)
    end
    return total
end

# Digamma ψ(x) at the DLMF-sum seed points, x > 0. At integer / half-odd-integer x (the only
# arguments the puncture reaches: a,b are (half-)integers there) the exact closed forms
#   ψ(j) = −γ + Σ_{i<j} 1/i,   ψ(j+1/2) = −γ − 2 ln 2 + 2 Σ_{i=1}^{j} 1/(2i−1)
# keep full BigFloat precision; anything else falls back to the asymptotic `digamma_big`.
function _psi0(x::T) where {T<:Real}
    γ = T(Base.MathConstants.eulergamma)
    if isinteger(x) && x > 0
        s = -γ
        for i in 1:(Int(x) - 1)
            s += one(T) / i
        end
        return s
    elseif isinteger(x - T(1) / 2) && x > 0
        s = -γ - 2 * log(T(2))
        for i in 1:Int(x - T(1) / 2)
            s += T(2) / (2i - 1)
        end
        return s
    end
    return real(digamma_big(complex(x)))
end
