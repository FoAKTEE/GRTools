# Static (ω=0) Teukolsky radial solutions + amplitudes — closed-form 2F1 path. The MST machinery
# degenerates at ε=2ω=0; static perturbations (no radiation) use a Hypergeometric2F1 closed form.
# λ = l(l+1) − s(s+1) (SWSH eigenvalue at aω=0); ν passed = λ. Self-contained — no MST ν/f_n.
# Pochhammer delegates to BigSpecial.pochhammer. Type-generic {T} (design §9 inheritance chain:
# "MST/CF/renorm-ν/STATIC, all type-generic {T}"): the working precision T = float(a) flows
# through τ/κ, the amplitudes, the norm and the radial value/deriv — a BigFloat `a` gives
# Complex{BigFloat} end-to-end. Γ: Float64 stays on the committed `SpecialFunctions.gamma`
# complex path (bit-identical regression, test_radial_static_precision.jl bake); any other
# precision routes through the canonical `BigSpecial.gamma_big` (reflection on Re z < 1/2, the
# validated Complex{BigFloat} half-plane). ₂F₁ goes through the local `_safe_2f1` (hypergeom.jl):
# byte-identical upstream `_₂F₁` for Float64, BigSpecial Pfaff/z→∞ fallback for BigFloat.
# The static branch of the public `TeukolskyRadial` dispatcher (Radial/teukolsky_radial.jl).

import SpecialFunctions

# Γ for the static path: ComplexF64 keeps the exact committed upstream call (Float64 bit-identity
# with the baked pre-{T} references); every other Complex{T} uses BigSpecial.gamma_big.
_static_gamma(z::ComplexF64) = SpecialFunctions.gamma(z)
_static_gamma(z::Complex) = gamma_big(z)

"""
    _static_amplitudes(s, l, m, a) -> (In=NamedTuple, Up=NamedTuple)

Static asymptotic amplitudes. κ=√(1−a²), τ=−ma/√(1−a²). MMA keys ℋ (horizon), ℐ (infinity) map
to fields H, I; the τ≠0 split keys ℐ±,ℋ± map to Im,Ip,Hm,Hp. τ==0: In=(H=1, I=…), Up=(H=…, I=1).
τ≠0: In=(H=1, Im=…, Ip=…), Up=(Hm=…, Hp=…, I=1).
"""
function _static_amplitudes(s::Integer, l::Integer, m::Integer, a::Real)
    af = float(a)                                # working precision T = float(a) (§9 inheritance)
    tau = -(m * af) / sqrt(1 - af^2)
    kappa = sqrt(1 - af^2)
    CT = Complex{typeof(kappa)}                  # Complex{T}
    g(z) = _static_gamma(convert(CT, z))
    if tau == 0
        ampInI = ((2 * kappa)^(-l + abs(s)) * g(2l + 1) * g(1 + abs(s))) /
                 (g(l + 1) * g(1 + l + abs(s)))
        ampUpH = s == 0 ?
            -((2 * kappa)^(-1 - l) * g(3 + 2l)) / (2 * g(l + 2) * g(1 + l)) :
            ((2 * kappa)^(-1 - l + abs(s)) * g(2 + 2l) * g(abs(s))) /
                (g(1 + l + abs(s)) * g(1 + l))
        return (In=(H=one(CT), I=convert(CT, ampInI)),
                Up=(H=convert(CT, ampUpH), I=one(CT)))
    else
        ampInIm = -((2.0^(l - s - im * tau) * complex(kappa)^(1 - s) * (-1)^(l + s) *
                     g(1 + l + s) * g(1 - s - im * tau)) / (g(2 + 2l) * g(-l - im * tau)))
        ampInIp = ((2 * kappa)^(-l - s - im * tau) * g(2l + 1) * g(1 - s - im * tau)) /
                  (g(1 + l - s) * g(1 + l - im * tau))
        ampUpHm = ((2 * kappa)^(-1 - l + s + im * tau) * g(2 + 2l) * g(s + im * tau)) /
                  (g(1 + l + s) * g(1 + l + im * tau))
        ampUpHp = ((2 * kappa)^(-1 - l - s - im * tau) * g(2 + 2l) * g(-s - im * tau)) /
                  (g(1 + l - s) * g(1 + l - im * tau))
        return (In=(H=one(CT), Im=convert(CT, ampInIm), Ip=convert(CT, ampInIp)),
                Up=(Hm=convert(CT, ampUpHm), Hp=convert(CT, ampUpHp), I=one(CT)))
    end
end

"Regularized ₂F₁: `₂F₁(a,b,c,z)/Γ(c)` (MMA Hypergeometric2F1Regularized; handles c≤0 integer)."
function _hyp2f1_regularized(a, b, c, z)
    cc = complex(c)
    if isreal(cc) && real(cc) <= 0 && isinteger(real(cc))
        # c = −n (non-positive integer): Γ(c)=∞ and the first n+1 series terms vanish, leaving the
        # finite DLMF 15.2.6 limit. Arises for the static In at s>0, τ=0 (c=1−s ≤ 0).
        n = Int(-real(cc)); za = complex(a); zb = complex(b)
        return pochhammer(za, n + 1) * pochhammer(zb, n + 1) / factorial(n + 1) *
               complex(z)^(n + 1) *
               _safe_2f1(za + (n + 1), zb + (n + 1), complex(n + 2), complex(z))
    end
    return _safe_2f1(complex(a), complex(b), cc, complex(z)) / _static_gamma(cc)
end

"""
    _static_norm(s, l, m, a, bc) -> Complex{T}   (T = typeof(float(a)))

Static radial normalisation. normUp=(2κ)^(−s−l−1); normIn=(2κ)^(−2s−iτ)Γ(1−s−iτ), or
Γ(s+1)·Pochhammer(l+s+1,−2s) when τ==0 && s>0.
"""
function _static_norm(s::Integer, l::Integer, m::Integer, a::Real, bc::Symbol)
    af = float(a)                                # working precision T = float(a) (§9 inheritance)
    tau = -(m * af) / sqrt(1 - af^2)
    kappa = sqrt(1 - af^2)
    CT = Complex{typeof(kappa)}                  # Complex{T}
    if bc === :Up
        return convert(CT, (2 * kappa)^(-s - l - 1))
    end
    if tau == 0 && s > 0
        return convert(CT, SpecialFunctions.gamma(one(kappa) * (s + 1)) *
                           pochhammer(one(kappa) * (l + s + 1), -2s))
    end
    return convert(CT, (2 * kappa)^(-2s - im * tau) *
                       _static_gamma(convert(CT, 1 - s - im * tau)))
end

"""
    _static_radial(s, l, m, a, bc, r; norm) -> Complex (precision of promote(float(a), float(r)))

Static radial value R^In(r)/R^Up(r). x=(1+κ−r)/(2κ); for r>rp=1+κ, x<0 (so −x>0) and 1−x>0.
In: ₂F₁ᵣₑ𝓰(−l−iτ, l+1−iτ, 1−s−iτ, x); Up: ₂F₁(l+1−iτ, l+1−s, 2l+2, 1/(1−x)). `norm` is the
precomputed `_static_norm`.
"""
function _static_radial(s::Integer, l::Integer, m::Integer, a::Real, bc::Symbol, r::Real; norm)
    tau = -(m * a) / sqrt(1 - a^2)
    kappa = sqrt(1 - a^2)
    x = (1 + kappa - r) / (2 * kappa)
    mx = complex(-x)            # −x > 0 for r > rp
    omx = complex(1 - x)        # 1 − x > 0 for r > rp
    if bc === :In
        return norm * mx^(-s - im * tau / 2) * omx^(-im * tau / 2) *
               _hyp2f1_regularized(-l - im * tau, l + 1 - im * tau, 1 - s - im * tau, x)
    else
        return norm * mx^(-s - im * tau / 2) * omx^(im * tau / 2 - l - 1) *
               _safe_2f1(complex(l + 1 - im * tau), complex(l + 1 - s),
                         complex(2l + 2), complex(1 / (1 - x)))
    end
end

"""
    _static_radial_deriv(s, l, m, a, bc, r; norm) -> Complex (precision of promote(float(a), float(r)))

Analytic dR/dr of the static (ω=0) radial solution [`_static_radial`] — a NEW closed-form primitive
(legacy `teukolsky_radial_deriv` throws for ω=0, so there is no legacy oracle). Product rule on the
x=(1+κ−r)/(2κ) power prefactor plus the regularized/ordinary ₂F₁ derivative
(d/dx ₂F̃₁(a,b;c;x)=ab ₂F̃₁(a+1,b+1;c+1;x); d/dx ₂F₁(a,b;c;y)=(ab/c) ₂F₁(a+1,b+1;c+1;y)). `norm` is
the same `_static_norm` used by the value, so R and dR share the one normalization.
"""
function _static_radial_deriv(s::Integer, l::Integer, m::Integer, a::Real, bc::Symbol, r::Real; norm)
    tau = -(m * a) / sqrt(1 - a^2)
    kappa = sqrt(1 - a^2)
    x = (1 + kappa - r) / (2 * kappa)
    xp = -1 / (2 * kappa)                         # dx/dr
    mx = complex(-x)                              # −x > 0 for r > rp
    omx = complex(1 - x)                          # 1 − x > 0 for r > rp
    p1 = -s - im * tau / 2                         # exponent of mx (both BCs)
    if bc === :In
        p2 = -im * tau / 2                         # exponent of omx
        A = -l - im * tau; B = l + 1 - im * tau; C = 1 - s - im * tau
        F  = _hyp2f1_regularized(A, B, C, x)
        Fp = A * B * _hyp2f1_regularized(A + 1, B + 1, C + 1, x)   # dF/dx (regularized identity)
        P = mx^p1 * omx^p2
        # dR/dr = norm·P·[ −xp·(p1/mx + p2/omx)·F + Fp·xp ]  (chain rule; dmx/dr=domx/dr=−xp)
        return norm * P * (-xp * (p1 / mx + p2 / omx) * F + Fp * xp)
    else
        q2 = im * tau / 2 - l - 1                  # exponent of omx (Up)
        A = l + 1 - im * tau; B = l + 1 - s; C = 2l + 2
        y = 1 / (1 - x)
        G  = _safe_2f1(complex(A), complex(B), complex(C), complex(y))
        Gp = (A * B / C) *
             _safe_2f1(complex(A + 1), complex(B + 1), complex(C + 1), complex(y))
        Q = mx^p1 * omx^q2
        dydr = xp / omx^2                          # d/dr [1/(1−x)]
        return norm * Q * (-xp * (p1 / mx + q2 / omx) * G + Gp * dydr)
    end
end

"""
    _teukolsky_radial_static(s, l, m, a) -> (In=TeukolskyRadialFunction, Up=...)

Static (ω=0) TeukolskyRadial. λ=l(l+1)−s(s+1); ν passed = λ; amplitudes = `_static_amplitudes`.
The carrier's `omega==0` / `method=:Static` routes its `(R)(r)` call to `_static_radial`; `norm`
holds the static prefactor.
"""
function _teukolsky_radial_static(s::Integer, l::Integer, m::Integer, a::Real)
    af = float(a)                                    # working precision T = float(a) (§9)
    lambda = oftype(af, l * (l + 1) - s * (s + 1))   # = SWSH eigenvalue at aω=0 (exact integer)
    amps = _static_amplitudes(s, l, m, af)
    normIn = _static_norm(s, l, m, af, :In)
    normUp = _static_norm(s, l, m, af, :Up)
    In = TeukolskyRadialFunction(s, l, m, af, zero(af), lambda, lambda,
                                 :In, normIn, amps.In, :Static)
    Up = TeukolskyRadialFunction(s, l, m, af, zero(af), lambda, lambda,
                                 :Up, normUp, amps.Up, :Static)
    return (In=In, Up=Up)
end
