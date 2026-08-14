# Martel-Poisson tensor-harmonic angular basis and the coordinate reconstruction of the regular
# metric perturbation hR. Port of Punctures/tensor_basis.jl (hR.m L63-153, Wardell et al.
# arXiv:2508.00087), REDESIGNED:
#
#   • the scalar harmonic seed is CONSUMED from `Harmonics.spin_weighted_spherical_Y` (s=0) — no
#     half-angle harmonic is re-implemented here.
#   • the legacy SWSH `_y_dtheta` (the exact k-th θ-derivative of the seed) is NOT exposed by the
#     consolidated Harmonics surface, so the θ-derivative is supplied LOCALLY by ForwardDiff over the
#     consumed seed (`_dthY`, nested automatic differentiation up to order 4 — the max the AngTerm
#     `_dtheta` machinery requests: basis k≤2 plus nθ≤2). AD is exact to machine precision, so this
#     reproduces the analytic θ-derivative the legacy term list carried.
#
# The eleven basis functions are fixed (l,m,θ,φ)-derivative combinations of the seed, stored as a
# half-angle-free term list `coeff·sin(θ)^sp·cos(θ)^cp·csc(θ)^isp·∂θ^k Y` supporting exact analytic
# θ-differentiation. The coordinate components (htt..hφφ) + their deriv[nt,nr,nθ,nφ] reconstruct the
# metric from the per-(l,m) Lorenz amplitudes hR[i][l,m] (+ dhR,ddhR) on the worldline (θ=π/2,φ=0,r=r0):
# time/azimuthal derivatives are the algebraic factors (-imΩ)^nt, (im)^nφ; r,θ derivatives analytic.

import ForwardDiff
import Harmonics

# Real θ-part of the consumed seed {}_0Y_lm(θ,0) (the cis(mφ) phase is applied in _eval_angular).
_seedY(l::Integer, m::Integer, θ) = real(Harmonics.spin_weighted_spherical_Y(0, l, m, θ, zero(θ)))

"""
∂θ^k of the real θ-part of `{}_0Y_lm(θ,0)`, by nested ForwardDiff over the consumed Harmonics seed
(`_y_dtheta` is not exposed). `k ∈ 0..4` — the max the basis (k≤2) and nθ-differentiation (nθ≤2)
request. AD is exact, so this equals the analytic θ-derivative.
"""
function _dthY(l::Integer, m::Integer, θ::Real, k::Integer)
    θf = float(θ)
    f0 = t -> _seedY(l, m, t)
    k == 0 && return f0(θf)
    f1 = t -> ForwardDiff.derivative(f0, t)
    k == 1 && return f1(θf)
    f2 = t -> ForwardDiff.derivative(f1, t)
    k == 2 && return f2(θf)
    f3 = t -> ForwardDiff.derivative(f2, t)
    k == 3 && return f3(θf)
    f4 = t -> ForwardDiff.derivative(f3, t)
    k == 4 && return f4(θf)
    throw(ArgumentError("_dthY: θ-derivative order k=$k > 4 not supported"))
end

_pow0(b, e::Integer) = e == 0 ? one(b) : b^e

_Lscalar(l) = l * (l + 1)
_L2scalar(l) = (l - 1) * l * (l + 1) * (l + 2)
_im(m) = ComplexF64(0, m)

"A single angular term `coeff·sin(θ)^sp·cos(θ)^cp·csc(θ)^isp·∂θ^k Y_lm`."
struct AngTerm
    c::ComplexF64
    sp::Int
    cp::Int
    isp::Int
    k::Int
end

# d/dθ of a term list (product rule; d(sin)=cos, d(cos)=-sin, d(csc)=-cos·csc²).
function _dtheta(terms::Vector{AngTerm})
    out = AngTerm[]
    for t in terms
        t.sp != 0 && push!(out, AngTerm(t.c * t.sp, t.sp - 1, t.cp + 1, t.isp, t.k))
        t.cp != 0 && push!(out, AngTerm(-t.c * t.cp, t.sp + 1, t.cp - 1, t.isp, t.k))
        t.isp != 0 && push!(out, AngTerm(-t.c * t.isp, t.sp, t.cp + 1, t.isp + 1, t.k))
        push!(out, AngTerm(t.c, t.sp, t.cp, t.isp, t.k + 1))
    end
    return out
end

function _eval_angular(terms::Vector{AngTerm}, l, m, θ, φ)
    s = sin(θ); c = cos(θ)
    acc = zero(ComplexF64)
    for t in terms
        trig = _pow0(s, t.sp) * _pow0(c, t.cp) / _pow0(s, t.isp)
        acc += t.c * trig * _dthY(l, m, θ, t.k)
    end
    return acc * cis(m * φ)
end

# Basis registry: key -> (l,m) -> (scalar prefactor, bare term list).
function _basis(key::Symbol, l::Integer, m::Integer)
    L = _Lscalar(l); L2 = _L2scalar(l)
    pZ = 1 / sqrt(complex(L)); p2 = sqrt(2) / sqrt(complex(L2))
    key === :Y      && return (1.0 + 0im, AngTerm[AngTerm(1, 0, 0, 0, 0)])
    key === :Ysin2  && return (1.0 + 0im, AngTerm[AngTerm(1, 2, 0, 0, 0)])
    key === :Ythth  && return (1.0 + 0im, AngTerm[AngTerm(1, 0, 0, 0, 2)])
    key === :Zth    && return (pZ, AngTerm[AngTerm(1, 0, 0, 0, 1)])
    key === :Zph    && return (pZ, AngTerm[AngTerm(_im(m), 0, 0, 0, 0)])
    key === :Xth    && return (-pZ, AngTerm[AngTerm(_im(m), 0, 0, 1, 0)])
    key === :Xph    && return (pZ, AngTerm[AngTerm(1, 1, 0, 0, 1)])
    key === :Zthth  && return (p2, AngTerm[AngTerm(1, 0, 0, 0, 2), AngTerm(L / 2, 0, 0, 0, 0)])
    key === :Zthph  && return (p2, AngTerm[AngTerm(_im(m), 0, 0, 0, 1), AngTerm(-_im(m), 0, 1, 1, 0)])
    key === :Zphph  && return (p2, AngTerm[AngTerm(-m^2, 0, 0, 0, 0), AngTerm(1, 1, 1, 0, 1), AngTerm(L / 2, 2, 0, 0, 0)])
    key === :Xthth  && return (-p2, AngTerm[AngTerm(_im(m), 0, 0, 1, 1), AngTerm(-_im(m), 0, 1, 2, 0)])
    key === :Xthph  && return (-p2 / 2, AngTerm[AngTerm(-m^2, 0, 0, 1, 0), AngTerm(-1, 1, 0, 0, 2), AngTerm(1, 0, 1, 0, 1)])
    key === :Xphph  && return (p2, AngTerm[AngTerm(_im(m), 1, 0, 0, 1), AngTerm(-_im(m), 0, 1, 0, 0)])
    throw(ArgumentError("unknown basis key $key"))
end

"""
    tensor_harmonic(key, l, m, θ, φ) -> Complex

Martel-Poisson angular basis function `key` ∈ {:Ythth, :Zth, :Zph, :Xth, :Xph, :Zthth, :Zthph,
:Zphph, :Xthth, :Xthph, :Xphph} at (l,m,θ,φ) (hR.m L63-93).
"""
function tensor_harmonic(key::Symbol, l::Integer, m::Integer, θ::Real, φ::Real)
    abs(m) <= l || throw(ArgumentError("need |m| <= l"))
    pref, terms = _basis(key, l, m)
    return pref * _eval_angular(terms, l, m, float(θ), float(φ))
end

# ---- coordinate reconstruction ------------------------------------------------

# A reconstruction channel: mode index i, l-dependent radial scalar Cextra(l), r^a f^b exponents,
# ± sign on the basis combination, basis key, and the minimum l at which the channel is active
# (parity/degeneracy guard, hR.m L105-123).
struct Channel
    i::Int
    Cextra::Function
    a::Int
    b::Int
    sign::Int
    basis::Symbol
    lmin::Int
end

const _HALF = 0.5
_invsqrtL(l) = 1 / sqrt(float(_Lscalar(l)))
_p2half(l) = sqrt(2) / (2 * sqrt(float(_L2scalar(l))))

# Coordinate components keyed by symbol (hR.m L96-123).
const _COMPONENTS = Dict{Symbol,Vector{Channel}}(
    :htt => [Channel(1, l -> _HALF, -1, 0, 1, :Y, 0),
             Channel(6, l -> _HALF, -1, 1, 1, :Y, 0)],
    :htr => [Channel(2, l -> _HALF, -1, -1, 1, :Y, 0)],
    :hrr => [Channel(1, l -> _HALF, -1, -2, 1, :Y, 0),
             Channel(6, l -> -_HALF, -1, -1, 1, :Y, 0)],
    :htth => [Channel(4, l -> _HALF * _invsqrtL(l), 0, 0, 1, :Zth, 1),
              Channel(8, l -> _HALF * _invsqrtL(l), 0, 0, -1, :Xth, 1)],
    :htph => [Channel(4, l -> _HALF * _invsqrtL(l), 0, 0, 1, :Zph, 1),
              Channel(8, l -> _HALF * _invsqrtL(l), 0, 0, -1, :Xph, 1)],
    :hrth => [Channel(5, l -> _HALF * _invsqrtL(l), 0, -1, 1, :Zth, 1),
              Channel(9, l -> _HALF * _invsqrtL(l), 0, -1, -1, :Xth, 1)],
    :hrph => [Channel(5, l -> _HALF * _invsqrtL(l), 0, -1, 1, :Zph, 1),
              Channel(9, l -> _HALF * _invsqrtL(l), 0, -1, -1, :Xph, 1)],
    :hthth => [Channel(3, l -> _HALF, 1, 0, 1, :Y, 0),
               Channel(7, l -> _p2half(l), 1, 0, 1, :Zthth, 2),
               Channel(10, l -> _p2half(l), 1, 0, -1, :Xthth, 2)],
    :hthph => [Channel(7, l -> _p2half(l), 1, 0, 1, :Zthph, 2),
               Channel(10, l -> _p2half(l), 1, 0, -1, :Xthph, 2)],
    :hphph => [Channel(3, l -> _HALF, 1, 0, 1, :Ysin2, 0),
               Channel(7, l -> _p2half(l), 1, 0, 1, :Zphph, 2),
               Channel(10, l -> _p2half(l), 1, 0, -1, :Xphph, 2)],
)

# Parity selection (hR.m L35): mode is zero when it violates the i-parity rule.
_parity_zero(i, l, m) =
    (i in (1, 2, 3, 4, 5, 6, 7) && isodd(l + m)) ||
    (i in (8, 9, 10) && iseven(l + m))

# Mode access with the reality rule hR[i][l,-m] = (-1)^m conj(hR[i][l,m]) (hR.m L30) + parity
# selection. `base(i,l,m)` returns the (h, dh, ddh) triple for valid m >= 0.
function _mode_triple(base, i, l, m)
    _parity_zero(i, l, m) && return (0.0 + 0im, 0.0 + 0im, 0.0 + 0im)
    if m < 0
        h, dh, ddh = base(i, l, -m)
        s = iseven(m) ? 1 : -1
        return (s * conj(h), s * conj(dh), s * conj(ddh))
    end
    return base(i, l, m)
end

# (g, g', g'') of g(r) = C·r^a·(1-2M/r)^b at r.
function _radial_derivs(C, a::Int, b::Int, r, M)
    f = 1 - 2M / r
    fp = 2M / r^2
    fpp = -4M / r^3
    fb = f^b; fbm1 = f^(b - 1); fbm2 = f^(b - 2)
    g = C * r^a * fb
    gp = C * (a * r^(a - 1) * fb + r^a * b * fbm1 * fp)
    gpp = C * (a * (a - 1) * r^(a - 2) * fb +
               2 * a * r^(a - 1) * b * fbm1 * fp +
               r^a * (b * (b - 1) * fbm2 * fp^2 + b * fbm1 * fpp))
    return (g, gp, gpp)
end

# ∂r^nr of g(r)·mode(r) given (g,g',g'') and (h,dh,ddh); nr ∈ 0..2.
function _radial_mode(g::NTuple{3}, h::NTuple{3}, nr::Int)
    nr == 0 && return g[1] * h[1]
    nr == 1 && return g[2] * h[1] + g[1] * h[2]
    nr == 2 && return g[3] * h[1] + 2 * g[2] * h[2] + g[1] * h[3]
    throw(ArgumentError("nr=$nr > 2 needs mode derivatives beyond ddhR"))
end

"""
    coord_deriv(comp, l, base, r0, M; nt=0, nr=0, nθ=0, nφ=0, θ=π/2, φ=0) -> Real

Reconstructed coordinate component `comp` ∈ {:htt,:htr,:hrr,:htth,:htph,:hrth,:hrph,:hthth,:hthph,
:hphph} at multipole `l`, optionally differentiated `deriv[nt,nr,nθ,nφ]` (hR.m L96-153). `base(i,l,m)`
supplies the (hR, dhR, ddhR) amplitudes for valid m ≥ 0; the reality + parity rules are applied
internally. Returns the real component value (hR.m wraps each sum in `Re`).
"""
function coord_deriv(comp::Symbol, l::Integer, base, r0::Real, M::Real;
                     nt::Int = 0, nr::Int = 0, nθ::Int = 0, nφ::Int = 0,
                     θ::Real = π / 2, φ::Real = 0.0)
    channels = _COMPONENTS[comp]
    Ω = sqrt(M / r0^3)
    total = zero(ComplexF64)
    for ch in channels
        l < ch.lmin && continue
        C0 = ch.Cextra(l)
        for m in -l:l
            ft = nt > 0 ? (ComplexF64(0, -m * Ω))^nt : (1.0 + 0im)
            fp = nφ > 0 ? (ComplexF64(0, m))^nφ : (1.0 + 0im)
            (ft == 0 || fp == 0) && continue
            bp, bare = _basis(ch.basis, l, m)
            C = C0 * ch.sign * bp
            g = _radial_derivs(C, ch.a, ch.b, float(r0), float(M))
            h = _mode_triple(base, ch.i, l, m)
            rad = _radial_mode(g, h, nr)
            terms = bare
            for _ in 1:nθ
                terms = _dtheta(terms)
            end
            ang = _eval_angular(terms, l, m, float(θ), float(φ))
            total += ft * fp * rad * ang
        end
    end
    return real(total)
end

"""
    coord_component(comp, l, base, r0, M; θ=π/2, φ=0) -> Real

Undifferentiated coordinate component (`coord_deriv` with all derivative orders zero).
"""
coord_component(comp::Symbol, l::Integer, base, r0::Real, M::Real;
                θ::Real = π / 2, φ::Real = 0.0) =
    coord_deriv(comp, l, base, r0, M; θ = θ, φ = φ)
