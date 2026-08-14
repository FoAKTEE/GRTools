# Geometry — canonical Kerr background metric, inverse, scalar bundle, horizons, Kinnersley
# tetrad (M = 1, signature −+++). Consolidates the MetricReconstruction / DiscModes /
# KerrLorenzCirc forks: the three coordinate charts (Boyer–Lindquist, ingoing, outgoing) are
# kept as `coords`-parametrized METHODS, NOT flattened. Type-generic {T}: every array is built
# at the input element type so BigFloat propagates with no Float64 down-cast.

# ---- background scalars ----------------------------------------------------------------------
"Δ = r² − 2 M r + a²  (M = 1)."
@inline delta_kerr(a, r) = r^2 - 2*r + a^2
"Σ = r² + a² cos²θ."
@inline sigma_kerr(a, r, θ) = r^2 + a^2*cos(θ)^2
"ζ = r − i a cosθ (the Kinnersley spin coefficient ρ⁻¹ scale)."
@inline zeta_kerr(a, r, θ) = complex(r, -a*cos(θ))
"K[ω,m] = (r²+a²) ω − a m (radial Teukolsky potential factor)."
@inline kfun(a, r, ω, m) = (r^2 + a^2)*ω - a*m
"Q[ω,m] = −a ω sinθ + m cscθ (angular Teukolsky potential factor)."
@inline qfun(a, ω, m, θ) = -a*ω*sin(θ) + m/sin(θ)

"""
    bg_scalars(a, r, θ) -> (Δ, Σ, ζ)

Geometric Kerr background-scalar bundle at (r, θ) (M = 1). The frequency-dependent factors
`K[ω,m]` and `Q[ω,m]` are `kfun`/`qfun` (separate signatures).
"""
@inline bg_scalars(a, r, θ) = (Δ = delta_kerr(a, r), Σ = sigma_kerr(a, r, θ), ζ = zeta_kerr(a, r, θ))

# ---- horizons --------------------------------------------------------------------------------
"""
    horizons(a) -> (r₊, r₋)

Outer/inner Kerr horizons r± = M ± √(M²−a²) (M = 1): the roots of Δ(a,r)=0.
"""
@inline horizons(a) = (one(a) + sqrt(one(a) - a^2), one(a) - sqrt(one(a) - a^2))

"Horizon angular velocity Ω₊ = a / (2 r₊) (M = 1)."
@inline omega_plus(a) = a / (2 * horizons(a)[1])

# ---- tortoise coordinate (canonical Kerr r* + inverse; DG-P2) ----------------------------------
"""
    tortoise_kerr(r, a) -> r*

Canonical Kerr tortoise coordinate (M = 1), the antiderivative of d(r*)/dr = (r²+a²)/Δ on the
exterior r > r₊, with the standard integration constant:

    r* = r + (2 r₊/(r₊−r₋)) ln((r−r₊)/2) − (2 r₋/(r₊−r₋)) ln((r−r₋)/2).

Convention: M = 1 like every neighbor in this file (the DG-Teukolsky reference parametrizes M
explicitly; here it is fixed to Spacetime's M = 1). At a = 0 this DELEGATES to `tortoise_schw`,
so the Schwarzschild reduction is bit-exact. Type-generic: `promote(r, a)` so BigFloat is
inherited with no Float64 down-cast.
"""
function tortoise_kerr(r, a)
    iszero(a) && return tortoise_schw(r)
    rr, aa = promote(r, a)
    rp, rm = horizons(aa)
    rr + (2 * rp / (rp - rm)) * log((rr - rp) / 2) - (2 * rm / (rp - rm)) * log((rr - rm) / 2)
end

"""
    tortoise_kerr_inverse(rs, a; tol, itmax) -> r

Invert the Kerr tortoise coordinate: solve `tortoise_kerr(r, a) = rs` for r ∈ (r₊, ∞) (M = 1),
rs ∈ (−∞, ∞). Newton iteration with dr/d(r*) = Δ/(r²+a²) from a two-branch initial guess
(linear r ≈ rs for large rs, exponential horizon approach r₊ + 2 e^{(rs−r₊)/2−1} otherwise) and
a bisection-style safeguard whenever a Newton step leaves the exterior — the strategy of the
DG-Teukolsky reference `inverse_tortoise`. Converges on the r*-residual |r*(r) − rs| <
tol·(1+|rs|) (so the round-trip error is directly controlled) with a stagnation fallback on the
step size. Type-generic: `tol` inherits precision through `eps`, so BigFloat is fully supported.
"""
function tortoise_kerr_inverse(rs, a; kwargs...)
    rsp, ap = promote(float(rs), float(a))
    _tortoise_kerr_inverse(rsp, ap; kwargs...)
end

function _tortoise_kerr_inverse(rs::T, a::T; tol::Real = 8 * eps(T), itmax::Int = 200) where {T<:AbstractFloat}
    rp, _ = horizons(a)
    # initial guess (DG strategy): large-r* linear branch, else exponential horizon approach
    r = rs > 4 ? rs : rp + 2 * exp((rs - rp) / 2 - 1)
    r = max(r, rp + eps(rp))                               # keep strictly outside the horizon
    for _ in 1:itmax
        f = tortoise_kerr(r, a) - rs
        abs(f) < tol * (1 + abs(rs)) && return r           # r*-residual converged
        rnew = r - f * delta_kerr(a, r) / (r^2 + a^2)      # Newton: r -= f / (dr*/dr)
        if !(rnew > rp) || !isfinite(rnew)                 # bisection-style safeguard (DG)
            rnew = (r + rp) / 2 + (r - rp) / 4 * sign(-f)
            rnew = max(rnew, rp + eps(rp))
        end
        abs(rnew - r) < tol * (1 + abs(r)) && return rnew  # step stagnation (ulp floor near r₊)
        r = rnew
    end
    return r
end

# ---- metric (coords-parametrized methods) ----------------------------------------------------
@inline function _geometry_promote(a::Real, r::Real, θ::Real)
    T = float(promote_type(typeof(a), typeof(r), typeof(θ)))
    return T(a), T(r), T(θ)
end

function _kerr_metric_BL(a::T, r, θ) where {T}
    s2 = sin(θ)^2; Σ = sigma_kerr(a, r, θ); Δ = delta_kerr(a, r)
    g = zeros(T, 4, 4)
    g[1,1] = -(1 - 2*r/Σ)
    g[1,4] = g[4,1] = -2*a*r*s2/Σ
    g[2,2] = Σ/Δ
    g[3,3] = Σ
    g[4,4] = (r^2 + a^2 + 2*a^2*r*s2/Σ)*s2
    g
end

function _kerr_metric_In(a::T, r, θ) where {T}
    s2 = sin(θ)^2; Σ = sigma_kerr(a, r, θ); Δ = delta_kerr(a, r)
    g = zeros(T, 4, 4)
    g[1,1] = -(1 - 2*r/Σ)
    g[1,2] = g[2,1] = one(T)
    g[1,4] = g[4,1] = -a*s2*(2*r)/Σ
    g[2,4] = g[4,2] = -a*s2
    g[3,3] = Σ
    g[4,4] = ((r^2+a^2)^2 - Δ*a^2*s2)/Σ * s2
    g
end

function _kerr_metric_Out(a::T, r, θ) where {T}
    s2 = sin(θ)^2; Σ = sigma_kerr(a, r, θ); Δ = delta_kerr(a, r)
    g = zeros(T, 4, 4)
    g[1,1] = -(1 - 2*r/Σ)
    g[1,2] = g[2,1] = -one(T)
    g[1,4] = g[4,1] = -a*s2*(2*r)/Σ
    g[2,4] = g[4,2] = a*s2
    g[3,3] = Σ
    g[4,4] = ((r^2+a^2)^2 - Δ*a^2*s2)/Σ * s2
    g
end

"""
    kerr_metric(a, r, θ; coords=:BL) -> 4×4 g_{μν}

Kerr background metric (M = 1, signature −+++) in coordinate chart `coords ∈ (:BL, :In, :Out)`:
Boyer–Lindquist, ingoing-Kerr, outgoing-Kerr. Order is (t, r, θ, φ).
"""
function kerr_metric(a::Real, r::Real, θ::Real; coords::Symbol = :BL)
    aa, rr, tt = _geometry_promote(a, r, θ)
    coords === :BL  && return _kerr_metric_BL(aa, rr, tt)
    coords === :In  && return _kerr_metric_In(aa, rr, tt)
    coords === :Out && return _kerr_metric_Out(aa, rr, tt)
    throw(ArgumentError("kerr_metric: unknown coords $(repr(coords)) (expected :BL, :In, :Out)"))
end

"""
    kerr_metric_inv(a, r, θ; coords=:BL) -> 4×4 g^{μν}

Inverse Kerr metric. For `:BL` the closed-form (t,φ) 2×2 block + diagonal (r,θ) inversion;
the non-orthogonal `:In`/`:Out` charts fall back to a dense type-generic inverse.
"""
function kerr_metric_inv(a::Real, r::Real, θ::Real; coords::Symbol = :BL)
    aa, rr, tt = _geometry_promote(a, r, θ)
    T = typeof(aa)
    if coords === :BL
        g = _kerr_metric_BL(aa, rr, tt); gi = zeros(T, 4, 4)
        det2 = g[1,1]*g[4,4] - g[1,4]^2
        gi[1,1] = g[4,4]/det2; gi[4,4] = g[1,1]/det2; gi[1,4] = gi[4,1] = -g[1,4]/det2
        gi[2,2] = 1/g[2,2]; gi[3,3] = 1/g[3,3]
        return gi
    end
    inv(kerr_metric(aa, rr, tt; coords = coords))
end

# ---- Kinnersley null tetrad (contravariant legs, coords-parametrized) -------------------------
function _kinnersley_BL(a::T, r, θ) where {T}
    Δ = delta_kerr(a, r); Σ = sigma_kerr(a, r, θ); sθ = sin(θ)
    ζ = zeta_kerr(a, r, θ); ζb = conj(ζ); rt2 = sqrt(T(2))
    lup    = Complex{T}[(r^2+a^2)/Δ, one(T), zero(T), a/Δ]
    nup    = Complex{T}[r^2+a^2, -Δ, zero(T), a] ./ (2*Σ)
    mup    = Complex{T}[im*a*sθ, zero(T), one(T), im/sθ] ./ (rt2*ζb)
    mbarup = Complex{T}[-im*a*sθ, zero(T), one(T), -im/sθ] ./ (rt2*ζ)
    (lup, nup, mup, mbarup)
end

function _kinnersley_In(a::T, r, θ) where {T}
    Δ = delta_kerr(a, r); Σ = sigma_kerr(a, r, θ); sθ = sin(θ)
    ζ = zeta_kerr(a, r, θ); ζb = conj(ζ); rt2 = sqrt(T(2))
    lup    = Complex{T}[2*(a^2+r^2)/Δ, one(T), zero(T), 2*a/Δ]
    nup    = Complex{T}[zero(T), -Δ/(2*Σ), zero(T), zero(T)]
    mup    = Complex{T}[im*a*sθ, zero(T), one(T), im/sθ] ./ (rt2*ζb)
    mbarup = Complex{T}[-im*a*sθ, zero(T), one(T), -im/sθ] ./ (rt2*ζ)
    (lup, nup, mup, mbarup)
end

function _kinnersley_Out(a::T, r, θ) where {T}
    Δ = delta_kerr(a, r); Σ = sigma_kerr(a, r, θ); sθ = sin(θ)
    ζ = zeta_kerr(a, r, θ); ζb = conj(ζ); rt2 = sqrt(T(2))
    lup    = Complex{T}[zero(T), one(T), zero(T), zero(T)]
    nup    = Complex{T}[r^2+a^2, -Δ/2, zero(T), a] ./ Σ
    mup    = Complex{T}[im*a*sθ, zero(T), one(T), im/sθ] ./ (rt2*ζb)
    mbarup = Complex{T}[-im*a*sθ, zero(T), one(T), -im/sθ] ./ (rt2*ζ)
    (lup, nup, mup, mbarup)
end

"""
    kinnersley_tetrad(a, r, θ; coords=:BL) -> (l, n, m, m̄)

Contravariant Kinnersley null tetrad legs in chart `coords ∈ (:BL, :In, :Out)` (M = 1).
Normalization (signature −+++): l·l = n·n = m·m = 0, l·n = −1, m·m̄ = 1.
"""
function kinnersley_tetrad(a::Real, r::Real, θ::Real; coords::Symbol = :BL)
    aa, rr, tt = _geometry_promote(a, r, θ)
    coords === :BL  && return _kinnersley_BL(aa, rr, tt)
    coords === :In  && return _kinnersley_In(aa, rr, tt)
    coords === :Out && return _kinnersley_Out(aa, rr, tt)
    throw(ArgumentError("kinnersley_tetrad: unknown coords $(repr(coords)) (expected :BL, :In, :Out)"))
end
