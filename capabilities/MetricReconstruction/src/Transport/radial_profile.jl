# radial_profile.jl — the HOLE-2 distributional radial-profile algebra.
#
# The Schwarzschild static thin-shell specialization of the Li–Yunes transport chain
# (arXiv:2605.11080v1, Supplemental Material) is entirely CLOSED-FORM: the sources are
# distributional (δ, δ′) on the shell r=R and the fields are piecewise-Laurent × Θ. This file
# provides the two {T}-generic value types the chain needs — NO DP5 ODE solve, NO symbolic CAS:
#
#   Laurent{T}      — a finite Laurent polynomial  Σ_k c_k r^k , k ∈ ℤ  (exact +,−,×,d/dr,∫,eval).
#   ShellProfile{T} — g⁻(r)Θ(R−r) + g⁺(r)Θ(r−R) + a·δ(r−R) + b·δ′(r−R)  at fixed shell radius R,
#                     with the algebra the chain uses: +, −, scale, multiply by a smooth Laurent
#                     (applying the distribution identity  f(r)δ′ = f(R)δ′ − f′(R)δ , LY:831),
#                     differentiate (Θ-jumps generate δ), and jump/continuity extraction.
#
# {T}-generic throughout (element type promoted from arguments); BigFloat (256…2048-bit) inherited.

# ================================================================================================
# Laurent{T}
# ================================================================================================
"""
    Laurent{T}

A finite Laurent polynomial `Σ_k c_k r^k` over `k ∈ ℤ`, coefficients of element type `T` (kept in a
sparse `Dict{Int,T}`, zero coefficients pruned). Supports exact `+`, `-`, `*` (Laurent×Laurent),
scalar scaling ([`lscale`](@ref)), differentiation ([`lderiv`](@ref)), antidifferentiation
([`lintegrate`](@ref) — errors on an `r^{-1}` term, which would integrate to a log), and point
evaluation ([`leval`](@ref)). {T}-generic (BigFloat inherited).
"""
struct Laurent{T}
    c::Dict{Int,T}
end

Laurent{T}() where {T} = Laurent{T}(Dict{Int,T}())

"""
    laurent(T, pairs::Pair{Int}...) -> Laurent{T}

Build a `Laurent{T}` from `power => coefficient` pairs (coefficients converted to `T`, duplicate
powers summed). Example `laurent(Float64, 0=>1, -1=>-2)` is `1 − 2/r`.
"""
function laurent(::Type{T}, pairs::Pair{Int,<:Any}...) where {T}
    d = Dict{Int,T}()
    for p in pairs
        d[p.first] = get(d, p.first, zero(T)) + T(p.second)
    end
    _prune(Laurent{T}(d))
end
laurent(pairs::Pair{Int,<:Any}...) =
    laurent(float(promote_type(map(p -> typeof(p.second), pairs)...)), pairs...)

lzero(::Type{T}) where {T} = Laurent{T}(Dict{Int,T}())

function _prune(a::Laurent{T}) where {T}
    for (k, v) in collect(a.c)
        v == 0 && delete!(a.c, k)
    end
    a
end

Base.eltype(::Laurent{T}) where {T} = T

function Base.:+(a::Laurent{T}, b::Laurent{U}) where {T,U}
    V = promote_type(T, U)
    d = Dict{Int,V}()
    for (k, v) in a.c; d[k] = get(d, k, zero(V)) + V(v); end
    for (k, v) in b.c; d[k] = get(d, k, zero(V)) + V(v); end
    _prune(Laurent{V}(d))
end
function Base.:-(a::Laurent{T}, b::Laurent{U}) where {T,U}
    V = promote_type(T, U)
    d = Dict{Int,V}()
    for (k, v) in a.c; d[k] = get(d, k, zero(V)) + V(v); end
    for (k, v) in b.c; d[k] = get(d, k, zero(V)) - V(v); end
    _prune(Laurent{V}(d))
end
Base.:-(a::Laurent{T}) where {T} = lscale(-one(T), a)

"""
    lscale(s, a::Laurent) -> Laurent

Scalar multiple `s·a`, element type promoted to `promote_type(T, typeof(s))`.
"""
function lscale(s, a::Laurent{T}) where {T}
    V = promote_type(T, typeof(s))
    d = Dict{Int,V}()
    for (k, v) in a.c; d[k] = V(s) * V(v); end
    _prune(Laurent{V}(d))
end

function Base.:*(a::Laurent{T}, b::Laurent{U}) where {T,U}
    V = promote_type(T, U)
    d = Dict{Int,V}()
    for (ka, va) in a.c, (kb, vb) in b.c
        d[ka + kb] = get(d, ka + kb, zero(V)) + V(va) * V(vb)
    end
    _prune(Laurent{V}(d))
end

"""
    leval(a::Laurent, r) -> value

Evaluate `Σ_k c_k r^k` at `r` (element type promoted with `typeof(r)`).
"""
function leval(a::Laurent{T}, r) where {T}
    V = promote_type(T, typeof(r))
    s = zero(V)
    for (k, v) in a.c
        s += V(v) * V(r)^k
    end
    s
end

"""
    lderiv(a::Laurent) -> Laurent

Exact term-by-term derivative `d/dr Σ c_k r^k = Σ k c_k r^{k-1}`.
"""
function lderiv(a::Laurent{T}) where {T}
    d = Dict{Int,T}()
    for (k, v) in a.c
        k != 0 && (d[k - 1] = get(d, k - 1, zero(T)) + T(k) * v)
    end
    _prune(Laurent{T}(d))
end

"""
    lintegrate(a::Laurent) -> Laurent

An antiderivative `∫ Σ c_k r^k dr = Σ c_k r^{k+1}/(k+1)` (integration constant = 0). Errors if an
`r^{-1}` term is present (its antiderivative is `log r`, outside the Laurent class).
"""
function lintegrate(a::Laurent{T}) where {T}
    d = Dict{Int,T}()
    for (k, v) in a.c
        k == -1 && error("lintegrate: r^{-1} term present (antiderivative is a log, not Laurent)")
        d[k + 1] = get(d, k + 1, zero(T)) + v / T(k + 1)
    end
    _prune(Laurent{T}(d))
end

"""
    lapprox(a::Laurent, b::Laurent; atol, rtol) -> Bool

Coefficient-wise approximate equality over the union of powers (the symbolic-equality predicate the
closed-form gates assert).
"""
function lapprox(a::Laurent{T}, b::Laurent{U}; atol = 0, rtol = atol > 0 ? 0 : sqrt(eps(float(real(promote_type(T, U)))))) where {T,U}
    V = promote_type(T, U)
    ks = union(keys(a.c), keys(b.c))
    all(isapprox(V(get(a.c, k, zero(T))), V(get(b.c, k, zero(U))); atol = atol, rtol = rtol) for k in ks)
end

# ================================================================================================
# ShellProfile{T}
# ================================================================================================
"""
    ShellProfile{T}

A distributional radial profile supported around a shell radius `R`:

    p(r) = inside(r)·Θ(R−r) + outside(r)·Θ(r−R) + dc·δ(r−R) + dpc·δ′(r−R),

with `inside, outside :: Laurent{T}` the interior/exterior smooth (Laurent) branches and `dc, dpc ::
T` the δ and δ′ coefficients at `r=R`. This is the value type the LY thin-shell transport chain
manipulates. See [`psmooth`](@ref) (multiply by a smooth Laurent, LY:831 identity), [`pderiv`](@ref)
(differentiate, Θ-jumps generate δ), and [`pjump`](@ref)/[`pcont`](@ref) (jump/continuity).
"""
struct ShellProfile{T}
    inside::Laurent{T}
    outside::Laurent{T}
    dc::T
    dpc::T
    R::T
end

"""
    shellprofile(inside, outside, R; dc=0, dpc=0) -> ShellProfile

Construct a [`ShellProfile`](@ref) from interior/exterior Laurent branches, shell radius `R`, and
optional δ (`dc`)/δ′ (`dpc`) coefficients.
"""
function shellprofile(inside::Laurent{T}, outside::Laurent{U}, R; dc = 0, dpc = 0) where {T,U}
    V = float(promote_type(T, U, typeof(R), typeof(dc), typeof(dpc)))
    ShellProfile{V}(lscale(one(V), inside), lscale(one(V), outside), V(dc), V(dpc), V(R))
end

"δ(r−R) with unit coefficient (a pure δ profile)."
pdelta(R::T) where {T} = ShellProfile{T}(lzero(T), lzero(T), one(T), zero(T), R)

Base.eltype(::ShellProfile{T}) where {T} = T

function Base.:+(x::ShellProfile{T}, y::ShellProfile{U}) where {T,U}
    V = promote_type(T, U)
    ShellProfile{V}(x.inside + y.inside, x.outside + y.outside, V(x.dc) + V(y.dc),
                    V(x.dpc) + V(y.dpc), V(x.R))
end
function Base.:-(x::ShellProfile{T}, y::ShellProfile{U}) where {T,U}
    V = promote_type(T, U)
    ShellProfile{V}(x.inside - y.inside, x.outside - y.outside, V(x.dc) - V(y.dc),
                    V(x.dpc) - V(y.dpc), V(x.R))
end

"""
    pscale(s, x::ShellProfile) -> ShellProfile

Scalar multiple `s·x`.
"""
function pscale(s, x::ShellProfile{T}) where {T}
    V = promote_type(T, typeof(s))
    ShellProfile{V}(lscale(s, x.inside), lscale(s, x.outside), V(s) * V(x.dc), V(s) * V(x.dpc), V(x.R))
end

"""
    psmooth(φ::Laurent, x::ShellProfile) -> ShellProfile

Multiply the profile by a smooth (Laurent) function `φ(r)`. The Θ branches multiply directly; the
distributional terms use the identities

    φ(r)·δ(r−R)  = φ(R)·δ(r−R),
    φ(r)·δ′(r−R) = φ(R)·δ′(r−R) − φ′(R)·δ(r−R)      (LY:831)

so the new δ coefficient is `φ(R)·dc − φ′(R)·dpc` and the new δ′ coefficient is `φ(R)·dpc`.
"""
function psmooth(φ::Laurent{S}, x::ShellProfile{T}) where {S,T}
    V = promote_type(S, T)
    R = x.R
    φR = leval(φ, R)
    φpR = leval(lderiv(φ), R)
    ShellProfile{V}(φ * x.inside, φ * x.outside,
                    V(φR) * V(x.dc) - V(φpR) * V(x.dpc),
                    V(φR) * V(x.dpc),
                    V(R))
end

"""
    psmooth_vd(valR, derR, x::ShellProfile) -> ShellProfile

Multiply by a smooth function specified only by its value `valR = φ(R)` and derivative `derR = φ′(R)`
at the shell (used when `φ` is rational-with-a-pole, e.g. `r³/f(r)`, and so not a Laurent — the
δ/δ′ identities need only the two boundary numbers). The caller must pass a profile whose Θ branches
are already handled (asserted empty here).
"""
function psmooth_vd(valR, derR, x::ShellProfile{T}) where {T}
    (isempty(x.inside.c) && isempty(x.outside.c)) ||
        error("psmooth_vd expects a pure distributional profile (empty Θ branches)")
    V = promote_type(T, typeof(valR), typeof(derR))
    ShellProfile{V}(lzero(V), lzero(V),
                    V(valR) * V(x.dc) - V(derR) * V(x.dpc),
                    V(valR) * V(x.dpc),
                    V(x.R))
end

"""
    pderiv(x::ShellProfile) -> ShellProfile

Distributional derivative `d/dr`. On the Θ branches it is the term-by-term Laurent derivative; the
Θ-jumps generate a δ:

    d/dr[g⁻Θ(R−r)] = g⁻′Θ(R−r) − g⁻(R)δ(r−R),
    d/dr[g⁺Θ(r−R)] = g⁺′Θ(r−R) + g⁺(R)δ(r−R),
    d/dr[a·δ]       = a·δ′.

A δ′ term would differentiate to δ′′ (outside this type's algebra) — asserted absent.
"""
function pderiv(x::ShellProfile{T}) where {T}
    x.dpc == 0 || error("pderiv of a δ′ term would need δ′′ — not supported")
    R = x.R
    gm = leval(x.inside, R)
    gp = leval(x.outside, R)
    ShellProfile{T}(lderiv(x.inside), lderiv(x.outside), T(gp) - T(gm), x.dc, R)
end

"""
    pjump(x::ShellProfile) -> value

The jump `outside(R) − inside(R)` of the smooth part across the shell (the `[Q]` bracket of LY:699).
"""
pjump(x::ShellProfile{T}) where {T} = leval(x.outside, x.R) - leval(x.inside, x.R)

"""
    pcont(x::ShellProfile) -> Bool

`true` if the smooth part is continuous across the shell (`inside(R) == outside(R)`, no δ/δ′).
"""
function pcont(x::ShellProfile{T}; atol = 0, rtol = atol > 0 ? 0 : sqrt(eps(float(real(T))))) where {T}
    isapprox(leval(x.inside, x.R), leval(x.outside, x.R); atol = atol, rtol = rtol) &&
        x.dc == 0 && x.dpc == 0
end

"""
    peval(x::ShellProfile, r) -> value

Evaluate the smooth (non-distributional) part at `r` (interior branch for `r<R`, exterior for
`r>R`, and the average of the two boundary limits exactly at `r=R`). The δ/δ′ pieces are not
sampled — they enter only the jump/continuity gates.
"""
function peval(x::ShellProfile{T}, r) where {T}
    r < x.R && return leval(x.inside, r)
    r > x.R && return leval(x.outside, r)
    (leval(x.inside, r) + leval(x.outside, r)) / 2
end
