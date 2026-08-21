# jets.jl — first-order λ-jets: the bookkeeping device that turns "linearise the fluid equations"
# into ordinary arithmetic.
#
# A field of the perturbed problem is written  F = F₀(r) + λ F₁(r) e^{−iωt+imφ}  (45a).  A `Jet`
# stores the pair (F₀, F₁) — the phase factor e^{−iωt+imφ} is IMPLICIT because every O(λ) term
# carries exactly one of it, so it cancels out of every linear equation.  Products, quotients and
# powers then obey ordinary dual-number rules, and the partial derivatives act as
#
#     ∂_t (F₀ , F₁) = (0 , −iω F₁),   ∂_φ (F₀ , F₁) = (0 , +im F₁),
#     ∂_r (F₀ , F₁) = (F₀′ , F₁′),    ∂_z (F₀ , F₁) = (0 , 0),
#
# with ′ = `ddr` of symbols.jl.  This is the Julia counterpart of the route-R2 "G = g + λ h, brute
# force to O(λ)" set-up: no δΓ, δq, δu-identity is ever hand-derived; they all come out of the
# arithmetic.

const SNum = Union{Num,Complex{Num}}

"""
    Jet(b, a)

First-order jet in the perturbation parameter λ: background value `b` = F₀ and O(λ) Fourier
amplitude `a` = F₁ (phase factor implicit).
"""
struct Jet
    b::SNum
    a::SNum
end

Jet(x::Number) = Jet(Num(x), Num(0))
Jet(x::Num) = Jet(x, Num(0))
Jet(x::Complex{Num}) = Jet(x, Num(0))
Jet(x::Jet) = x

Base.zero(::Type{Jet}) = Jet(Num(0), Num(0))
Base.zero(::Jet) = zero(Jet)
Base.one(::Type{Jet}) = Jet(Num(1), Num(0))
Base.one(::Jet) = one(Jet)

Base.:+(x::Jet, y::Jet) = Jet(x.b + y.b, x.a + y.a)
Base.:-(x::Jet, y::Jet) = Jet(x.b - y.b, x.a - y.a)
Base.:-(x::Jet) = Jet(-x.b, -x.a)
Base.:*(x::Jet, y::Jet) = Jet(x.b * y.b, x.b * y.a + x.a * y.b)
Base.:/(x::Jet, y::Jet) = Jet(x.b / y.b, (x.a * y.b - x.b * y.a) / y.b^2)
Base.inv(x::Jet) = Jet(1 / x.b, -x.a / x.b^2)
Base.:^(x::Jet, n::Integer) = Jet(x.b^n, n * x.b^(n - 1) * x.a)

for op in (:+, :-, :*, :/)
    @eval Base.$op(x::Jet, y::Union{Number,Num,Complex{Num}}) = $op(x, Jet(y))
    @eval Base.$op(x::Union{Number,Num,Complex{Num}}, y::Jet) = $op(Jet(x), y)
end

iszero_entry(x::Jet) = iszero_entry(x.b) && iszero_entry(x.a)

"Background (O(λ⁰)) part of a jet."
bg(x::Jet) = x.b
bg(x) = x
"First-order (O(λ)) Fourier amplitude of a jet."
pert(x::Jet) = x.a
pert(x) = zero(x)

"""
    ∂jet(x::Jet, α) -> Jet

Partial derivative ∂_α with α ∈ {1,2,3,4} = (t, r, φ, z) for the Fourier ansatz (45a):
∂_t → −iω, ∂_φ → +im, ∂_r → `ddr`, ∂_z → 0.  Background fields are stationary, axisymmetric and
z-independent (Sec. II.A), so only their r-derivative survives.
"""
function ∂jet(x::Jet, α::Int)
    α == IT && return Jet(Num(0), -im * omega * x.a)
    α == IPHI && return Jet(Num(0), im * m * x.a)
    α == IR && return Jet(ddr(x.b), ddr(x.a))
    α == IZ && return Jet(Num(0), Num(0))
    error("∂jet: bad index $α")
end

"`∂_α` for a plain background scalar (only ∂_r survives)."
∂bg(x, α::Int) = α == IR ? ddr(x) : Num(0)
