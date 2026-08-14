"""
    MetricReconstruction.Radiation — separable jet algebra + GHP/Teukolsky ladder operators

The 2D separable-mode representation underlying the CCK reconstruction. Every field flowing through
the GHP/Teukolsky operator chains is a SUM of separable terms `Σ_k ρ_k(r) σ_k(θ)`, with `ρ_k` a
*radial jet* `[ρ, ρ', ρ'', ρ''', ρ'''']` and `σ_k` an *angular jet* `[σ, σ', …, σ'''']` (length
`NJ=5`, derivatives to 4th order — the deepest chain composes four derivatives on one factor, so the
representation is EXACT, no finite differences, no CAS). A radial ladder `Dn` raises r-derivative
order on the radial jets; an angular ladder `Ln` on the angular jets. Type-generic via `Complex{T}`.

`RTJet` is the companion analytic `(r,θ)` 2-jet `(v, ∂_r, ∂_θ, ∂_rr, ∂_θθ, ∂_rθ)` used by the
linearized-Einstein residual: background scalars (ζ, Δ, Σ, …) and H-components are built as RTJets,
Leibniz-combined, giving exact analytic derivatives.

GHP/Teukolsky ladder (`M=1`):
  `Dn[n,ω,m] f = f' − iK/Δ f + 2n(r−1)/Δ f`   (`sgn=−1`; `sgn=+1` flips the K sign ⇒ adjoint),
  `Ln[n,ω,m] f = f_θ + Q f + n cotθ f`         (`sgn=+1`; `sgn=−1` flips the Q sign ⇒ adjoint),
with `K=(r²+a²)ω−am`, `Q=−aω sinθ + m cscθ`; `bar=true` is the GHP overbar (i → −i, ω → ω̄).
"""
const NJ = 5   # jet length: value + derivatives up to 4th order

# Canonical 1D derivative-tower jet algebra now lives in NumericsCore.jets (Move A). The underscore
# `_jet_*` names below are thin length-fixed adapters over these length-generic kernels.
using NumericsCore: jmul, jadd, jderiv, jinv

# ───────────────────────────── separable mode field ──────────────────────────────────────────────
"""
    SepField{C}

A separable mode field `Σ_k radial[k] ⊗ angular[k]`, each entry a length-`NJ` jet of `Complex{C}`.
Operators map `SepField → SepField`.
"""
struct SepField{C}
    radial::Vector{Vector{C}}    # each: [ρ, ρ', ρ'', ρ''', ρ'''']  (∂_r)
    angular::Vector{Vector{C}}   # each: [σ, σ', σ'', σ''', σ'''']  (∂_θ)
end

"Build a one-term separable field from a radial jet `rj` and angular jet `aj`."
sepfield(rj::Vector{C}, aj::Vector{C}) where {C} = SepField{C}([copy(rj)], [copy(aj)])

"Field value `f(r,θ) = Σ_k ρ_k σ_k` (the jet [1] entries)."
function value(F::SepField{C}) where {C}
    acc = zero(C)
    @inbounds for k in eachindex(F.radial)
        acc += F.radial[k][1] * F.angular[k][1]
    end
    acc
end

"Sum two separable fields (concatenate terms)."
Base.:+(A::SepField{C}, B::SepField{C}) where {C} =
    SepField{C}(vcat(A.radial, B.radial), vcat(A.angular, B.angular))

"Scale a separable field by a pure number (applied to the radial jets)."
Base.:*(c::Number, F::SepField{C}) where {C} =
    SepField{C}([c .* rj for rj in F.radial], [copy(aj) for aj in F.angular])

# ───────────────────────────── Taylor (value+derivs) jet primitives ──────────────────────────────
# CANONICAL HOME = NumericsCore.jets (Move A). Every Radiation call site passes length-NJ(=5) jets,
# so these thin adapters over the length-generic NumericsCore kernels are equivalent to the former
# fixed-5 copy: jmul/jadd/jderiv agree BIT-EXACT for all eltypes; jinv agrees bit-exact for real
# eltype, and for complex eltype the reciprocal's leading term is one(C)/u[1] (Lorenz convention) vs
# the former inv(u[1]) — Julia inv(::Complex) differs by ≤1 ulp, absorbed by the tolerance-based
# Radiation gates. Characterization oracle: core/NumericsCore/test/test_jets.jl.
_jet_mul(u::Vector{C}, v::Vector{C}) where {C} = jmul(u, v)
_jet_add(u::Vector{C}, v::Vector{C}) where {C} = jadd(u, v)
_jet_deriv(u::Vector{C}) where {C} = jderiv(u)
_jet_inv(u::Vector{C}) where {C} = jinv(u)

# ───────────────────────────── analytic coefficient jets ─────────────────────────────────────────
"Jet of `Q(θ) = −aω sinθ + m cscθ` (analytic derivatives to 4th order)."
function _Q_jet(::Type{C}, a, θ, ω, m) where {C}
    sθ = sin(θ); cθ = cos(θ)
    csc = 1/sθ; cot = cθ/sθ
    c0 = csc
    c1 = -csc*cot
    c2 =  csc*(cot^2 + csc^2)
    c3 = -csc*cot*(cot^2 + 5*csc^2)
    c4 =  csc*(cot^4 + 18*cot^2*csc^2 + 5*csc^4)
    s0 = sθ; s1 = cθ; s2 = -sθ; s3 = -cθ; s4 = sθ
    Qj = zeros(C, NJ)
    Qj[1] = -a*ω*s0 + m*c0
    Qj[2] = -a*ω*s1 + m*c1
    Qj[3] = -a*ω*s2 + m*c2
    Qj[4] = -a*ω*s3 + m*c3
    Qj[5] = -a*ω*s4 + m*c4
    Qj
end

"Jet of `cotθ` (analytic derivatives to 4th order)."
function _cot_jet(::Type{C}, θ) where {C}
    sθ = sin(θ); cθ = cos(θ); csc2 = 1/sθ^2; cot = cθ/sθ
    j = zeros(C, NJ)
    j[1] = cot
    j[2] = -csc2
    j[3] = 2*csc2*cot
    j[4] = -2*csc2*(csc2 + 2*cot^2)
    j[5] = 8*csc2*cot*(2*csc2 + cot^2)
    j
end

"""
    _radcoef_jet(C, sgn, n, a, r, ω, m; bar=false)

Radial-ladder coefficient jet `c(r) = sgn·i·K/Δ + 2n(r−1)/Δ`, `K=(r²+a²)ω−am`. `bar=true` is the
GHP overbar (i → −i); the caller passes ω̄ for complex ω (so `K → K̄`).
"""
function _radcoef_jet(::Type{C}, sgn::Int, n::Int, a, r, ω, m; bar::Bool=false) where {C}
    Kj  = zeros(C, NJ); Kj[1] = (r^2+a^2)*ω - a*m; Kj[2] = 2*r*ω; Kj[3] = 2*ω
    Δj  = zeros(C, NJ); Δj[1] = r^2 - 2r + a^2; Δj[2] = 2r - 2; Δj[3] = C(2)
    invΔ = _jet_inv(Δj)
    rm1 = zeros(C, NJ); rm1[1] = C(r - 1); rm1[2] = one(C)
    isign = bar ? -1 : 1
    term1 = (sgn * isign * im) .* _jet_mul(Kj, invΔ)
    term2 = (2*n) .* _jet_mul(rm1, invΔ)
    term1 .+ term2
end

"Angular-ladder coefficient jet `cθ(θ) = sgn·Q + n cotθ`."
function _angcoef_jet(::Type{C}, sgn::Int, n::Int, a, θ, ω, m) where {C}
    Qj = _Q_jet(C, a, θ, ω, m)
    cotj = (n == 0) ? zeros(C, NJ) : (C(n) .* _cot_jet(C, θ))
    (sgn .* Qj) .+ cotj
end

"Δ jet `[Δ, Δ', Δ'', 0, 0]`."
function _Delta_jet(::Type{C}, a, r) where {C}
    j = zeros(C, NJ); j[1] = r^2 - 2r + a^2; j[2] = 2r - 2; j[3] = C(2); j
end

"`ζ̄⁻¹ = 1/(r + i a cosθ)` as an r-jet at fixed θ."
function _zbarinv_rjet(::Type{C}, a, r, θ) where {C}
    z = zeros(C, NJ); z[1] = complex(r, a*cos(θ)); z[2] = one(C); _jet_inv(z)
end

"`ζ⁻¹ = 1/(r − i a cosθ)` as an r-jet at fixed θ."
function _zinv_rjet(::Type{C}, a, r, θ) where {C}
    z = zeros(C, NJ); z[1] = complex(r, -a*cos(θ)); z[2] = one(C); _jet_inv(z)
end

# ───────────────────────────── ladder applications ───────────────────────────────────────────────
"""
    Dn(sgn, n, F; a, r, ω, m, bar=false) -> SepField

Radial ladder. `sgn=−1` → `Dn`, `sgn=+1` → `Dnadj`; `bar=true` → barred variant. Acts on each radial
jet `ρ ↦ ρ' + c·ρ` (Leibniz to order 4); angular jets untouched.
"""
function Dn(sgn::Int, n::Int, F::SepField{C}; a, r, ω, m, bar::Bool=false) where {C}
    cj = _radcoef_jet(C, sgn, n, a, r, ω, m; bar=bar)
    newr = Vector{Vector{C}}(undef, length(F.radial))
    @inbounds for k in eachindex(F.radial)
        newr[k] = _jet_deriv(F.radial[k]) .+ _jet_mul(cj, F.radial[k])
    end
    SepField{C}(newr, [copy(aj) for aj in F.angular])
end

"""
    Ln(sgn, n, F; a, θ, ω, m) -> SepField

Angular ladder. `sgn=+1` → `Ln`, `sgn=−1` → `Lnadj`. Acts on each angular jet `σ ↦ σ' + cθ·σ`.
"""
function Ln(sgn::Int, n::Int, F::SepField{C}; a, θ, ω, m) where {C}
    cj = _angcoef_jet(C, sgn, n, a, θ, ω, m)
    newa = Vector{Vector{C}}(undef, length(F.angular))
    @inbounds for k in eachindex(F.angular)
        newa[k] = _jet_deriv(F.angular[k]) .+ _jet_mul(cj, F.angular[k])
    end
    SepField{C}([copy(rj) for rj in F.radial], newa)
end

"Multiply each radial jet by an r-coefficient jet (Leibniz)."
mul_rjet(F::SepField{C}, cj::Vector{C}) where {C} =
    SepField{C}([_jet_mul(cj, rj) for rj in F.radial], [copy(aj) for aj in F.angular])

"Multiply each angular jet by a θ-coefficient jet (Leibniz)."
mul_ajet(F::SepField{C}, cj::Vector{C}) where {C} =
    SepField{C}([copy(rj) for rj in F.radial], [_jet_mul(cj, aj) for aj in F.angular])

# ───────────────────────────── RTJet: analytic (r,θ) 2-jet ───────────────────────────────────────
"Analytic `(r,θ)` 2-jet `(v, ∂_r, ∂_θ, ∂_rr, ∂_θθ, ∂_rθ)`, exact to 2nd order."
struct RTJet{C}
    v::C; r::C; θ::C; rr::C; θθ::C; rθ::C
end
RTJet{C}() where {C} = RTJet{C}(zero(C),zero(C),zero(C),zero(C),zero(C),zero(C))
"Constant `(r,θ)` 2-jet: value `c`, all spatial derivatives zero (the `RTJet` embedding of a scalar)."
rtconst(c::C) where {C} = RTJet{C}(c,zero(C),zero(C),zero(C),zero(C),zero(C))

Base.:+(a::RTJet{C}, b::RTJet{C}) where {C} = RTJet{C}(a.v+b.v, a.r+b.r, a.θ+b.θ, a.rr+b.rr, a.θθ+b.θθ, a.rθ+b.rθ)
Base.:-(a::RTJet{C}, b::RTJet{C}) where {C} = RTJet{C}(a.v-b.v, a.r-b.r, a.θ-b.θ, a.rr-b.rr, a.θθ-b.θθ, a.rθ-b.rθ)
Base.:-(a::RTJet{C}) where {C} = RTJet{C}(-a.v, -a.r, -a.θ, -a.rr, -a.θθ, -a.rθ)
Base.:*(c::Number, a::RTJet{C}) where {C} = RTJet{C}(c*a.v, c*a.r, c*a.θ, c*a.rr, c*a.θθ, c*a.rθ)

"Leibniz product of two `(r,θ)` 2-jets (exact to 2nd order)."
function Base.:*(a::RTJet{C}, b::RTJet{C}) where {C}
    v  = a.v*b.v
    r  = a.r*b.v + a.v*b.r
    θ  = a.θ*b.v + a.v*b.θ
    rr = a.rr*b.v + 2*a.r*b.r + a.v*b.rr
    θθ = a.θθ*b.v + 2*a.θ*b.θ + a.v*b.θθ
    rθ = a.rθ*b.v + a.r*b.θ + a.θ*b.r + a.v*b.rθ
    RTJet{C}(v, r, θ, rr, θθ, rθ)
end

"Reciprocal `(r,θ)` 2-jet `1/a` (exact to 2nd order)."
function rtinv(a::RTJet{C}) where {C}
    iv = inv(a.v)
    v  = iv
    r  = -a.r*iv^2
    θ  = -a.θ*iv^2
    rr = (2*a.r^2*iv - a.rr)*iv^2
    θθ = (2*a.θ^2*iv - a.θθ)*iv^2
    rθ = (2*a.r*a.θ*iv - a.rθ)*iv^2
    RTJet{C}(v, r, θ, rr, θθ, rθ)
end
Base.:/(a::RTJet{C}, b::RTJet{C}) where {C} = a*rtinv(b)
"Integer power of a `(r,θ)` 2-jet by repeated Leibniz product (`n ≥ 1`)."
Base.:^(a::RTJet{C}, n::Integer) where {C} = (p = a; for _ in 2:n; p = p * a; end; p)

"Extract the `(r,θ)` 2-jet of a SepField `Σ_k ρ_k(r) σ_k(θ)` from its stored derivative-jets."
function sepfield_rtjet(F::SepField{C}) where {C}
    v=zero(C); r=zero(C); θ=zero(C); rr=zero(C); θθ=zero(C); rθ=zero(C)
    @inbounds for k in eachindex(F.radial)
        ρ = F.radial[k]; σ = F.angular[k]
        v  += ρ[1]*σ[1]
        r  += ρ[2]*σ[1]
        θ  += ρ[1]*σ[2]
        rr += ρ[3]*σ[1]
        θθ += ρ[1]*σ[3]
        rθ += ρ[2]*σ[2]
    end
    RTJet{C}(v, r, θ, rr, θθ, rθ)
end

# RTJet base builders for the background scalars (M=1).
_rt_r(::Type{C}, r) where {C}   = RTJet{C}(C(r), one(C), zero(C), zero(C), zero(C), zero(C))
_rt_sin(::Type{C}, θ) where {C} = RTJet{C}(C(sin(θ)), zero(C), C(cos(θ)), zero(C), C(-sin(θ)), zero(C))
_rt_cos(::Type{C}, θ) where {C} = RTJet{C}(C(cos(θ)), zero(C), C(-sin(θ)), zero(C), C(-cos(θ)), zero(C))
"`sin 2θ` as RTJet."
_rt_sin2(::Type{C}, θ) where {C} = RTJet{C}(C(sin(2θ)), zero(C), C(2cos(2θ)), zero(C), C(-4sin(2θ)), zero(C))
"`ζ = r − i a cosθ` as RTJet."
_rt_zeta(::Type{C}, a, r, θ) where {C}    = _rt_r(C, r) - rtconst(C(im*a))*_rt_cos(C, θ)
"`ζ̄ = r + i a cosθ` as RTJet."
_rt_zetabar(::Type{C}, a, r, θ) where {C} = _rt_r(C, r) + rtconst(C(im*a))*_rt_cos(C, θ)
"`Δ = r² − 2r + a²` as RTJet."
_rt_Delta(::Type{C}, a, r) where {C} = RTJet{C}(C(r^2-2r+a^2), C(2r-2), zero(C), C(2), zero(C), zero(C))
"`Σ = r² + a² cos²θ` as RTJet."
function _rt_Sigma(::Type{C}, a, r, θ) where {C}
    c = _rt_cos(C, θ); _rt_r(C,r)*_rt_r(C,r) + rtconst(C(a^2))*(c*c)
end
