# Per-(s,l,m,a,ω) radial jet: `[R, R', R'', R''', R'''']` at r, built by collapsing the radial
# Teukolsky ODE (no finite differences, precision-preserving). The canonical shared home:
# MR.Radiation re-exports it; ModifiedGravity calls it directly. R from TeukolskyRadial, R′ from
# teukolsky_radial_deriv, higher orders from the ODE and its r-derivatives, λ from Harmonics.
# Complex{T}; BigFloat inherited.

const _RADIAL_JET_NJ = 5   # jet length: value + derivatives up to 4th order

# Memoized λ + radial solves (the expensive steps, hit per ladder component per grid point).
# Keys carry the BigFloat precision so different setprecision regimes never collide.
# Thread-safe (consumers run the (l,m) sum + (r,θ) grid threaded).
const _radial_jet_lambda_cache = Dict{Any,Any}()
const _radial_jet_cache        = Dict{Any,Any}()
const _radial_jet_lock         = ReentrantLock()

"Clear the λ and radial-solution caches used by [`radial_jet`](@ref) (e.g. when changing precision regimes)."
function clear_radial_jet_caches!()
    lock(_radial_jet_lock) do
        empty!(_radial_jet_lambda_cache); empty!(_radial_jet_cache)
    end
    nothing
end

# Internal: the Teukolsky separation constant λ (= notebook `λ[s,ω,l,m]`), precision-keyed cache.
# The eqS separation constant is `A = λ + 2amω − a²ω²`, matching the SWSH convention
# `λ = A − 2mγ + γ²` (γ=aω); so this IS the notebook λ directly.
function _radial_jet_lambda(s::Integer, l::Integer, m::Integer, a::Real, ω::Real)
    key = (s, l, m, a, ω, typeof(a), precision(BigFloat))
    lock(_radial_jet_lock) do
        get!(_radial_jet_lambda_cache, key) do
            SpinWeightedSpheroidalEigenvalue(s, l, m, a*ω)
        end
    end
end

"""
    radial_jet(s, l, m, a, ω, r; bc=:Up, method=:MST, λ=nothing, order=4) -> Vector{Complex}

`[R, R', R'', R''', R'''']` at r. R from `TeukolskyRadial` (unit transmission); R′ from
`teukolsky_radial_deriv`; R″ and up by collapsing the radial Teukolsky ODE eqR
`Δ R'' + (s+1) Δ' R' + V R = 0`, `V = (K² − 2is(r−1)K)/Δ + 4isωr − λ`, `K=(r²+a²)ω − am`, and its
r-derivatives. Exact (no finite differences) and precision-preserving; `Complex{T}` follows the
input eltype so BigFloat flows through. Canonical home for MR.Radiation and ModifiedGravity.
"""
function radial_jet(s::Integer, l::Integer, m::Integer, a::Real, ω::Real, r::Real;
                    bc::Symbol=:Up, method::Symbol=:MST, λ=nothing, order::Int=4)
    λv = λ === nothing ? _radial_jet_lambda(s, l, m, a, ω) : λ
    rkey = (s, l, m, a, ω, method, λv, typeof(a), precision(BigFloat))
    sols = lock(_radial_jet_lock) do
        get!(_radial_jet_cache, rkey) do
            TeukolskyRadial(s, l, m, a, ω; method=method, lambda=λv)
        end
    end
    R = bc === :Up ? sols.Up : sols.In
    C = typeof(complex(float(a*ω*r)))
    R0 = C(R(r))
    R1 = C(teukolsky_radial_deriv(R, r))
    jet = zeros(C, _RADIAL_JET_NJ); jet[1] = R0; jet[2] = R1
    order < 2 && return jet
    aT = C(a); ωT = C(ω); mT = C(m); sT = C(s); λT = C(λv)
    Δ(x)  = x^2 - 2x + aT^2
    dΔ(x) = 2x - 2
    K(x)  = (x^2 + aT^2)*ωT - aT*mT
    dK(x) = 2x*ωT
    # eqR = Δ^{-s} ∂(Δ^{s+1} ∂R) + V R = Δ R'' + (s+1)Δ' R' + V R = 0
    V(x)  = (K(x)^2 - 2*im*sT*(x-1)*K(x))/Δ(x) + 4*im*sT*ωT*x - λT
    N(x)   = K(x)^2 - 2*im*sT*(x-1)*K(x)
    dN(x)  = 2*K(x)*dK(x) - 2*im*sT*(K(x) + (x-1)*dK(x))
    d2N(x) = 2*(dK(x)^2 + K(x)*(2*ωT)) - 2*im*sT*(dK(x) + (dK(x) + (x-1)*2*ωT))
    g(x)   = N(x)/Δ(x)
    dg(x)  = (dN(x)*Δ(x) - N(x)*dΔ(x))/Δ(x)^2
    d2g(x) = (d2N(x) - 2*dg(x)*dΔ(x) - g(x)*C(2))/Δ(x)
    dV2(x)  = dg(x) + 4*im*sT*ωT
    d2V2(x) = d2g(x)
    rx = C(r)
    R2 = -((sT+1)*dΔ(rx)*R1 + V(rx)*R0)/Δ(rx)
    jet[3] = R2
    order < 3 && return jet
    # ∂_r eqR: Δ R''' + (s+2)Δ' R'' + (s+1)Δ'' R' + V' R + V R' = 0
    R3 = -((sT+2)*dΔ(rx)*R2 + (sT+1)*C(2)*R1 + dV2(rx)*R0 + V(rx)*R1)/Δ(rx)
    jet[4] = R3
    order < 4 && return jet
    # ∂²_r eqR: Δ R'''' + (s+3)Δ' R''' + (2s+3)Δ'' R'' + V'' R + 2V' R' + V R'' = 0
    R4 = -((sT+3)*dΔ(rx)*R3 + (2*sT+3)*C(2)*R2 + d2V2(rx)*R0 + 2*dV2(rx)*R1 + V(rx)*R2)/Δ(rx)
    jet[5] = R4
    jet
end
