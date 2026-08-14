# Radiation mode ingredients: λ (teuk_lambda), angular jets, and the Perturbation-owned radial_jet
# (re-exported so `Radiation.radial_jet` works unchanged). No new ODE solve here.

using Harmonics: SpinWeightedSpheroidalEigenvalue, SpinWeightedSpheroidalHarmonicS, swsh_derivative
using Perturbation: radial_jet, clear_radial_jet_caches!

# λ memoization; the key carries the BigFloat precision so setprecision regimes never collide.
const _LAMBDA_CACHE = Dict{Any,Any}()
const _CACHE_LOCK = ReentrantLock()

"Clear the local λ cache AND the relocated radial_jet caches in Perturbation (e.g. when changing precision regimes)."
function clear_mode_caches!()
    lock(_CACHE_LOCK) do
        empty!(_LAMBDA_CACHE)
    end
    clear_radial_jet_caches!()
    nothing
end

"""
    teuk_lambda(s, l, m, a, ω) -> Complex

Teukolsky λ (= notebook `λ[s,ω,l,m]`). The eqS separation constant is `A = λ + 2amω − a²ω²`,
matching the SWSH convention `λ = A − 2mγ + γ²` (γ=aω); so this IS the notebook λ directly. For
BigFloat inputs the eigenvalue carries precision through the GrTool Harmonics surface.
"""
function teuk_lambda(s::Integer, l::Integer, m::Integer, a::Real, ω::Real)
    key = (s, l, m, a, ω, typeof(a), precision(BigFloat))
    lock(_CACHE_LOCK) do
        get!(_LAMBDA_CACHE, key) do
            SpinWeightedSpheroidalEigenvalue(s, l, m, a*ω)
        end
    end
end

"""
    angular_jet(s, l, m, a, ω, θ; order=4) -> Vector{Complex}

`[S, S_θ, S_θθ, S_θθθ, S_θθθθ]` at θ (φ=0; the e^{imφ} factor lives in the plus/minus assembly). S
and its θ-derivatives from `SpinWeightedSpheroidalHarmonicS`/`swsh_derivative` (SphericalExpansion).
"""
function angular_jet(s::Integer, l::Integer, m::Integer, a::Real, ω::Real, θ::Real; order::Int=4)
    γ = a*ω
    S = SpinWeightedSpheroidalHarmonicS(s, l, m, Float64(γ); method=:SphericalExpansion)
    C = typeof(complex(float(γ*θ)))
    jet = zeros(C, NJ)
    for d in 0:order
        jet[d+1] = C(swsh_derivative(S, d, 0, Float64(θ), 0.0))
    end
    jet
end
