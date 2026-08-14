# Spin-weighted spherical and spheroidal harmonics: eigenvalues, S functions, limits.
#
# Run from the repository root:
#     julia --project=. examples/02_spheroidal_harmonics.jl
#
# Conventions: spheroidicity gamma = a*omega; eigenvalue lambda with
# lambda(gamma=0) = l(l+1) - s(s+1); fields ~ e^{i(m phi - omega t)}.

using Harmonics
using Printf

s, l, m = -2, 2, 2

# Spherical limit: the eigenvalue collapses to the exact l(l+1) - s(s+1).
λ0 = SpinWeightedSpheroidalEigenvalue(s, l, m, 0.0)
@printf "lambda(s=%d, l=%d, m=%d, gamma=0)   = %.12f   (exact %d)\n" s l m λ0 (l * (l + 1) - s * (s + 1))

# Two independent backends agree at finite spheroidicity.
γ = 0.35
λse = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:SphericalExpansion)
λlv = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:Leaver)
@printf "lambda(gamma=%.2f): SphericalExpansion = %.12f, Leaver = %.12f, |diff| = %.2e\n" γ λse λlv abs(λse - λlv)

# The spheroidal S function evaluates at (theta, phi); gamma -> 0 recovers sYlm.
S = SpinWeightedSpheroidalHarmonicS(s, l, m, γ)
S0 = SpinWeightedSpheroidalHarmonicS(s, l, m, 0.0)
θ, φ = 1.1, 0.9
@printf "\nS_{%d,%d,%d}(theta=%.1f, phi=%.1f; gamma=%.2f) = %.10f%+.10fim\n" s l m θ φ γ real(S(θ, φ)) imag(S(θ, φ))
Y = SpinWeightedSphericalHarmonicY(s, l, m, θ, φ)
@printf "gamma->0 limit vs sYlm:  |S0 - Y| = %.2e\n" abs(S0(θ, φ) - Y)

# The phi-derivative identity d/dphi = i*m.
dS = swsh_derivative(S, 0, 1, θ, φ)
@printf "swsh_derivative (0,1) vs im*m*S: |diff| = %.2e\n" abs(dS - im * m * S(θ, φ))

# Complex spheroidicity (quasinormal-mode applications) is first-class.
γc = 0.5 + 0.2im
λc = SpinWeightedSpheroidalEigenvalue(s, l, m, γc)
println("\nComplex gamma = ", γc, ":  lambda = ", λc)

# BigFloat gamma inherits the working precision.
setprecision(BigFloat, 256) do
    λb = SpinWeightedSpheroidalEigenvalue(s, l, m, big"0.35")
    @printf "BigFloat-256 eigenvalue agrees with Float64 to %.2e\n" Float64(abs(λb - λse))
end
