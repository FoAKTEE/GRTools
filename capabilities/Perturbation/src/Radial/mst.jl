# MST minimal- (recessive-) solution coefficients f_n. Sasaki–Tagoshi Eq.124: the three-term
# recurrence α_n f_{n+1} + β_n f_n + γ_n f_{n-1} = 0, f_0 = 1, solved by the convergent continued
# fraction. The CF tail reuses the promoted NumericsCore.cf_convergent primitive (no local copy).
# Type-generic accumulator: ComplexF64 for Float64 inputs, Complex{BigFloat} for the high-ℓ small-ω
# extended-precision amplitude path (inferred from the recurrence-ratio type). ν, λ passed in.

"MST three-term recurrence coefficient closures (α, β, γ) — Sasaki–Tagoshi §2.2."
function _mst_recurrence_coeffs(s::Integer, m::Integer, a::Real, omega::Real, nu, lambda)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa
    alpha = n -> (im * eps * kappa * (n + nu + 1 + s + im * eps) * (n + nu + 1 + s - im * eps) *
                  (n + nu + 1 + im * tau)) / ((n + nu + 1) * (2n + 2nu + 3))
    beta = n -> -lambda - s * (s + 1) + (n + nu) * (n + nu + 1) + eps^2 + eps * (eps - m * q) +
                (eps * (eps - m * q) * (s^2 + eps^2)) / ((n + nu) * (n + nu + 1))
    gamma = n -> -((im * eps * kappa * (n + nu - s + im * eps) * (n + nu - s - im * eps) *
                    (n + nu - im * tau)) / ((n + nu) * (2n + 2nu - 1)))
    return alpha, beta, gamma
end

"""
    _mst_fn_array(s, m, a, omega, nu, lambda, N) -> (Dict n=>f_n, alpha, beta, gamma)

MST minimal-solution coefficients `f_n` for n = −N..N, `f_0 = 1`. Upward ratio
`f_n/f_{n-1} = cf_convergent(−α(k−1)γ(k), β(k), n)/α(n−1)`; downward via the index reflection
`i → 2n−i`. Returns the coefficient table + the recurrence closures (for residual checks).
"""
function _mst_fn_array(s::Integer, m::Integer, a::Real, omega::Real, nu, lambda, N::Integer)
    alpha, beta, gamma = _mst_recurrence_coeffs(s, m, a, omega, nu, lambda)
    # Type-generic accumulator inferred from the recurrence-ratio type: ComplexF64 for Float64
    # inputs, Complex{BigFloat} for extended precision (the high-ℓ small-ω amplitude path).
    T = typeof(complex(beta(0) / beta(0)) + complex(nu) + complex(lambda) + complex(2 * omega))
    f = Dict{Int,T}(0 => one(T))
    for nf in 1:N
        R = cf_convergent(k -> -alpha(k - 1) * gamma(k), k -> beta(k), nf) / alpha(nf - 1)
        f[nf] = f[nf-1] * R
    end
    for nf in -1:-1:-N
        R = cf_convergent(k -> -alpha(2nf - k) * gamma(2nf - k + 1), k -> beta(2nf - k), nf) / gamma(nf + 1)
        f[nf] = f[nf+1] * R
    end
    return f, alpha, beta, gamma
end

"Single MST coefficient `f_n`; for many `n` prefer `_mst_fn_array`."
_mst_fn(s::Integer, m::Integer, a::Real, omega::Real, nu, lambda, n::Integer) =
    _mst_fn_array(s, m, a, omega, nu, lambda, abs(n))[1][n]
