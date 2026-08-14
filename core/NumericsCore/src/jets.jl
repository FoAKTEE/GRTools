# 1D derivative-tower ("jet") algebra. A jet stores DERIVATIVE VALUES [f, f′, …, f^(N)]
# at the expansion point — NOT Taylor coefficients f^(k)/k! — so products carry the Leibniz
# binomial C(k,j) (a plain Cauchy convolution would omit it). Length- and eltype-generic; Base-only.

"""
    jadd(f, g)

Elementwise sum of two jets (derivative towers), `f .+ g`. Linear, so identical for the
derivative-value and Taylor-coefficient conventions. Eltype promotes from `f`, `g`.
"""
jadd(f, g) = f .+ g

"""
    jscale(s, f)

Scale a jet by a scalar `s` (any `Number`), `s .* f`. Eltype-generic; `s` may be `im`, a rational,
`-one(CT)`, etc.
"""
jscale(s, f) = s .* f

"""
    jmul(f, g) -> Vector

Leibniz product of two jets: `(f·g)^(k) = Σ_j C(k,j) f^(j) g^(k−j)`. Output length `= length(f)`
(requires `length(g) ≥ length(f)`). The binomial factor `C(k,j)` is exactly what distinguishes the
stored-derivative-value tower from a plain Cauchy convolution of Taylor coefficients.
"""
jmul(f, g) = (n = length(f) - 1; [sum(binomial(k, j) * f[j+1] * g[k-j+1] for j in 0:k) for k in 0:n])

"""
    jderiv(f) -> Vector

Differentiate a jet by one derivative order: `(f′)^(k) = f^(k+1)`. Implemented as a shift-down with
the top slot zero-padded, so the output has the SAME length as the input (the highest-order
derivative is truncated). The bare shift — no `(k+1)` factorial rescaling — is correct precisely
because the tower stores derivative VALUES (a Taylor-coefficient tower would rescale term k by k+1).
"""
jderiv(f) = vcat(f[2:end], [zero(eltype(f))])

"""
    jinv(f) -> Vector

Reciprocal jet `1/f` about the expansion point, solved order-by-order from the Leibniz relation
`Σ_j C(k,j) f^(j) (1/f)^(k−j) = δ_{k0}`: `g[1] = 1/f[1]`, and for `k ≥ 1`,
`g^(k) = −(1/f[1]) Σ_{j=1..k} C(k,j) f^(j) g^(k−j)`. Requires `f[1] ≠ 0`. Eltype and length are
taken from `f`, so `Complex{BigFloat}` is preserved at full precision.
"""
function jinv(f)
    n = length(f) - 1
    g = zeros(eltype(f), n + 1)
    g[1] = one(eltype(f)) / f[1]
    for k in 1:n
        g[k+1] = -sum(binomial(k, j) * f[j+1] * g[k-j+1] for j in 1:k) / f[1]
    end
    g
end

"""
    constjet(s, N::Integer, ::Type{CT}) where {CT} -> Vector{CT}

Constant-function jet of value `s`: `[CT(s), 0, 0, …, 0]` of length `N+1`. All derivatives above
order 0 vanish (correct for a constant in either convention). The only primitive whose output
element type is pinned explicitly by the `::Type{CT}` argument rather than inferred from an input.
"""
constjet(s, N::Integer, ::Type{CT}) where {CT} = vcat([CT(s)], zeros(CT, N))

# ── NTuple methods: @generated fully-unrolled (zero-allocation) for hot tower kernels. ──────────
# CONSTRAINT: the tuple `jinv` uses `inv(f[1])` while the Vector `jinv` uses `one(T)/f[1]` — they
# differ at the last ulp / signed-zero for complex eltype (test/test_jets.jl); the tower kernels'
# characterized output is pinned to the `inv` form. Do not unify.
_jet_term(c, a, b) = c == 1 ? :($a * $b) : :($c * $a * $b)

"""
    jmul(f::NTuple{N}, g::NTuple{M}) -> NTuple{N}

Allocation-free `NTuple` Leibniz product (requires `M ≥ N`), identical recurrence and accumulation
order to the `Vector` [`jmul`](@ref) — bit-exact to it on equal-length input. `@generated` so the whole
product is unrolled at compile time (register-resident, no heap).
"""
@generated function jmul(f::NTuple{N,T1}, g::NTuple{M,T2}) where {N,M,T1,T2}
    N <= M || return :(throw(ArgumentError("jmul(NTuple): need length(g) ≥ length(f)")))
    rows = map(0:N-1) do k
        s = _jet_term(binomial(k, 0), :(f[1]), :(g[$(k+1)]))
        for j in 1:k
            s = :($s + $(_jet_term(binomial(k, j), :(f[$(j+1)]), :(g[$(k-j+1)]))))
        end
        s
    end
    Expr(:tuple, rows...)
end

"""
    jinv(f::NTuple{N}) -> NTuple{N}

Allocation-free `NTuple` reciprocal jet `1/f`, same order-by-order Leibniz solve as the `Vector`
[`jinv`](@ref) but following the `inv(f[1])` reciprocal convention of the former MetricReconstruction
`j4inv` (bit-exact to it; see the section note re: the ulp/signed-zero difference vs the `Vector` form).
`@generated` (fully unrolled, no heap). Requires `f[1] ≠ 0`.
"""
@generated function jinv(f::NTuple{N,T}) where {N,T}
    ex = Expr(:block, :(g1 = inv(f[1])))
    gsyms = Any[:g1]                                     # gsyms[i] = g_i (derivative order i-1, 0-based)
    for k in 1:N-1                                       # produce g_{k+1} (0-based order k)
        s = _jet_term(binomial(k, 1), :(f[2]), gsyms[k]) # j=1 term: C(k,1) f^(1) g^(k-1)
        for j in 2:k
            s = :($s + $(_jet_term(binomial(k, j), :(f[$(j+1)]), gsyms[k-j+1])))
        end
        gk = Symbol(:g, k + 1)
        push!(ex.args, :($gk = -($s) * g1))
        push!(gsyms, gk)
    end
    push!(ex.args, Expr(:tuple, gsyms...))
    ex
end
