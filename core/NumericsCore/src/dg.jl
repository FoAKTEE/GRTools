# 1D nodal discontinuous-Galerkin spectral-element core — port of DG-Teulolsky/src/DGCore.jl
# (LGL reference element, differentiation/mass matrices, multi-domain mesh with forced interfaces,
# barycentric Lagrange interpolation) + DG-Teulolsky/src/TimeStepping.jl (RK4 + CFL helper).
# Conventions follow Hesthaven & Warburton, *Nodal DG Methods* (2008). Problem-agnostic: no
# Teukolsky/physics enters here.
#
# The original's FastGaussQuadrature dependency is ELIMINATED (§7.7 sink budget: LinearAlgebra
# ONLY): the Legendre-Gauss-Lobatto rule is computed self-contained — endpoints ±1, interior nodes
# by Newton iteration on the roots of P'_N (Chebyshev-Gauss-Lobatto seeds; P_N, P'_N by the
# three-term recurrence; d/dx P'_N from the Legendre ODE (1-x²)P″_N = 2xP'_N − N(N+1)P_N), weights
# w_i = 2/(N(N+1) P_N(x_i)²). Everything is type-generic {T}, so BigFloat(256+) inputs give
# genuine BigFloat quadrature/operators — a NEW capability the original (Float64-only
# `gausslobatto`) lacked. Characterized against the original FastGaussQuadrature-based
# RefElement(N), N = 2..16, baked in test/fixtures/dg_lgl_ref.jl (captured once by
# scratchpad/capture_dg_lgl_ref.jl).

# ----------------------------------------------------------------------------
# Self-contained Legendre-Gauss-Lobatto rule (replaces FastGaussQuadrature.gausslobatto).
# ----------------------------------------------------------------------------

"""
    _legendre_PdP(N, x) -> (P_N(x), P'_N(x))

Legendre polynomial `P_N` and its derivative at `x`, by the stable three-term recurrences
`(n+1)P_{n+1} = (2n+1)xP_n − nP_{n−1}` and `(n+1)P'_{n+1} = (2n+1)(P_n + xP'_n) − nP'_{n−1}`.
Type-generic in `x::T` (BigFloat-clean).
"""
function _legendre_PdP(N::Int, x::T) where {T<:AbstractFloat}
    P, dP = one(T), zero(T)
    N == 0 && return P, dP
    Pm, dPm = P, dP          # P_{n-1}, P'_{n-1}
    P, dP = x, one(T)        # P_1, P'_1
    for n in 1:(N - 1)
        Pn = ((2n + 1) * x * P - n * Pm) / (n + 1)
        dPn = ((2n + 1) * (P + x * dP) - n * dPm) / (n + 1)
        Pm, dPm = P, dP
        P, dP = Pn, dPn
    end
    return P, dP
end

"""
    _lgl_nodes_weights(T, N) -> (r, w)

The `Np = N+1`-point Legendre-Gauss-Lobatto quadrature rule on `[-1,1]` in precision `T`, exact
for polynomials of degree `≤ 2N−1`. Endpoints `±1`; interior nodes are the roots of `P'_N`, found
by Newton iteration from Chebyshev-Gauss-Lobatto seeds `−cos(πk/N)` (the Newton step uses
`P″_N = (2xP'_N − N(N+1)P_N)/(1−x²)` from the Legendre ODE); weights `w_i = 2/(N(N+1)P_N(x_i)²)`.
Self-contained substitute for `FastGaussQuadrature.gausslobatto(N+1)` that additionally inherits
BigFloat precision (converges to `eps(T)`, e.g. ~1e-77 residuals at `BigFloat(256)`).
"""
function _lgl_nodes_weights(::Type{T}, N::Int) where {T<:AbstractFloat}
    N >= 1 || throw(ArgumentError("polynomial degree N must be ≥ 1"))
    Np = N + 1
    r = Vector{T}(undef, Np)
    r[1] = -one(T)
    r[Np] = one(T)
    NN1 = T(N) * (N + 1)
    for k in 1:(N - 1)
        x = -cospi(T(k) / N)                       # CGL seed (exactly 0 at the symmetric node)
        for _ in 1:100                             # quadratic convergence: ~6 its at 256 bits
            P, dP = _legendre_PdP(N, x)
            δ = dP * (1 - x^2) / (2x * dP - NN1 * P)   # Newton on P'_N
            x -= δ
            abs(δ) <= 4 * eps(T) * max(one(T), abs(x)) && break
        end
        r[k + 1] = x
    end
    w = Vector{T}(undef, Np)
    for i in 1:Np
        P, _ = _legendre_PdP(N, r[i])
        w[i] = 2 / (NN1 * P^2)
    end
    return r, w
end

# ----------------------------------------------------------------------------
# Orthonormal Legendre basis on [-1,1] and its derivative, by recurrence.
# P̃_n(x) = sqrt((2n+1)/2) P_n(x)  ⇒  ∫_{-1}^{1} P̃_n P̃_m dx = δ_{nm}.
# ----------------------------------------------------------------------------
"orthonormal-Legendre Vandermonde `V[i,n+1] = P̃_n(x_i)` and its derivative `Vr` (DGCore
`legendre_vandermonde`, ported verbatim; type-generic)."
function _legendre_vandermonde(x::AbstractVector{T}, N::Int) where {T<:AbstractFloat}
    nx = length(x)
    Np = N + 1
    P = zeros(T, nx, Np)      # unnormalized P_n
    dP = zeros(T, nx, Np)     # P_n'
    @inbounds for i in 1:nx
        xi = x[i]
        P[i, 1] = one(T)
        dP[i, 1] = zero(T)
        if N >= 1
            P[i, 2] = xi
            dP[i, 2] = one(T)
        end
        for n in 1:(N - 1)
            # (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
            P[i, n + 2] = ((2n + 1) * xi * P[i, n + 1] - n * P[i, n]) / (n + 1)
            dP[i, n + 2] =
                ((2n + 1) * (P[i, n + 1] + xi * dP[i, n + 1]) - n * dP[i, n]) / (n + 1)
        end
    end
    V = similar(P)
    Vr = similar(dP)
    for n in 0:N
        s = sqrt(T(2n + 1) / 2)
        @views V[:, n + 1] .= s .* P[:, n + 1]
        @views Vr[:, n + 1] .= s .* dP[:, n + 1]
    end
    return V, Vr
end

"""
    RefElement{T}(N)

Degree-`N` nodal DG reference element on `[-1,1]` with `Np = N+1` LGL nodes (Hesthaven &
Warburton conventions; DGCore port).

Fields: `r` (nodes), `w` (LGL weights), `V` (orthonormal Legendre Vandermonde),
`D` (differentiation matrix, d/dr on `[-1,1]`), `M` (mass matrix `(VVᵀ)⁻¹`), `invM` (`= VVᵀ`).

Construct with `RefElement(N; T=Float64)`; the LGL rule and all operators are computed natively
in `T`, so `T=BigFloat` yields full-precision spectral operators.
"""
struct RefElement{T<:AbstractFloat}
    N::Int
    Np::Int
    r::Vector{T}
    w::Vector{T}
    V::Matrix{T}
    D::Matrix{T}
    M::Matrix{T}
    invM::Matrix{T}
end

function RefElement(N::Int; T::Type{<:AbstractFloat}=Float64)
    @assert N >= 1 "polynomial degree N must be ≥ 1"
    r, w = _lgl_nodes_weights(T, N)
    V, Vr = _legendre_vandermonde(r, N)
    D = Vr / V
    invM = V * transpose(V)
    M = inv(invM)
    return RefElement{T}(N, N + 1, r, w, V, D, M, invM)
end

"""
    lagrange_interp(el, u, ξ)

Evaluate the nodal field `u` (length `Np`) of reference element `el` at reference
coordinate `ξ ∈ [-1,1]` using barycentric Lagrange interpolation.
"""
function lagrange_interp(el::RefElement{T}, u::AbstractVector, ξ::Real) where {T}
    r = el.r
    Np = el.Np
    # barycentric weights for LGL nodes
    wbary = ones(T, Np)
    @inbounds for j in 1:Np, k in 1:Np
        if k != j
            wbary[j] /= (r[j] - r[k])
        end
    end
    num = zero(eltype(u))
    den = zero(T)
    @inbounds for j in 1:Np
        d = T(ξ) - r[j]
        if abs(d) < eps(T) * 10
            return u[j]
        end
        t = wbary[j] / d
        num += t * u[j]
        den += t
    end
    return num / den
end

# ----------------------------------------------------------------------------
# Multi-domain mesh: K subdomains partitioning [a,b], optionally with forced
# interface points (e.g. layer start R and particle r_*p) that MUST be nodes.
# ----------------------------------------------------------------------------
"""
    Mesh{T}

A 1D partition `ρ₀ < ρ₁ < … < ρ_K` into `K` subdomains, with per-element left
edge `xL`, width `h`, and the affine map to physical coordinates `ρ = xL + (ξ+1)/2 h`.
"""
struct Mesh{T<:AbstractFloat}
    K::Int
    edges::Vector{T}   # length K+1
    xL::Vector{T}      # length K
    h::Vector{T}       # length K
end

"""
    build_mesh(a, b, K; forced=Float64[])

Uniform `K`-subdomain mesh on `[a,b]`. Any coordinates in `forced` replace distinct nearest
interior interfaces so they land exactly on subdomain boundaries (required for the hyperboloidal
layer start `R` and the particle interface `r_*p`). The two endpoints are immutable. If there are
more distinct interior requests than interior slots, placement is impossible and an
`ArgumentError` is thrown.
"""
function build_mesh(a::Real, b::Real, K::Int; forced::AbstractVector=Float64[], T::Type{<:AbstractFloat}=Float64)
    K >= 1 || throw(ArgumentError("number of subdomains K must be at least 1 (got $K)"))
    aa, bb = T(a), T(b)
    aa < bb || throw(ArgumentError("mesh endpoints must satisfy a < b (got [$a,$b])"))
    edges = collect(range(aa, bb; length=K + 1))

    requested = T[]
    for raw in forced
        f = T(raw)
        aa <= f <= bb || throw(ArgumentError("forced interface $raw outside [$a,$b]"))
        push!(requested, f)
    end
    requested = sort!(unique(requested))
    interior = filter(f -> aa < f < bb, requested)
    length(interior) <= K - 1 ||
        throw(ArgumentError("cannot place $(length(interior)) distinct interior forced interfaces " *
                            "in a $K-subdomain mesh with only $(K - 1) interior slots"))

    available = collect(2:K) # endpoints 1 and K+1 are never candidates
    for f in interior
        pos = argmin(abs.(edges[available] .- f))
        slot = available[pos]
        edges[slot] = f
        deleteat!(available, pos)
    end
    sort!(edges)
    h = diff(edges)
    all(width -> width > zero(T), h) ||
        throw(ArgumentError("forced interfaces collapsed mesh; no injective placement exists"))
    edges[1] == aa && edges[end] == bb || error("internal error: mesh endpoints were not preserved")
    all(f -> f in edges, requested) || error("internal error: a forced interface was not placed")
    xL = edges[1:end-1]
    return Mesh{T}(K, edges, xL, h)
end

"""
    node_coords(mesh, el)

Physical coordinates of all DG nodes: returns the `(K × Np)` matrix `ρ[k,i]` of the affine images
of the reference nodes `el.r` in each subdomain of `mesh`.
"""
function node_coords(mesh::Mesh{T}, el::RefElement{T}) where {T}
    ρ = Matrix{T}(undef, mesh.K, el.Np)
    @inbounds for k in 1:mesh.K, i in 1:el.Np
        ρ[k, i] = mesh.xL[k] + (el.r[i] + 1) / 2 * mesh.h[k]
    end
    return ρ
end

"""
    element_of(mesh, ρ)

Index of the `mesh` element containing `ρ` (the first whose right edge is `≥ ρ`, with a `1e-12`
interface tolerance; clamps to the last element).
"""
function element_of(mesh::Mesh, ρ::Real)
    for k in 1:mesh.K
        if ρ <= mesh.edges[k+1] + 1e-12
            return k
        end
    end
    return mesh.K
end

# ----------------------------------------------------------------------------
# Explicit time stepping (DG-Teulolsky TimeStepping.jl, ported verbatim): classical RK4 with an
# in-place closure RHS `rhs!(du, u, t)`, plus the DG CFL time-step helper.
# ----------------------------------------------------------------------------

"""
    rk4_step!(rhs!, u, t, dt, work)

Advance state `u` (any mutable array) by one RK4 step in place. `work` is a tuple
of 5 scratch arrays `(k1,k2,k3,k4,utmp)` sized like `u`.
"""
function rk4_step!(rhs!, u, t, dt, work)
    k1, k2, k3, k4, utmp = work
    rhs!(k1, u, t)
    @. utmp = u + (dt / 2) * k1
    rhs!(k2, utmp, t + dt / 2)
    @. utmp = u + (dt / 2) * k2
    rhs!(k3, utmp, t + dt / 2)
    @. utmp = u + dt * k3
    rhs!(k4, utmp, t + dt)
    @. u = u + (dt / 6) * (k1 + 2k2 + 2k3 + k4)
    return u
end

"""
    rk4!(rhs!, u, tspan, dt; callback=nothing) -> (u, t)

Integrate `rhs!(du, u, t)` from `tspan[1]` to `tspan[2]` with fixed-step classical RK4 (the last
step is shortened to land exactly on the end time). `callback(u, t, step)` is invoked after every
step when supplied.
"""
function rk4!(rhs!, u, tspan, dt; callback=nothing)
    t0, t1 = tspan
    (dt isa Real && isfinite(dt) && dt > zero(dt)) ||
        throw(ArgumentError("rk4!: dt must be finite and positive (got $dt)"))
    t0 <= t1 || throw(ArgumentError("rk4!: tspan must be forward ordered (got $tspan)"))
    work = (similar(u), similar(u), similar(u), similar(u), similar(u))
    t = t0
    step = 0
    callback === nothing || callback(u, t, step)
    while t < t1
        remaining = t1 - t
        h = min(dt, remaining)
        rk4_step!(rhs!, u, t, h, work)
        # Assign the endpoint directly on the shortened final step: this both executes arbitrarily
        # short positive spans and avoids leaving t one rounding unit shy of t1.
        t = h == remaining ? t1 : t + h
        step += 1
        callback === nothing || callback(u, t, step)
    end
    return u, t
end

"""
    cfl_dt(C, Δx_min, λ_max=1.0)

DG CFL estimate `Δt = C · Δx_min / λ_max` (DG-Teulolsky Paper 1 §III.A convention, `λ_max = 1`).
`Δx_min` is the smallest physical spacing between neighboring LGL nodes.
"""
cfl_dt(C, Δx_min, λ_max=1.0) = C * Δx_min / λ_max
