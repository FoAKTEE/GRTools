# FlatWave — TEST-SUPPORT harness for Perturbation.TimeDomain (DG-P3). Port of
# DG-Teulolsky/src/FlatWave.jl, kept OUT of the runtime module by design (the repo used it as
# its §IV.A verification problem, not as physics): Paper 1 §IV.A / Paper 2 §IV.A — the ℓ=2
# flat-space wave equation `[-∂_tt + ∂_rr - ℓ(ℓ+1)/r²]ψ = 0` (the `a=M=0` limit of the scalar
# Teukolsky equation), solved with the **System-1** first-order hyperboloidal reduction and the
# nodal DG scheme. It has a closed-form outgoing solution, so it gives an unambiguous
# convergence test — and the place where the `2N+1` superconvergence at 𝓘⁺ is demonstrated.
#
# System (single mode, `f=μ=0`), with `H = H(ρ)`:
#     ψ̇ = -π
#     π̇ - H φ̇ = -(1-H) φ' + (ℓ(ℓ+1)/r²) ψ
#     φ̇ - H π̇ = -(1-H) π'
# Inverting the 2×2 `E = [[1,-H],[-H,1]]` gives flux matrix `A = 1/(1+H) [[H,1],[1,H]]`
# (finite at ρ=s) and a finite potential source `s_pot = [Ppot ψ, H·Ppot ψ]`,
# `Ppot = ℓ(ℓ+1)(Ω-ρΩ')/(ρ²(1+H))`.
#
# Only the pieces the repo's own testsets exercise are ported (profile_derivs, FlatWaveProblem,
# solve_flatwave, rel_L2_in_time + internals). DG core from NumericsCore.dg; Layer machinery
# from Perturbation.TimeDomain.
module TDFlatWave

using LinearAlgebra: mul!, dot
using NumericsCore: RefElement, Mesh, build_mesh, element_of, node_coords, rk4!
using Perturbation.TimeDomain: Layer, Hfun, denomΩ, rstar_of_ρ

export FlatWaveProblem, solve_flatwave, exact_flatwave, rel_L2_in_time, profile_derivs

# ----------------------------------------------------------------------------
# Sine-Gaussian profile f(u)=sin(f0(u-u0))·exp(-cc(u-u0)²) and derivatives f',f'',f'''.
# ----------------------------------------------------------------------------
"""
    profile_derivs(u; f0=2.0, cc=1.0, u0=-10.0) -> (f, f1, f2, f3)

Return the sine-Gaussian profile and its first three derivatives at `u`.
"""
@inline function profile_derivs(u::T; f0=2.0, cc=1.0, u0=-10.0) where {T<:Real}
    g = u - u0
    θ = f0 * g
    s, c = sin(θ), cos(θ)
    E = exp(-cc * g^2)
    f = s * E
    f1 = (f0 * c - 2cc * g * s) * E
    f2 = ((-f0^2 - 2cc + 4cc^2 * g^2) * s - 4cc * f0 * g * c) * E
    f3 = ((12cc^2 * g + 6cc * f0^2 * g - 8cc^3 * g^3) * s +
          (-f0^3 - 6cc * f0 + 12cc^2 * f0 * g^2) * c) * E
    return f, f1, f2, f3
end

"Exact outgoing ℓ=2 solution ψ at retarded time `u=τ-ρ` and areal radius `r`."
@inline function exact_flatwave(u::Real, r::Real; kw...)
    f, f1, f2, _ = profile_derivs(u; kw...)
    return f2 + 3 / r * f1 + 3 / r^2 * f          # r=Inf ⇒ ψ=f2
end

# ----------------------------------------------------------------------------
# Problem container with all precomputes.
# ----------------------------------------------------------------------------
struct FlatWaveProblem{T<:AbstractFloat}
    ℓ::Int
    el::RefElement{T}
    mesh::Mesh{T}
    layer::Layer{T}
    ρ::Matrix{T}          # (Np,K) node coords
    H::Matrix{T}          # (Np,K) layer function
    Ppot::Matrix{T}       # (Np,K) finite potential coefficient
    invMc1::Vector{T}     # invM[:,1]
    invMcN::Vector{T}     # invM[:,Np]
    twoOverh::Vector{T}   # 2/h_k
    HL::T                 # H at left physical boundary
    probes::Vector{T}
    probe_k::Vector{Int}
    probe_L::Matrix{T}    # (nprobe,Np) interpolation rows
    # scratch
    Dpi::Vector{T}
    Dphi::Vector{T}
end

function lagrange_row(el::RefElement{T}, ξ::Real) where {T}
    r = el.r; Np = el.Np
    wb = ones(T, Np)
    @inbounds for j in 1:Np, k in 1:Np
        k != j && (wb[j] /= (r[j] - r[k]))
    end
    L = zeros(T, Np)
    den = zero(T); exact = 0
    @inbounds for j in 1:Np
        d = T(ξ) - r[j]
        if abs(d) < eps(T) * 10
            exact = j
        else
            L[j] = wb[j] / d
            den += L[j]
        end
    end
    if exact != 0
        fill!(L, zero(T)); L[exact] = one(T); return L
    end
    L ./= den
    return L
end

function FlatWaveProblem(; ℓ=2, N=4, K=16, R=30.0, s=50.0, P=4, a=1.0,
                         probes=[15.0, 40.0, 50.0], T::Type{<:AbstractFloat}=Float64)
    el = RefElement(N; T=T)
    layer = Layer{T}(T(R), T(s), P)
    mesh = build_mesh(T(a), T(s), K; forced=[T(R)], T=T)
    coords = node_coords(mesh, el)             # (K,Np)
    ρ = permutedims(coords)                     # (Np,K)
    Hm = similar(ρ); Ppot = similar(ρ)
    ℓℓ = T(ℓ * (ℓ + 1))
    @inbounds for k in 1:mesh.K, i in 1:el.Np
        ρv = ρ[i, k]
        Hv = Hfun(layer, ρv)
        Hm[i, k] = Hv
        Ppot[i, k] = ℓℓ * denomΩ(layer, ρv) / (ρv^2 * (1 + Hv))
    end
    probe_k = [element_of(mesh, p) for p in probes]
    probe_L = Matrix{T}(undef, length(probes), el.Np)
    for (ip, p) in enumerate(probes)
        k = probe_k[ip]
        ξ = 2 * (T(p) - mesh.xL[k]) / mesh.h[k] - 1
        probe_L[ip, :] = lagrange_row(el, ξ)
    end
    return FlatWaveProblem{T}(ℓ, el, mesh, layer, ρ, Hm, Ppot,
        el.invM[:, 1], el.invM[:, el.Np], 2 ./ mesh.h, Hfun(layer, mesh.edges[1]),
        T.(probes), probe_k, probe_L, Vector{T}(undef, el.Np), Vector{T}(undef, el.Np))
end

# initial data from the exact solution at τ=0 (u = -ρ, r = ρ/Ω)
function initial_state(prob::FlatWaveProblem{T}; kw...) where {T}
    Np, K = prob.el.Np, prob.mesh.K
    u0 = zeros(T, Np * K * 3)
    U = reshape(u0, Np, K, 3)
    @inbounds for k in 1:K, i in 1:Np
        ρv = prob.ρ[i, k]
        r = rstar_of_ρ(prob.layer, ρv)         # = ρ/Ω (Inf at ρ=s)
        u = -ρv
        f, f1, f2, f3 = profile_derivs(u; kw...)
        ψ = f2 + 3 / r * f1 + 3 / r^2 * f
        π = -(f3 + 3 / r * f2 + 3 / r^2 * f1)
        φ = π - 3 / r^2 * f1 - 6 / r^3 * f
        U[i, k, 1] = ψ; U[i, k, 2] = π; U[i, k, 3] = φ
    end
    return u0
end

# A·U for the 2×2 flux matrix A = 1/(1+H)[[H,1],[1,H]] acting on (π,φ)
@inline flux_AU(H, π, φ) = (invH = 1 / (1 + H); (invH * (H * π + φ), invH * (π + H * φ)))

function flatwave_rhs!(du, u, ::Real, prob::FlatWaveProblem{T}) where {T}
    el = prob.el; Np = el.Np; K = prob.mesh.K
    U = reshape(u, Np, K, 3)
    dU = reshape(du, Np, K, 3)
    Dpi = prob.Dpi; Dphi = prob.Dphi
    # ---- volume terms ----
    @inbounds for k in 1:K
        invJ = prob.twoOverh[k]
        mul!(Dpi, el.D, view(U, :, k, 2))
        mul!(Dphi, el.D, view(U, :, k, 3))
        for i in 1:Np
            H = prob.H[i, k]
            invH = 1 / (1 + H)
            dπ = invJ * Dpi[i]; dφ = invJ * Dphi[i]
            AUπ = invH * (H * dπ + dφ)
            AUφ = invH * (dπ + H * dφ)
            ψ = U[i, k, 1]
            dU[i, k, 1] = -U[i, k, 2]
            dU[i, k, 2] = -AUπ + prob.Ppot[i, k] * ψ
            dU[i, k, 3] = -AUφ + H * prob.Ppot[i, k] * ψ
        end
    end
    # ---- interior interface fluxes (LLF, λ=1) ----
    @inbounds for k in 1:(K - 1)
        H = prob.H[Np, k]
        πL = U[Np, k, 2]; φL = U[Np, k, 3]
        πR = U[1, k + 1, 2]; φR = U[1, k + 1, 3]
        FLπ, FLφ = flux_AU(H, πL, φL)
        FRπ, FRφ = flux_AU(H, πR, φR)
        Fsπ = 0.5 * (FLπ + FRπ) - 0.5 * (πR - πL)
        Fsφ = 0.5 * (FLφ + FRφ) - 0.5 * (φR - φL)
        dRπ = FLπ - Fsπ; dRφ = FLφ - Fsφ        # (F-F*) on right face of elem k
        dLπ = FRπ - Fsπ; dLφ = FRφ - Fsφ        # (F-F*) on left  face of elem k+1
        ck = prob.twoOverh[k]; ck1 = prob.twoOverh[k + 1]
        for i in 1:Np
            dU[i, k, 2] += ck * prob.invMcN[i] * dRπ
            dU[i, k, 3] += ck * prob.invMcN[i] * dRφ
            dU[i, k + 1, 2] -= ck1 * prob.invMc1[i] * dLπ
            dU[i, k + 1, 3] -= ck1 * prob.invMc1[i] * dLφ
        end
    end
    # ---- left physical boundary: weak Sommerfeld (zero incoming char p=π+φ) ----
    @inbounds begin
        H = prob.HL
        πi = U[1, 1, 2]; φi = U[1, 1, 3]
        q = πi - φi                                # outgoing char (λ<0) kept
        πg = q / 2; φg = -q / 2                     # ghost: p=0
        FLπ, FLφ = flux_AU(H, πg, φg)
        FRπ, FRφ = flux_AU(H, πi, φi)
        Fsπ = 0.5 * (FLπ + FRπ) - 0.5 * (πi - πg)
        Fsφ = 0.5 * (FLφ + FRφ) - 0.5 * (φi - φg)
        dLπ = FRπ - Fsπ; dLφ = FRφ - Fsφ
        c1 = prob.twoOverh[1]
        for i in 1:Np
            dU[i, 1, 2] -= c1 * prob.invMc1[i] * dLπ
            dU[i, 1, 3] -= c1 * prob.invMc1[i] * dLφ
        end
    end
    # ---- right physical boundary ρ=s: pure outflow ⇒ no correction ----
    return du
end

"""
    solve_flatwave(prob; tfinal=50.0, dt=3e-4, record_every=50, profile...) -> (τs, Ψnum, u)

Evolve to `tfinal`; return recorded times `τs` and `Ψnum` (nprobe × nsamples) of the
interpolated numerical ψ at each probe location (+ the final state).
"""
function solve_flatwave(prob::FlatWaveProblem{T}; tfinal=50.0, dt=3e-4,
                        record_every=50, kw...) where {T}
    u = initial_state(prob; kw...)
    Np, K = prob.el.Np, prob.mesh.K
    nprb = length(prob.probes)
    τs = T[]; cols = Vector{T}[]
    record!(uu, t, step) = begin
        if step % record_every == 0
            U = reshape(uu, Np, K, 3)
            push!(τs, t)
            v = Vector{T}(undef, nprb)
            for ip in 1:nprb
                v[ip] = dot(view(prob.probe_L, ip, :), view(U, :, prob.probe_k[ip], 1))
            end
            push!(cols, v)
        end
    end
    rhs!(du, uu, t) = flatwave_rhs!(du, uu, t, prob)
    rk4!(rhs!, u, (0.0, tfinal), dt; callback=record!)
    Ψ = Matrix{T}(undef, nprb, length(τs))
    for j in eachindex(cols), ip in 1:nprb
        Ψ[ip, j] = cols[j][ip]
    end
    return τs, Ψ, u            # u is the final state (mutated in place by rk4!)
end

"relative L2-in-time error at probe `ip` against the exact solution."
function rel_L2_in_time(prob::FlatWaveProblem, τs, Ψ, ip::Int; kw...)
    ρp = prob.probes[ip]
    r = rstar_of_ρ(prob.layer, ρp)
    num = 0.0; den = 0.0
    for j in 2:length(τs)
        dτ = τs[j] - τs[j - 1]
        for (jj, w) in ((j - 1, 0.5), (j, 0.5))
            ψe = exact_flatwave(τs[jj] - ρp, r; kw...)
            ψn = Ψ[ip, jj]
            num += w * dτ * (ψn - ψe)^2
            den += w * dτ * ψe^2
        end
    end
    return num / den
end

end # module TDFlatWave
