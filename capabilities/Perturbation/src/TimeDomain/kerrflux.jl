# Kerr (a≠0) circular-orbit scalar **flux** in the coupled m≠0 sector (Paper 1 §IV.C Kerr) —
# port of DG-Teulolsky/src/KerrFlux.jl. Combines the multi-mode coupling of kerrmodes.jl
# (complex fields + `μ` friction) with the characteristic delta-flux source of flux1d.jl,
# generalized to the 2n×2n block system: the source injects `±(T P^± T⁻¹) G(t)` using the
# eigendecomposition of the flux matrix `A` at the particle interface. Energy flux at 𝓘⁺:
# `Ė_{ℓm}=|π_{ℓm}|²/(4π)`. Orbital quantities are the Kerr circular-equatorial ones
# (repo commit a643b5f — the Schwarzschild forms gave a ~5% flux error), via `orbital_params`
# (→ Spacetime).

"""
    KerrFluxProblem{T}

Precomputed Kerr coupled-mode circular-orbit flux problem (Paper 1 §IV.C Kerr, T1.3/F1.7):
the [`KerrTailProblem`](@ref) coupled-block machinery (complex fields, `μ` friction) plus the
particle interface `kp`, the eigen-decomposed characteristic injection matrices
`T P^± T⁻¹` and the per-mode source prefactors
`Gpre_ℓ = -4πq/(u^t √(r_p²+a²))·Y_{ℓm}(π/2,0)`. Build with
[`KerrFluxProblem(ℓs, m, rp; ...)`](@ref KerrFluxProblem(::Vector{Int}, ::Any, ::Any));
evolve with [`solve_kerr_flux`](@ref).
"""
struct KerrFluxProblem{T<:AbstractFloat}
    ℓs::Vector{Int}
    m::Int
    n::Int
    el::RefElement{T}
    mesh::Mesh{T}
    layer::Layer{T}
    M::T
    a::T
    rp::T
    H::Matrix{T}
    Aflux::Array{T,4}
    Sinv::Array{T,4}
    Vmode::Array{T,3}
    μmode::Array{Complex{T},3}
    invMc1::Vector{T}
    invMcN::Vector{T}
    twoOverh::Vector{T}
    kp::Int
    Minjp::Matrix{T}            # T P^+ T^-1 at the particle (2n×2n)
    Minjm::Matrix{T}            # T P^- T^-1
    Gpre::Vector{Complex{T}}    # per-mode source prefactor
    Ω::T
    ramp_τ::T
    ramp_δ::T
    Dπ::Matrix{Complex{T}}
    Dφ::Matrix{Complex{T}}
end

"""
    KerrFluxProblem(ℓs, m, rp; N, K, ρL, R, s, P=4, M=1, a=0.9, q=1, ramp_τ=200, ramp_δ=0.001, T=Float64)

Build the Kerr coupled circular-orbit flux problem for the mode set `ℓs` (fixed `m`, same
ℓ-parity) and orbital radius `rp` (particle `r_*,p` and layer start `R` forced onto the mesh).
"""
function KerrFluxProblem(ℓs::Vector{Int}, m, rp; N=18, K=320, ρL=-1600.0, R=100.0,
                         s=1600.0, P=4, M=1.0, a=0.9, q=1.0, ramp_τ=200.0,
                         ramp_δ=0.001, T::Type{<:AbstractFloat}=Float64)
    n = length(ℓs)
    rsp = _tortoise(rp; M=M, a=a)
    el = RefElement(N; T=T)
    layer = Layer{T}(T(R), T(s), P)
    mesh = build_mesh(T(ρL), T(s), K; forced=[T(rsp), T(R)], T=T)
    jp = findmin(abs.(mesh.edges .- T(rsp)))[2]
    kp = jp - 1
    ρ = permutedims(node_coords(mesh, el))
    Np, Kk = el.Np, mesh.K
    Hm = similar(ρ)
    Aflux = Array{T,4}(undef, 2n, 2n, Np, Kk)
    Sinv = Array{T,4}(undef, n, n, Np, Kk)
    Vmode = Array{T,3}(undef, n, Np, Kk)
    μmode = Array{Complex{T},3}(undef, n, Np, Kk)
    @inbounds for k in 1:Kk, i in 1:Np
        ρv = ρ[i, k]; Hv = Hfun(layer, ρv); Hm[i, k] = Hv
        if k == Kk && i == Np
            r = Inf; f = 0.0
        else
            r = _inverse_tortoise(rstar_of_ρ(layer, ρv); M=M, a=a)
            f = coupling_f(r; M=M, a=a)
        end
        A, Si = _Amat(_Eππ(ℓs, m, f, n), Hv, n)
        Aflux[:, :, i, k] = A; Sinv[:, :, i, k] = Si
        for q2 in 1:n
            Vmode[q2, i, k] = isinf(r) ? 0.0 : potential_V(r, ℓs[q2], m; M=M, a=a)
            μmode[q2, i, k] = isinf(r) ? 0.0 + 0im : friction_mu(r, m; M=M, a=a)
        end
    end
    # source injection matrices from eigendecomposition of A at the particle
    Ap = Aflux[:, :, Np, kp]
    F = eigen(Ap)
    vals = real.(F.values); V = real.(F.vectors); Vi = inv(V)
    Pp = V * (Float64.(vals .> 1e-12) .* Vi)
    Pm = V * (Float64.(vals .< -1e-12) .* Vi)
    op = orbital_params(rp; M=M, a=a)        # Kerr orbital quantities (repo a643b5f fix)
    Gpre = Complex{T}[(-4π * q / (op.ut * sqrt(rp^2 + a^2))) * Ylm_equator(ℓ, m) for ℓ in ℓs]
    return KerrFluxProblem{T}(ℓs, m, n, el, mesh, layer, T(M), T(a), T(rp), Hm, Aflux,
        Sinv, Vmode, μmode, el.invM[:, 1], el.invM[:, el.Np], 2 ./ mesh.h, kp,
        T.(Pp), T.(Pm), Gpre, T(op.Ω), T(ramp_τ), T(ramp_δ),
        Matrix{Complex{T}}(undef, Np, n), Matrix{Complex{T}}(undef, Np, n))
end

function _rhs!(du, u, t, p::KerrFluxProblem{T}, sc) where {T}
    el = p.el; Np = el.Np; K = p.mesh.K; n = p.n
    U = reshape(u, Np, K, 3, n); dU = reshape(du, Np, K, 3, n)
    Dπ, Dφ, Fπ, Fφ, FπR, FφR, tπ, tφ = sc
    @inbounds for k in 1:K
        invJ = p.twoOverh[k]
        for q in 1:n
            mul!(view(Dπ, :, q), el.D, view(U, :, k, 2, q))
            mul!(view(Dφ, :, q), el.D, view(U, :, k, 3, q))
        end
        for i in 1:Np
            A = view(p.Aflux, :, :, i, k); H = p.H[i, k]; Si = view(p.Sinv, :, :, i, k)
            for q in 1:n
                tπ[q] = invJ * Dπ[i, q]; tφ[q] = invJ * Dφ[i, q]
            end
            _AUvec!(Fπ, Fφ, A, tπ, tφ, n)
            for q in 1:n
                w = zero(Complex{T})
                for c in 1:n
                    w += Si[q, c] * (p.μmode[c, i, k] * U[i, k, 2, c] + p.Vmode[c, i, k] * U[i, k, 1, c])
                end
                dU[i, k, 1, q] = -U[i, k, 2, q]
                dU[i, k, 2, q] = -Fπ[q] - w
                dU[i, k, 3, q] = -Fφ[q] - H * w
            end
        end
    end
    @inbounds for k in 1:(K-1)
        A = view(p.Aflux, :, :, Np, k)
        for q in 1:n; tπ[q] = U[Np, k, 2, q]; tφ[q] = U[Np, k, 3, q]; end
        _AUvec!(Fπ, Fφ, A, tπ, tφ, n)
        for q in 1:n; tπ[q] = U[1, k+1, 2, q]; tφ[q] = U[1, k+1, 3, q]; end
        _AUvec!(FπR, FφR, A, tπ, tφ, n)
        ck = p.twoOverh[k]; ck1 = p.twoOverh[k+1]
        for q in 1:n
            uLπ = U[Np, k, 2, q]; uLφ = U[Np, k, 3, q]; uRπ = U[1, k+1, 2, q]; uRφ = U[1, k+1, 3, q]
            Fsπ = 0.5 * (Fπ[q] + FπR[q]) - 0.5 * (uRπ - uLπ)
            Fsφ = 0.5 * (Fφ[q] + FφR[q]) - 0.5 * (uRφ - uLφ)
            dRπ = Fπ[q] - Fsπ; dRφ = Fφ[q] - Fsφ; dLπ = FπR[q] - Fsπ; dLφ = FφR[q] - Fsφ
            for i in 1:Np
                dU[i, k, 2, q] += ck * p.invMcN[i] * dRπ
                dU[i, k, 3, q] += ck * p.invMcN[i] * dRφ
                dU[i, k+1, 2, q] -= ck1 * p.invMc1[i] * dLπ
                dU[i, k+1, 3, q] -= ck1 * p.invMc1[i] * dLφ
            end
        end
    end
    # delta source injection at the particle (kp | kp+1)
    @inbounds begin
        ramp = _ramp(t, p.ramp_τ, p.ramp_δ); ph = cis(-p.m * p.Ω * t)
        # 2n source vector G = [Gpre_q·e^{-imΩt}·ramp ; 0]
        Gp = zeros(Complex{T}, 2p.n); Gm = zeros(Complex{T}, 2p.n)
        @views for q in 1:n
            g = p.Gpre[q] * ph * ramp
            for r in 1:2n
                Gp[r] += p.Minjp[r, q] * g
                Gm[r] += p.Minjm[r, q] * g
            end
        end
        ckp = p.twoOverh[p.kp]; ckp1 = p.twoOverh[p.kp+1]
        for q in 1:n
            for i in 1:Np
                dU[i, p.kp, 2, q] += ckp * p.invMcN[i] * (-Gm[q])
                dU[i, p.kp, 3, q] += ckp * p.invMcN[i] * (-Gm[n+q])
                dU[i, p.kp+1, 2, q] -= ckp1 * p.invMc1[i] * Gp[q]
                dU[i, p.kp+1, 3, q] -= ckp1 * p.invMc1[i] * Gp[n+q]
            end
        end
    end
    # left boundary: weak Sommerfeld per mode
    @inbounds begin
        A = view(p.Aflux, :, :, 1, 1); c1 = p.twoOverh[1]
        for q in 1:n
            qq = U[1, 1, 2, q] - U[1, 1, 3, q]; tπ[q] = qq / 2; tφ[q] = -qq / 2
        end
        _AUvec!(Fπ, Fφ, A, tπ, tφ, n)
        for q in 1:n; tπ[q] = U[1, 1, 2, q]; tφ[q] = U[1, 1, 3, q]; end
        _AUvec!(FπR, FφR, A, tπ, tφ, n)
        for q in 1:n
            πi = U[1, 1, 2, q]; φi = U[1, 1, 3, q]; qq = πi - φi
            Fsπ = 0.5 * (Fπ[q] + FπR[q]) - 0.5 * (πi - qq / 2)
            Fsφ = 0.5 * (Fφ[q] + FφR[q]) - 0.5 * (φi + qq / 2)
            dLπ = FπR[q] - Fsπ; dLφ = FφR[q] - Fsφ
            for i in 1:Np
                dU[i, 1, 2, q] -= c1 * p.invMc1[i] * dLπ
                dU[i, 1, 3, q] -= c1 * p.invMc1[i] * dLφ
            end
        end
    end
    return du
end

"""
    solve_kerr_flux(p; tfinal, dt, record_every, target=1) -> (τs, Edot)

Trivial IC; evolve; record 𝓘⁺ energy flux `|π(scri)|²/(4π)` of mode `ℓs[target]`.
"""
function solve_kerr_flux(p::KerrFluxProblem{T}; tfinal=4000.0, dt=0.017, record_every=200,
                         target=1) where {T}
    Np, K, n = p.el.Np, p.mesh.K, p.n
    u = zeros(Complex{T}, Np * K * 3 * n)
    sc = (Matrix{Complex{T}}(undef, Np, n), Matrix{Complex{T}}(undef, Np, n),
        zeros(Complex{T}, n), zeros(Complex{T}, n), zeros(Complex{T}, n),
        zeros(Complex{T}, n), zeros(Complex{T}, n), zeros(Complex{T}, n))
    τs = T[]; Es = T[]
    rec(uu, t, step) = begin
        if step % record_every == 0
            UU = reshape(uu, Np, K, 3, n)
            push!(τs, t); push!(Es, abs2(UU[Np, K, 2, target]) / (4π))
        end
    end
    rhs!(du, uu, t) = _rhs!(du, uu, t, p, sc)
    rk4!(rhs!, u, (0.0, tfinal), dt; callback=rec)
    return τs, Es
end
