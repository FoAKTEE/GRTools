# Kerr multi-mode coupled hyperboloidal DG evolution (Paper 1 §II.E–F) — port of
# DG-Teulolsky/src/KerrModes.jl. For a fixed `m` and ℓ-parity, the modes
# `ℓ = ℓ_min, ℓ_min+2, …, ℓ_max` couple through `sin²θ Y_{ℓm}` (coefficients `c^L_{ℓm}`). The
# block system is `E U̇ + Â U' + B̂ + V̂ = Ĝ δ` with `E=[[E_ππ,-HI],[-HI,I]]`, `E_ππ` symmetric
# tridiagonal (`diag 1-f c^ℓ_{ℓm}`, `offdiag -f c^ℓ_{ℓ±2,m}`). Inverting,
# `A=(1-H)[[H S⁻¹, S⁻¹],[I+H²S⁻¹, H S⁻¹]]`, `S=E_ππ-H²I` (`A=½[[I,I],[I,I]]` at 𝓘⁺).
#
# This covers the **m=0** (real, frictionless `μ=0`) sector used by the Kerr Price-tail test
# (F1.4) and the m=0 rows of the scattering matrix (T1.1). The `E_ππ`/`A`-matrix builders live
# in common.jl (shared with kerrflux.jl); the tortoise inverse delegates to Spacetime.

"""
    cminus(l, m)

Mode-coupling ladder coefficient `c⁻_{lm} = √((l²-m²)/((2l-1)(2l+1)))` (0 for `l ≤ |m|`;
Paper 1 App A).
"""
cminus(l, m) = l <= abs(m) ? 0.0 : sqrt((l^2 - m^2) / ((2l - 1) * (2l + 1)))

"""
    mode_coupling(l, m) -> (C₋₋, C₀, C₊₊)

`(c^ℓ_{ℓ-2,m}, c^ℓ_{ℓ,m}, c^ℓ_{ℓ+2,m})` — the `sin²θ Y_{ℓm}` inter-mode coupling coefficients
(Paper 1 App A closed form).
"""
function mode_coupling(l, m)
    cl = cminus(l, m); clp1 = cminus(l + 1, m); clp2 = cminus(l + 2, m); clm1 = cminus(l - 1, m)
    return (-clm1 * cl, 1 - cl^2 - clp1^2, -clp1 * clp2)
end

"""
    mode_list(m, ℓmax; parity=:even)

Even/odd-ℓ coupled-mode list for fixed `m` up to `ℓmax` (step 2 from `ℓ_min`).
"""
function mode_list(m, ℓmax; parity=:even)
    base = abs(m)
    ℓmin = (iseven(base) == (parity == :even)) ? base : base + 1
    return collect(ℓmin:2:ℓmax)
end

"""
    KerrTailProblem{T}

Precomputed Kerr m=0 coupled-mode tail problem (Paper 1 §II.E–F, F1.4): per-node coupled flux
matrices `Aflux` (2n×2n), lower-order `S⁻¹`, per-mode potentials `Vmode`, interface/boundary
flux matrices and probe interpolation rows. Build with
[`KerrTailProblem(ℓs, m; ...)`](@ref KerrTailProblem(::Vector{Int}, ::Any)); evolve with
[`solve_kerr_tail`](@ref).
"""
struct KerrTailProblem{T<:AbstractFloat}
    ℓs::Vector{Int}
    m::Int
    n::Int
    el::RefElement{T}
    mesh::Mesh{T}
    layer::Layer{T}
    M::T
    a::T
    ρ::Matrix{T}
    H::Matrix{T}
    Aflux::Array{T,4}     # (2n,2n,Np,K) real flux matrix per node
    Sinv::Array{T,4}      # (n,n,Np,K) S⁻¹ for the lower-order term
    Vmode::Array{T,3}     # (n,Np,K) potential V_{ℓm}(r) per mode
    invMc1::Vector{T}
    invMcN::Vector{T}
    twoOverh::Vector{T}
    Aint::Array{T,3}      # (2n,2n,K-1) flux matrix at interior interfaces
    AL::Matrix{T}         # flux matrix at left boundary
    probe_ρ::Vector{T}
    probe_k::Vector{Int}
    probe_L::Matrix{T}
end

"""
    KerrTailProblem(ℓs, m; N, K, ρL, R, s, P=4, M=1, a=0.99995, extract_rstar=[100.0], T=Float64)

Build the Kerr m=0 coupled-mode problem for the mode set `ℓs` (same parity, step 2).
`extract_rstar` are tortoise-coordinate extraction locations.
"""
function KerrTailProblem(ℓs::Vector{Int}, m; N=15, K=400, ρL=-2100.0, R=149.35,
                         s=400.0, P=4, M=1.0, a=0.99995, extract_rstar=[100.0],
                         T::Type{<:AbstractFloat}=Float64)
    @assert m == 0 "KerrTailProblem currently supports m=0 (real sector)"
    n = length(ℓs)
    el = RefElement(N; T=T)
    layer = Layer{T}(T(R), T(s), P)
    mesh = build_mesh(T(ρL), T(s), K; forced=[T(R)], T=T)
    ρ = permutedims(node_coords(mesh, el))
    Np, Kk = el.Np, mesh.K
    Hm = similar(ρ)
    Aflux = Array{T,4}(undef, 2n, 2n, Np, Kk)
    Sinv = Array{T,4}(undef, n, n, Np, Kk)
    Vmode = Array{T,3}(undef, n, Np, Kk)
    @inbounds for k in 1:Kk, i in 1:Np
        ρv = ρ[i, k]; Hv = Hfun(layer, ρv); Hm[i, k] = Hv
        if k == Kk && i == Np
            r = Inf; f = 0.0
        else
            r = _inverse_tortoise(rstar_of_ρ(layer, ρv); M=M, a=a)
            f = coupling_f(r; M=M, a=a)
        end
        A, Si = _Amat(_Eππ(ℓs, m, f, n), Hv, n)   # scri branch returns Si = 0 (lower-order ~0)
        Aflux[:, :, i, k] = A; Sinv[:, :, i, k] = Si
        for q in 1:n
            Vmode[q, i, k] = isinf(r) ? 0.0 : potential_V(r, ℓs[q], m; M=M, a=a)
        end
    end
    Aint = Array{T,3}(undef, 2n, 2n, Kk - 1)
    for k in 1:(Kk-1)
        Aint[:, :, k] = Aflux[:, :, Np, k]
    end
    AL = Aflux[:, :, 1, 1]
    probes = T.(extract_rstar)
    pρ = [isinf(p) ? layer.s : _ρ_of_rstar(layer, p, T(ρL)) for p in probes]
    pk = [element_of(mesh, x) for x in pρ]
    pL = Matrix{T}(undef, length(pρ), Np)
    for (ip, x) in enumerate(pρ)
        kk = pk[ip]; ξ = 2 * (x - mesh.xL[kk]) / mesh.h[kk] - 1
        pL[ip, :] = _lagrange_row(el, ξ)
    end
    return KerrTailProblem{T}(ℓs, m, n, el, mesh, layer, T(M), T(a), ρ, Hm, Aflux,
        Sinv, Vmode, el.invM[:, 1], el.invM[:, el.Np], 2 ./ mesh.h, Aint, AL,
        pρ, pk, pL)
end

# KerrModes-private bisection (KEPT VERBATIM alongside the public ρ_of_rstar: the repo's two
# bisections stop at different tolerances — 1e-9/1e-11 here vs 1e-10/1e-12 in Mode1D — so the
# probe locations of the two problem families must not be cross-contaminated).
function _ρ_of_rstar(L::Layer, rstar, ρlo)
    a, b = ρlo, L.s - 1e-9; fa = rstar_of_ρ(L, a) - rstar
    for _ in 1:200
        mm = (a + b) / 2; fm = rstar_of_ρ(L, mm) - rstar
        (abs(fm) < 1e-9 || (b - a) < 1e-11) && return mm
        sign(fm) == sign(fa) ? (a = mm; fa = fm) : (b = mm)
    end
    return (a + b) / 2
end

function _rhs!(du, u, ::Real, p::KerrTailProblem{T}, sc) where {T}
    el = p.el; Np = el.Np; K = p.mesh.K; n = p.n
    U = reshape(u, Np, K, 3, n); dU = reshape(du, Np, K, 3, n)
    Dπ, Dφ, Fπ, Fφ, FπR, FφR, tmpπ, tmpφ = sc
    @inbounds for k in 1:K
        invJ = p.twoOverh[k]
        for q in 1:n
            mul!(view(Dπ, :, q), el.D, view(U, :, k, 2, q))
            mul!(view(Dφ, :, q), el.D, view(U, :, k, 3, q))
        end
        for i in 1:Np
            A = view(p.Aflux, :, :, i, k)
            for q in 1:n
                tmpπ[q] = invJ * Dπ[i, q]; tmpφ[q] = invJ * Dφ[i, q]
            end
            _AUvec!(Fπ, Fφ, A, tmpπ, tmpφ, n)
            Si = view(p.Sinv, :, :, i, k); H = p.H[i, k]
            for q in 1:n
                # lower-order  -E⁻¹[V⊙ψ;0] : w=V⊙ψ ; contributes -[S⁻¹w ; H S⁻¹w]
                w = 0.0
                for c in 1:n
                    w += Si[q, c] * (p.Vmode[c, i, k] * real(U[i, k, 1, c]))
                end
                dU[i, k, 1, q] = -U[i, k, 2, q]
                dU[i, k, 2, q] = -Fπ[q] - w
                dU[i, k, 3, q] = -Fφ[q] - H * w
            end
        end
    end
    # interior interface LLF
    @inbounds for k in 1:(K-1)
        A = view(p.Aint, :, :, k)
        for q in 1:n
            tmpπ[q] = U[Np, k, 2, q]; tmpφ[q] = U[Np, k, 3, q]
        end
        _AUvec!(Fπ, Fφ, A, tmpπ, tmpφ, n)          # F(U^L)
        for q in 1:n
            tmpπ[q] = U[1, k+1, 2, q]; tmpφ[q] = U[1, k+1, 3, q]
        end
        _AUvec!(FπR, FφR, A, tmpπ, tmpφ, n)        # F(U^R)
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
    # left boundary: weak Sommerfeld per mode (p=π+φ incoming → 0)
    @inbounds begin
        c1 = p.twoOverh[1]; AL = p.AL
        for q in 1:n
            πi = U[1, 1, 2, q]; φi = U[1, 1, 3, q]
            qq = πi - φi; πg = qq / 2; φg = -qq / 2
            tmpπ[q] = πg; tmpφ[q] = φg
        end
        _AUvec!(Fπ, Fφ, AL, tmpπ, tmpφ, n)         # F(ghost)
        for q in 1:n
            tmpπ[q] = U[1, 1, 2, q]; tmpφ[q] = U[1, 1, 3, q]
        end
        _AUvec!(FπR, FφR, AL, tmpπ, tmpφ, n)       # F(interior)
        for q in 1:n
            πi = U[1, 1, 2, q]; φi = U[1, 1, 3, q]; qq = πi - φi
            πg = qq / 2; φg = -qq / 2
            Fsπ = 0.5 * (Fπ[q] + FπR[q]) - 0.5 * (πi - πg)
            Fsφ = 0.5 * (Fφ[q] + FφR[q]) - 0.5 * (φi - φg)
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
    solve_kerr_tail(p; init_mode, tfinal, dt, record_every, ic) -> (τs, Ψ)

Momentum-Gaussian IC in `ℓs[init_mode]`; evolve; return per-mode ψ time series at
the extraction points: `Ψ[probe, mode, sample]`.
"""
function solve_kerr_tail(p::KerrTailProblem{T}; init_mode=1, tfinal=20000.0, dt=0.04,
                         record_every=100, ic=(amp=1 / sqrt(2π * 6^2), ρ0=25.0, w=6.0)) where {T}
    Np, K, n = p.el.Np, p.mesh.K, p.n
    u = zeros(T, Np * K * 3 * n)
    U = reshape(u, Np, K, 3, n)
    @inbounds for k in 1:K, i in 1:Np
        U[i, k, 2, init_mode] = ic.amp * exp(-(p.ρ[i, k] - ic.ρ0)^2 / (2π * ic.w^2))
    end
    sc = (zeros(T, Np, n), zeros(T, Np, n), zeros(T, n), zeros(T, n),
        zeros(T, n), zeros(T, n), zeros(T, n), zeros(T, n))
    nprb = length(p.probe_ρ)
    τs = T[]; cols = Array{T,2}[]
    rec(uu, t, step) = begin
        if step % record_every == 0
            UU = reshape(uu, Np, K, 3, n)
            push!(τs, t)
            v = Matrix{T}(undef, nprb, n)
            for ip in 1:nprb, q in 1:n
                v[ip, q] = dot(view(p.probe_L, ip, :), view(UU, :, p.probe_k[ip], 1, q))
            end
            push!(cols, v)
        end
    end
    rhs!(du, uu, t) = _rhs!(du, uu, t, p, sc)
    rk4!(rhs!, u, (0.0, tfinal), dt; callback=rec)
    Ψ = Array{T,3}(undef, nprb, n, length(τs))
    for j in eachindex(cols), ip in 1:nprb, q in 1:n
        Ψ[ip, q, j] = cols[j][ip, q]
    end
    return τs, Ψ
end
