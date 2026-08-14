# Single-mode (`a=0`, no inter-mode coupling) hyperboloidal DG evolution of the rescaled scalar
# wave/Teukolsky equation with a real background potential `V_ℓ(r)` (Paper 1 §II, System-1
# reduction) — port of DG-Teulolsky/src/Mode1D.jl. Reuses `NumericsCore.dg` + `hyperboloidal.jl`.
# Used for the Schwarzschild **Price-tail** test (Paper 1 §IV.B) and, with the source added in
# flux1d.jl, the Schwarzschild circular-orbit fluxes (§IV.C).
#
# System (single mode):  ψ̇=-π;  π̇-Hφ̇=-(1-H)φ'+C ψ;  φ̇-Hπ̇=-(1-H)π',
# with flux matrix `A=1/(1+H)[[H,1],[1,H]]` and potential coefficient
# `C(ρ) = -V_ℓ(r)/(1-H²)` (finite at ρ=s by the analytic limit).

"""
    ρ_of_rstar(L::Layer, rstar; ρlo=-1e6, ρhi=nothing) -> ρ

ρ such that `r_*(ρ) = ρ/Ω(ρ)` equals the target tortoise value (bisection on `[ρlo, s)`).
"""
function ρ_of_rstar(L::Layer, rstar; ρlo=-1e6, ρhi=nothing)
    hi = ρhi === nothing ? L.s - 1e-9 : ρhi
    rstar < ρlo && return ρlo
    f(ρ) = rstar_of_ρ(L, ρ) - rstar
    a, b = ρlo, hi
    fa = f(a)
    for _ in 1:200
        m = (a + b) / 2
        fm = f(m)
        (abs(fm) < 1e-10 || (b - a) < 1e-12) && return m
        if sign(fm) == sign(fa)
            a, fa = m, fm
        else
            b = m
        end
    end
    return (a + b) / 2
end

"""
    Mode1DProblem{T}

Precomputed single-mode Schwarzschild hyperboloidal DG problem (Paper 1 §IV.B): reference
element, mesh (layer start `R` forced), layer, node arrays `ρ`/`H`/potential coefficient `C`,
lifting rows, extraction-probe interpolation rows and RHS scratch. Build with
[`Mode1DProblem(ℓ; ...)`](@ref Mode1DProblem(::Any)); evolve with [`solve_mode1d`](@ref).
"""
struct Mode1DProblem{T<:AbstractFloat}
    ℓ::Int
    el::RefElement{T}
    mesh::Mesh{T}
    layer::Layer{T}
    M::T
    ρ::Matrix{T}
    H::Matrix{T}
    C::Matrix{T}          # potential coefficient -V_ℓ(r)/(1-H²) (with scri limit)
    invMc1::Vector{T}
    invMcN::Vector{T}
    twoOverh::Vector{T}
    HL::T
    probe_ρ::Vector{T}
    probe_k::Vector{Int}
    probe_L::Matrix{T}
    Dpi::Vector{T}
    Dphi::Vector{T}
end

"""
    Mode1DProblem(ℓ; N, K, ρL, R, s, P=4, M=1, extract_rstar=[500.0, Inf], T=Float64)

Build the Schwarzschild single-mode problem. `extract_rstar` are tortoise-coordinate extraction
locations (e.g. `500.0` inside the layer, and `Inf` for 𝓘⁺=ρ=s).
"""
function Mode1DProblem(ℓ; N=10, K=400, ρL=-200.0, R=150.0, s=1200.0, P=4, M=1.0,
                       extract_rstar=[500.0, Inf], T::Type{<:AbstractFloat}=Float64)
    el = RefElement(N; T=T)
    layer = Layer{T}(T(R), T(s), P)
    mesh = build_mesh(T(ρL), T(s), K; forced=[T(R)], T=T)
    ρ = permutedims(node_coords(mesh, el))
    Hm = similar(ρ); C = similar(ρ)
    ℓℓ = T(ℓ * (ℓ + 1))
    @inbounds for k in 1:mesh.K, i in 1:el.Np
        ρv = ρ[i, k]
        Hv = Hfun(layer, ρv); Hm[i, k] = Hv
        if k == mesh.K && i == el.Np      # ρ=s node: analytic limit (2M terms vanish)
            C[i, k] = ℓℓ * denomΩ(layer, ρv) / (ρv^2 * (1 + Hv))
        else
            rs = rstar_of_ρ(layer, ρv)
            r = _inverse_tortoise(rs; M=M, a=0.0)
            om = 1 - Hv                    # 1-H = omega_factor
            C[i, k] = -potential_V(r, ℓ, 0; M=M, a=0.0) / (om * (1 + Hv))
        end
    end
    probes = T.(extract_rstar)
    pρ = [isinf(p) ? layer.s : ρ_of_rstar(layer, p; ρlo=T(ρL)) for p in probes]
    pk = [element_of(mesh, x) for x in pρ]
    pL = Matrix{T}(undef, length(pρ), el.Np)
    for (ip, x) in enumerate(pρ)
        k = pk[ip]
        ξ = 2 * (x - mesh.xL[k]) / mesh.h[k] - 1
        pL[ip, :] = _lagrange_row(el, ξ)
    end
    return Mode1DProblem{T}(ℓ, el, mesh, layer, T(M), ρ, Hm, C,
        el.invM[:, 1], el.invM[:, el.Np], 2 ./ mesh.h, Hfun(layer, mesh.edges[1]),
        pρ, pk, pL, Vector{T}(undef, el.Np), Vector{T}(undef, el.Np))
end

function _rhs!(du, u, ::Real, p::Mode1DProblem{T}) where {T}
    el = p.el; Np = el.Np; K = p.mesh.K
    U = reshape(u, Np, K, 3); dU = reshape(du, Np, K, 3)
    Dpi = p.Dpi; Dphi = p.Dphi
    @inbounds for k in 1:K
        invJ = p.twoOverh[k]
        mul!(Dpi, el.D, view(U, :, k, 2)); mul!(Dphi, el.D, view(U, :, k, 3))
        for i in 1:Np
            H = p.H[i, k]; iH = 1 / (1 + H)
            dπ = invJ * Dpi[i]; dφ = invJ * Dphi[i]
            AUπ = iH * (H * dπ + dφ); AUφ = iH * (dπ + H * dφ)
            ψ = U[i, k, 1]
            dU[i, k, 1] = -U[i, k, 2]
            dU[i, k, 2] = -AUπ + p.C[i, k] * ψ
            dU[i, k, 3] = -AUφ + H * p.C[i, k] * ψ
        end
    end
    @inbounds for k in 1:(K - 1)
        H = p.H[Np, k]
        πL = U[Np, k, 2]; φL = U[Np, k, 3]; πR = U[1, k + 1, 2]; φR = U[1, k + 1, 3]
        FLπ, FLφ = _AU(H, πL, φL); FRπ, FRφ = _AU(H, πR, φR)
        Fsπ = 0.5 * (FLπ + FRπ) - 0.5 * (πR - πL)
        Fsφ = 0.5 * (FLφ + FRφ) - 0.5 * (φR - φL)
        dRπ = FLπ - Fsπ; dRφ = FLφ - Fsφ; dLπ = FRπ - Fsπ; dLφ = FRφ - Fsφ
        ck = p.twoOverh[k]; ck1 = p.twoOverh[k + 1]
        for i in 1:Np
            dU[i, k, 2] += ck * p.invMcN[i] * dRπ
            dU[i, k, 3] += ck * p.invMcN[i] * dRφ
            dU[i, k + 1, 2] -= ck1 * p.invMc1[i] * dLπ
            dU[i, k + 1, 3] -= ck1 * p.invMc1[i] * dLφ
        end
    end
    @inbounds begin                      # left boundary: weak Sommerfeld (zero incoming p=π+φ)
        H = p.HL; πi = U[1, 1, 2]; φi = U[1, 1, 3]
        q = πi - φi; πg = q / 2; φg = -q / 2
        FLπ, FLφ = _AU(H, πg, φg); FRπ, FRφ = _AU(H, πi, φi)
        Fsπ = 0.5 * (FLπ + FRπ) - 0.5 * (πi - πg)
        Fsφ = 0.5 * (FLφ + FRφ) - 0.5 * (φi - φg)
        dLπ = FRπ - Fsπ; dLφ = FRφ - Fsφ; c1 = p.twoOverh[1]
        for i in 1:Np
            dU[i, 1, 2] -= c1 * p.invMc1[i] * dLπ
            dU[i, 1, 3] -= c1 * p.invMc1[i] * dLφ
        end
    end
    return du
end

"""
    solve_mode1d(p; tfinal, dt, momentum_ic, record_every) -> (τs, Ψ)

Evolve with non-zero-momentum Gaussian initial data `π = amp·exp(-(ρ-ρ0)²/(2πw²))`,
`ψ=φ=0` (Paper 1 eq. non-zero-momentum-initial-data). Returns probe time series Ψ
(nprobe × nsamples) of ψ at the extraction locations.
"""
function solve_mode1d(p::Mode1DProblem{T}; tfinal=10000.0, dt=0.0866,
                      momentum_ic=(amp=1 / sqrt(2π * 10^2), ρ0=30.0, w=10.0),
                      record_every=200) where {T}
    Np, K = p.el.Np, p.mesh.K
    u = zeros(T, Np * K * 3)
    U = reshape(u, Np, K, 3)
    @inbounds for k in 1:K, i in 1:Np
        ρv = p.ρ[i, k]
        U[i, k, 2] = momentum_ic.amp * exp(-(ρv - momentum_ic.ρ0)^2 / (2π * momentum_ic.w^2))
    end
    nprb = length(p.probe_ρ)
    τs = T[]; cols = Vector{T}[]
    rec(uu, t, step) = begin
        if step % record_every == 0
            UU = reshape(uu, Np, K, 3)
            push!(τs, t)
            v = [dot(view(p.probe_L, ip, :), view(UU, :, p.probe_k[ip], 1)) for ip in 1:nprb]
            push!(cols, v)
        end
    end
    rhs!(du, uu, t) = _rhs!(du, uu, t, p)
    rk4!(rhs!, u, (0.0, tfinal), dt; callback=rec)
    Ψ = Matrix{T}(undef, nprb, length(τs))
    for j in eachindex(cols), ip in 1:nprb
        Ψ[ip, j] = cols[j][ip]
    end
    return τs, Ψ
end
