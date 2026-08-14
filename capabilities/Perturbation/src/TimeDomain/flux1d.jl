# Schwarzschild (`a=0`) **circular-orbit scalar flux** solver (Paper 1 §IV.C, Paper 2 §IV.B) —
# port of DG-Teulolsky/src/Flux1D.jl. Single complex mode `(ℓ,m)` driven by a point particle
# through the **characteristic-projected delta-flux injection** at the particle interface
# (Paper 1 §III.C, eq. NumericalFluxModified):
#
#     -(F*)^{kp+1}_{left,mod}  = -(F*)_left  + T P⁺ T⁻¹ G(t)
#     -(F*)^{kp}_{right,mod}   = -(F*)_right - T P⁻ T⁻¹ G(t)
#
# At the particle (bulk, H=0) `A=[[0,1],[1,0]]`, giving `T P⁺T⁻¹=½[[1,1],[1,1]]`,
# `T P⁻T⁻¹=½[[1,-1],[-1,1]]`. Source on the π-row only: `G=[G_π,0]`.
# Energy flux at 𝓘⁺: `Ė_{ℓm}=|π_{ℓm}|²/(4π)`.
#
# REUSE SUBSTITUTIONS (repo → GrTool; DG-P3 directive):
#   Flux1D.Ylm_equator (local closed form)      → Harmonics.Ylm_real(ℓ, m, π/2) + exact parity 0
#   Flux1D.orbital_params (repo commit a643b5f) → Spacetime geodesic_energy /
#     geodesic_angular_momentum / kerr_frequencies (circular equatorial e=0, x=1; M=1 canon)
# Both characterized against baked repo values in test/fixtures/td_substitution.jl (≤1e-12 rel).

"""
    Ylm_equator(ℓ, m) -> Float64

Real value `Y_{ℓm}(π/2, 0)` of the scalar spherical harmonic at the equator (the `e^{imφ}`
phase is applied by the caller). DELEGATES to `Harmonics.Ylm_real` — same convention
(`|m|`, no Condon-Shortley sign flip for `m<0`, `√((2ℓ+1)/4π·(ℓ-|m|)!/(ℓ+|m|)!)` norm) — with
the repo's exact-zero parity shortcut for odd `ℓ+|m|` restored (at `θ=π/2` the recurrence
would otherwise leave an `O(1e-17)` residue from `cos(π/2) ≠ 0` in floating point).
"""
Ylm_equator(ℓ, m) = isodd(ℓ + abs(m)) ? 0.0 : Ylm_real(Int(ℓ), Int(m), π / 2)

"""
    orbital_params(rp; M=1, a=0) -> (E, Lz, ut, Ω)

Kerr circular equatorial (prograde) geodesic quantities at Boyer–Lindquist radius `rp`
(repo commit a643b5f closed forms: `v²=M/rp`, `E=(1−2v²+ãv³)/√(1−3v²+2ãv³)`,
`u^t=(1+ãv³)/√(1−3v²+2ãv³)`, `Ω=v³/(M(1+ãv³))`). DELEGATES to the `Spacetime` circular
`(e=0, x=1)` geodesic surface — `geodesic_energy`, `geodesic_angular_momentum` and the Mino
`kerr_frequencies` with `u^t = Υ_t/Σ` (`Σ = r_p²` on the equator), `Ω = Υ_φ/Υ_t` — normalized
from Spacetime's `M = 1` convention by the exact homogeneity scalings (`E`, `u^t`
dimensionless; `Lz ∝ M`; `Ω ∝ 1/M`).
"""
function orbital_params(rp; M=1.0, a=0.0)
    p = rp / M; ã = a / M
    E = geodesic_energy(ã, p, 0.0, 1.0)
    Lz = M * geodesic_angular_momentum(ã, p, 0.0, 1.0)
    f = kerr_frequencies(ã, p, 0.0, 1.0; Time="Mino")
    ut = f.Ut / p^2                     # u^t = Υ_t/Σ, Σ = r² (equatorial circular)
    Ω = f.Uph / (M * f.Ut)              # Ω_φ = Υ_φ/Υ_t (BL), rescaled by 1/M
    return (E=E, Lz=Lz, ut=ut, Ω=Ω)
end

"""
    Flux1DProblem{T}

Precomputed Schwarzschild circular-orbit flux problem (Paper 1 §IV.C): the [`Mode1DProblem`](@ref)
machinery plus the particle interface `kp` (`r_*,p` forced onto the mesh), the per-mode source
prefactor `Gpre = -4πq/(u^t r_p)·Y_{ℓm}(π/2,0)`, orbital frequency `Ω` and the smooth source
ramp. Build with [`Flux1DProblem(ℓ, m, rp; ...)`](@ref Flux1DProblem(::Any, ::Any, ::Any));
evolve with [`solve_flux`](@ref) / [`selfforce_ft`](@ref).
"""
struct Flux1DProblem{T<:AbstractFloat}
    ℓ::Int
    m::Int
    el::RefElement{T}
    mesh::Mesh{T}
    layer::Layer{T}
    M::T
    rp::T
    rsp::T                 # r_*,p
    H::Matrix{T}
    C::Matrix{T}
    invMc1::Vector{T}
    invMcN::Vector{T}
    twoOverh::Vector{T}
    HL::T
    kp::Int                # left element at the particle (right face = r_*,p)
    Gpre::Complex{T}       # -4πq/(ut rp) · Y_{ℓm}(π/2,0)
    Ω::T
    ramp_τ::T
    ramp_δ::T
    scri_k::Int
    Dpi::Vector{Complex{T}}
    Dphi::Vector{Complex{T}}
end

"""
    Flux1DProblem(ℓ, m, rp; N, K, ρL, R, s, P=4, M=1, q=1, ramp_τ=400, ramp_δ=0.00025, T=Float64)

Build the Schwarzschild circular-orbit flux problem for mode `(ℓ,m)` and orbital radius `rp`
(the particle tortoise radius `r_*,p` and the layer start `R` are forced mesh interfaces).
"""
function Flux1DProblem(ℓ, m, rp; N=10, K=200, ρL=-100.0, R=150.0, s=200.0, P=4,
                       M=1.0, q=1.0, ramp_τ=400.0, ramp_δ=0.00025,
                       T::Type{<:AbstractFloat}=Float64)
    rsp = _tortoise(rp; M=M, a=0.0)
    el = RefElement(N; T=T)
    layer = Layer{T}(T(R), T(s), P)
    mesh = build_mesh(T(ρL), T(s), K; forced=[T(rsp), T(R)], T=T)
    # locate the particle interface edge
    jp = findmin(abs.(mesh.edges .- T(rsp)))[2]
    kp = jp - 1
    ρ = permutedims(node_coords(mesh, el))
    Hm = similar(ρ); C = similar(ρ)
    ℓℓ = T(ℓ * (ℓ + 1))
    @inbounds for k in 1:mesh.K, i in 1:el.Np
        ρv = ρ[i, k]; Hv = Hfun(layer, ρv); Hm[i, k] = Hv
        if k == mesh.K && i == el.Np
            C[i, k] = ℓℓ * denomΩ(layer, ρv) / (ρv^2 * (1 + Hv))
        else
            r = _inverse_tortoise(rstar_of_ρ(layer, ρv); M=M, a=0.0)
            C[i, k] = -potential_V(r, ℓ, m; M=M, a=0.0) / ((1 - Hv) * (1 + Hv))
        end
    end
    op = orbital_params(rp; M=M)
    Gpre = Complex{T}(-4π * q / (op.ut * rp) * Ylm_equator(ℓ, m))
    return Flux1DProblem{T}(ℓ, m, el, mesh, layer, T(M), T(rp), T(rsp), Hm, C,
        el.invM[:, 1], el.invM[:, el.Np], 2 ./ mesh.h, Hfun(layer, mesh.edges[1]),
        kp, Gpre, T(op.Ω), T(ramp_τ), T(ramp_δ), mesh.K,
        Vector{Complex{T}}(undef, el.Np), Vector{Complex{T}}(undef, el.Np))
end

function _rhs!(du, u, t, p::Flux1DProblem{T}) where {T}
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
    # delta-flux injection at the particle interface (kp | kp+1)
    @inbounds begin
        Gπ = p.Gpre * cis(-p.m * p.Ω * t) * _ramp(t, p.ramp_τ, p.ramp_δ)   # G(t)=Gπ·e^{-imΩt}·ramp
        kp = p.kp
        # left element kp, right face: += (2/h)·invMcN·(-T P⁻ T⁻¹ G) = (2/h)invMcN·½[-Gπ, +Gπ]
        ckp = p.twoOverh[kp]
        for i in 1:Np
            dU[i, kp, 2] += ckp * p.invMcN[i] * (-0.5 * Gπ)
            dU[i, kp, 3] += ckp * p.invMcN[i] * (0.5 * Gπ)
        end
        # right element kp+1, left face: += -(2/h)·invMc1·(T P⁺ T⁻¹ G) = -(2/h)invMc1·½[Gπ, Gπ]
        ckp1 = p.twoOverh[kp + 1]
        for i in 1:Np
            dU[i, kp + 1, 2] -= ckp1 * p.invMc1[i] * (0.5 * Gπ)
            dU[i, kp + 1, 3] -= ckp1 * p.invMc1[i] * (0.5 * Gπ)
        end
    end
    @inbounds begin   # left boundary: weak Sommerfeld (zero incoming p=π+φ)
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
    solve_flux(p; tfinal, dt, record_every) -> (τs, Escri, Ehor)

Trivial initial data; evolve and record the energy flux at 𝓘⁺ `Ė^{𝓘⁺}=|π(scri)|²/(4π)`
and at the horizon side `Ė^{𝓗}=|π(left boundary)|²/(4π)`. The late-time values are the
circular-orbit fluxes `Ė_{ℓm}/q²` to infinity / into the horizon (Paper 2 §IV.B).
"""
function solve_flux(p::Flux1DProblem{T}; tfinal=8000.0, dt=0.0018, record_every=200) where {T}
    Np, K = p.el.Np, p.mesh.K
    u = zeros(Complex{T}, Np * K * 3)
    τs = T[]; Es = T[]; Eh = T[]
    rec(uu, t, step) = begin
        if step % record_every == 0
            UU = reshape(uu, Np, K, 3)
            push!(τs, t)
            push!(Es, abs2(UU[Np, K, 2]) / (4π))     # 𝓘⁺
            push!(Eh, abs2(UU[1, 1, 2]) / (4π))       # horizon side (left boundary)
        end
    end
    rhs!(du, uu, t) = _rhs!(du, uu, t, p)
    rk4!(rhs!, u, (0.0, tfinal), dt; callback=rec)
    return τs, Es, Eh
end

"""
    selfforce_ft(p; tfinal, dt, twin) -> (ft_outflow, ft_inflow)

Dissipative temporal self-force `F_t^{ℓm}` extracted from `∂_tψ=-π` AT the particle,
from the **left (outflow/Radau)** node and the **right (inflow)** node (Paper 2 §IV.B.2).
`F_t = (1/(r_p u^t)) Y_{ℓm}(π/2,0) · ⟨e^{imΩτ} π_p⟩` (steady-state, late-time mean over
`t > twin`). The outflow value is far more accurate (Radau outflow point).
"""
function selfforce_ft(p::Flux1DProblem{T}; tfinal=6000.0, dt=0.004, twin=4000.0) where {T}
    Np, K = p.el.Np, p.mesh.K
    u = zeros(Complex{T}, Np * K * 3)
    pref = Ylm_equator(p.ℓ, p.m) / p.rp
    accL = ComplexF64(0); accR = ComplexF64(0); cnt = 0
    rec(uu, t, step) = begin
        if step % 50 == 0 && t > twin
            UU = reshape(uu, Np, K, 3)
            ph = cis(p.m * p.Ω * t)
            accL += ph * UU[Np, p.kp, 2]        # left of particle (outflow, Radau)
            accR += ph * UU[1, p.kp + 1, 2]     # right of particle (inflow)
            cnt += 1
        end
    end
    rhs!(du, uu, t) = _rhs!(du, uu, t, p)
    rk4!(rhs!, u, (0.0, tfinal), dt; callback=rec)
    ftL = pref * real(accL / cnt)
    ftR = pref * real(accR / cnt)
    return ftL, ftR
end
