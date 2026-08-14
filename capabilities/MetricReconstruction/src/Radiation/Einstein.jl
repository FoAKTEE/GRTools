# =============================================================================================================
# Einstein — the linearized (Lichnerowicz) Einstein tensor G⁽¹⁾_μν[h] of the reconstructed coordinate h, and the
# vacuum residual gate (legacy MetricReconstruction/einstein.jl, the LinEinsteinab oracle):
#   ∇_μ h_νρ = ∂_μ h_νρ − Γ^σ_μν h_σρ − Γ^σ_μρ h_νσ ;  ∇_α∇_μ h_νρ = ∂_α(∇_μ h_νρ) − Γ·∇h terms
#   G⁽¹⁾_μν = −½□h_μν − ½∇_μ∇_ν h + ∇^ρ∇_(μ h_ν)ρ + ½ g_μν(□h − ∇^ρ∇^σ h_ρσ)
# For a CCK radiation-gauge reconstruction of a VACUUM perturbation about Kerr, G⁽¹⁾_μν[h] ≡ 0 — the strongest
# correctness gate (G-M4). The background g/g⁻¹/Γ and analytic ∂g are DELEGATED to Spacetime; the Kerr ∂Γ and the
# exact h-derivatives use the analytic RTJet (r,θ) 2-jet machinery (∂_t=∓iω, ∂_φ=±im exact on the ±m phase split,
# the Projection h_coordinate_jet). A finite-difference path (deriv=:fd) is kept as a diagnostic cross-check.
# Type-generic {T}.
# =============================================================================================================

import Spacetime: kerr_metric, kerr_metric_inv, christoffel
import Spacetime

export lin_einstein_residual

# ∂_r g, ∂_θ g (BL, M=1): DELEGATE to Spacetime._dg_BL. Returns (gr, gθ).
_dg_BL(a, r, θ) = Spacetime._dg_BL(a, r, θ)

# Christoffels Γ^λ_μν (BL): DELEGATE to Spacetime.christoffel (analytic ∂g, ∂_t=∂_φ=0).
christoffel_BL(a, r, θ) = christoffel(a, r, θ; coords=:BL)

# ─── ANALYTIC metric (r,θ) 2-jet + ∂Γ (RTJet; the Kerr ∂Γ Spacetime does not expose closed-form) ─────────────
"BL metric g_{μν} as a 4×4 of RTJets (r,θ 2-jets). M=1, signature −+++."
function _metric_BL_rtjet(::Type{C}, a, r, θ) where {C}
    rj = _rt_r(C,r); sθ = _rt_sin(C,θ); cθ = _rt_cos(C,θ)
    s2 = sθ*sθ
    Σ = _rt_Sigma(C,a,r,θ); Δ = _rt_Delta(C,a,r)
    a2 = rtconst(C(a^2)); r2a2 = rj*rj + a2
    Z = rtconst(zero(C))
    g = Matrix{RTJet{C}}(undef,4,4); for i in 1:4,j in 1:4; g[i,j]=Z; end
    g[1,1] = rtconst(C(-1)) + rtconst(C(2))*rj*rtinv(Σ)
    g[1,4] = g[4,1] = rtconst(C(-2*a))*rj*s2*rtinv(Σ)
    g[2,2] = Σ*rtinv(Δ)
    g[3,3] = Σ
    g[4,4] = (r2a2 + rtconst(C(2*a^2))*rj*s2*rtinv(Σ))*s2
    g
end

"Inverse BL metric as 4×4 RTJets (analytic (t,φ) 2×2 block + diagonal r,θ)."
function _metric_inv_BL_rtjet(::Type{C}, a, r, θ) where {C}
    g = _metric_BL_rtjet(C, a, r, θ)
    det2 = g[1,1]*g[4,4] - g[1,4]*g[1,4]
    idet2 = rtinv(det2)
    Z = rtconst(zero(C))
    gi = Matrix{RTJet{C}}(undef,4,4); for i in 1:4,j in 1:4; gi[i,j]=Z; end
    gi[1,1] = g[4,4]*idet2
    gi[4,4] = g[1,1]*idet2
    gi[1,4] = gi[4,1] = rtconst(C(-1))*g[1,4]*idet2
    gi[2,2] = rtinv(g[2,2])
    gi[3,3] = rtinv(g[3,3])
    gi
end

"Christoffel Γ^λ_μν + ∂_α Γ (RTJet metric). Returns (Γ[λ,μ,ν], dΓ[α,λ,μ,ν]) with α∈{2,3} analytic."
function _christoffel_and_d_BL_analytic(::Type{C}, a, r, θ) where {C}
    g  = _metric_BL_rtjet(C, a, r, θ)
    gi = _metric_inv_BL_rtjet(C, a, r, θ)
    dg  = zeros(C,4,4,4)
    ddg = zeros(C,4,4,4,4)
    dgi = zeros(C,4,4,4)
    giv = zeros(C,4,4)
    @inbounds for μ in 1:4, ν in 1:4
        dg[2,μ,ν]=g[μ,ν].r; dg[3,μ,ν]=g[μ,ν].θ
        ddg[2,2,μ,ν]=g[μ,ν].rr; ddg[3,3,μ,ν]=g[μ,ν].θθ
        ddg[2,3,μ,ν]=ddg[3,2,μ,ν]=g[μ,ν].rθ
        dgi[2,μ,ν]=gi[μ,ν].r; dgi[3,μ,ν]=gi[μ,ν].θ
        giv[μ,ν]=gi[μ,ν].v
    end
    Γ = zeros(C,4,4,4)
    @inbounds for λ in 1:4, μ in 1:4, ν in 1:4
        acc = zero(C)
        for σ in 1:4
            giv[λ,σ]==0 && continue
            acc += giv[λ,σ]*(dg[μ,σ,ν]+dg[ν,σ,μ]-dg[σ,μ,ν])
        end
        Γ[λ,μ,ν]=acc/2
    end
    dΓ = zeros(C,4,4,4,4)
    @inbounds for α in (2,3), λ in 1:4, μ in 1:4, ν in 1:4
        acc = zero(C)
        for σ in 1:4
            acc += dgi[α,λ,σ]*(dg[μ,σ,ν]+dg[ν,σ,μ]-dg[σ,μ,ν])
            acc += giv[λ,σ]*(ddg[α,μ,σ,ν]+ddg[α,ν,σ,μ]-ddg[α,σ,μ,ν])
        end
        dΓ[α,λ,μ,ν]=acc/2
    end
    (Γ, dΓ)
end

"∂_α Γ^λ_μν analytically (RTJet); thin wrapper used by lin_einstein_residual (deriv=:analytic)."
function _dchristoffel_BL_analytic(a::T, r, θ) where {T}
    C = typeof(complex(float(a*r*θ)))
    _, dΓ = _christoffel_and_d_BL_analytic(C, a, r, θ)
    dΓ
end

"∂_α Γ^λ_μν via central finite difference in (r,θ) (diagnostic; deriv=:fd)."
function _dchristoffel_BL(a::T, r, θ; hr=nothing, hθ=nothing) where {T}
    hr = hr === nothing ? (eps(T))^(T(1)/3)*max(abs(r), one(T)) : T(hr)
    hθ = hθ === nothing ? (eps(T))^(T(1)/3) : T(hθ)
    Γrp = christoffel_BL(a, r+hr, θ); Γrm = christoffel_BL(a, r-hr, θ)
    Γθp = christoffel_BL(a, r, θ+hθ); Γθm = christoffel_BL(a, r, θ-hθ)
    dΓ = zeros(T, 4, 4, 4, 4)
    dΓ[2,:,:,:] = (Γrp .- Γrm) ./ (2hr)
    dΓ[3,:,:,:] = (Γθp .- Γθm) ./ (2hθ)
    dΓ
end

# ─── h-field jets: h, ∂_α h, ∂_α∂_β h at a BL point ──────────────────────────────────────────────────────────
"""
    _h_jets_analytic(gauge, a, ω, l, m, r, θ, t, φ, weyl; bc=:Up) -> (h, dh, ddh)

EXACT analytic h[μ,ν], ∂_α h, ∂_α∂_β h at the BL point. Spatial (r,θ) derivatives from the RTJet projection
(`h_coordinate_jet`); temporal/azimuthal from the exact phases ∂_t=∓iω, ∂_φ=±im on the h_plus/h_minus split.
"""
function _h_jets_analytic(gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ, t, φ, weyl; bc::Symbol=:Up)
    hp, hm = h_coordinate_jet(gauge, :BL, a, ω, l, m, r, θ, weyl; bc=bc)
    C = typeof(hp[1,1].v)
    Ep = exp(-im*ω*t + im*m*φ); Em = exp(im*ω*t - im*m*φ)
    iωp = -im*ω; iωm = im*ω
    imp = im*m;  imm = -im*m
    h   = zeros(C,4,4)
    dh  = zeros(C,4,4,4)
    ddh = zeros(C,4,4,4,4)
    @inbounds for μ in 1:4, ν in 1:4
        P = hp[μ,ν]; Mn = hm[μ,ν]
        h[μ,ν] = P.v*Ep + Mn.v*Em
        dh[1,μ,ν] = iωp*P.v*Ep + iωm*Mn.v*Em
        dh[2,μ,ν] = P.r*Ep + Mn.r*Em
        dh[3,μ,ν] = P.θ*Ep + Mn.θ*Em
        dh[4,μ,ν] = imp*P.v*Ep + imm*Mn.v*Em
        ddh[1,1,μ,ν] = iωp^2*P.v*Ep + iωm^2*Mn.v*Em
        ddh[4,4,μ,ν] = imp^2*P.v*Ep + imm^2*Mn.v*Em
        ddh[2,2,μ,ν] = P.rr*Ep + Mn.rr*Em
        ddh[3,3,μ,ν] = P.θθ*Ep + Mn.θθ*Em
        ddh[1,4,μ,ν] = ddh[4,1,μ,ν] = iωp*imp*P.v*Ep + iωm*imm*Mn.v*Em
        ddh[1,2,μ,ν] = ddh[2,1,μ,ν] = iωp*P.r*Ep + iωm*Mn.r*Em
        ddh[1,3,μ,ν] = ddh[3,1,μ,ν] = iωp*P.θ*Ep + iωm*Mn.θ*Em
        ddh[4,2,μ,ν] = ddh[2,4,μ,ν] = imp*P.r*Ep + imm*Mn.r*Em
        ddh[4,3,μ,ν] = ddh[3,4,μ,ν] = imp*P.θ*Ep + imm*Mn.θ*Em
        ddh[2,3,μ,ν] = ddh[3,2,μ,ν] = P.rθ*Ep + Mn.rθ*Em
    end
    (h, dh, ddh)
end

"""
    _h_jets(gauge, a, ω, l, m, r, θ, t, φ, weyl; bc=:Up, hr, hθ) -> (h, dh, ddh)

Finite-difference h-jets (diagnostic; deriv=:fd). 4th-order central FD in r,θ,t,φ of `h_coordinate`.
"""
function _h_jets(gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ, t, φ, weyl;
                 bc::Symbol=:Up, hr=nothing, hθ=nothing)
    T = typeof(float(real(a*ω*r*θ)))
    hr = hr === nothing ? T(1e-3)*max(abs(r), one(T)) : T(hr)
    hθ = hθ === nothing ? T(1e-3) : T(hθ)
    Hc(rr, tt, tcoord, fcoord) = h_coordinate(gauge, :BL, a, ω, l, m, rr, tt, tcoord, fcoord, weyl; bc=bc)
    h0 = Hc(r, θ, t, φ)
    C = eltype(h0)
    Hrp1=Hc(r+hr,θ,t,φ); Hrm1=Hc(r-hr,θ,t,φ); Hrp2=Hc(r+2hr,θ,t,φ); Hrm2=Hc(r-2hr,θ,t,φ)
    hr1 = (-Hrp2 .+ 8 .*Hrp1 .- 8 .*Hrm1 .+ Hrm2) ./ (12*hr)
    hr2 = (-Hrp2 .+ 16 .*Hrp1 .- 30 .*h0 .+ 16 .*Hrm1 .- Hrm2) ./ (12*hr^2)
    Hθp1=Hc(r,θ+hθ,t,φ); Hθm1=Hc(r,θ-hθ,t,φ); Hθp2=Hc(r,θ+2hθ,t,φ); Hθm2=Hc(r,θ-2hθ,t,φ)
    hθ1 = (-Hθp2 .+ 8 .*Hθp1 .- 8 .*Hθm1 .+ Hθm2) ./ (12*hθ)
    hθ2 = (-Hθp2 .+ 16 .*Hθp1 .- 30 .*h0 .+ 16 .*Hθm1 .- Hθm2) ./ (12*hθ^2)
    Hpp=Hc(r+hr,θ+hθ,t,φ); Hpm=Hc(r+hr,θ-hθ,t,φ); Hmp=Hc(r-hr,θ+hθ,t,φ); Hmm=Hc(r-hr,θ-hθ,t,φ)
    hrθ = (Hpp .- Hpm .- Hmp .+ Hmm) ./ (4*hr*hθ)
    ht = T(1e-3); hpf = T(1e-3)
    Ht1=Hc(r,θ,t+ht,φ); Ht_1=Hc(r,θ,t-ht,φ); Ht2=Hc(r,θ,t+2ht,φ); Htm2=Hc(r,θ,t-2ht,φ)
    hdt = (-Ht2 .+ 8 .*Ht1 .- 8 .*Ht_1 .+ Htm2) ./ (12*ht)
    hdtt = (-Ht2 .+ 16 .*Ht1 .- 30 .*h0 .+ 16 .*Ht_1 .- Htm2) ./ (12*ht^2)
    Hf1=Hc(r,θ,t,φ+hpf); Hf_1=Hc(r,θ,t,φ-hpf); Hf2=Hc(r,θ,t,φ+2hpf); Hfm2=Hc(r,θ,t,φ-2hpf)
    hdf = (-Hf2 .+ 8 .*Hf1 .- 8 .*Hf_1 .+ Hfm2) ./ (12*hpf)
    hdff = (-Hf2 .+ 16 .*Hf1 .- 30 .*h0 .+ 16 .*Hf_1 .- Hfm2) ./ (12*hpf^2)
    dh = zeros(C, 4, 4, 4)
    dh[1,:,:] = hdt; dh[2,:,:] = hr1; dh[3,:,:] = hθ1; dh[4,:,:] = hdf
    ddh = zeros(C, 4, 4, 4, 4)
    ddh[1,1,:,:] = hdtt; ddh[2,2,:,:] = hr2; ddh[3,3,:,:] = hθ2; ddh[4,4,:,:] = hdff
    ddh[2,3,:,:] = hrθ; ddh[3,2,:,:] = hrθ
    function mixedFD(d1::Symbol, d2::Symbol)
        f(s1, s2) = Hc(r + (d1===:r ? s1 : (d2===:r ? s2 : zero(T))),
                       θ + (d1===:θ ? s1 : (d2===:θ ? s2 : zero(T))),
                       t + (d1===:t ? s1 : (d2===:t ? s2 : zero(T))),
                       φ + (d1===:φ ? s1 : (d2===:φ ? s2 : zero(T))))
        step(d) = d===:r ? hr : d===:θ ? hθ : d===:t ? ht : hpf
        d1a = step(d1); d2a = step(d2)
        (f(d1a,d2a) .- f(d1a,-d2a) .- f(-d1a,d2a) .+ f(-d1a,-d2a)) ./ (4*d1a*d2a)
    end
    for (i, di) in ((1,:t),(2,:r),(3,:θ),(4,:φ)), (j, dj) in ((1,:t),(2,:r),(3,:θ),(4,:φ))
        i >= j && continue
        (i,j) == (2,3) && continue
        ((i==1&&j==1)||(i==2&&j==2)||(i==3&&j==3)||(i==4&&j==4)) && continue
        mij = mixedFD(di, dj)
        ddh[i,j,:,:] = mij; ddh[j,i,:,:] = mij
    end
    (h0, dh, ddh)
end

# ─── the linearized Einstein tensor at a BL point ────────────────────────────────────────────────────────────
"""
    lin_einstein_residual(gauge, a, ω, l, m, r, θ, t, φ, weyl; bc=:Up, deriv=:analytic) -> 4×4 Complex G⁽¹⁾_μν[h]

Linearized Einstein tensor of the reconstructed h at the BL point. For a vacuum CCK reconstruction this is ≈ 0
(G-M4). Built from delegated Spacetime Christoffels (analytic ∂g) + EXACT analytic second derivatives of h
(deriv=:analytic, default; the RTJet projection). deriv=:fd uses the finite-difference sampler (diagnostic).
"""
function lin_einstein_residual(gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ, t, φ, weyl;
                               bc::Symbol=:Up, deriv::Symbol=:analytic)
    g  = kerr_metric(a, r, θ; coords=:BL)
    gi = kerr_metric_inv(a, r, θ; coords=:BL)
    Γ  = christoffel_BL(a, r, θ)
    dΓ = deriv === :analytic ? _dchristoffel_BL_analytic(a, r, θ) : _dchristoffel_BL(a, r, θ)
    h, dh, ddh = deriv === :analytic ? _h_jets_analytic(gauge, a, ω, l, m, r, θ, t, φ, weyl; bc=bc) :
                                       _h_jets(gauge, a, ω, l, m, r, θ, t, φ, weyl; bc=bc)
    C = eltype(h)
    gC = C.(g)
    Ch = zeros(C, 4, 4, 4)
    @inbounds for μ in 1:4, ν in 1:4, ρ in 1:4
        acc = dh[μ, ν, ρ]
        for σ in 1:4
            acc -= Γ[σ,μ,ν]*h[σ,ρ] + Γ[σ,μ,ρ]*h[ν,σ]
        end
        Ch[μ,ν,ρ] = acc
    end
    CCh = zeros(C, 4, 4, 4, 4)
    @inbounds for α in 1:4, μ in 1:4, ν in 1:4, ρ in 1:4
        acc = ddh[α, μ, ν, ρ]
        for σ in 1:4
            acc -= dΓ[α,σ,μ,ν]*h[σ,ρ] + Γ[σ,μ,ν]*dh[α,σ,ρ]
            acc -= dΓ[α,σ,μ,ρ]*h[ν,σ] + Γ[σ,μ,ρ]*dh[α,ν,σ]
        end
        for σ in 1:4
            acc -= Γ[σ,α,μ]*Ch[σ,ν,ρ] + Γ[σ,α,ν]*Ch[μ,σ,ρ] + Γ[σ,α,ρ]*Ch[μ,ν,σ]
        end
        CCh[α,μ,ν,ρ] = acc
    end
    box(μ, ν) = sum(gi[α,β]*CCh[α,β,μ,ν] for α in 1:4, β in 1:4)
    trCC(μ, ν) = sum(gi[ρ,σ]*CCh[μ,ν,ρ,σ] for ρ in 1:4, σ in 1:4)
    divdiv(μ, ν) = sum(gi[ρ,σ]*CCh[σ,μ,ν,ρ] for ρ in 1:4, σ in 1:4)
    boxh = sum(gi[α,β]*trCC(α,β) for α in 1:4, β in 1:4)
    divdivh = sum(gi[μ,ν]*divdiv(μ,ν) for μ in 1:4, ν in 1:4)
    G = zeros(C, 4, 4)
    @inbounds for μ in 1:4, ν in 1:4
        G[μ,ν] = -box(μ,ν)/2 - trCC(μ,ν)/2 +
                 (divdiv(μ,ν) + divdiv(ν,μ))/2 +
                 gC[μ,ν]*(boxh - divdivh)/2
    end
    G
end
