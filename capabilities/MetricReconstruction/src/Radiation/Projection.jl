# =============================================================================================================
# Projection — tetrad→coordinate projection of the CCK reconstruction (legacy MetricReconstruction/projection.jl
# + projection_jet.jl). Build the 4×4 tetrad-frame matrix H_ab (index order l,n,m,m̄) from `component_value`,
# form the covariant DUAL co-tetrad legs (e^a)_μ from the Spacetime Kinnersley tetrad + Kerr metric (BHPToolkit
# CoordProjBL map e^l=−n_down, e^n=−l_down, e^m=m̄_down, e^m̄=m_down), and contract h_μν = Σ_{ab} (e^a)_μ(e^b)_ν H_ab.
#
# The covariant legs DELEGATE the Kinnersley tetrad + metric to Spacetime (kinnersley_tetrad/kerr_metric). The
# analytic-jet projection (`h_coordinate_jet`, RTJet r,θ 2-jets feeding the Einstein residual) builds its own
# RTJet metric/tetrad legs (the Spacetime closed forms transcribed onto the Radiation RTJet algebra) and the
# _jet H-blocks (HIRG_*_jet/HORG_*_jet, green TetradComponents). Type-generic {T}.
# =============================================================================================================

import Spacetime: kinnersley_tetrad, kerr_metric

export tetrad_matrix, covariant_legs, h_coordinate, h_coordinate_jet

# the nonzero tetrad components per gauge, with the index pair (i,j) in {l=1,n=2,m=3,m̄=4}
const _IRG_COMPS = ((:nn, 2, 2), (:nm, 2, 3), (:nmbar, 2, 4), (:mm, 3, 3), (:mbarmbar, 4, 4))
const _ORG_COMPS = ((:ll, 1, 1), (:lm, 1, 3), (:lmbar, 1, 4), (:mm, 3, 3), (:mbarmbar, 4, 4))

# legacy projection.jl used coords ∈ :BL,:Ingoing,:Outgoing; Spacetime uses :BL,:In,:Out.
_coords_st(c::Symbol) = c === :Ingoing ? :In : c === :Outgoing ? :Out : c

"""
    tetrad_matrix(gauge, a, ω, l, m, r, θ, t, φ, weyl; bc=:Up) -> 4×4 Complex tetrad-frame H_ab

Symmetric 4×4 matrix in tetrad-frame index order (l,n,m,m̄). Only the gauge's nonzero components are filled
(IRG: nn,nm,nm̄,mm,m̄m̄; ORG: ll,lm,lm̄,mm,m̄m̄), each from `component_value`.
"""
function tetrad_matrix(gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ, t, φ, weyl; bc::Symbol=:Up)
    C = typeof(complex(float(a*real(ω)*r*θ)))
    H = zeros(C, 4, 4)
    comps = gauge === :IRG ? _IRG_COMPS : gauge === :ORG ? _ORG_COMPS :
            error("unknown gauge $gauge")
    for (comp, i, j) in comps
        v = component_value(gauge, comp, a, ω, l, m, r, θ, t, φ, weyl; bc=bc)
        H[i, j] = v
        i != j && (H[j, i] = v)
    end
    H
end

"""
    covariant_legs(coords, a, r, θ) -> D (4×4 Complex)

Columns are the DUAL CO-TETRAD one-forms (e^a)_μ in tetrad-index order (l,n,m,m̄), so that h_μν =
Σ_{a,b}(e^a)_μ(e^b)_ν H_ab. For a null tetrad with l·n=−1, m·m̄=+1 the dual basis is
e^l=−n_down, e^n=−l_down, e^m=m̄_down, e^m̄=m_down (the BHPToolkit CoordProjBL map). Kinnersley tetrad +
metric are DELEGATED to Spacetime.
"""
function covariant_legs(coords::Symbol, a::T, r, θ) where {T}
    cs = _coords_st(coords)
    lup, nup, mup, mbarup = kinnersley_tetrad(a, r, θ; coords=cs)
    g = kerr_metric(a, r, θ; coords=cs)
    lo(v) = [sum(g[μ,ν]*v[ν] for ν in 1:4) for μ in 1:4]
    ldn = lo(lup); ndn = lo(nup); mdn = lo(mup); mbardn = lo(mbarup)
    D = zeros(Complex{T}, 4, 4)
    D[:,1] = -ndn          # e^l   = −n_down
    D[:,2] = -ldn          # e^n   = −l_down
    D[:,3] =  mbardn       # e^m   =  m̄_down
    D[:,4] =  mdn          # e^m̄  =  m_down
    D
end

"""
    h_coordinate(gauge, coords, a, ω, l, m, r, θ, t, φ, weyl; bc=:Up) -> 4×4 Complex coordinate h_μν

Project the tetrad-frame H_ab into coordinate components (coords ∈ :BL,:Ingoing,:Outgoing) at (t,r,θ,φ).
"""
function h_coordinate(gauge::Symbol, coords::Symbol, a, ω, l::Integer, m::Integer, r, θ, t, φ, weyl; bc::Symbol=:Up)
    H = tetrad_matrix(gauge, a, ω, l, m, r, θ, t, φ, weyl; bc=bc)
    Z = covariant_legs(coords, a, r, θ)
    C = eltype(H)
    h = zeros(C, 4, 4)
    @inbounds for μ in 1:4, ν in 1:4
        acc = zero(C)
        for ai in 1:4, bj in 1:4
            H[ai, bj] == 0 && continue
            acc += Z[μ, ai] * Z[ν, bj] * H[ai, bj]
        end
        h[μ, ν] = acc
    end
    h
end

# ─── ANALYTIC-JET projection (RTJet r,θ 2-jets); feeds the analytic Einstein residual ────────────────────────

# RTJet of the dual co-tetrad legs D[μ,a] (tetrad order a=l,n,m,m̄) in coords ∈ :BL,:Ingoing,:Outgoing.
function _legs_rtjet(coords::Symbol, ::Type{C}, a, r, θ) where {C}
    sθ = _rt_sin(C,θ); cθ = _rt_cos(C,θ); rj = _rt_r(C,r)
    Δ = _rt_Delta(C,a,r); Σ = _rt_Sigma(C,a,r,θ); ζ = _rt_zeta(C,a,r,θ); ζb = _rt_zetabar(C,a,r,θ)
    rt2 = sqrt(C(2))
    Z = rtconst(zero(C))
    a2 = rtconst(C(a^2)); rr = rj*rj
    r2a2 = rr + a2
    g = Matrix{RTJet{C}}(undef,4,4); for i in 1:4, j in 1:4; g[i,j]=Z; end
    s2 = sθ*sθ
    if coords === :BL
        g[1,1] = rtconst(C(-1)) + rtconst(C(2))*rj*rtinv(Σ)
        g[1,4] = g[4,1] = rtconst(C(-2*a))*rj*s2*rtinv(Σ)
        g[2,2] = Σ*rtinv(Δ)
        g[3,3] = Σ
        g[4,4] = (r2a2 + rtconst(C(2*a^2))*rj*s2*rtinv(Σ))*s2
        lup = [r2a2*rtinv(Δ), rtconst(one(C)), Z, rtconst(C(a))*rtinv(Δ)]
        nup = [r2a2, rtconst(C(-1))*Δ, Z, rtconst(C(a))]
        nup = [rtconst(C(1)/(2))*rtinv(Σ)*x for x in nup]
    elseif coords === :Ingoing
        g[1,1] = rtconst(C(-1)) + rtconst(C(2))*rj*rtinv(Σ)
        g[1,2] = g[2,1] = rtconst(one(C))
        g[1,4] = g[4,1] = rtconst(C(-2*a))*rj*s2*rtinv(Σ)
        g[2,4] = g[4,2] = rtconst(C(-a))*s2
        g[3,3] = Σ
        g[4,4] = ((r2a2*r2a2 - Δ*a2*s2)*rtinv(Σ))*s2
        lup = [rtconst(C(2))*r2a2*rtinv(Δ), rtconst(one(C)), Z, rtconst(C(2*a))*rtinv(Δ)]
        nup = [Z, rtconst(C(-1))*Δ*rtinv(rtconst(C(2))*Σ), Z, Z]
    elseif coords === :Outgoing
        g[1,1] = rtconst(C(-1)) + rtconst(C(2))*rj*rtinv(Σ)
        g[1,2] = g[2,1] = rtconst(C(-1))
        g[1,4] = g[4,1] = rtconst(C(-2*a))*rj*s2*rtinv(Σ)
        g[2,4] = g[4,2] = rtconst(C(a))*s2
        g[3,3] = Σ
        g[4,4] = ((r2a2*r2a2 - Δ*a2*s2)*rtinv(Σ))*s2
        lup = [Z, rtconst(one(C)), Z, Z]
        nup = [r2a2, rtconst(C(-1)/2)*Δ, Z, rtconst(C(a))]
        nup = [rtinv(Σ)*x for x in nup]
    else
        error("unknown coords $coords")
    end
    ias = rtconst(C(im*a))*sθ
    icsc = rtconst(C(im))*rtinv(sθ)
    mup    = [ias, Z, rtconst(one(C)), icsc];  mup    = [rtinv(rt2*ζb)*x for x in mup]
    mbarup = [rtconst(C(-1))*ias, Z, rtconst(one(C)), rtconst(C(-1))*icsc]; mbarup = [rtinv(rt2*ζ)*x for x in mbarup]
    lo(v) = [reduce(+, (g[μ,ν]*v[ν] for ν in 1:4)) for μ in 1:4]
    ldn=lo(lup); ndn=lo(nup); mdn=lo(mup); mbardn=lo(mbarup)
    D = Matrix{RTJet{C}}(undef,4,4)
    for μ in 1:4
        D[μ,1] = rtconst(C(-1))*ndn[μ]
        D[μ,2] = rtconst(C(-1))*ldn[μ]
        D[μ,3] = mbardn[μ]
        D[μ,4] = mdn[μ]
    end
    D
end

# RTJet of the four H-blocks (mirrors _hblocks) for one component.
function _hblocks_jet(gauge::Symbol, comp::Symbol, a, ω, l::Integer, m::Integer, r, θ; bc::Symbol=:Up)
    mω = -ω
    C = typeof(complex(float(a*real(ω)*r*θ)))
    Z = rtconst(zero(C))
    if gauge === :IRG
        if comp === :nn
            pA = HIRG_nn_jet(a,ω,l,m,r,θ;barred=false,bc=bc); pB = HIRG_nn_jet(a,mω,l,-m,r,θ;barred=true,bc=bc)
            mA = HIRG_nn_jet(a,ω,l,m,r,θ;barred=true,bc=bc);  mB = HIRG_nn_jet(a,mω,l,-m,r,θ;barred=false,bc=bc)
        elseif comp === :nm
            pA = Z; pB = HIRG_nm_jet(a,mω,l,-m,r,θ;barred=true,bc=bc); mA = HIRG_nm_jet(a,ω,l,m,r,θ;barred=true,bc=bc); mB = Z
        elseif comp === :nmbar
            pA = HIRG_nm_jet(a,ω,l,m,r,θ;barred=false,bc=bc); pB = Z; mA = Z; mB = HIRG_nm_jet(a,mω,l,-m,r,θ;barred=false,bc=bc)
        elseif comp === :mm
            pA = Z; pB = HIRG_mm_jet(a,mω,l,-m,r,θ;barred=true,bc=bc); mA = HIRG_mm_jet(a,ω,l,m,r,θ;barred=true,bc=bc); mB = Z
        elseif comp === :mbarmbar
            pA = HIRG_mm_jet(a,ω,l,m,r,θ;barred=false,bc=bc); pB = Z; mA = Z; mB = HIRG_mm_jet(a,mω,l,-m,r,θ;barred=false,bc=bc)
        else error("IRG component $comp zero/unknown") end
    elseif gauge === :ORG
        if comp === :ll
            pA = HORG_ll_jet(a,ω,l,m,r,θ;barred=false,bc=bc); pB = HORG_ll_jet(a,mω,l,-m,r,θ;barred=true,bc=bc)
            mA = HORG_ll_jet(a,ω,l,m,r,θ;barred=true,bc=bc);  mB = HORG_ll_jet(a,mω,l,-m,r,θ;barred=false,bc=bc)
        elseif comp === :lm
            pA = HORG_lm_jet(a,ω,l,m,r,θ;barred=false,bc=bc); pB = Z; mA = Z; mB = HORG_lm_jet(a,mω,l,-m,r,θ;barred=false,bc=bc)
        elseif comp === :lmbar
            pA = Z; pB = HORG_lm_jet(a,mω,l,-m,r,θ;barred=true,bc=bc); mA = HORG_lm_jet(a,ω,l,m,r,θ;barred=true,bc=bc); mB = Z
        elseif comp === :mm
            pA = HORG_mm_jet(a,ω,l,m,r,θ;barred=false,bc=bc); pB = Z; mA = Z; mB = HORG_mm_jet(a,mω,l,-m,r,θ;barred=false,bc=bc)
        elseif comp === :mbarmbar
            pA = Z; pB = HORG_mm_jet(a,mω,l,-m,r,θ;barred=true,bc=bc); mA = HORG_mm_jet(a,ω,l,m,r,θ;barred=true,bc=bc); mB = Z
        else error("ORG component $comp zero/unknown") end
    else error("unknown gauge $gauge") end
    (plusA=pA, plusB=pB, minusA=mA, minusB=mB)
end

# tetrad-frame H matrices (plus,minus) as RTJets (index order l,n,m,m̄).
function _tetrad_matrix_jet(gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ, weyl; bc::Symbol=:Up)
    C = typeof(complex(float(a*real(ω)*r*θ)))
    Z = rtconst(zero(C))
    Hp = Matrix{RTJet{C}}(undef,4,4); Hm = Matrix{RTJet{C}}(undef,4,4)
    for i in 1:4, j in 1:4; Hp[i,j]=Z; Hm[i,j]=Z; end
    comps = gauge === :IRG ? _IRG_COMPS : _ORG_COMPS
    for (comp,i,j) in comps
        b = _hblocks_jet(gauge, comp, a, ω, l, m, r, θ; bc=bc)
        hp = (weyl.A)*b.plusA + (weyl.Bbar)*b.plusB
        hm = (weyl.Abar)*b.minusA + (weyl.B)*b.minusB
        Hp[i,j]=hp; Hm[i,j]=hm
        if i != j; Hp[j,i]=hp; Hm[j,i]=hm; end
    end
    (Hp, Hm)
end

"""
    h_coordinate_jet(gauge, coords, a, ω, l, m, r, θ, weyl; bc=:Up) -> (hp, hm)

The coordinate h_plus/h_minus 4×4 RTJets (analytic r,θ 2-jets): h_{plus,μν}=Σ_{ab} D[μ,a]D[ν,b] Hp[a,b], same
for minus, with D the RTJet dual co-tetrad legs and Hp/Hm the ±m RTJet H-matrices. .v matches `h_coordinate`'s
plus/minus split; feeds the analytic linearized-Einstein residual.
"""
function h_coordinate_jet(gauge::Symbol, coords::Symbol, a, ω, l::Integer, m::Integer, r, θ, weyl; bc::Symbol=:Up)
    C = typeof(complex(float(a*real(ω)*r*θ)))
    D = _legs_rtjet(coords, C, a, r, θ)
    Hp, Hm = _tetrad_matrix_jet(gauge, a, ω, l, m, r, θ, weyl; bc=bc)
    Z = rtconst(zero(C))
    hp = Matrix{RTJet{C}}(undef,4,4); hm = Matrix{RTJet{C}}(undef,4,4)
    for μ in 1:4, ν in 1:4
        accp = Z; accm = Z
        for ai in 1:4, bj in 1:4
            DD = D[μ,ai]*D[ν,bj]
            accp = accp + DD*Hp[ai,bj]
            accm = accm + DD*Hm[ai,bj]
        end
        hp[μ,ν]=accp; hm[μ,ν]=accm
    end
    (hp, hm)
end
