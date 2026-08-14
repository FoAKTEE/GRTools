# =============================================================================================================
# PlusMinus — assembly of a full tetrad component h_ab(t,r,θ,φ) from the six Hertz H-blocks (CCK reconstruction,
# legacy MetricReconstruction/plus_minus.jl). The ±m phase split (MetricPerturbations.nb cells 7/15/24/28):
#   h_ab(t,r,θ,φ) = h_ab^plus(r,θ) e^{−iωt+imφ} + h_ab^minus(r,θ) e^{+iω̄t−imφ},  ω̄ = ω (real-ω modes),
#   h_ab^plus  = A·h^{plus,A}  + B̄·h^{plus,B},   h_ab^minus = Ā·h^{minus,A} + B·h^{minus,B},
# with the A/B blocks mapped to the H-functions (HIRG_nn/nm/mm, HORG_ll/lm/mm — green TetradComponents) by the
# verified cell-24 (IRG) / cell-28 (ORG) substitution tables. `weyl` is the (A,Abar,B,Bbar) Weyl amplitude tuple
# (WeylAmplitudes). Leaf of the Radiation reconstruction tail; consumed by Projection. Type-generic {T}.
# =============================================================================================================

export component_value

"""
    _hblocks(gauge, comp, a, ω, l, m, r, θ; bc=:Up, λ=nothing) -> (plusA, plusB, minusA, minusB)

The four A/B building blocks of one tetrad component (gauge ∈ :IRG/:ORG; comp the nonzero component symbol),
per the cell-24 (IRG) / cell-28 (ORG) tables. The "−ω̄" arguments are taken as "−ω" (real-ω modes). Convention
(both gauges, diagonal block): plusA=H(ω,l,m), plusB=H̄(−ω,l,−m), minusA=H̄(ω,l,m), minusB=H(−ω,l,−m); the
off-diagonal blocks zero the corresponding entries per the verified tables.
"""
function _hblocks(gauge::Symbol, comp::Symbol, a, ω, l::Integer, m::Integer, r, θ; bc::Symbol=:Up, λ=nothing)
    mω = -ω
    Z = zero(complex(float(a*real(ω))))
    if gauge === :IRG
        if comp === :nn
            pA = HIRG_nn(a, ω, l, m, r, θ; barred=false, bc=bc)
            pB = HIRG_nn(a, mω, l, -m, r, θ; barred=true, bc=bc)
            mA = HIRG_nn(a, ω, l, m, r, θ; barred=true, bc=bc)
            mB = HIRG_nn(a, mω, l, -m, r, θ; barred=false, bc=bc)
        elseif comp === :nm
            pA = Z
            pB = HIRG_nm(a, mω, l, -m, r, θ; barred=true, bc=bc)
            mA = HIRG_nm(a, ω, l, m, r, θ; barred=true, bc=bc)
            mB = Z
        elseif comp === :nmbar
            pA = HIRG_nm(a, ω, l, m, r, θ; barred=false, bc=bc)
            pB = Z
            mA = Z
            mB = HIRG_nm(a, mω, l, -m, r, θ; barred=false, bc=bc)
        elseif comp === :mm
            pA = Z
            pB = HIRG_mm(a, mω, l, -m, r, θ; barred=true, bc=bc)
            mA = HIRG_mm(a, ω, l, m, r, θ; barred=true, bc=bc)
            mB = Z
        elseif comp === :mbarmbar
            pA = HIRG_mm(a, ω, l, m, r, θ; barred=false, bc=bc)
            pB = Z
            mA = Z
            mB = HIRG_mm(a, mω, l, -m, r, θ; barred=false, bc=bc)
        else
            error("IRG component $comp is identically zero or unknown")
        end
    elseif gauge === :ORG
        if comp === :ll
            pA = HORG_ll(a, ω, l, m, r, θ; barred=false, bc=bc)
            pB = HORG_ll(a, mω, l, -m, r, θ; barred=true, bc=bc)
            mA = HORG_ll(a, ω, l, m, r, θ; barred=true, bc=bc)
            mB = HORG_ll(a, mω, l, -m, r, θ; barred=false, bc=bc)
        elseif comp === :lm
            pA = HORG_lm(a, ω, l, m, r, θ; barred=false, bc=bc)
            pB = Z
            mA = Z
            mB = HORG_lm(a, mω, l, -m, r, θ; barred=false, bc=bc)
        elseif comp === :lmbar
            pA = Z
            pB = HORG_lm(a, mω, l, -m, r, θ; barred=true, bc=bc)
            mA = HORG_lm(a, ω, l, m, r, θ; barred=true, bc=bc)
            mB = Z
        elseif comp === :mm
            pA = HORG_mm(a, ω, l, m, r, θ; barred=false, bc=bc)
            pB = Z
            mA = Z
            mB = HORG_mm(a, mω, l, -m, r, θ; barred=false, bc=bc)
        elseif comp === :mbarmbar
            pA = Z
            pB = HORG_mm(a, mω, l, -m, r, θ; barred=true, bc=bc)
            mA = HORG_mm(a, ω, l, m, r, θ; barred=true, bc=bc)
            mB = Z
        else
            error("ORG component $comp is identically zero or unknown")
        end
    else
        error("unknown gauge $gauge (expected :IRG or :ORG)")
    end
    (plusA=pA, plusB=pB, minusA=mA, minusB=mB)
end

"""
    component_value(gauge, comp, a, ω, l, m, r, θ, t, φ, weyl; bc=:Up) -> Complex

The full tetrad component h_comp(t,r,θ,φ) (cells 7/15/24/28). `weyl` is a NamedTuple (A,Abar,B,Bbar).
"""
function component_value(gauge::Symbol, comp::Symbol, a, ω, l::Integer, m::Integer, r, θ, t, φ, weyl; bc::Symbol=:Up)
    b = _hblocks(gauge, comp, a, ω, l, m, r, θ; bc=bc)
    hplus  = weyl.A * b.plusA + weyl.Bbar * b.plusB
    hminus = weyl.Abar * b.minusA + weyl.B * b.minusB
    phplus  = exp(-im*ω*t + im*m*φ)
    phminus = exp( im*ω*t - im*m*φ)        # ω̄ = ω
    hplus*phplus + hminus*phminus
end
