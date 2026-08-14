# ASSEMBLE — parity-aware BLS field → metric perturbation h_μν(r,θ): the DDKW dual-leg assembly
# (the iter-245 normalization-fix form). The BLS fields h̄¹..h̄¹⁰ (h1Lorenz corpus / KerrLorenzCirc
# reconstruction) are in the DDKW convention (arXiv:2306.16459) and ALREADY embed the ladder factors
# λ̂=√(l(l+1)) and Λ̂=√((l−1)l(l+1)(l+2)); the assembly (1) inverts the BLS fields back to the 10
# projected tetrad radial modes (the exact algebraic inverse of KerrLorenzCirc.to_BLS_kerr), (2)
# attaches the spin-weighted SPHERICAL harmonic ₛY_lm(θ), (3) contracts with the dual legs of the
# Kinnersley tetrad, h_μν = Σ_ab h_ab w^a_μ w^b_ν. Reuses Spacetime.kerr_metric (collapses the legacy
# DiscModes.kerr_metric fork — AHA), NumericsCore.lambdahat/Lambdahat (the ladder seed — AHA),
# Spacetime.schwarzschild_lapse (the Schwarzschild lapse seed), Harmonics SWSH spherical-Y.
# Type-generic in the radial fields; the angular ₛY_lm are Float64 (the SWSH angular cap) and promoted.
# From DiscEMRIForce/assemble.jl + kerr_source.jl::kerr_parity.

import Spacetime: kerr_metric, schwarzschild_lapse
import Harmonics: SpinWeightedSphericalHarmonicY
import NumericsCore: lambdahat, Lambdahat

"Parity of a Kerr (l,m) mode: `:even` if (l+m) even (7 BLS fields), `:odd` otherwise (3 fields)."
kerr_parity(l::Integer, m::Integer) = iseven(l + m) ? :even : :odd

"Schwarzschild lapse seed f(r)=1−2/r (reuses Spacetime.schwarzschild_lapse)."
_f(r) = schwarzschild_lapse(r)

"Kerr radial factor f_K(r,a)=(r²−2r+a²)/r² = Δ/r² (reduces to the Schwarzschild lapse at a=0)."
_fK(r::T, a::T) where {T<:Real} = (r^2 - 2r + a^2)/r^2

"Spin-1 angular ladder λ̂=√(l(l+1)) (reuses NumericsCore.lambdahat)."
_λ̂(l) = lambdahat(l)
"Spin-2 angular ladder Λ̂=√((l−1)l(l+1)(l+2)) (reuses NumericsCore.Lambdahat)."
_Λ̂(l) = Lambdahat(l)

"Spin-weighted spherical harmonic ₛY_lm(θ,0) (Float64 angular cap); 0 when l<|s|."
_sY(s::Integer, l::Integer, m::Integer, θ) =
    l < abs(s) ? 0.0 + 0.0im : Complex{Float64}(SpinWeightedSphericalHarmonicY(s, l, m, Float64(θ), 0.0))

"""
    tetmodes_from_bls(mode, k, l; a=0) -> NamedTuple

The 10 projected tetrad radial modes at radial index `k`, recovered from the stored BLS fields by the
EXACT algebraic inverse of `KerrLorenzCirc.to_BLS_kerr` (a=0 limit = DDKW eq:hBL). Parity-aware (even:
h̄¹..h̄⁷; odd: h̄⁸,h̄⁹,h̄¹⁰). Guards `l<2` (no spin-2 content; Λ̂=0).
"""
function tetmodes_from_bls(mode::H1Mode{T}, k::Integer, l::Integer; a=zero(T)) where {T}
    z = zero(Complex{T})
    if mode.parity === :even
        h1 = mode.val_right[1][k]; h3 = mode.val_right[2][k]; h5 = mode.val_right[3][k]
        h6 = mode.val_right[4][k]; h7 = mode.val_right[5][k]; h2 = mode.val_right[6][k]
        h4 = mode.val_right[7][k]; h8 = z; h9 = z; h10 = z
    else
        h9 = mode.val_right[1][k]; h10 = mode.val_right[2][k]; h8 = mode.val_right[3][k]
        h1 = z; h2 = z; h3 = z; h4 = z; h5 = z; h6 = z; h7 = z
    end
    r = mode.r_right[k]; aT = T(a)
    Δ = r^2 - 2r + aT^2; Σc = r^2 + aT^2; f = Δ/r^2
    lplp = (h1 + h2)/(r*f^2); lmlm = (h1 - h2)/(r*f^2)
    if l ≥ 2
        mpmp = (r/T(_Λ̂(l)))*(h7 - im*h10); mmmm = (r/T(_Λ̂(l)))*(h7 + im*h10)
    else
        mpmp = z; mmmm = z
    end
    if l >= 1
        lf = T(_λ̂(l))*f
        Aq = (-(h4 + h5) + im*(h8 + h9))/(2lf); Bq = ((h4 + h5) + im*(h8 + h9))/(2lf)
        Cq = ((h4 - h5) - im*(h8 - h9))/(2lf); Dq = (-(h4 - h5) - im*(h8 - h9))/(2lf)
        ρlpmp = r*Aq; ρclpmm = r*Bq; ρclmmp = r*Cq; ρlmmm = r*Dq
    else
        # The monopole has no spin-1 sector. Its compact file omits h̄⁴/h̄⁵, so the formal
        # inverse is 0/λ̂(0)=0/0; use the exact low-l contraction instead.
        ρlpmp = z; ρclpmm = z; ρclmmp = z; ρlmmm = z
    end
    hlplm = -h6/(r*f); ΣΔlplm = Σc*Δ*hlplm
    htr = (r*h3 + Δ*hlplm)/Σc
    (lplp=lplp, lmlm=lmlm, mpmp=mpmp, mmmm=mmmm, ρlpmp=ρlpmp, ρclpmm=ρclpmm,
     ρclmmp=ρclmmp, ρlmmm=ρlmmm, ΣΔlplm=ΣΔlplm, htr=htr, r=r, Δ=Δ)
end

"""
    hmunu(mode, k, l, m, θ; a=0) -> 4×4 Complex{T}

`h_μν^{lm}(r,θ)` via the DDKW dual-leg assembly (the iter-245 normalization fix): invert the BLS
fields to the 10 tetrad modes (`tetmodes_from_bls`), attach ₛY_lm, divide the mixed l±m± entries by
ρ=r+iacosθ, recover h_{m+m−} from the trace identity, and contract with the dual Kinnersley legs
`h_μν = Σ_ab h_ab w^a_μ w^b_ν`, using `Spacetime.kerr_metric` for the lowered legs. Parity-aware.
"""
function hmunu(mode::H1Mode{T}, k::Integer, l::Integer, m::Integer, θ; a=zero(T)) where {T}
    p = tetmodes_from_bls(mode, k, l; a=a)
    r = p.r; aT = T(a); Δ = p.Δ
    ct = cos(θ); st = sin(θ)
    Σθ = r^2 + aT^2*ct^2
    ρ = r + im*aT*ct; ρc = r - im*aT*ct
    CT = Complex{T}
    Y0 = CT(_sY(0, l, m, θ)); Y2 = CT(_sY(2, l, m, θ)); Ym2 = CT(_sY(-2, l, m, θ))
    Y1 = CT(_sY(1, l, m, θ)); Ym1 = CT(_sY(-1, l, m, θ))
    H_lplp = p.lplp*Y0;  H_lmlm = p.lmlm*Y0
    H_mpmp = p.mpmp*Y2;  H_mmmm = p.mmmm*Ym2
    H_lpmp = p.ρlpmp*Y1/ρ;   H_lpmm = p.ρclpmm*Ym1/ρc
    H_lmmp = p.ρclmmp*Y1/ρc; H_lmmm = p.ρlmmm*Ym1/ρ
    H_lplm = p.ΣΔlplm*Y0/(Σθ*Δ); H_tr = p.htr*Y0
    H_mpmm = Σθ*H_tr - Δ*H_lplm
    # unnormalized Kinnersley legs (BL components) and their dual lowered forms
    lp  = CT[(r^2 + aT^2)/Δ, 1, 0, aT/Δ]
    lm_ = CT[-(r^2 + aT^2)/Δ, 1, 0, -aT/Δ]
    mp  = CT[im*aT*st, 0, 1, im/st]
    mm_ = CT[-im*aT*st, 0, 1, -im/st]
    g = kerr_metric(aT, r, θ)
    lo(v) = CT[sum(g[μ,ν]*v[ν] for ν in 1:4) for μ in 1:4]
    w = ((Δ/(2Σθ)).*lo(lm_), (Δ/(2Σθ)).*lo(lp), (one(T)/(2Σθ)).*lo(mm_), (one(T)/(2Σθ)).*lo(mp))
    hab = zeros(CT, 4, 4)
    hab[1,1] = H_lplp; hab[2,2] = H_lmlm; hab[3,3] = H_mpmp; hab[4,4] = H_mmmm
    hab[1,2] = hab[2,1] = H_lplm; hab[3,4] = hab[4,3] = H_mpmm
    hab[1,3] = hab[3,1] = H_lpmp; hab[1,4] = hab[4,1] = H_lpmm
    hab[2,3] = hab[3,2] = H_lmmp; hab[2,4] = hab[4,2] = H_lmmm
    H = zeros(CT, 4, 4)
    @inbounds for μ in 1:4, ν in 1:4
        H[μ,ν] = sum(hab[Ai,Bi]*w[Ai][μ]*w[Bi][ν] for Ai in 1:4, Bi in 1:4)
    end
    H
end

"""
    hmunu_eq(mode, k, l, m) -> 4×4 Complex{T}

Equatorial (θ=π/2) `h_μν` via the DDKW assembly.
"""
hmunu_eq(mode::H1Mode{T}, k::Integer, l::Integer, m::Integer) where {T} = hmunu(mode, k, l, m, T(π)/2)
