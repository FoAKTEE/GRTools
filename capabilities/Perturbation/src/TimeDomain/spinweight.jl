# Spin-weighted mode-coupling coefficients for the generic-s (s ∈ {0, ±1, ±2}) m-mode
# ℓ-coupled 1+1 reduction (DG-P4). Closed forms from Barack & Giudice [arXiv:1702.04204],
# transcribed via the dossier progress/emri-disk/DG_SPIN_FORMULATION.md §2.2 (commit 5f6bd0e):
#
#   ₛYₗₘ cosθ  = c₋^{ℓ+1} ₛY_{ℓ+1,m} + c₀^ℓ ₛYₗₘ + c₊^{ℓ−1} ₛY_{ℓ−1,m}   [arXiv:1702.04204, Eq. (32)]
#   ₛYₗₘ sin²θ = C₋₋^{ℓ+2} ₛY_{ℓ+2,m} + C₋^{ℓ+1} ₛY_{ℓ+1,m} + C₀^ℓ ₛYₗₘ
#                + C₊^{ℓ−1} ₛY_{ℓ−1,m} + C₊₊^{ℓ−2} ₛY_{ℓ−2,m}            [arXiv:1702.04204, Eq. (33)]
#
# Row-ℓ index placement per [arXiv:1702.04204, Eq. (37)] (dossier §3.2/§3.3): in the evolution
# equation FOR mode ℓ, the coupling TO mode ℓ±1 / ℓ±2 carries c±^ℓ / C±±^ℓ; the resulting
# matrices are symmetric because c₊^ℓ = c₋^{ℓ+1} ⇒ C₊^ℓ = C₋^{ℓ+1}, C₊₊^ℓ = C₋₋^{ℓ+2}.
#
# s = 0 REDUCTION (gate A): c₀ ≡ 0 kills the ℓ±1 band and c₋^ℓ → √((ℓ²−m²)/((2ℓ−1)(2ℓ+1))),
# so `spin_cminus(0,ℓ,m) == cminus(ℓ,m)` and `spin_mode_coupling(0,ℓ,m)[[1,3,5]] ==
# mode_coupling(ℓ,m)` BIT-EXACTLY (identical integer arithmetic + operation order; the two
# integer rationals (ℓ²·(ℓ²−m²))/(ℓ²·(2ℓ−1)(2ℓ+1)) and (ℓ²−m²)/((2ℓ−1)(2ℓ+1)) are the same
# exact rational, hence the same correctly-rounded Float64 quotient).
#
# CG-EQUIVALENCE gate (gate B, dossier §2.2 DERIVED/CHECK3):
#   ⟨ₛY_{ℓ′m}|cosθ|ₛYₗₘ⟩ = √((2ℓ+1)/(2ℓ′+1)) · ⟨ℓ,m;1,0|ℓ′,m⟩ · ⟨ℓ,−s;1,0|ℓ′,−s⟩
# baked one-shot from GrTools core/Wigner.clebsch_gordan (scratchpad/capture_spin_cg_fixture.jl)
# into test/fixtures/spin_cg_reference.jl — Wigner is NOT a runtime dep (firewall preserved).

"""
    spin_cminus(s, l, m) -> c₋^ℓ

Spin-weighted ladder coefficient
`c₋^ℓ = √((ℓ²−s²)(ℓ²−m²) / (ℓ²(2ℓ−1)(2ℓ+1)))` [arXiv:1702.04204, Eq. (34)] — the cosθ
coupling of mode ℓ to mode ℓ−1 (row-ℓ convention). Returns 0 for `ℓ ≤ max(|s|,|m|)` (the
formula's numerator vanishes at ℓ = max(|s|,|m|); below that the mode pair is unphysical).
`spin_cminus(0, l, m) ≡ cminus(l, m)` bit-exactly (s=0 gate).
"""
function spin_cminus(s::Integer, l::Integer, m::Integer)
    l <= max(abs(s), abs(m)) && return 0.0
    return sqrt(((l^2 - s^2) * (l^2 - m^2)) / (l^2 * (2l - 1) * (2l + 1)))
end

"""
    spin_czero(s, l, m) -> c₀^ℓ

Diagonal cosθ coefficient `c₀^ℓ = −ms/(ℓ(ℓ+1))` [arXiv:1702.04204, Eq. (34)] (0 for ℓ ≤ 0:
ℓ=0 only occurs at s=m=0 where the numerator vanishes anyway, and ℓ=−1 only appears through
the C₋^{ℓ=0} combination c₀^{ℓ−1}, which multiplies c₋^0 = 0).
"""
spin_czero(s::Integer, l::Integer, m::Integer) = l <= 0 ? 0.0 : -(m * s) / (l * (l + 1))

"""
    spin_cplus(s, l, m) -> c₊^ℓ

Raising cosθ coefficient `c₊^ℓ = c₋^{ℓ+1}` [arXiv:1702.04204, Eq. (34)].
"""
spin_cplus(s::Integer, l::Integer, m::Integer) = spin_cminus(s, l + 1, m)

"""
    spin_cos_coupling(s, l, m) -> (c₋, c₀, c₊)

The cosθ inter-mode coupling row of mode ℓ: `ₛYₗₘcosθ`-projection coefficients onto
ℓ−1, ℓ, ℓ+1 [arXiv:1702.04204, Eqs. (32),(34)].
"""
spin_cos_coupling(s::Integer, l::Integer, m::Integer) =
    (spin_cminus(s, l, m), spin_czero(s, l, m), spin_cplus(s, l, m))

"""
    spin_mode_coupling(s, l, m) -> (C₋₋, C₋, C₀, C₊, C₊₊)

The sin²θ inter-mode coupling row of mode ℓ onto ℓ−2 … ℓ+2 [arXiv:1702.04204, Eq. (35)]:

    C₊₊^ℓ = −c₊^ℓ c₊^{ℓ+1},  C₊^ℓ = −c₊^ℓ(c₀^{ℓ+1}+c₀^ℓ),  C₀^ℓ = 1−(c₋^ℓ)²−(c₊^ℓ)²−(c₀^ℓ)²,
    C₋^ℓ = −c₋^ℓ(c₀^ℓ+c₀^{ℓ−1}),  C₋₋^ℓ = −c₋^{ℓ−1}c₋^ℓ.

For s=0: C₋ = C₊ = 0 (parity decoupling restored) and `(C₋₋, C₀, C₊₊) == mode_coupling(l, m)`
bit-exactly [s=0 limit = arXiv:2307.01349, Eqs. (52)–(53); gate A]. For s·m ≠ 0 the ℓ±1 band
is nonzero — the coupling matrix becomes PENTADIAGONAL and the even/odd-ℓ split is broken
(dossier §3.3).
"""
function spin_mode_coupling(s::Integer, l::Integer, m::Integer)
    cm = spin_cminus(s, l, m)
    cp = spin_cplus(s, l, m)
    c0 = spin_czero(s, l, m)
    c0p = spin_czero(s, l + 1, m)
    c0m = spin_czero(s, l - 1, m)
    cmm1 = spin_cminus(s, l - 1, m)   # c₋^{ℓ−1}
    cpp1 = spin_cplus(s, l + 1, m)    # c₊^{ℓ+1} = c₋^{ℓ+2}
    Cmm = -cmm1 * cm                  # C₋₋^ℓ  (== repo C₋₋ order: -clm1*cl)
    Cm = -cm * (c0 + c0m)             # C₋^ℓ
    C0 = 1 - cm^2 - cp^2 - c0^2       # C₀^ℓ   (assoc. ((1−cm²)−cp²)−c0², bit-= repo at s=0)
    Cp = -cp * (c0p + c0)             # C₊^ℓ
    Cpp = -cp * cpp1                  # C₊₊^ℓ  (== repo C₊₊ order: -clp1*clp2)
    return (Cmm, Cm, C0, Cp, Cpp)
end

"""
    spin_mode_list(s, m, ℓmax) -> Vector{Int}

Joint (no parity split) coupled-mode list `ℓ = max(|s|,|m|), …, ℓmax` for the generic-s
sector. For s·a ≠ 0 the ℓ±1 couplings (β± ∝ a·s, dossier §3.1 item 5) break the even/odd-ℓ
decoupling for EVERY m, so the joint consecutive list is required (contrast `mode_list`,
whose step-2 parity split is exact only at s=0).
"""
spin_mode_list(s::Integer, m::Integer, ℓmax::Integer) = collect(max(abs(s), abs(m)):ℓmax)

"cosθ coupling matrix 𝒦 over the sorted mode list `ℓs` (row ℓᵢ → col ℓⱼ, [1702.04204 Eq. (32)]
placement; symmetric since c₊^ℓ = c₋^{ℓ+1}). Entries beyond |Δℓ|>1 vanish identically."
function _spin_cos_matrix(s, ℓs, m)
    n = length(ℓs)
    Km = zeros(n, n)
    for i in 1:n, j in 1:n
        dl = ℓs[j] - ℓs[i]
        if dl == 0
            Km[i, j] = spin_czero(s, ℓs[i], m)
        elseif dl == 1
            Km[i, j] = spin_cplus(s, ℓs[i], m)
        elseif dl == -1
            Km[i, j] = spin_cminus(s, ℓs[i], m)
        end
    end
    return Km
end

"sin²θ coupling matrix 𝒮 over the sorted mode list `ℓs` (row ℓᵢ → col ℓⱼ, [1702.04204
Eq. (33)/(35)] placement; symmetric pentadiagonal for consecutive ℓs)."
function _spin_sin2_matrix(s, ℓs, m)
    n = length(ℓs)
    Sm = zeros(n, n)
    for i in 1:n
        Cmm, Cm, C0, Cp, Cpp = spin_mode_coupling(s, ℓs[i], m)
        for j in 1:n
            dl = ℓs[j] - ℓs[i]
            if dl == 0
                Sm[i, j] = C0
            elseif dl == -1
                Sm[i, j] = Cm
            elseif dl == 1
                Sm[i, j] = Cp
            elseif dl == -2
                Sm[i, j] = Cmm
            elseif dl == 2
                Sm[i, j] = Cpp
            end
        end
    end
    return Sm
end
