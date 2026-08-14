# sources.jl — the stress-energy → NP-Ricci source layer (HOLE-6: this plumbing is new, D4 is
# vacuum). Two pieces:
#   (1) the on-shell projection  T_{ab}^{(1)} → (Φ00, Φ22, Φ11, Λ)  via  R_{μν}=8π(T_{μν}−½g_{μν}T)
#       and Φ00=R_ll/2, Φ22=R_nn/2, Φ11=(R_mm̄+R_ln)/4, Λ=(R_mm̄−R_ln)/12  (LY:757);
#   (2) the operator-valued Bianchi sources R_g, R_h (LY:505–515), type-D reduced.
# Plus the thin-shell stress-energy itself, DERIVED from the linearized Israel-junction EOS
# (ANCHOR-1) through the Kinnersley tetrad contractions (so G-E1 is EOS → T_ab → Φ,Λ end-to-end).
#
# NP normalization (mostly-plus, Kinnersley): l·n = −1, m·m̄ = +1  (verified from LY:735–738), so
#   g_ll = g_nn = 0 ,  g_ln = −1 ,  g_mm̄ = +1 ,  T = g^{μν}T_{μν} = −2T_ln + 2T_mm̄ .

"""
    shell_eos(T, M, R) -> (σ, p)

The O(ε) surface energy density and pressure of the static thin shell (ANCHOR-1, LY:718–721),
`σ = M/(4πR²√f(R))`, `p = M²/(8πR³ f^{3/2}(R))` with `f(R)=1−2M/R` (the ε=1 coefficient of the
linearized perfect-fluid junction data). {T}-generic.
"""
function shell_eos(::Type{T}, M, R) where {T}
    Mt, Rt = T(M), T(R)
    fR = 1 - 2Mt / Rt
    π_ = T(pi)
    σ = Mt / (4π_ * Rt^2 * sqrt(fR))
    p = Mt^2 / (8π_ * Rt^3 * fR^(T(3) / 2))
    (σ = σ, p = p)
end

"""
    shell_T_np(T, M, R) -> (Tll, Tnn, Tln, Tmm̄)

The δ-coefficients of the shell stress-energy in the Kinnersley NP basis (ANCHOR-3, LY:749–756),
DERIVED from the EOS `(σ,p)` of [`shell_eos`](@ref) via the tetrad contractions of the perfect-fluid
`T^{μν}=[σu^μu^ν+p(ê_θê_θ+ê_φê_φ)]δ(ℓ)` (LY:723–730), `ℓ=(r−R)/√f(R)` so `δ(ℓ)=√f(R)δ(r−R)`, with
`u·l=−f^{-1/2}`, `u·n=−f^{1/2}/2`, `u·m=0`, and the angular unit-basis contractions `|m·ê_θ|²=|m·ê_φ|²=½`:

    Tll = σ·f(R)^{-1/2},   Tnn = σ·f(R)^{3/2}/4,   Tln = σ·f(R)^{1/2}/2,   Tmm̄ = p·f(R)^{1/2}.

Returns the four coefficients of `δ(r−R)` (values, M=area-mass). Reproduces ANCHOR-3 exactly (G-E1a).
"""
function shell_T_np(::Type{T}, M, R) where {T}
    eos = shell_eos(T, M, R)
    fR = 1 - 2 * T(M) / T(R)
    sqfR = sqrt(fR)
    (Tll  = eos.σ / sqfR,
     Tnn  = eos.σ * fR * sqfR / 4,
     Tln  = eos.σ * sqfR / 2,
     Tmm̄ = eos.p * sqfR)
end

"""
    np_ricci_from_T(Tll, Tnn, Tln, Tmm̄) -> (Φ00, Φ22, Φ11, Λ)

The on-shell NP-Ricci projection (LY:757). With `R_{μν}=8π(T_{μν}−½g_{μν}T)`, the NP normalization
`g_ln=−1, g_mm̄=+1, g_ll=g_nn=0`, and the trace `T=−2Tln+2Tmm̄`:

    R_ll = 8π Tll ,  R_nn = 8π Tnn ,  R_ln = 8π(Tln + ½T) ,  R_mm̄ = 8π(Tmm̄ − ½T) ,

then `Φ00=R_ll/2, Φ22=R_nn/2, Φ11=(R_mm̄+R_ln)/4, Λ=(R_mm̄−R_ln)/12`. Operates on the δ-coefficients
(or on any scalar). Reproduces ANCHOR-4 from ANCHOR-3 (G-E1c; hand-verified in the dossier).
"""
function np_ricci_from_T(Tll, Tnn, Tln, Tmm̄)
    π_ = oftype(float(real(Tll)), pi)
    Tr = -2Tln + 2Tmm̄
    R_ll = 8π_ * Tll
    R_nn = 8π_ * Tnn
    R_ln = 8π_ * (Tln + Tr / 2)
    R_mm̄ = 8π_ * (Tmm̄ - Tr / 2)
    (Φ00 = R_ll / 2, Φ22 = R_nn / 2, Φ11 = (R_mm̄ + R_ln) / 4, Λ = (R_mm̄ - R_ln) / 12)
end

"""
    shell_ricci_profiles(T, M, R) -> (Φ00, Φ22, Φ11, Λ)

The four shell NP-Ricci scalars as pure-δ [`ShellProfile`](@ref)s (each `= coeff · δ(r−R)`), built
from [`shell_T_np`](@ref) through [`np_ricci_from_T`](@ref). These are the ε-order sources that drive
the transport chain.
"""
function shell_ricci_profiles(::Type{T}, M, R) where {T}
    Tn = shell_T_np(T, M, R)
    φ = np_ricci_from_T(Tn.Tll, Tn.Tnn, Tn.Tln, Tn.Tmm̄)
    Rt = T(R)
    mk(c) = ShellProfile{T}(lzero(T), lzero(T), T(c), zero(T), Rt)
    (Φ00 = mk(φ.Φ00), Φ22 = mk(φ.Φ22), Φ11 = mk(φ.Φ11), Λ = mk(φ.Λ))
end

"""
    Rg_shell(T, M, R) -> ShellProfile

The Bianchi source `R_g` (LY:505–510, `eq:Bianchi_Delta_Psi2_source`) reduced to the static
Schwarzschild shell. Type-D (`ν=λ=0`) and Schwarzschild (`π=0`, `Φ21=Φ12=Φ20=0`) collapse it to

    R_g = − D_{[2,2,0,-1]}Φ22 − 2ΔΛ − 2μΦ11 ,

with `D=∂_r` (l^r=1, static), `Δ=−(f/2)∂_r` (n^r=−f/2), and background `ρ̄=−1/r, μ=−f/2r`. Built
here from the [`shell_ricci_profiles`](@ref) sources through the distributional algebra — it
reproduces LY:810–814 exactly (`−M(2R−3M)/(12R³)δ′ + M²/(12R⁴f(R))δ`, G-E4 sub-check).
"""
function Rg_shell(::Type{T}, M, R) where {T}
    src = shell_ricci_profiles(T, M, R)
    Φ22 = src.Φ22; Φ11 = src.Φ11; Λ = src.Λ
    np = schw_np_laurent(T, M)
    f = schw_f(T, M)
    # D_{[2,2,0,-1]}Φ22 = (D + 2ε + 2ε̄ − ρ̄)Φ22 = ∂_rΦ22 − ρ̄·Φ22  (ε=0, ρ̄=ρ=−1/r)
    DΦ22 = pderiv(Φ22) - psmooth(np.ρ, Φ22)
    # 2ΔΛ = 2·(−f/2 ∂_r)Λ = −f·∂_rΛ
    twoΔΛ = pscale(-one(T), psmooth(f, pderiv(Λ)))
    # 2μΦ11
    twoμΦ11 = psmooth(lscale(2, np.μ), Φ11)
    pscale(-one(T), DΦ22) - twoΔΛ - twoμΦ11
end

"""
    Rh_shell(T, M, R) -> ShellProfile

The Bianchi source `R_h` (LY:511–515) for the static Schwarzschild shell. Type-D reduces it to
`R_h = Δ_{[0,2,2,0]}Φ21 − δ̄_{[2,2,0,-1]}Φ22`; with `Φ21=0` and a spherically symmetric (ℓ=0) `Φ22`
the angular operator `δ̄_{[2,2,0,-1]}Φ22 = (δ̄ + 2α + 2β̄)Φ22 = 0` (background `α+β̄=0`, `δ̄` annihilates
the monopole), so `R_h = 0` — the ANCHOR-5 vanishing that forces `Ψ3^{(1)}=0`. Monopole scope
(the angular δ̄ on ℓ≥1 sectors is [FUTURE]).
"""
Rh_shell(::Type{T}, M, R) where {T} = ShellProfile{T}(lzero(T), lzero(T), zero(T), zero(T), T(R))
