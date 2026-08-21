# diagnostics.jl — route R3 stage 3: the angular-momentum diagnostics of arXiv:2608.17003
# (Sec. II.B and Sec. III) in the master variables.
#
#   • the Lagrangian displacement ξ^α reconstructed from (16) — ξ^r = iδû^r/(u^tω_d) falls out of the
#     r-component, u_αζ^α = 0 fixes ζ^t in terms of ζ^φ and the φ-component fixes ζ^φ, while the
#     t-component is left over as a CONSISTENCY CHECK; an arbitrary trivial displacement F u^α is
#     carried through so that gauge independence is certified, not assumed;
#   • U^{αβγδ} (32), V^{αβγδ} (33), W^{αβγδ} (42), Π^{αβ} (36);
#   • δ_ξT^{αβ} = δT^{αβ} − δ_{ξ=0}T^{αβ} (40),(41) straight from the λ-jets;
#   • the per-mode torque-density integrand ½δ_ξT̂^{αβ}(imĥ_{αβ})* of (43b)/(72b) and the per-mode AMF
#     integrand u_φδp̂δû^{r*} + (ε+p)δû^rδû_φ* of (72a), with the Fourier bilinear identity (73)
#     ⟨XY⟩ = 2πX̂₀Ŷ₀ + 4πΣ_{m≥1}Re[X̂_mŶ_m*].
#
# The gate is the Noether identity (69) — the angular-momentum current of Lagrangian perturbation
# theory equals the ADVECTED angular-momentum flux — on-shell and φ-averaged (ERRATUM E-M3: it is NOT
# a pointwise off-shell identity), for EVERY gauge function.
#
# ORDER OF OPERATIONS (important): the algebraic eliminations and the change to master variables are
# applied BEFORE any d/dr, because d/dr of the SYMBOL δû^φ is not the derivative of the SOLVED δû^φ.
# After that the only derivatives produced are y0′, y1′ (replaced by C·y + S) and ĥ′.
#
# Conjugation: ω = mΩ• (46) is real, hence so is ω_d = ω − mΩ; the amplitudes y0, y1, ĥ_{μν} and the
# gauge function are complex and each gets an independent barred symbol, with `cbar` performing the
# formal conjugation (i → −i plus the symbol swap).  Treating the bars as INDEPENDENT makes the
# certificates stronger: an identity that survives independent bars holds a fortiori on the conjugate
# locus.

# --- conjugate alphabet ---------------------------------------------------------------------------
@variables y0c y1c Fg Fg_r Fgc Fgc_r
@variables httc htrc htphic hrrc hrphic hphiphic hzzc
@variables httc_r htrc_r htphic_r hrrc_r hrphic_r hphiphic_r hzzc_r

_rule!(Fg, Fg_r)
_rule!(Fgc, Fgc_r)
for x in (Fg_r, Fgc_r, y0c, y1c, httc, htrc, htphic, hrrc, hrphic, hphiphic, hzzc,
          httc_r, htrc_r, htphic_r, hrrc_r, hrphic_r, hphiphic_r, hzzc_r)
    _rule!(x, 0)          # barred symbols are never differentiated in these expressions
end
append!(ALL_SYMBOLS, [y0c, y1c, Fg, Fg_r, Fgc, Fgc_r,
                      httc, htrc, htphic, hrrc, hrphic, hphiphic, hzzc,
                      httc_r, htrc_r, htphic_r, hrrc_r, hrphic_r, hphiphic_r, hzzc_r])

"Symbol ↔ barred-symbol pairing used by [`cbar`](@ref)."
const CONJ_PAIRS = [(y0, y0c), (y1, y1c), (Fg, Fgc), (Fg_r, Fgc_r),
                    (htt, httc), (htr, htrc), (htphi, htphic), (hrr, hrrc),
                    (hrphi, hrphic), (hphiphi, hphiphic), (hzz, hzzc),
                    (htt_r, httc_r), (htr_r, htrc_r), (htphi_r, htphic_r), (hrr_r, hrrc_r),
                    (hrphi_r, hrphic_r), (hphiphi_r, hphiphic_r), (hzz_r, hzzc_r)]

const _CONJ_SUB = Dict{Num,Num}(vcat([[a => b, b => a] for (a, b) in CONJ_PAIRS]...))

"""
    cbar(e)

Formal complex conjugation of a Fourier amplitude: i → −i together with the amplitude ↔ barred-symbol
swap.  Background quantities, ω and ω_d are real and are left alone.
"""
cbar(e::Num) = Symbolics.substitute(e, _CONJ_SUB)
cbar(e::Complex{Num}) = cbar(real(e)) - im * cbar(imag(e))

"Real part of a Fourier bilinear, ½(z + z̄)."
re_part(z) = (z + cbar(z)) / 2
"Imaginary part of a Fourier bilinear, (z − z̄)/(2i)."
im_part(z) = (z - cbar(z)) / (2 * im)

# --- reduction to master variables ------------------------------------------------------------------

"Algebraic eliminations + change to the master variables (49).  No d/dr may precede this."
to_master(e) = lin_subs(eliminate(e), MASTER_SUB)

"""
    onshell_yp(e)

Impose the derived master system: replace y0′, y1′ by C·y + S.  Only d/dr of a master-variable
expression produces them.
"""
onshell_yp(e) = lin_subs(e, Dict(y0_r => YPRIME[1], y1_r => YPRIME[2]))

"δû^α in the master variables: (δû^t, i y1, δû^φ, 0)."
const DUUP_M = [to_master(DUT_SOL), to_master(dur), to_master(DUPHI_SOL), Num(0)]
"δp̂ = (ε+p)y0."
const DP_M = (eps + p) * y0
"δε̂ in the master variables."
const DEPS_M = to_master(DEPS_SOL)
"The COVARIANT δû_φ = g_{φμ}δû^μ + ĥ_{φμ}u^μ in the master variables (ERRATUM E-M1)."
const DUPHI_LOW_M = to_master(DUPHI_LOW)

# --- Lagrangian displacement from (16) ---------------------------------------------------------------

"q^α{}_β = δ^α{}_β + u^αu_β (2)."
const QMIX = Num[(α == β ? 1 : 0) + UUP[α] * ULOW[β] for α in 1:4, β in 1:4]

"""
    lie_u_vec(z)

𝓛_u z^β for a first-order vector amplitude with the ansatz (45a):
u^μ∂_μz^β − z^μ∂_μu^β = −i u^tω_d z^β − z^r∂_ru^β  (u^r = 0; the background is stationary,
axisymmetric and z-independent).
"""
lie_u_vec(z) = [-im * UT * (omega - m * Ω_expr) * z[β] - z[IR] * ddr(UUP[β]) for β in 1:4]

"Right-hand side of (16), δû^α = q^α{}_β𝓛_uζ^β + ½u^αu^βu^γĥ_{βγ}, for a trial ζ."
function rhs16(z)
    lie = lie_u_vec(z)
    [sum(QMIX[α, β] * lie[β] for β in 1:4) +
     UUP[α] * sum(UUP[β] * UUP[γ] * HMODE[β, γ] for β in 1:4, γ in 1:4) / 2 for α in 1:4]
end

"ξ^r = iδû^r/(u^tω_d) in the master variables (= −y1/(u^tω_d))."
const XI_R_M = im * DUUP_M[IR] / (UT * (omega - m * Ω_expr))

"""
    ZETA_UP

ζ^α — the u-orthogonal part of the Lagrangian displacement, solved from (16) in the master variables.
"""
const ZETA_UP = let
    @variables zphi
    ζ = [-ULOW[IPHI] * zphi / ULOW[IT], XI_R_M, Num(zphi), Num(0)]     # u_αζ^α = 0, ζ^r from (16)|_r
    sol = solve_linear(rhs16(ζ)[IPHI] - DUUP_M[IPHI], zphi)
    [lin_subs(ζ[α], Dict(zphi => sol)) for α in 1:4]
end

"""
    XI_UP

ξ^α = ζ^α + F u^α — the displacement with an ARBITRARY trivial (u-parallel) part `Fg`.  Physical
results must not depend on it; the (69) gate certifies exactly that.
"""
const XI_UP = [ZETA_UP[α] + Fg * UUP[α] for α in 1:4]
"ξ_α = g_{αμ}ξ^μ."
const XI_LOW = [sum(METRIC[α, μ] * XI_UP[μ] for μ in 1:4) for α in 1:4]

"∂_α of a FIRST-ORDER amplitude with the ansatz (45a): ∂_t → −iω, ∂_φ → +im, ∂_r → d/dr, ∂_z → 0."
∂amp(x, α::Int) = α == IT ? -im * omega * x :
                  α == IPHI ? im * m * x :
                  α == IR ? ddr(x) : Num(0)

"∇_γξ_δ = ∂_γξ_δ − Γ^σ{}_{γδ}ξ_σ (background connection: ξ is already first order), on-shell."
const NABLA_XI = [onshell_yp(∂amp(XI_LOW[δ], γ) - sum(GAMMA[σ, γ, δ] * XI_LOW[σ] for σ in 1:4))
                  for γ in 1:4, δ in 1:4]

# --- the Lagrangian tensors (32), (33), (42) ----------------------------------------------------------

"U^{αβγδ} (32): (ε+p)u^αu^γq^{βδ} + p(g^{αβ}g^{γδ} − g^{αδ}g^{βγ}) − Γ₁p q^{αβ}q^{γδ}."
const UTENS = Num[(eps + p) * UUP[α] * UUP[γ] * QUP[β, δ] +
                  p * (METRIC_INV[α, β] * METRIC_INV[γ, δ] - METRIC_INV[α, δ] * METRIC_INV[β, γ]) -
                  Gamma1 * p * QUP[α, β] * QUP[γ, δ]
                  for α in 1:4, β in 1:4, γ in 1:4, δ in 1:4]

"V^{αβγδ} (33): 2V = (ε+p)(u^αu^γq^{βδ} + u^βu^γq^{αδ} − u^αu^βq^{γδ}) − Γ₁p q^{αβ}q^{γδ}."
const VTENS = Num[((eps + p) * (UUP[α] * UUP[γ] * QUP[β, δ] + UUP[β] * UUP[γ] * QUP[α, δ] -
                                UUP[α] * UUP[β] * QUP[γ, δ]) - Gamma1 * p * QUP[α, β] * QUP[γ, δ]) / 2
                  for α in 1:4, β in 1:4, γ in 1:4, δ in 1:4]

"W^{αβγδ} (42)."
const WTENS = Num[((eps + p) * UUP[α] * UUP[β] * UUP[γ] * UUP[δ] +
                   p * (METRIC_INV[α, β] * METRIC_INV[γ, δ] - METRIC_INV[α, γ] * METRIC_INV[β, δ] -
                        METRIC_INV[α, δ] * METRIC_INV[β, γ]) -
                   Gamma1 * p * QUP[α, β] * QUP[γ, δ]) / 2
                  for α in 1:4, β in 1:4, γ in 1:4, δ in 1:4]

"T^{αβ} = (ε+p)u^αu^β + p g^{αβ} (1),(2) — background."
const TUP_BG = Num[(eps + p) * UUP[α] * UUP[β] + p * METRIC_INV[α, β] for α in 1:4, β in 1:4]

"∇_γT^{αβ} of the background."
const NABLA_T = [[∂bg(TUP_BG[α, β], γ) +
                  sum(GAMMA[α, γ, σ] * TUP_BG[σ, β] + GAMMA[β, γ, σ] * TUP_BG[α, σ] for σ in 1:4)
                  for α in 1:4, β in 1:4] for γ in 1:4]

"T^{αβ} as a λ-jet; its O(λ) part is δT^{αβ}."
const TUP_JET = [(EJET + PJET) * UJET[α] * UJET[β] + PJET * GINVJET[α, β] for α in 1:4, β in 1:4]

"δ_{ξ=0}T^{αβ} = [W^{αβγδ} − ½T^{αβ}g^{γδ}]ĥ_{γδ} (41)."
const DT_XI0 = [sum((WTENS[α, β, γ, δ] - TUP_BG[α, β] * METRIC_INV[γ, δ] / 2) * HMODE[γ, δ]
                    for γ in 1:4, δ in 1:4) for α in 1:4, β in 1:4]

"δ_ξT^{αβ} = δT^{αβ} − δ_{ξ=0}T^{αβ} (40), in the master variables."
const DXI_T = [to_master(pert(TUP_JET[α, β]) - DT_XI0[α, β]) for α in 1:4, β in 1:4]

"Π^{αβ} = U^{αβγδ}∇_γξ_δ + V^{γδαβ}ĥ_{γδ} (36)."
const PI_UP = [sum(UTENS[α, β, γ, δ] * NABLA_XI[γ, δ] for γ in 1:4, δ in 1:4) +
               sum(VTENS[γ, δ, α, β] * HMODE[γ, δ] for γ in 1:4, δ in 1:4) for α in 1:4, β in 1:4]

"Δg_{αβ} = ĥ_{αβ} + 2∇_{(α}ξ_{β)} (18)."
const DELTA_G = [HMODE[α, β] + NABLA_XI[α, β] + NABLA_XI[β, α] for α in 1:4, β in 1:4]
"q^{αβ}Δg_{αβ} — the trace that drives (19) and (21)."
const TRACE_DG = sum(QUP[α, β] * DELTA_G[α, β] for α in 1:4, β in 1:4)

# --- per-mode diagnostics integrands ------------------------------------------------------------------

"""
    AMF_MODE

Per-mode AMF integrand of (72a): u_φδp̂δû^{r*} + (ε+p)δû^rδû_φ*, in the master variables (δp̂ = (ε+p)y0,
δû^r = iy1, δû_φ the COVARIANT component — ERRATUM E-M1).  The azimuthally averaged flux is
F_J = √−g·4πΣ_{m≥1}Re[AMF_MODE] by the bilinear identity (73), with √−g = e^{ν+ψ+μ}.
"""
const AMF_MODE = ULOW[IPHI] * DP_M * cbar(DUUP_M[IR]) +
                 (eps + p) * DUUP_M[IR] * cbar(DUPHI_LOW_M)

"""
    TORQUE_MODE

Per-mode torque-density integrand of (43b)/(72b): ½Σ_{αβ}δ_ξT̂^{αβ}(imĥ_{αβ})*.  The azimuthally
averaged torque density is ⟨𝒯⟩ = 4πΣ_{m≥1}Re[TORQUE_MODE].
"""
const TORQUE_MODE = sum(DXI_T[α, β] * cbar(im * m * HMODE[α, β]) for α in 1:4, β in 1:4) / 2

"""
    JR_NOETHER_MODE

Per-mode Noether angular-momentum current j^r_J = −Π^{rβ}𝓛_φξ_β (35b; φ^r = 0 so the 𝓛_eff term drops):
the bilinear whose real part times 4π is ⟨j^r_J⟩.
"""
const JR_NOETHER_MODE = -sum(PI_UP[IR, β] * cbar(im * m * XI_LOW[β]) for β in 1:4)

"""
    stage3_certificates(; npts=50, seed=2608, verbose=true)

Stage-3 gates: consistency of the reconstructed ξ against the unused t-component of (16); the tensor
identity (39) that turns (38) into (43); the Lagrangian relations (19),(21) as an independent
cross-validation of the energy-equation row; and the Noether identity (69) — on-shell, φ-averaged,
for an arbitrary trivial displacement.
"""
function stage3_certificates(; npts::Int = 50, seed::Int = 2608, verbose::Bool = true)
    C = Certificate[]

    res16 = [rhs16(ZETA_UP)[α] - DUUP_M[α] for α in 1:4]
    push!(C, certify_exact("(16) reconstruction of ξ — all four components", "(16)", res16;
                           npts, seed, onshell = true, verbose))
    push!(C, certify_bool("ξ^r = iδû^r/(u^tω_d)", "(16)",
                          _is_zero_expr(ZETA_UP[IR] - XI_R_M),
                          "the r-component of (16) inverts to the recipe used by the closure"; verbose))

    res39 = [DXI_T[α, β] / 2 - sum(VTENS[α, β, γ, δ] * NABLA_XI[γ, δ] for γ in 1:4, δ in 1:4) +
             sum(XI_UP[γ] * NABLA_T[γ][α, β] for γ in 1:4) / 2 for α in 1:4, β in 1:4]
    push!(C, certify_exact("(39)  ½δ_ξT^{αβ} = V^{αβγδ}∇_γξ_δ − ½ξ^γ∇_γT^{αβ}", "(39)-(42)",
                           vec(res39); npts, seed, onshell = true, verbose))

    push!(C, certify_exact("(19)  Δε/(ε+p) = −½q^{αβ}Δg_{αβ}", "(19)",
                           (DEPS_M + XI_UP[IR] * eps_r) / (eps + p) + TRACE_DG / 2;
                           npts, seed, onshell = true, verbose))
    push!(C, certify_exact("(21)  Δp = −½Γ₁p q^{αβ}Δg_{αβ}", "(21)",
                           (DP_M + XI_UP[IR] * p_r) + Gamma1 * p * TRACE_DG / 2;
                           npts, seed, onshell = true, verbose))

    res69 = re_part(JR_NOETHER_MODE) - re_part(AMF_MODE)
    push!(C, certify_exact("(69)  ⟨j^r_J⟩ = u_φ⟨δpδu^r⟩ + (ε+p)⟨δu^rδu_φ⟩", "(69),(72a)",
                           res69; npts, seed, onshell = true, verbose))
    push!(C, certify_exact("(69) is gauge invariant  (∂/∂F of the residual)", "(69)",
                           [dsym(res69, Fg), dsym(res69, Fgc)]; npts, seed, onshell = true, verbose))
    C
end
