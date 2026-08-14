"""
    SelfForce.Puncture

Puncture / effective-source machinery (legacy `Punctures`): static h0 AGM elliptic integrals,
auto-generated coefficients, terminating-`2F1` α-integration via `BigSpecial`, Wigner-D β=π/2 mode
mixing, digamma regularisation sums, the singular/retarded field subtraction, and the worldline
self-field. Type-generic `{T}`.
"""
module Puncture

# ── Static h0 first-order field: l-mode static coefficients over the auto-gen table + consumed
#    NumericsCore elliptic K/E + BigSpecial Gamma + Spacetime worldline connection ─────────────────
include("Puncture/Static.jl")

export h0_static, ellipK, ellipE, fact_big
export christoffel, christoffel_dr, christoffel_dtheta

# ── Terminating-2F1 polar-angle α-integration: per-(p,l,m) angular integrals (legacy
#    Punctures/alpha_integration.jl), with _poch→BigSpecial.pochhammer + _f21_term→f21_terminating ──
include("Puncture/AlphaIntegration.jl")

export exactint, modeint, alphaint
export alphaintv1, alphaintv2, alphaintv3, alphaintv4, alphaintt1, alphaintt2

# ── Regularized l-sum (tail fit) → value ± uncertainty (legacy Punctures/regsum.jl), with the
#    analytic ∞-tail digamma collapsed onto BigSpecial.digamma_big ─────────────────────────────────
include("Puncture/RegSum.jl")

export RegResult, reg_sum, noreg_sum

# ── Lorenz wave-operator second radial derivative ddhret[i,l,m] from the Barack-Lousto field
#    equations + even-static m=0 routing predicate (legacy Punctures/subtraction_fieldeq.jl) ─────────
include("Puncture/FieldEqDeriv.jl")

export ddhret, is_evenstatic_ddh, EVENSTATIC_M0_COMPONENTS

# ── Martel-Poisson tensor-harmonic angular basis (11 fns) + coordinate reconstruction of hR (legacy
#    Punctures/tensor_basis.jl); seed CONSUMED from Harmonics, θ-derivative via local ForwardDiff ────
include("Puncture/TensorBasis.jl")

export tensor_harmonic, coord_component, coord_deriv

# ── First-order puncture-mode amplitudes hP[i][l,m,r0,M,Δr] + one-sided radial derivatives dhP:
#    prefactor × Gauss-Kronrod β-integral over the vendored machine-generated rational integrand,
#    contracted against the AlphaIntegration α-moments; dhP via exact ForwardDiff over the integrand ──
include("Puncture/PunctureModes.jl")

export hP, dhP, rat

# ── Wigner-D β=π/2 mode mixing h[i][l,m] = Σ_{mp} D^l_{m,mp}(π,π/2,π/2) hP[i][l,mp]: the fixed-angle
#    rotation CONSUMED from Wigner.wignerD_fixed (AHA), per-component _MIX_PATTERN selecting the mp set ─
include("Puncture/ModeMixing.jl")

export h_mode

# ── Singular-field l-modes hS[i,l,m] on the worldtube Δr-grid (= mode-mixed h_mode) + the two
#    identically-zero source modes short-circuited + worldline left/right split helpers ─────────────
include("Puncture/SingularField.jl")

export hS_mode, worldline_index, split_worldline, SINGULAR_ZERO_MODES

# ── Regular-field subtraction hR = hret − hS (+ ddhR): retarded-mode HDF5 read CONSUMED from
#    Data.read_h1, physical-component→block map CONSUMED from Data.hret_block ─────────────────────────
include("Puncture/Subtraction.jl")

export hret_block, read_hret_mode, compute_hR, compute_ddhR

# ── Worldline regular self-field: assemble hR1/dhR1/d2hR1 (parity-routed component map, reusing
#    coord_deriv + reg_sum/noreg_sum) + covariant 1st/2nd derivs on the a=0 background (Christoffel
#    CONSUMED from Spacetime). CLOSES the Puncture namespace ────────────────────────────────────────
include("Puncture/WorldlineSelfField.jl")

export assemble_hR1Data, assemble_dhR1Data, assemble_d2hR1Data
export covariant_first, covariant_second

end # module Puncture
