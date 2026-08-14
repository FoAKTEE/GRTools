"""
    SelfForce.SecondOrder

Second-order Ricci source over the h1 corpus (legacy `SecondOrderRicci`). Reads the h1 corpus via
`Data`, evaluates the 10 quadratic combinations `r2_1..r2_10` on the `Kernels.r2` GPU host, couples
3 harmonics through `Wigner` — the coupling-evaluating `r2_term` HOST WRAPPER lives HERE (§5.5
Kernels purity: Kernels only takes precomputed coupling arrays per `Kernels.R2_COUPLING_SPEC`) —
supplies the analytic `h1_mass`/`h1_spin` and `d2hdr2`, and assembles `second_order_ricci` /
`ricci_source` then `write_src`. Type-generic `{T}`.

`a_il` (Barack-Lousto tensor-harmonic normalization) is owned by `NumericsCore` and re-exported here
as the SecondOrder front door — it is not re-owned.
"""
module SecondOrder

using NumericsCore: a_il
using Wigner: coupling
using Data: h1_parity_fields

# ── 3-harmonic spin-weighted angular coupling front door (delegates to Wigner.coupling) ───────────
include("SecondOrder/Coupling.jl")

# ── Mode machinery: a_il re-export + 2nd-order Ricci mode selection ───────────────────────────────
include("SecondOrder/Modes.jl")

# ── Analytic first-order fields (mass/spin) + Lorenz field-equation 2nd radial derivative ─────────
include("SecondOrder/AnalyticFields.jl")

# ── Directory-level h1-corpus reader + zero-pad field shim (consumes Data.h1_parity_fields) ───────
include("SecondOrder/ReadCorpus.jl")

# ── End-to-end Ricci-source assembly + writer (consumes Kernels.r2_term GPU host) ─────────────────
include("SecondOrder/Assembly.jl")

# ── Q0 schema only: typed normalization/modes/bandwidth/corpus fingerprints (no compiler) ─────────
include("SecondOrder/PlanSchema.jl")

export a_il, selection_modes, coupling
export h1_mass, h1_spin, d2hdr2
export read_h1, h1_parity_fields, H1RadialJetCorpus, ZeroPadField
export r2_term   # coupling-evaluating host wrapper over Kernels.r2_term (legacy-parity export)
export second_order_ricci, ricci_source, ricci_source_mass, ricci_source_spin, write_src
export R2Normalization, OrderedRicciBilinear, SecondVariationRicci
export ModeID, R2TargetID, EffectiveBandwidth
export H1CorpusConventions, BilinearCorpusSchema, QuadraticPlanSchema
export R2_FINGERPRINT_SCHEMA_VERSION, R2_COUPLING_SPEC_VERSION
export legacy_effective_bandwidth, quadratic_target_order, legacy_corpus_schema
export assert_bilinear_compatible, assert_corpus_matches, quadratic_plan_schema, validate_schema

end # module SecondOrder
