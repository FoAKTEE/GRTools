# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_incl_correction.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from incl_correction.json via gen/json_to_jl.py. BHPToolkit (SpinWeightedSpheroidalHarmonics`) oracle.
const REF_INCL_CORRECTION = [
    (l=2, m=2, n=0, k=0, a=0.9, p=12.0, e=0.2, x=0.5, omega=0.045986687799471085, dE1I=2.5284257324000807e-06, dE1H=-6.017846240342998e-09, dL1I=7.625199227798406e-05, dL1H=-2.0907349312035986e-07),
    (l=2, m=2, n=0, k=-10, a=0.9, p=12.0, e=0.2, x=0.5, omega=-0.17502479936699417, dE1I=2.1611657811262876e-22, dE1H=1.4004492784501084e-23, dL1I=-1.6135590182641215e-21, dL1H=-1.3637094191535584e-22),
]
