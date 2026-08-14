# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_incl_spin_mode_k10.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from incl_spin_mode_k10.json via gen/json_to_jl.py. BHPToolkit (SpinWeightedSpheroidalHarmonics`) oracle.
const REF_INCL_SPIN_MODE_K10 = [
    (l=2, m=2, n=0, k=-10, spar=0.0, a=0.9, p=12.0, e=0.2, x=0.5, omega=-0.17502479936699417, CP=(-8.495358344630758e-12 + -2.261024129152578e-12im), efI=2.0075984963783626e-22, efH=5.548550117726206e-24, jfI=-2.294073186930276e-21, jfH=-6.340301646159225e-23),
    (l=2, m=2, n=0, k=-10, spar=0.1, a=0.9, p=12.0, e=0.2, x=0.5, omega=-0.18180901771237007, CP=(-9.359774421646009e-12 + -2.3369079746236284e-12im), efI=2.240540258620118e-22, efH=7.17595140926819e-24, jfI=-2.4647185126589835e-21, jfH=-7.893944425375933e-23),
]
