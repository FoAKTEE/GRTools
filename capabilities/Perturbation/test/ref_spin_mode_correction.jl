# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_spin_mode_correction.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from spin_mode_correction.json via gen/json_to_jl.py. BHPToolkit (SpinWeightedSpheroidalHarmonics`) oracle.
const REF_SPIN_MODE_CORRECTION = [
    (l=2, m=2, n=0, k=0, a=0.9, p=12.0, e=0.2, x=0.5, nmax=4, kmax=4, omega=0.045986687799471085, efI=2.0772765758418024e-06, efH=-3.2441239797290215e-09, jfI=9.034251759552442e-05, jfH=-1.4108969943107462e-07, dE1I=-1.8040904769593304e-07, dE1H=3.5449987647647265e-10, dL1I=-4.378696682042671e-06, dL1H=1.000232138653752e-08, stepsr=64, stepsth=64),
    (l=2, m=2, n=1, k=0, a=0.9, p=12.0, e=0.2, x=0.5, nmax=4, kmax=4, omega=0.06333747392664302, efI=1.2700614733724982e-06, efH=-2.9601057353096825e-09, jfI=4.010458247335452e-05, jfH=-9.347091229871448e-08, dE1I=-3.446337505743826e-07, dE1H=7.692675685664882e-10, dL1I=-1.0313887943206348e-05, dL1H=2.2965913885341684e-08, stepsr=96, stepsth=64),
    (l=2, m=1, n=0, k=0, a=0.9, p=12.0, e=0.2, x=0.5, nmax=4, kmax=4, omega=0.022993343899735542, efI=1.5234981858700702e-09, efH=-1.2925210782766077e-12, jfI=6.62582264029719e-08, jfH=-5.6212836371810785e-11, dE1I=6.922614164136149e-10, dE1H=1.6410214386081034e-12, dL1I=3.265009920086308e-08, dL1H=6.921190145930226e-11, stepsr=32, stepsth=32),
]
