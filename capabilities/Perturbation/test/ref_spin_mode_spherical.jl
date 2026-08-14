# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_spin_mode_spherical.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from spin_mode_spherical.json via gen/json_to_jl.py. BHPToolkit (SpinWeightedSpheroidalHarmonics`) oracle.
const REF_SPIN_MODE_SPHERICAL = [
    (l=2, m=2, k=0, a=0.9, p=12.0, e=0.0, x=0.5, kmax=4, spar=0.1, omega=0.04817787653632065, CP=(-0.00030909567201558224 + 7.32359832587575e-05im), CM=(-0.0004084518509192486 + -0.00010919969302333102im), efI=3.459407132534999e-06, efH=-5.902207467853887e-09, jfI=0.00014360978030764798, jfH=-2.450173354321373e-07, stepsth=64),
    (l=2, m=2, k=2, a=0.9, p=12.0, e=0.0, x=0.5, kmax=4, spar=0.1, omega=0.09456847110486204, CP=(5.602607734640552e-06 + -1.8233299615828524e-06im), CM=(2.7899129934008454e-06 + 1.1439382709245709e-06im), efI=3.0888646354205217e-10, efH=-4.555490773096454e-13, jfI=6.532546417072644e-09, jfH=-9.634269688139735e-12, stepsth=96),
    (l=2, m=1, k=0, a=0.9, p=12.0, e=0.0, x=0.5, kmax=4, spar=0.1, omega=0.024088938268160325, CP=(-6.443876388395458e-07 + -3.903433891758475e-06im), CM=(3.7002745323876123e-06 + -3.739252982755331e-05im), efI=2.1464735142247688e-09, efH=-1.313270316146306e-12, jfI=8.910619016620922e-08, jfH=-5.451756742148024e-11, stepsth=32),
]
