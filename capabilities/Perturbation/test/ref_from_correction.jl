# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_from_correction.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from from_correction.json via gen/json_to_jl.py. BHPToolkit (SpinWeightedSpheroidalHarmonics`) oracle.
const REF_FROM_CORRECTION = [
    (l=2, m=2, n=0, k=0, a=0.9, p=12.0, e=0.2, x=0.5, nmax=4, kmax=4, spar=0.1, omega=0.04581018581308623, CP=(-0.0002270885364309784 + 5.225378650833227e-05im), CM=(-0.0003073677457312951 + -8.030309348078012e-05im), efI=2.05903226562483e-06, efH=-3.2083234377147435e-09, jfI=8.989408050105061e-05, jfH=-1.4007030885250186e-07, cfKI=0.0005452747815357657, cfKH=-8.496311062250928e-07, stepsr=64, stepsth=64),
    (l=2, m=2, n=1, k=0, a=0.9, p=12.0, e=0.2, x=0.5, nmax=4, kmax=4, spar=0.1, omega=0.06324767896995062, CP=(-0.00024035299350130682 + 6.6018904803043e-05im), CM=(-0.00025673182966304525 + -7.885637633271531e-05im), efI=1.2359153446393776e-06, efH=-2.883555890346618e-09, jfI=3.908176125250601e-05, jfH=-9.118297895853581e-08, cfKI=0.00023620223283260286, cfKH=-5.510914179936903e-07, stepsr=96, stepsth=64),
    (l=2, m=1, n=0, k=0, a=0.9, p=12.0, e=0.2, x=0.5, nmax=4, kmax=4, spar=0.1, omega=0.022905092906543113, CP=(-5.096673368656627e-07 + -3.199891144574127e-06im), CM=(3.455101738501221e-06 + -3.5403199299545535e-05im), efI=1.5924896159749645e-09, efH=-1.1337059764630643e-12, jfI=6.952556894104764e-08, jfH=-4.94958034481147e-11, cfKI=4.2172453630065777e-07, cfKH=-3.0022903912780465e-10, stepsr=32, stepsth=32),
]
