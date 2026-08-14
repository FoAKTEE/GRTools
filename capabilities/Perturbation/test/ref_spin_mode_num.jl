# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_spin_mode_num.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from spin_mode_num.json via gen/json_to_jl.py. BHPToolkit (SpinWeightedSpheroidalHarmonics`) oracle.
const REF_SPIN_MODE_NUM = [
    (l=2, m=2, n=0, k=0, a=0.9, p=12.0, e=0.2, x=0.5, nmax=4, kmax=4, spar=0.1, omega=0.04581016095690081, CP=(-0.00022709117833139028 + 5.225095712351999e-05im), CM=(-0.0003073622172648629 + -8.029992621254215e-05im), efI=2.0590687876799386e-06, efH=-3.208197908260715e-09, jfI=8.989572377259949e-05, jfH=-1.4006490443371536e-07, stepsr=64, stepsth=64),
    (l=2, m=2, n=1, k=0, a=0.9, p=12.0, e=0.2, x=0.5, nmax=4, kmax=4, spar=0.1, omega=0.06324766632442667, CP=(-0.00024034410869636867 + 6.601835355219153e-05im), CM=(-0.0002567321843565122 + -7.885805328356416e-05im), efI=1.23582942974299e-06, efH=-2.8835733666231465e-09, jfI=3.9079052289576873e-05, jfH=-9.118354981927582e-08, stepsr=96, stepsth=64),
    (l=2, m=1, n=0, k=0, a=0.9, p=12.0, e=0.2, x=0.5, nmax=4, kmax=4, spar=0.1, omega=0.022905080478450406, CP=(-5.097263073014421e-07 + -3.2001682586832127e-06im), CM=(3.4571203355681695e-06 + -3.541195516321694e-05im), efI=1.592769472623349e-09, efH=-1.134273567694625e-12, jfI=6.953782476869535e-08, jfH=-4.9520610449798416e-11, stepsr=32, stepsth=32),
]
