# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_incl_spin_mode.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from incl_spin_mode.json via gen/json_to_jl.py. BHPToolkit (SpinWeightedSpheroidalHarmonics`) oracle.
const REF_INCL_SPIN_MODE = [
    (l=2, m=2, n=0, k=0, spar=0.0, a=0.9, p=12.0, e=0.2, x=0.5, omega=0.045986687799471085, CP=(-0.00022894681956752257 + 5.2797617981444476e-05im), efI=2.0772954417177845e-06, efH=-3.2441238012192917e-09, jfI=9.034333808844901e-05, jfH=-1.4108969166753533e-07),
    (l=2, m=2, n=0, k=0, spar=0.1, a=0.9, p=12.0, e=0.2, x=0.5, omega=0.04776919937577677, CP=(-0.00025162180206300856 + 5.961535739516104e-05im), efI=2.331897415292784e-06, efH=-3.910615335230742e-09, jfI=9.763184000422094e-05, jfH=-1.6372957413281547e-07),
    (l=2, m=2, n=1, k=0, spar=0.0, a=0.9, p=12.0, e=0.2, x=0.5, omega=0.06333747392664302, CP=(-0.0002439847349477071 + 6.706541113670015e-05im), efI=1.2700686658152388e-06, efH=-2.9601063884513457e-09, jfI=4.010480958827701e-05, jfH=-9.347093292289233e-08),
    (l=2, m=2, n=0, k=2, spar=0.0, a=0.9, p=12.0, e=0.2, x=0.5, omega=0.09018898523276414, CP=(1.7417673883457098e-06 + -5.590533004977698e-07im), efI=0.0, efH=0.0, jfI=7.259784922810668e-10, jfH=0.0),
    (l=2, m=2, n=0, k=2, spar=0.1, a=0.9, p=12.0, e=0.2, x=0.5, omega=0.09368484279340614, CP=(1.923543787241171e-06 + -6.290802859278029e-07im), efI=0.0, efH=0.0, jfI=7.927694792256463e-10, jfH=0.0),
]
