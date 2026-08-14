# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_spin_mode_spherical_correction_analytic.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from spin_mode_spherical_correction_analytic.json
# (gen/gen_spin_mode_spherical_correction_analytic.wl). BHPToolkit (TeukolskySpinFluxes`) oracle:
# TeukolskySpinModeSphericalCorrectionAnalytical[l,m,k,orbit,{angparNew,TeukolskySolverHS1spin}],
# geodesic spherical orbit a=0.9,p=12,e=0,x=0.5, hh=1e-7 (hand-built callbacks' ω-deriv FD step).
# Geodesic Fluxes (efI/efH/jfI/jfH + CarterConstantK kfI/kfH) + FluxesCorrection (dE1*/dL1*/dK1*).
# NB: the .wl's Υτz (.wl L4500) calls the absent KerrGeoMinoFrequencyτz; the oracle overrides it with the
# proper-time polar frequency KerrGeoProperTime[go]["Υτz"] (the genuinely-intended quantity), so dK1 is valid.
const REF_SPIN_MODE_SPHERICAL_CORR_ANALYTIC = [
    (l=2, m=2, k=0, a=0.9, p=12.0, e=0.0, x=0.5, hh=1.0e-7, omega=0.04835306488146746,
     efI=3.5191978105466852e-6, efH=-6.012029336798213e-9, jfI=1.455625540665782e-4, jfH=-2.486721100942032e-7,
     kfI=8.61300898104153e-4, kfH=-1.4714052877887018e-6,
     dE1I=4.165590274513807e-6, dE1H=-1.0971586506299655e-8, dL1I=1.181763729628846e-4, dL1H=-3.6135105791051517e-7,
     dK1I=7.783425470828508e-4, dK1H=-2.273240781745112e-6, stepsth=64),
    (l=2, m=1, k=0, a=0.9, p=12.0, e=0.0, x=0.5, hh=1.0e-7, omega=0.02417653244073373,
     efI=2.055274284819029e-9, efH=-1.496404952198038e-12, jfI=8.501112762375337e-8, jfH=-6.189493699587895e-11,
     kfI=5.030150854435799e-7, kfH=-3.6623543166374515e-10,
     dE1I=4.627466314335466e-9, dE1H=-2.2750572813564757e-12, dL1I=1.597946974890873e-7, dL1H=-7.108833436843269e-11,
     dK1I=9.91701448542366e-7, dK1H=-4.542619811026582e-10, stepsth=32),
]
