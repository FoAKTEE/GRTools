# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_spin_mode_spherical_correction.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from spin_mode_spherical_correction.json (gen/gen_spin_mode_spherical_correction.wl).
# BHPToolkit (TeukolskySpinFluxes`) oracle: TeukolskySpinModeSphericalCorrection[l,m,k,oc,odv,
# {angparNew,TeukolskySolverHS1spin}], a=0.9,p=12,e=0,x=0.5,kmax=4, hh=1e-6 (ω-deriv FD), h=1e-5
# (orbit r/x-deriv FD). Geodesic Fluxes (efI/efH/jfI/jfH) + FluxesCorrection (dE1*/dL1*) +
# FluxesDerivatives r/x (dEdr*/dLdr*/dEdx*/dLdx*).
const REF_SPIN_MODE_SPHERICAL_CORR = [
    (l=2, m=2, k=0, a=0.9, p=12.0, e=0.0, x=0.5, kmax=4, hh=1.0e-6, h=1.0e-5,
     omega=0.04835306488146746,
     efI=3.519197810546684e-06, efH=-6.012029336798213e-09, jfI=0.0001455625540665781, jfH=-2.486721100942032e-07,
     dE1I=-6.018986709506123e-07, dE1H=1.1009910802461443e-09, dL1I=-1.9622100783154188e-05, dL1H=3.653000301649091e-08,
     dEdrI=-1.465952140344973e-06, dEdrH=3.7211552688983727e-09, dLdrI=-4.2390908133605674e-05, dLdrH=1.2274809759219085e-07,
     dEdxI=-1.6952363359434444e-06, dEdxH=7.894797827637291e-11, dLdxI=-6.159824647295683e-05, dLdxH=-1.1291113674057449e-08,
     stepsth=64),
    (l=2, m=1, k=0, a=0.9, p=12.0, e=0.0, x=0.5, kmax=4, hh=1.0e-6, h=1.0e-5,
     omega=0.02417653244073373,
     efI=2.055274284819029e-09, efH=-1.4964049521980386e-12, jfI=8.501112762375337e-08, jfH=-6.189493699587898e-11,
     dE1I=9.149225545066096e-10, dE1H=1.892786487903043e-12, dL1I=4.0923459589142516e-08, dL1H=7.604771836009714e-11,
     dEdrI=-1.214071981412752e-09, dEdrH=1.1435958761152815e-12, dLdrI=-3.956189019037498e-08, dLdrH=3.95441499832605e-11,
     dEdxI=-2.815879130094579e-09, dEdxH=-4.336207398832886e-12, dLdxI=-1.1149526645019393e-07, dLdxH=-1.8297921394022818e-10,
     stepsth=32),
]
