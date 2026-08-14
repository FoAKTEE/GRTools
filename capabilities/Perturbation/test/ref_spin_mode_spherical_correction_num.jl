# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_spin_mode_spherical_correction_num.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# AUTO-DERIVED from spin_mode_spherical_correction_num.json (gen/gen_spin_mode_spherical_correction_num.wl).
# BHPToolkit (TeukolskySpinFluxes`) oracle: TeukolskySpinModeSphericalCorrectionNum[l,m,k,oc,δω],
# a=0.9,p=12,e=0,x=0.5,kmax=4,δω=1e-6. Geodesic Fluxes (efI/efH/jfI/jfH) + FluxesCorrection (dE1*/dL1*).
const REF_SPIN_MODE_SPHERICAL_CORR_NUM = [
    (l=2, m=2, k=0, a=0.9, p=12.0, e=0.0, x=0.5, kmax=4, dw=1.0e-6, omega=0.04835306488146746, efI=3.5192343863388845e-06, efH=-6.012029847027697e-09, jfI=0.00014556406693002496, jfH=-2.486721311985317e-07, dE1I=-6.019637484977878e-07, dE1H=1.1009896033085827e-09, dL1I=-1.9624737735625075e-05, dL1H=3.652994116213991e-08, stepsth=64),
    (l=2, m=2, k=2, a=0.9, p=12.0, e=0.0, x=0.5, kmax=4, dw=1.0e-6, omega=0.09483988428389654, efI=3.49842994203846e-10, efH=-5.331271331295798e-13, jfI=7.377550001149633e-09, jfH=-1.1242677849198993e-11, dE1I=-4.235612940170271e-10, dE1H=8.061317721905654e-13, dL1I=-8.721003312964943e-09, dL1H=1.6678103805854036e-11, stepsth=96),
    (l=2, m=1, k=0, a=0.9, p=12.0, e=0.0, x=0.5, kmax=4, dw=1.0e-6, omega=0.02417653244073373, efI=2.0552724051659005e-09, efH=-1.4964049518850133e-12, jfI=8.501104987674259e-08, jfH=-6.18949369829315e-11, dE1I=9.14921250657754e-10, dE1H=1.89278648757813e-12, dL1I=4.092340284193488e-08, dL1H=7.604771834712704e-11, stepsth=32),
]
