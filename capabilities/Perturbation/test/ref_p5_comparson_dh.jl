# BAKED FIXTURE (Phase-7B rule): copied 2026-07-02 from packages/TeukolskySpinFluxes/test/ref_p5_comparson_dh.jl.
# The values are BHPToolkit MATHEMATICA-ORACLE outputs (not legacy-Julia outputs) baked as
# literals; tolerances in test_spinfluxes_modes.jl are FROZEN verbatim from the legacy suite
# packages/TeukolskySpinFluxes/test/runtests.jl. No live packages/* dep at test time.
# EXTRACTED from comparson_DH.nb (Supplement for arXiv:2603.13482), cached Out[50] cell:
# Table[ TeukolskySpinModeAnalytical[2,2,0,-3, orbits[[i]]], {i,1,4} ] with
# orbits[[i]] = KerrSpinOrbit[0.9,12,0.2,0.5, ({-2,-1,1,2}*0.05)[[i]], "Gauge"->"pex"].
# These are the published TeukolskySpinModeAnalytical golden values for the DH P5 orbit.
# The headline 3.81833507e-15 is the cached modesk[[1]] (k=-10, DH NUMERIC trajectory,
# arXiv:2303.16798) Energy.ℐ = |Amp_ℐ|²/(4πω²) — a DIFFERENT solver, NOT reproducible
# from this mirror (see ref-code .nb line 60; reconciliation in p5_reproduction.md).
const REF_P5_COMPARSON_DH = [
    (l=2, m=2, n=0, k=-3, spar=-0.1,  gauge="pex",
     omega=-0.02030718076784457812582922965898259842,
     AmpI=2.2827420197712649240031550727012e-7 - 1.56834494002737166934715262532249e-6im,
     EnergyI=4.847061269764846532736116492e-10, EnergyH=1.2215508949353275957887431111e-10,
     AngMomI=-4.77374119546907233490227562761e-8, stepsr=32, stepsz=160),
    (l=2, m=2, n=0, k=-3, spar=-0.05, gauge="pex",
     omega=-0.02030791725092721650073867608325950968,
     AmpI=2.3333049632520984595190624660612e-7 - 1.60337277314643727814512908640056e-6im,
     EnergyI=5.0655848914589795161350242192e-10, EnergyH=1.1846678659425390460192730507e-10,
     AngMomI=-4.988778345773194083798612205003e-8, stepsr=32, stepsz=160),
    (l=2, m=2, n=0, k=-3, spar=0.05,  gauge="pex",
     omega=-0.0203337608617653238959187501529895308,
     AmpI=2.4478777090470734867928317425699e-7 - 1.6805836844942657872611619444005e-6im,
     EnergyI=5.5512721447135403012358198426e-10, EnergyH=1.1200854343257828518167057706e-10,
     AngMomI=-5.460152878213394493193461819398e-8, stepsr=32, stepsz=192),
    (l=2, m=2, n=0, k=-3, spar=0.1,   gauge="pex",
     omega=-0.0203589999299665112454693856198783244,
     AmpI=2.5125195304209137626310008510799e-7 - 1.72305250243270416015652436354069e-6im,
     EnergyI=5.8211905848181787128018708188e-10, EnergyH=1.0919844853020969718829939902e-10,
     AngMomI=-5.718542762260085195200232307056e-8, stepsr=32, stepsz=192),
]
