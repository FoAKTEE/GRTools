# Bindings of the SpinFluxes legacy names onto the SINGLE owner of the linear-in-spin
# orbit-correction engine, `Spacetime.SpinningSecondary` (numeric `Correction` engine, the 11
# `_oc_*` four-velocity/metric helpers, and the `spin_orbit_correction`/`spin_orbit_numeric`
# assemblies). The former local copy `spinning_orbit_correction_local.jl` was code-identical
# (docstrings aside) and is gone; the BHPToolkit-facing names below are kept because the
# TeukolskySpinMode* entrypoints, docs and tests speak them.
using Spacetime.SpinningSecondary: spin_orbit_correction, spin_orbit_numeric,
                                    _numeric_correction, _chith_hat, _d2z_dλ2,
                                    _oc_h1, _oc_h2, _oc_h3,
                                    _oc_Sgtt, _oc_Sgphph, _oc_Sgtph,
                                    _oc_dSgttdr, _oc_dSgttdz, _oc_dSgtphdr,
                                    _oc_dSgphphdr, _oc_dSgphphdz

"""
    KerrSpinOrbitCorrection(a, p, e, x, nmax, kmax) -> NamedTuple

BHPToolkit name of `Spacetime.SpinningSecondary.spin_orbit_correction` (SpinningTrajectory.wl
L907-911): numeric `Correction` engine + orbit-correction assembly. `kmax` must be even. Unexported.
"""
const KerrSpinOrbitCorrection = spin_orbit_correction

"""
    KerrSpinOrbitNum(orbitCorrection, spar) -> NamedTuple

BHPToolkit name of `Spacetime.SpinningSecondary.spin_orbit_numeric` (SpinningTrajectory.wl
L914-938): the spin-baked frequency assembly of a `KerrSpinOrbitCorrection` object. Unexported.
"""
const KerrSpinOrbitNum = spin_orbit_numeric
