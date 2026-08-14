# 3-harmonic spin-weighted angular coupling — the SecondOrder front door onto Wigner.coupling.
# Wigner is the canonical HOME of this primitive (it dropped SecondOrderRicci's WignerSymbols dep and
# is GSL-grounded vs gsl_sf_coupling_3j); SecondOrder does NOT re-own it (AHA). `coupling` is brought
# into this namespace by the `using Wigner: coupling` in SecondOrder.jl and re-exported, so there is
# exactly one implementation. This file documents that delegation and exists as the explicit
# SecondOrder/Coupling.jl crosswalk slot for SecondOrderRicci/coupling.jl.

"""
    coupling([T=Float64,] lc, mc, sc, la, ma, sa, lb, mb, sb) -> T

Angular 3-harmonic coupling coefficient (the SecondOrder front door — identical object to
[`Wigner.coupling`](@ref)): couples two spin-weighted harmonics `(la,ma,sa)` and `(lb,mb,sb)` onto a
third `(lc,mc,sc)` for the second-order Ricci summation, returning 0 on any 3j selection-rule
violation. Type-generic `{T}` (feeds both the √-normalisation and the 3j evaluation) with a Float64
convenience method. The legacy `SecondOrderRicci/coupling.jl/coupling` lives canonically in `Wigner`.
"""
coupling
