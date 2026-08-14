"""
    angular_precision_floor([T=BigFloat]) -> (value_tier, eigenvalue_tier)

Report the angular-solve precision policy as a `NamedTuple`:

  - `value_tier = eps(Float64)` — the harmonic-VALUE CEILING. `S_lm(θ)` and its θ-derivative go
    through Float64 expansion coefficients and a LAPACK eigensolve, so they are Float64-capped
    irrespective of BigFloat input (design §4.2; the accepted KerrLorenzCirc cap).
  - `eigenvalue_tier = eps(T)` — the spheroidal-EIGENVALUE FLOOR. The type-generic Leaver
    continued-fraction path inherits the precision of `T`; with the default `T=BigFloat` this
    TRACKS the active `setprecision(BigFloat, bits)` (so calling under `setprecision(BigFloat,512)`
    reports the 512-bit floor).

The asymmetry `eigenvalue_tier ≤ value_tier` is the contract the BigFloat radial chain relies on:
the eigenvalue is BigFloat-native while the angular value carrier is not.
"""
angular_precision_floor(::Type{T}=BigFloat) where {T<:AbstractFloat} =
    (value_tier = eps(Float64), eigenvalue_tier = eps(T))
