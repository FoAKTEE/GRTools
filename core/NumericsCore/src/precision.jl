# BigFloat precision policy + precision-inheriting tolerance/step selectors.
#
# Canonicalizes the inline BigFloat plumbing scattered across SWSH/Teukolsky/H1Lorenz/KLC into one
# type-generic policy. Every selector threads the element type T so BigFloat precision (256/512/1024/
# 2048 bits) is inherited; a hard Float64 floor would CAP the result at ~14 digits no matter the
# mantissa (the bug the H1Lorenz tolerances were written to avoid).

"""`bc_tol(::Type{T})` — series/boundary-condition truncation tolerance, precision-inherited.
Float64 keeps the validated C value 1e-14; BigFloat targets ~2^(10-bits) with an 8*eps floor so the
asymptotic series actually reaches full mantissa precision. (H1Lorenz `_bc_tol`.)"""
function bc_tol(::Type{T}) where {T}
    if real(T) <: BigFloat
        bits = precision(real(T))
        return max(big(2.0)^(10 - bits), 8*eps(real(T)))
    end
    return T(1e-14)
end

"""`rk_tol(::Type{T})` — adaptive-RK step/solve tolerance, precision-inherited. Float64 keeps the C
value 1e-12; BigFloat targets 2^(20-bits). (H1Lorenz `_rk_tol`.)"""
function rk_tol(::Type{T}) where {T}
    if real(T) <: BigFloat
        bits = precision(real(T))
        return big(2.0)^(20 - bits)
    end
    return T(1e-12)
end

"""`fd_step(::Type{T}; kind)` — finite-difference step for background-field derivatives, precision-
inherited. `kind=:radial` -> 1e-6, `kind=:theta` -> 1e-4 (the validated DiscEMRIForce `_fdr`/`_fdθh`
scales). Written as `T(1)/T(10)^k` so BigFloat inherits precision rather than a Float64 literal."""
function fd_step(::Type{T}; kind::Symbol=:radial) where {T}
    kind === :radial && return T(1)/T(10)^6
    kind === :theta  && return T(1)/T(10)^4
    throw(ArgumentError("fd_step: kind must be :radial or :theta, got $kind"))
end

"""`bigprec_for_mag(magpow)` — BigFloat mantissa bits sized for a working magnitude of ~10^magpow
digits: `clamp(128 + 6*ceil(magpow), 160, 320)` (KerrLorenzCirc `_bigprec_for_mag`)."""
bigprec_for_mag(magpow) = clamp(128 + 6*ceil(Int, magpow), 160, 320)

"""`with_precision(f, bits)` — run `f()` under global BigFloat precision `bits`, restoring the prior
precision afterwards (via `setprecision`'s do-block, so it is restored on exception too)."""
with_precision(f, bits::Integer) = setprecision(f, BigFloat, bits)
