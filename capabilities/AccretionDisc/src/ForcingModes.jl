"""
Shared forcing vocabulary for Newtonian eccentric discs and relativistic disc response.
Owned by ECC::shared-forcing-api; it contains no orbit, potential or response solver.
"""
module ForcingModes

export HarmonicMode, NewtonianMode, mode_power
export ClosureTag, Adiabatic, LocallyIsothermalOptionA, Cooling
export SofteningPolicy, FixedSoftening, SemimajorScaleHeight, InstantaneousScaleHeight,
       LocalScaleHeight
export ModeMaskTag, FR22Mask, FR24Mask, FD25Mask

"""
    HarmonicMode(m, k, omega)
    HarmonicMode{T}(m, k, omega)

A harmonic with phase `m*phi - omega*t`. `omega` is the finite, real TOTAL
forcing frequency, in inverse units of the chosen coordinate time. In GR the
caller supplies `omega = m*Omega_phi + k*Omega_r`; its sign may be negative or zero.
The pattern speed is `omega/m` only for `m > 0`.

Indices are stored as `Int`, with `m >= 0` and signed `k`. The constructor admits
axisymmetric, static and secular modes; mode-selection/feedback policies decide
which to use. Dynamic independent axisymmetric modes use `k > 0`, and `(0,0)`
is excluded from dissipative feedback by those policies, not by this container.
The untyped constructor floats the frequency (integer inputs give Float64;
BigFloat inputs retain their type). Explicit `T <: Real` is also supported.
"""
struct HarmonicMode{T<:Real}
    m::Int
    k::Int
    omega::T
    function HarmonicMode{T}(m::Integer, k::Integer, omega::Real) where {T<:Real}
        m >= 0 || throw(DomainError(m, "HarmonicMode requires m >= 0"))
        w = convert(T, omega)
        isfinite(w) || throw(DomainError(w, "total forcing frequency must be finite"))
        new{T}(Int(m), Int(k), w)
    end
end

HarmonicMode(m::Integer, k::Integer, omega::Real) =
    HarmonicMode{typeof(float(omega))}(m, k, omega)

"""
    NewtonianMode(m, l, n)

Adapter returning `HarmonicMode(m, l-m, l*n)` for the Newtonian phase
`m*phi - l*n*t`. Mean motion `n` must be finite and positive. The pattern speed
is `(l/m)*n` for `m > 0`; it is NOT the stored frequency. Signed and zero `l`
are retained so paper-specific exclusions remain the mode-selection owner's job.
"""
function NewtonianMode(m::Integer, l::Integer, n::Real)
    nn = float(n)
    isfinite(nn) && nn > 0 || throw(DomainError(n, "mean motion must be finite and positive"))
    mi, li = Int(m), Int(l)
    HarmonicMode(mi, Base.Checked.checked_sub(li, mi), li * nn)
end

"""
    mode_power(mode::HarmonicMode, Tdisc::Real)

Gravitational power into the disc, `(mode.omega / mode.m) * Tdisc`, with the
same reaction convention as the supplied torque into the disc. Requires `m > 0`;
axisymmetric power needs a separate integrand even when its torque is zero.
This is gravitational work, with no cooling/heat-reservoir contribution.
Frequency and torque are promoted before division so a BigFloat torque does not
inherit a Float64-rounded pattern speed.
"""
function mode_power(mode::HarmonicMode, Tdisc::Real)
    mode.m > 0 || throw(DomainError(mode.m, "axisymmetric power cannot be inferred from torque"))
    omega, torque = promote(mode.omega, Tdisc)
    (omega / mode.m) * torque
end

"""
Thermodynamic model tag: `Adiabatic`, `LocallyIsothermalOptionA`, or `Cooling`.
OPTION A means the prescribed Eulerian `p = Theta(r)*epsilon` law with conserved
stress energy and no bath four-force. These are tags only; cooling parameters,
closure equations and their evaluation belong to the downstream disc models.
"""
@enum ClosureTag Adiabatic LocallyIsothermalOptionA Cooling

"""
Abstract softening policy with scalar type `T`. Policies specify where a later
potential evaluator must evaluate the disc scale height; they do not evaluate it.
"""
abstract type SofteningPolicy{T<:Real} end

# The same finite, nonnegative invariant applies to lengths and all scale factors.
function _softening_parameter(x::Real)
    isfinite(x) && x >= 0 || throw(DomainError(x, "softening parameter must be finite and nonnegative"))
    x
end

"""
    FixedSoftening(eps)
    FixedSoftening{T}(eps)

A constant length `eps`, independent of time and field radius. FD25 supplies
`eps = (3//10)*h_p*a_p`. Zero is allowed as an explicitly unsoftened policy.
"""
struct FixedSoftening{T<:Real} <: SofteningPolicy{T}
    eps::T
    FixedSoftening{T}(eps::Real) where {T<:Real} =
        new{T}(_softening_parameter(convert(T, eps)))
end
FixedSoftening(eps::Real) = FixedSoftening{typeof(float(eps))}(eps)

"""
    SemimajorScaleHeight(b)
    SemimajorScaleHeight{T}(b)

Length `b*H(a)`, frozen at the orbit's semimajor axis. FR22's provisional fixed-H
reading uses `b = 0.3`; the competing reading is `InstantaneousScaleHeight`.
"""
struct SemimajorScaleHeight{T<:Real} <: SofteningPolicy{T}
    b::T
    SemimajorScaleHeight{T}(b::Real) where {T<:Real} =
        new{T}(_softening_parameter(convert(T, b)))
end
SemimajorScaleHeight(b::Real) = SemimajorScaleHeight{typeof(float(b))}(b)

"""
    InstantaneousScaleHeight(b)
    InstantaneousScaleHeight{T}(b)

Length `b*H(r_p(t))`, evaluated at the perturber's instantaneous radius. This is
FR24's literal policy and FR22's alternate reading; it varies along an eccentric
orbit even for constant aspect ratio (`qTemp = 1`).
"""
struct InstantaneousScaleHeight{T<:Real} <: SofteningPolicy{T}
    b::T
    InstantaneousScaleHeight{T}(b::Real) where {T<:Real} =
        new{T}(_softening_parameter(convert(T, b)))
end
InstantaneousScaleHeight(b::Real) = InstantaneousScaleHeight{typeof(float(b))}(b)

"""
    LocalScaleHeight(b)
    LocalScaleHeight{T}(b)

Length `b*H(r)` at the field radius: a separately labelled FR22 sensitivity case.
Potential radial derivatives must include derivatives of this length; the shared
policy does not implement those derivatives. All scale factors are finite and
nonnegative, floated by default, with BigFloat inherited or explicit `T`.
"""
struct LocalScaleHeight{T<:Real} <: SofteningPolicy{T}
    b::T
    LocalScaleHeight{T}(b::Real) where {T<:Real} =
        new{T}(_softening_parameter(convert(T, b)))
end
LocalScaleHeight(b::Real) = LocalScaleHeight{typeof(float(b))}(b)

"""
Paper-mask identity tags `FR22Mask`, `FR24Mask`, `FD25Mask`. Enumeration, cutoffs,
domain/resonance exclusions and exclusion reasons belong to downstream
`newtonian_mode_set`; these tags alone do not accept or reject any mode.
"""
@enum ModeMaskTag FR22Mask FR24Mask FD25Mask

end # module ForcingModes
