abstract type WeylScalarContent end

"Radial differentiation at fixed outgoing Schwarzschild time `u`."
struct FixedOutgoingURadialDerivative end

"Finite-radius free Hertz IVP data; this policy is neither retarded nor physical."
struct FiniteHertzRadialIVP end

"The reduced Weyl scalar `psi0L = T0[h]`, linear in the metric perturbation at order `N`."
struct LinearizedPsi0 <: WeylScalarContent end

"The full perturbative-order-two Weyl coefficient, including quadratic/tetrad-rotation pieces."
struct FullSecondOrderWeyl <: WeylScalarContent end

"Canonical outgoing-Schwarzschild IRG convention `psi0L = (1/4) d_r^4 bar(Phi)`."
struct QuarterThorn4 end

"""
    ReducedWeylMode{N}(psi0; M, l, m, omega, source_fingerprint, phase,
                        content=LinearizedPsi0(), tetrad=KinnersleyMinusPlus(),
                        slicing=SchwarzschildOutgoing())

A separated radial Weyl mode with enough metadata to prevent cross-mode Hertz inversion.  At
second order, `content=LinearizedPsi0()` means the part linear in `h2`; a
`FullSecondOrderWeyl()` carrier can be represented for an explicit diagnostic rejection, but it
cannot be inverted because it also contains first-order quadratic and tetrad-rotation terms.
`phase` is a required producer assertion of the canonical
`exp(-i*omega*u + i*m*phi)` separated-mode convention.
"""
struct ReducedWeylMode{N,T<:AbstractFloat,W<:WeylScalarContent,F}
    psi0::F
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    source_fingerprint::UInt128
    content::W
    tetrad::KinnersleyMinusPlus
    slicing::SchwarzschildOutgoing
    phase::ExpMinusIomegaUPlusImPhi
end

function ReducedWeylMode{N}(psi0; M, l::Integer, m::Integer, omega=zero(M),
        source_fingerprint, content::WeylScalarContent=LinearizedPsi0(),
        tetrad=KinnersleyMinusPlus(), slicing=SchwarzschildOutgoing(), phase) where {N}
    N isa Integer && !(N isa Bool) && N > 0 ||
        throw(ArgumentError("perturbative order N must be a positive integer"))
    tetrad isa KinnersleyMinusPlus ||
        throw(ArgumentError("Hertz inversion requires the Kinnersley (-,+,+,+) tetrad"))
    slicing isa SchwarzschildOutgoing ||
        throw(ArgumentError("Hertz inversion requires outgoing Schwarzschild slicing"))
    phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "reduced Weyl phase must be exp(-i*omega*u + i*m*phi)"))
    lI = _checked_mode_index(l, "l")
    mI = _checked_mode_index(m, "m")
    lI >= 2 || throw(ArgumentError("a spin-2 reduced Weyl mode requires l >= 2"))
    -lI <= mI <= lI || throw(ArgumentError("mode requires |m| <= l"))
    _require_nonbool_number(M, "M")
    _require_nonbool_number(omega, "omega")
    T = promote_type(typeof(float(M)), typeof(float(real(omega))), typeof(float(imag(omega))))
    T <: AbstractFloat || throw(ArgumentError("M and omega must promote to a floating type"))
    MT, omegaT = T(M), Complex{T}(omega)
    _finite(MT) && MT > zero(T) || throw(ArgumentError("M must be finite and positive"))
    _finite(omegaT) || throw(ArgumentError("omega must be finite"))
    fp = _checked_fingerprint(source_fingerprint, "source_fingerprint")
    ReducedWeylMode{N,T,typeof(content),typeof(psi0)}(
        psi0, MT, lI, mI, omegaT, fp, content, tetrad, slicing, phase)
end

perturbative_order(::ReducedWeylMode{N}) where {N} = N

function _validate_reduced_weyl_mode(src::ReducedWeylMode{N,T}) where {N,T}
    N isa Integer && !(N isa Bool) && N > 0 ||
        throw(ArgumentError("reduced Weyl perturbative order must be positive"))
    _finite(src.M) && src.M > zero(T) ||
        throw(ArgumentError("reduced Weyl M must be finite and positive"))
    src.l >= 2 && -src.l <= src.m <= src.l ||
        throw(ArgumentError("invalid spin-2 reduced Weyl mode key"))
    _finite(src.omega) || throw(ArgumentError("reduced Weyl omega must be finite"))
    src.source_fingerprint != 0 ||
        throw(ArgumentError("reduced Weyl source fingerprint must be nonzero"))
    src.content isa WeylScalarContent ||
        throw(ArgumentError("invalid reduced Weyl content marker"))
    src.tetrad isa KinnersleyMinusPlus || throw(ArgumentError(
        "reduced Weyl tetrad is not Kinnersley (-,+,+,+)"))
    src.slicing isa SchwarzschildOutgoing ||
        throw(ArgumentError("reduced Weyl mode is not on outgoing Schwarzschild slices"))
    src.phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "reduced Weyl phase is not exp(-i*omega*u + i*m*phi)"))
    nothing
end

"""
    HertzRadialBoundary(mode; rb, b0=0, b1=0, b2=0, b3=0,
        boundary_fingerprint, radial_derivative, policy,
        normalization=QuarterThorn4())

Keyed inner-ray data for `bar(Phi)` and its first three radial derivatives.  The independent,
nonzero boundary fingerprint records the provenance of these four integration constants.
`radial_derivative` and `policy` are mandatory producer assertions: the data are derivatives at
fixed outgoing `u` for a finite-radius free IVP, not fixed-`t`, retarded, or physical data.
"""
struct HertzRadialBoundary{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    source_fingerprint::UInt128
    tetrad::KinnersleyMinusPlus
    slicing::SchwarzschildOutgoing
    phase::ExpMinusIomegaUPlusImPhi
    radial_derivative::FixedOutgoingURadialDerivative
    policy::FiniteHertzRadialIVP
    rb::T
    b0::Complex{T}
    b1::Complex{T}
    b2::Complex{T}
    b3::Complex{T}
    boundary_fingerprint::UInt128
    normalization::QuarterThorn4
    gauge::HertzIRG
end

function HertzRadialBoundary(src::ReducedWeylMode{N,T}; rb, b0=zero(T), b1=zero(T),
        b2=zero(T), b3=zero(T), boundary_fingerprint,
        normalization=QuarterThorn4(), gauge=HertzIRG(), radial_derivative,
        policy) where {N,T}
    normalization isa QuarterThorn4 || throw(ArgumentError(
        "Hertz boundary normalization must be canonical positive QuarterThorn4"))
    gauge isa HertzIRG || throw(ArgumentError("Hertz inversion requires IRG"))
    radial_derivative isa FixedOutgoingURadialDerivative || throw(ArgumentError(
        "Hertz boundary data must be radial derivatives at fixed outgoing u"))
    policy isa FiniteHertzRadialIVP || throw(ArgumentError(
        "Hertz boundary data require the finite radial IVP policy"))
    rbT = _checked_real(T, rb, "Hertz boundary rb")
    rbT >= 2src.M && _finite(rbT) ||
        throw(ArgumentError("Hertz boundary rb must be finite and satisfy rb >= 2M"))
    data = (_checked_complex(T, b0, "Hertz boundary b0"),
        _checked_complex(T, b1, "Hertz boundary b1"),
        _checked_complex(T, b2, "Hertz boundary b2"),
        _checked_complex(T, b3, "Hertz boundary b3"))
    all(_finite, data) || throw(ArgumentError("Hertz boundary derivatives must be finite"))
    fp = _checked_fingerprint(boundary_fingerprint, "boundary_fingerprint")
    HertzRadialBoundary{N,T}(src.M, src.l, src.m, src.omega, src.source_fingerprint,
        src.tetrad, src.slicing, src.phase, radial_derivative, policy,
        rbT, data..., fp, normalization, gauge)
end

perturbative_order(::HertzRadialBoundary{N}) where {N} = N

function _validate_hertz_boundary(bc::HertzRadialBoundary{N,T}) where {N,T}
    N isa Integer && !(N isa Bool) && N > 0 ||
        throw(ArgumentError("Hertz boundary perturbative order must be positive"))
    _finite(bc.M) && bc.M > zero(T) ||
        throw(ArgumentError("Hertz boundary M must be finite and positive"))
    bc.l >= 2 && -bc.l <= bc.m <= bc.l ||
        throw(ArgumentError("invalid spin-2 Hertz boundary mode key"))
    _finite(bc.omega) || throw(ArgumentError("Hertz boundary omega must be finite"))
    bc.source_fingerprint != 0 ||
        throw(ArgumentError("Hertz boundary source fingerprint must be nonzero"))
    bc.boundary_fingerprint != 0 ||
        throw(ArgumentError("Hertz boundary fingerprint must be nonzero"))
    bc.tetrad isa KinnersleyMinusPlus || throw(ArgumentError(
        "Hertz boundary tetrad is not Kinnersley (-,+,+,+)"))
    bc.slicing isa SchwarzschildOutgoing ||
        throw(ArgumentError("Hertz boundary is not on outgoing Schwarzschild slices"))
    bc.phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "Hertz boundary phase is not exp(-i*omega*u + i*m*phi)"))
    bc.radial_derivative isa FixedOutgoingURadialDerivative || throw(ArgumentError(
        "Hertz boundary data are not radial derivatives at fixed outgoing u"))
    bc.policy isa FiniteHertzRadialIVP ||
        throw(ArgumentError("Hertz boundary is not a finite-radius free IVP"))
    _finite(bc.rb) && bc.rb >= T(2) * bc.M ||
        throw(ArgumentError("Hertz boundary rb must be finite and satisfy rb >= 2M"))
    all(_finite, (bc.b0, bc.b1, bc.b2, bc.b3)) ||
        throw(ArgumentError("Hertz boundary derivatives must be finite"))
    bc.normalization isa QuarterThorn4 ||
        throw(ArgumentError("Hertz boundary normalization is not QuarterThorn4"))
    bc.gauge isa HertzIRG || throw(ArgumentError("Hertz boundary is not IRG"))
    nothing
end

"A solved radial mode for the conjugate IRG Hertz potential in the `QuarterThorn4` convention."
struct ConjugateHertzMode{N,T<:AbstractFloat,F}
    source::ReducedWeylMode{N,T,LinearizedPsi0,F}
    boundary::HertzRadialBoundary{N,T}
    normalization::QuarterThorn4
    gauge::HertzIRG
    phase::ExpMinusIomegaUPlusImPhi
    radial_derivative::FixedOutgoingURadialDerivative
    policy::FiniteHertzRadialIVP
end

perturbative_order(::ConjugateHertzMode{N}) where {N} = N

function _same_hertz_key(src::ReducedWeylMode, bc::HertzRadialBoundary)
    src.M == bc.M && src.l == bc.l && src.m == bc.m && src.omega == bc.omega &&
        src.source_fingerprint == bc.source_fingerprint &&
        typeof(src.tetrad) === typeof(bc.tetrad) && typeof(src.slicing) === typeof(bc.slicing) &&
        typeof(src.phase) === typeof(bc.phase)
end

function _validate_conjugate_hertz_mode(mode::ConjugateHertzMode)
    src, bc = mode.source, mode.boundary
    _validate_reduced_weyl_mode(src)
    _validate_hertz_boundary(bc)
    _same_hertz_key(src, bc) ||
        throw(ArgumentError("conjugate Hertz mode source and boundary keys differ"))
    mode.normalization isa QuarterThorn4 && bc.normalization isa QuarterThorn4 ||
        throw(ArgumentError("conjugate Hertz mode normalization is not QuarterThorn4"))
    mode.gauge isa HertzIRG && bc.gauge isa HertzIRG ||
        throw(ArgumentError("conjugate Hertz mode is not IRG"))
    src.tetrad isa KinnersleyMinusPlus && bc.tetrad isa KinnersleyMinusPlus ||
        throw(ArgumentError("conjugate Hertz mode tetrad is not Kinnersley (-,+,+,+)"))
    src.slicing isa SchwarzschildOutgoing && bc.slicing isa SchwarzschildOutgoing ||
        throw(ArgumentError("conjugate Hertz mode is not on outgoing Schwarzschild slices"))
    src.phase isa ExpMinusIomegaUPlusImPhi &&
        bc.phase isa ExpMinusIomegaUPlusImPhi &&
        mode.phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "conjugate Hertz mode phase is not exp(-i*omega*u + i*m*phi)"))
    bc.radial_derivative isa FixedOutgoingURadialDerivative &&
        mode.radial_derivative isa FixedOutgoingURadialDerivative || throw(ArgumentError(
        "conjugate Hertz mode derivatives are not radial derivatives at fixed outgoing u"))
    bc.policy isa FiniteHertzRadialIVP && mode.policy isa FiniteHertzRadialIVP ||
        throw(ArgumentError("conjugate Hertz mode is not a finite-radius free IVP"))
    nothing
end

function _psi0(src::ReducedWeylMode{N,T}, r::T) where {N,T}
    value = _checked_complex(T, src.psi0(r), "linearized psi0")
    _finite(value) || throw(ArgumentError("linearized psi0 must be finite along the integration ray"))
    value
end

"""
    invert_hertz(mode, boundary) -> ConjugateHertzMode

Validate an explicitly linearized `psi0` carrier and the complete boundary key.  This installs
the positive mostly-plus IRG identity `(1/4) d_r^4 bar(Phi) = psi0`; it performs no metric
reconstruction or `m <-> -m` partner assembly.
"""
function invert_hertz(src::ReducedWeylMode{N,T,LinearizedPsi0,F},
        bc::HertzRadialBoundary{N,T}) where {N,T<:AbstractFloat,F}
    _validate_reduced_weyl_mode(src)
    _validate_hertz_boundary(bc)
    _same_hertz_key(src, bc) ||
        throw(ArgumentError("Weyl mode and Hertz boundary keys differ"))
    bc.normalization isa QuarterThorn4 ||
        throw(ArgumentError("Hertz inversion requires positive QuarterThorn4 normalization"))
    bc.gauge isa HertzIRG || throw(ArgumentError("Hertz inversion requires IRG"))
    bc.phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "Hertz boundary phase is not exp(-i*omega*u + i*m*phi)"))
    bc.radial_derivative isa FixedOutgoingURadialDerivative || throw(ArgumentError(
        "Hertz boundary data are not radial derivatives at fixed outgoing u"))
    bc.policy isa FiniteHertzRadialIVP || throw(ArgumentError(
        "Hertz boundary is not a finite-radius free IVP"))
    _psi0(src, bc.rb)
    ConjugateHertzMode{N,T,F}(src, bc, bc.normalization, bc.gauge, bc.phase,
        bc.radial_derivative, bc.policy)
end

function invert_hertz(src::ReducedWeylMode, bc::HertzRadialBoundary)
    src.content isa FullSecondOrderWeyl && throw(ArgumentError(
        "full second-order Weyl data are not invertible; supply the linearized psi0L = T0[hN]"))
    src.content isa LinearizedPsi0 ||
        throw(ArgumentError("unsupported Weyl content; Hertz inversion requires LinearizedPsi0"))
    throw(ArgumentError("Weyl mode and Hertz boundary order, precision, or keys differ"))
end

function _hertz_radius(mode::ConjugateHertzMode{N,T}, r) where {N,T}
    _validate_conjugate_hertz_mode(mode)
    rT = _checked_real(T, r, "Hertz evaluation radius")
    _finite(rT) && rT >= mode.boundary.rb && rT >= 2mode.source.M || throw(ArgumentError(
        "Hertz evaluation radius must be finite, outside the horizon, and r >= rb"))
    rT
end

function _taylor_term(value::Complex{T}, delta::T, order::Int) where {T}
    term = value
    for k in 1:order
        term *= delta / T(k)
    end
    term
end

function _fused_scale_add(scale::T, value::Complex{T}, offset::Complex{T}) where {T}
    Complex{T}(fma(scale, real(value), real(offset)),
               fma(scale, imag(value), imag(offset)))
end

function _cancellation_risk(product::T, offset::T, result::T) where {T}
    isfinite(product) || return true
    magnitude = max(abs(product), abs(offset))
    magnitude != zero(T) && abs(result) <= sqrt(eps(T)) * magnitude
end

function _horner_component(delta::T, denominator::T, value::T, offset::T) where {T}
    scale = delta / denominator
    result = fma(scale, value, offset)
    product = fma(scale, value, zero(T))
    if product == zero(T) && delta != zero(T) && value != zero(T)
        scaled_offset = denominator * offset
        if isfinite(scaled_offset)
            numerator = fma(delta, value, scaled_offset)
            if isfinite(numerator)
                alternate = numerator / denominator
                isfinite(alternate) && return alternate, true
            end
        end
    end
    result, _cancellation_risk(product, offset, result)
end

function _horner_step(delta::T, denominator::T, value::Complex{T},
        offset::Complex{T}) where {T}
    real_result, real_risk = _horner_component(
        delta, denominator, real(value), real(offset))
    imag_result, imag_risk = _horner_component(
        delta, denominator, imag(value), imag(offset))
    Complex{T}(real_result, imag_result), real_risk || imag_risk
end

function _homogeneous_taylor(data::NTuple{4,Complex{T}}, delta::T, n::Int) where {T}
    value = data[4]
    risk = false
    for j in 2:-1:n
        denominator = T(j + 1 - n)
        value, step_risk = _horner_step(delta, denominator, value, data[j + 1])
        risk |= step_risk
        _finite(value) || return value, true, risk
    end
    value, false, risk
end

_fallback_precision(::Type{Float16}) = 256
_fallback_precision(::Type{Float32}) = 1024
_fallback_precision(::Type{Float64}) = 8192
_fallback_precision(::Type{BigFloat}) = max(256, precision(BigFloat))
_fallback_precision(::Type{T}) where {T} = max(256, 4precision(T))

_big_complex(z) = Complex{BigFloat}(BigFloat(real(z)), BigFloat(imag(z)))

function _wide_homogeneous(data, delta_big::BigFloat, n::Int)
    homogeneous = zero(Complex{BigFloat})
    for j in n:3
        k = j - n
        homogeneous += _big_complex(data[j + 1]) * delta_big^k / BigFloat(factorial(k))
    end
    homogeneous
end

function _widened_sourced_composition(mode::ConjugateHertzMode{N,T}, data,
        delta::T, r::T, n::Int, atol::T, rtol::T) where {N,T}
    bits = max(_fallback_precision(T), precision(BigFloat))
    setprecision(BigFloat, bits) do
        power = 3 - n
        function wide_kernel(q)
            qT = T(q)
            term = _big_complex(_psi0(mode.source, qT))
            ray_delta = BigFloat(r - qT)
            for k in 1:power
                term *= ray_delta / BigFloat(k)
            end
            term
        end
        integral = mode.boundary.rb == r ? zero(Complex{BigFloat}) : first(quadgk(
            wide_kernel, BigFloat(mode.boundary.rb), BigFloat(r);
            atol=BigFloat(atol), rtol=BigFloat(rtol)))
        total = 4integral + _wide_homogeneous(data, BigFloat(delta), n)
        Complex{T}(T(real(total)), T(imag(total)))
    end
end

_quad_nonfinite_error(err) = err isa DomainError &&
    occursin("integrand produced", err.msg)

function _bar_hertz_derivative(mode::ConjugateHertzMode{N,T}, r::T, n::Int,
        atol::T, rtol::T) where {N,T}
    n == 4 && return T(4) * _psi0(mode.source, r)
    bc = mode.boundary
    delta = r - bc.rb
    data = (bc.b0, bc.b1, bc.b2, bc.b3)
    homogeneous, failed, risk = _homogeneous_taylor(data, delta, n)
    (failed || risk) &&
        return _widened_sourced_composition(mode, data, delta, r, n, atol, rtol)
    power = 3 - n
    integral = try
        _quad(q -> _taylor_term(_psi0(mode.source, q), r - q, power),
              bc.rb, r; atol, rtol)
    catch err
        _quad_nonfinite_error(err) || rethrow()
        return _widened_sourced_composition(mode, data, delta, r, n, atol, rtol)
    end
    _finite(integral) ||
        return _widened_sourced_composition(mode, data, delta, r, n, atol, rtol)
    fast = _fused_scale_add(T(4), integral, homogeneous)
    real_risk = _cancellation_risk(T(4) * real(integral), real(homogeneous), real(fast))
    imag_risk = _cancellation_risk(T(4) * imag(integral), imag(homogeneous), imag(fast))
    (real_risk || imag_risk || !_finite(fast)) &&
        return _widened_sourced_composition(mode, data, delta, r, n, atol, rtol)
    fast
end

"""
    bar_hertz_derivative(mode, r, n; atol=0, rtol=sqrt(eps(T)))

Evaluate `d_r^n bar(Phi)` for `0 <= n <= 4`.  For `n <= 3`,

`sum(b_j (r-rb)^(j-n)/(j-n)!, j=n:3) + 4 integral((r-q)^(3-n)/(3-n)! psi0(q), q=rb:r)`;

for `n=4`, the defining identity gives `4psi0(r)` exactly.
"""
function bar_hertz_derivative(mode::ConjugateHertzMode{N,T}, r, n::Integer;
        atol=zero(T), rtol=sqrt(eps(T))) where {N,T<:AbstractFloat}
    n isa Bool && throw(ArgumentError("Hertz derivative order must be an integer from 0 through 4"))
    0 <= n <= 4 ||
        throw(ArgumentError("Hertz derivative order must be an integer from 0 through 4"))
    rT = _hertz_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    _bar_hertz_derivative(mode, rT, Int(n), at, rt)
end

function bar_hertz_derivative(mode::ConjugateHertzMode{N,T}, r, n;
        atol=zero(T), rtol=sqrt(eps(T))) where {N,T<:AbstractFloat}
    throw(ArgumentError("Hertz derivative order must be an integer from 0 through 4"))
end

"Return the immutable radial jet `(bar(Phi), d bar(Phi), ..., d^4 bar(Phi))`."
function bar_hertz_jet(mode::ConjugateHertzMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    rT = _hertz_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    ntuple(i -> _bar_hertz_derivative(mode, rT, i - 1, at, rt), 5)
end

export WeylScalarContent, LinearizedPsi0, FullSecondOrderWeyl, QuarterThorn4
export FixedOutgoingURadialDerivative, FiniteHertzRadialIVP
export ReducedWeylMode, HertzRadialBoundary, ConjugateHertzMode
export invert_hertz, bar_hertz_derivative, bar_hertz_jet
