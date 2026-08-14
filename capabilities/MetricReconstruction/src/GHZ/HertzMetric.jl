"""
    ConjugateHertzRadialJet{N}(derivatives; M, l, m, omega, r,
        source_fingerprint, boundary_fingerprint, phase, radial_derivative, policy, ...)

Solver-independent radial data for the conjugate Hertz potential `bar(Phi)`.  The five entries
of `derivatives` are `(bar(Phi), d_r bar(Phi), ..., d_r^4 bar(Phi))`, where every derivative
is taken at fixed outgoing Schwarzschild time `u`.  These are not Boyer--Lindquist fixed-`t`
derivatives.  All convention, finite-IVP, and provenance metadata needed by the
metric-injection kernel are retained explicitly.
"""
struct ConjugateHertzRadialJet{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    source_fingerprint::UInt128
    boundary_fingerprint::UInt128
    r::T
    derivatives::NTuple{5,Complex{T}}
    normalization::QuarterThorn4
    gauge::HertzIRG
    tetrad::KinnersleyMinusPlus
    slicing::SchwarzschildOutgoing
    phase::ExpMinusIomegaUPlusImPhi
    radial_derivative::FixedOutgoingURadialDerivative
    policy::FiniteHertzRadialIVP
end

_hertz_fingerprint(value, label) = _checked_fingerprint(value, label)

function _hertz_jet_float_type(M, omega, r, derivatives)
    types = (typeof(float(M)), typeof(float(real(omega))), typeof(float(imag(omega))),
        typeof(float(r)),
        (typeof(float(part)) for z in derivatives for part in (real(z), imag(z)))...)
    T = promote_type(types...)
    T <: AbstractFloat ||
        throw(ArgumentError("Hertz jet data must promote to an AbstractFloat precision"))
    T
end

function ConjugateHertzRadialJet{N}(derivatives::NTuple{5}; M, l::Integer, m::Integer,
        omega=zero(M), r, source_fingerprint, boundary_fingerprint,
        normalization=QuarterThorn4(), gauge=HertzIRG(),
        tetrad=KinnersleyMinusPlus(), slicing=SchwarzschildOutgoing(),
        phase, radial_derivative, policy) where {N}
    N isa Integer && !(N isa Bool) && N > 0 ||
        throw(ArgumentError("perturbative order N must be a positive integer"))
    lI, mI = _checked_mode_index(l, "l"), _checked_mode_index(m, "m")
    lI >= 2 || throw(ArgumentError("a spin-2 Hertz jet requires l >= 2"))
    -lI <= mI <= lI || throw(ArgumentError("Hertz jet requires |m| <= l"))
    M isa Bool && throw(ArgumentError("M must be a non-Boolean real scalar"))
    r isa Bool && throw(ArgumentError("r must be a non-Boolean real scalar"))
    (real(omega) isa Bool || imag(omega) isa Bool) &&
        throw(ArgumentError("omega must have non-Boolean scalar components"))
    any(z -> real(z) isa Bool || imag(z) isa Bool, derivatives) &&
        throw(ArgumentError("Hertz radial derivatives must have non-Boolean scalar components"))
    normalization isa QuarterThorn4 || throw(ArgumentError(
        "Hertz jet requires the positive QuarterThorn4 normalization"))
    gauge isa HertzIRG || throw(ArgumentError("Hertz jet requires IRG"))
    tetrad isa KinnersleyMinusPlus || throw(ArgumentError(
        "Hertz jet requires the Kinnersley (-,+,+,+) tetrad"))
    slicing isa SchwarzschildOutgoing || throw(ArgumentError(
        "Hertz jet requires outgoing Schwarzschild slicing"))
    phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "Hertz jet phase must be exp(-i*omega*u + i*m*phi)"))
    radial_derivative isa FixedOutgoingURadialDerivative || throw(ArgumentError(
        "Hertz jet derivatives must be radial derivatives at fixed outgoing u"))
    policy isa FiniteHertzRadialIVP || throw(ArgumentError(
        "Hertz jet data require the finite radial IVP policy, not a retarded/physical claim"))

    T = _hertz_jet_float_type(M, omega, r, derivatives)
    MT, omegaT, rT = T(M), Complex{T}(omega), T(r)
    _finite(MT) && MT > zero(T) || throw(ArgumentError("M must be finite and positive"))
    _finite(omegaT) || throw(ArgumentError("omega must be finite"))
    _finite(rT) && rT >= T(2) * MT ||
        throw(ArgumentError("Hertz jet radius must be finite and satisfy r >= 2M"))
    values = ntuple(i -> Complex{T}(derivatives[i]), 5)
    all(_finite, values) || throw(ArgumentError("Hertz radial derivatives must be finite"))
    sfp = _hertz_fingerprint(source_fingerprint, "source_fingerprint")
    bfp = _hertz_fingerprint(boundary_fingerprint, "boundary_fingerprint")

    ConjugateHertzRadialJet{N,T}(MT, lI, mI, omegaT, sfp, bfp, rT, values,
        normalization, gauge, tetrad, slicing, phase, radial_derivative, policy)
end

function ConjugateHertzRadialJet{N}(derivatives; kwargs...) where {N}
    throw(ArgumentError("a conjugate Hertz radial jet requires exactly five derivatives"))
end

perturbative_order(::ConjugateHertzRadialJet{N}) where {N} = N

"""
    conjugate_hertz_radial_jet(mode, r; atol=0, rtol=sqrt(eps(T)))

Materialize a validated `ConjugateHertzMode` as pure radial data.  This is the only adapter from
the inversion solver into the metric kernel; the kernel itself has no dependency on a radial
solver object.  The inversion is a finite-radius free IVP, and its radial derivatives are at
fixed outgoing Schwarzschild `u`; the adapter propagates the source and boundary assertions
without inventing retarded or physical boundary data.
"""
function conjugate_hertz_radial_jet(mode::ConjugateHertzMode{N,T}, r;
        atol=zero(T), rtol=sqrt(eps(T))) where {N,T}
    _validate_conjugate_hertz_mode(mode)
    values = bar_hertz_jet(mode, r; atol, rtol)
    src, bc = mode.source, mode.boundary
    ConjugateHertzRadialJet{N}(values; M=src.M, l=src.l, m=src.m, omega=src.omega,
        r=T(r), source_fingerprint=src.source_fingerprint,
        boundary_fingerprint=bc.boundary_fingerprint, normalization=mode.normalization,
        gauge=mode.gauge, tetrad=src.tetrad, slicing=src.slicing,
        phase=mode.phase, radial_derivative=mode.radial_derivative, policy=mode.policy)
end

"""
    PhysicalHertzPair(mode, partner; pair_fingerprint)

An explicit, caller-asserted pairing of a `bar(Phi)_(l,m,omega)` jet with the independent partner
jet whose key is `(l,-m,-conj(omega))`.  Each constituent retains its own content and boundary
fingerprints; the independent, nonzero `pair_fingerprint` records that the producer deliberately
formed this keyed pair.  The carrier validates keys and metadata, not functional reality or
physical boundary data.  Supplying the partner explicitly is mandatory, including for the
degenerate static `m=0` key; no same-mode conjugation is performed implicitly.
"""
struct PhysicalHertzPair{N,T<:AbstractFloat}
    mode::ConjugateHertzRadialJet{N,T}
    partner::ConjugateHertzRadialJet{N,T}
    pair_fingerprint::UInt128

    function PhysicalHertzPair{N,T}(mode::ConjugateHertzRadialJet{N,T},
            partner::ConjugateHertzRadialJet{N,T}, pair_fingerprint::UInt128,
            ::Val{:validated}) where {N,T<:AbstractFloat}
        new{N,T}(mode, partner, pair_fingerprint)
    end
end

function _validate_hertz_jet(jet::ConjugateHertzRadialJet{N,T}) where {N,T}
    N isa Integer && !(N isa Bool) && N > 0 ||
        throw(ArgumentError("Hertz jet perturbative order must be positive"))
    _finite(jet.M) && jet.M > zero(T) ||
        throw(ArgumentError("Hertz jet M must be finite and positive"))
    jet.l >= 2 && -jet.l <= jet.m <= jet.l ||
        throw(ArgumentError("invalid spin-2 Hertz mode key"))
    _finite(jet.omega) || throw(ArgumentError("Hertz jet omega must be finite"))
    jet.source_fingerprint != 0 ||
        throw(ArgumentError("Hertz jet source fingerprint must be nonzero"))
    jet.boundary_fingerprint != 0 ||
        throw(ArgumentError("Hertz jet boundary fingerprint must be nonzero"))
    _finite(jet.r) && jet.r >= T(2) * jet.M ||
        throw(ArgumentError("Hertz jet radius must be finite and satisfy r >= 2M"))
    all(_finite, jet.derivatives) ||
        throw(ArgumentError("Hertz radial derivatives must be finite"))
    jet.normalization isa QuarterThorn4 ||
        throw(ArgumentError("Hertz jet normalization is not QuarterThorn4"))
    jet.gauge isa HertzIRG || throw(ArgumentError("Hertz jet is not IRG"))
    jet.tetrad isa KinnersleyMinusPlus ||
        throw(ArgumentError("Hertz jet tetrad is not Kinnersley (-,+,+,+)"))
    jet.slicing isa SchwarzschildOutgoing ||
        throw(ArgumentError("Hertz jet is not on outgoing Schwarzschild slices"))
    jet.phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "Hertz jet phase is not exp(-i*omega*u + i*m*phi)"))
    jet.radial_derivative isa FixedOutgoingURadialDerivative || throw(ArgumentError(
        "Hertz jet derivatives are not radial derivatives at fixed outgoing u"))
    jet.policy isa FiniteHertzRadialIVP || throw(ArgumentError(
        "Hertz jet is not finite-radius free-IVP data"))
    nothing
end

function PhysicalHertzPair(mode::ConjugateHertzRadialJet{N,T},
        partner::ConjugateHertzRadialJet{N,T}; pair_fingerprint=nothing) where {N,T<:AbstractFloat}
    pair_fp = _hertz_fingerprint(pair_fingerprint, "pair_fingerprint")
    pair = PhysicalHertzPair{N,T}(mode, partner, pair_fp, Val(:validated))
    _validate_physical_hertz_pair(pair)
    pair
end

function PhysicalHertzPair(mode::ConjugateHertzRadialJet,
        partner::ConjugateHertzRadialJet; pair_fingerprint=nothing)
    _hertz_fingerprint(pair_fingerprint, "pair_fingerprint")
    throw(ArgumentError("Hertz partners have different perturbative orders or precisions"))
end

PhysicalHertzPair(::ConjugateHertzRadialJet; pair_fingerprint=nothing) =
    throw(ArgumentError("an explicit (l,-m,-conj(omega)) Hertz partner is required"))
PhysicalHertzPair(::ConjugateHertzRadialJet, ::Nothing; pair_fingerprint=nothing) =
    throw(ArgumentError("an explicit (l,-m,-conj(omega)) Hertz partner is required"))

perturbative_order(::PhysicalHertzPair{N}) where {N} = N

function _validate_physical_hertz_pair(pair::PhysicalHertzPair)
    mode, partner = pair.mode, pair.partner
    _validate_hertz_jet(mode)
    _validate_hertz_jet(partner)
    pair.pair_fingerprint != 0 ||
        throw(ArgumentError("pair_fingerprint must be nonzero"))
    mode.M == partner.M || throw(ArgumentError("Hertz partners have different M"))
    mode.l == partner.l || throw(ArgumentError("Hertz partners have different l"))
    partner.m == -mode.m || throw(ArgumentError(
        "Hertz partner must have azimuthal key -m"))
    partner.omega == -conj(mode.omega) || throw(ArgumentError(
        "Hertz partner must have frequency key -conj(omega)"))
    mode.r == partner.r || throw(ArgumentError("Hertz partners must be evaluated at the same r"))
    typeof(mode.normalization) === typeof(partner.normalization) ||
        throw(ArgumentError("Hertz partner normalizations differ"))
    typeof(mode.gauge) === typeof(partner.gauge) ||
        throw(ArgumentError("Hertz partner gauges differ"))
    typeof(mode.tetrad) === typeof(partner.tetrad) ||
        throw(ArgumentError("Hertz partner tetrads differ"))
    typeof(mode.slicing) === typeof(partner.slicing) ||
        throw(ArgumentError("Hertz partner slicings differ"))
    typeof(mode.phase) === typeof(partner.phase) ||
        throw(ArgumentError("Hertz partner phase conventions differ"))
    typeof(mode.radial_derivative) === typeof(partner.radial_derivative) ||
        throw(ArgumentError("Hertz partner radial-derivative conventions differ"))
    typeof(mode.policy) === typeof(partner.policy) ||
        throw(ArgumentError("Hertz partner radial-IVP policies differ"))
    nothing
end

function _physical_hertz_jet_unchecked(pair::PhysicalHertzPair{N,T}) where {N,T}
    phase = isodd(pair.mode.m) ? -one(T) : one(T)
    ntuple(i -> phase * conj(pair.partner.derivatives[i]), 5)
end

"Return `(Phi, d_r Phi, ..., d_r^4 Phi)` using only the explicitly supplied partner jet."
function physical_hertz_jet(pair::PhysicalHertzPair{N,T}) where {N,T}
    _validate_physical_hertz_pair(pair)
    _physical_hertz_jet_unchecked(pair)
end

"The five nonzero tetrad-mode coefficients of an outgoing-Schwarzschild IRG Hertz metric."
struct IRGHertzMetricMode{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    mode_source_fingerprint::UInt128
    partner_source_fingerprint::UInt128
    mode_boundary_fingerprint::UInt128
    partner_boundary_fingerprint::UInt128
    pair_fingerprint::UInt128
    r::T
    k_nn::Complex{T}
    k_nbar::Complex{T}
    k_nm::Complex{T}
    k_barbar::Complex{T}
    k_mm::Complex{T}
    normalization::QuarterThorn4
    gauge::HertzIRG
    tetrad::KinnersleyMinusPlus
    slicing::SchwarzschildOutgoing
    phase::ExpMinusIomegaUPlusImPhi
    radial_derivative::FixedOutgoingURadialDerivative
    policy::FiniteHertzRadialIVP
end

perturbative_order(::IRGHertzMetricMode{N}) where {N} = N

function _hertz_couplings(::Type{T}, l::Int) where {T<:AbstractFloat}
    # Convert before shifting, so `l+1` and `l+2` cannot overflow Int.  Reject precisions too
    # coarse to distinguish neighboring harmonic labels rather than silently using wrong factors.
    ell = T(l)
    _finite(ell) || throw(ArgumentError("l is not finite when represented in $T"))
    oneT = one(T)
    ellm1, ellp1, ellp2 = ell - oneT, ell + oneT, ell + T(2)
    ellp1 > ell && ellp2 > ellp1 || throw(ArgumentError(
        "precision $T cannot resolve the angular labels l, l+1, and l+2 for l=$l"))
    lambda2 = sqrt(ellm1) * sqrt(ell) * sqrt(ellp1) * sqrt(ellp2)
    mu = sqrt(ellm1) * sqrt(ellp2)
    d = mu / sqrt(T(2))
    all(_finite, (lambda2, mu, d)) ||
        throw(ArgumentError("Hertz angular couplings are not finite in $T for l=$l"))
    lambda2, mu, d
end

function _real_complex_muladd(a::T, z::Complex{T}, w::Complex{T}) where {T<:AbstractFloat}
    Complex{T}(muladd(a, real(z), real(w)), muladd(a, imag(z), imag(w)))
end

_real_complex_scale(a::T, z::Complex{T}) where {T<:AbstractFloat} =
    Complex{T}(a * real(z), a * imag(z))

function _hertz_scaled_component_lost(original::Complex{T},
        scaled::Complex{T}) where {T<:AbstractFloat}
    (!iszero(real(original)) && iszero(real(scaled))) ||
        (!iszero(imag(original)) && iszero(imag(scaled)))
end


_hertz_has_subnormal(z::Complex{T}) where {T<:Union{Float16,Float32,Float64}} =
    (!iszero(real(z)) && issubnormal(real(z))) ||
    (!iszero(imag(z)) && issubnormal(imag(z)))
_hertz_has_subnormal(::Complex{T}) where {T<:AbstractFloat} = false

_hertz_metric_fallback_precision(::Type{Float16}) = 256
_hertz_metric_fallback_precision(::Type{Float32}) = 1024
_hertz_metric_fallback_precision(::Type{Float64}) = 8192
_hertz_metric_fallback_precision(::Type{BigFloat}) = max(256, 2precision(BigFloat))
_hertz_metric_fallback_precision(::Type{T}) where {T} = max(256, 4precision(T))

_hertz_big_complex(z) = Complex{BigFloat}(BigFloat(real(z)), BigFloat(imag(z)))

function _hertz_round_wide(::Type{T}, z::Complex{BigFloat},
        target_precision::Int) where {T<:AbstractFloat}
    T === BigFloat && return setprecision(BigFloat, target_precision) do
        Complex{BigFloat}(BigFloat(real(z)), BigFloat(imag(z)))
    end
    Complex{T}(T(real(z)), T(imag(z)))
end

function _hertz_metric_coefficients_wide(pair::PhysicalHertzPair{N,T}) where
        {N,T<:AbstractFloat}
    target_precision = T === BigFloat ? precision(BigFloat) : precision(T)
    bits = _hertz_metric_fallback_precision(T)
    setprecision(BigFloat, bits) do
        mode = pair.mode
        P = _physical_hertz_jet_unchecked(pair)
        B = mode.derivatives
        r = BigFloat(mode.r)
        ell = BigFloat(mode.l)
        lambda2 = sqrt((ell - 1) * ell * (ell + 1) * (ell + 2))
        c = sqrt((ell - 1) * (ell + 2)) / (2sqrt(BigFloat(2)))
        PB, BB = _hertz_big_complex.(P), _hertz_big_complex.(B)
        values = (
            -lambda2 * (PB[1] + BB[1]) / (4r^2),
            c * (PB[2] / r - 2PB[1] / r^2),
            -c * (BB[2] / r - 2BB[1] / r^2),
            -PB[3] / 2 + PB[2] / r,
            -BB[3] / 2 + BB[2] / r,
        )
        ntuple(i -> _hertz_round_wide(T, values[i], target_precision), 5)
    end
end

function _hertz_metric_unchecked(pair::PhysicalHertzPair{N,T}) where {N,T<:AbstractFloat}
    mode = pair.mode
    lambda2, _, d = _hertz_couplings(T, mode.l)
    P = _physical_hertz_jet_unchecked(pair)
    B = mode.derivatives
    r = mode.r
    two = T(2)

    P0_over_r2 = (P[1] / r) / r
    B0_over_r2 = (B[1] / r) / r
    P1_over_r = P[2] / r
    B1_over_r = B[2] / r
    c = d / two
    sum0 = _real_complex_muladd(one(T), P0_over_r2, B0_over_r2)
    k_nn = -_real_complex_scale(lambda2 / T(4), sum0)
    scaled_P0 = _real_complex_scale(d, P[1])
    scaled_P1 = _real_complex_scale(c, P[2])
    scaled_B0 = _real_complex_scale(d, B[1])
    scaled_B1 = _real_complex_scale(c, B[2])
    P0_term = (scaled_P0 / r) / r
    P1_term = scaled_P1 / r
    B0_term = (scaled_B0 / r) / r
    B1_term = scaled_B1 / r
    k_nbar = _real_complex_muladd(-one(T), P0_term, P1_term)
    k_nm = _real_complex_muladd(one(T), B0_term, -B1_term)

    if r >= one(T)
        k_barbar = _real_complex_muladd(-inv(two), P[3], P[2] / r)
        k_mm = _real_complex_muladd(-inv(two), B[3], B[2] / r)
    else
        half_r = r / two
        k_barbar = _real_complex_muladd(-half_r, P[3], P[2]) / r
        k_mm = _real_complex_muladd(-half_r, B[3], B[2]) / r
    end
    coefficients = (k_nn, k_nbar, k_nm, k_barbar, k_mm)
    division_underflow =
        _hertz_scaled_component_lost(P[1], P0_over_r2) ||
        _hertz_scaled_component_lost(B[1], B0_over_r2) ||
        _hertz_scaled_component_lost(P[2], P1_over_r) ||
        _hertz_scaled_component_lost(B[2], B1_over_r) ||
        _hertz_scaled_component_lost(P[1], scaled_P0) ||
        _hertz_scaled_component_lost(P[2], scaled_P1) ||
        _hertz_scaled_component_lost(B[1], scaled_B0) ||
        _hertz_scaled_component_lost(B[2], scaled_B1) ||
        _hertz_scaled_component_lost(scaled_P0, P0_term) ||
        _hertz_scaled_component_lost(scaled_P1, P1_term) ||
        _hertz_scaled_component_lost(scaled_B0, B0_term) ||
        _hertz_scaled_component_lost(scaled_B1, B1_term)
    subnormal_risk = any(_hertz_has_subnormal,
        (P0_over_r2, B0_over_r2, P1_over_r, B1_over_r,
         scaled_P0, scaled_P1, scaled_B0, scaled_B1,
         P0_term, P1_term, B0_term, B1_term)) ||
        any(_hertz_has_subnormal, coefficients)
    zero_cancellation_risk =
        (iszero(k_nn) && (!iszero(P[1]) || !iszero(B[1]))) ||
        (iszero(k_nbar) && (!iszero(P[1]) || !iszero(P[2]))) ||
        (iszero(k_nm) && (!iszero(B[1]) || !iszero(B[2]))) ||
        (iszero(k_barbar) && (!iszero(P[2]) || !iszero(P[3]))) ||
        (iszero(k_mm) && (!iszero(B[2]) || !iszero(B[3])))
    if !all(_finite, coefficients) || division_underflow || subnormal_risk ||
            zero_cancellation_risk
        coefficients = _hertz_metric_coefficients_wide(pair)
    end
    all(_finite, coefficients) ||
        throw(ArgumentError("Hertz metric coefficients are not finite"))

    IRGHertzMetricMode{N,T}(mode.M, mode.l, mode.m, mode.omega,
        mode.source_fingerprint, pair.partner.source_fingerprint,
        mode.boundary_fingerprint, pair.partner.boundary_fingerprint,
        pair.pair_fingerprint, r, coefficients...,
        mode.normalization, mode.gauge, mode.tetrad, mode.slicing, mode.phase,
        mode.radial_derivative, mode.policy)
end

"""
    hertz_metric(pair) -> IRGHertzMetricMode

Inject a caller-asserted explicit Hertz partner pair into the pure Schwarzschild `(-,+,+,+)` IRG
tetrad kernel.  Pair validation establishes matching keys and metadata, not functional reality
or physical boundary data.  The supplied radial derivatives are at fixed outgoing `u`, so the
Schwarzschild `a=0` formulas have no explicit `i*omega` terms.  They must not be interpreted as
fixed-`t` Boyer--Lindquist radial derivatives, whose conversion would introduce `-i*omega/f`.
The finite-IVP policy makes no retarded or physical-data claim.
"""
function hertz_metric(pair::PhysicalHertzPair{N,T}) where {N,T<:AbstractFloat}
    _validate_physical_hertz_pair(pair)
    _hertz_metric_unchecked(pair)
end

"""
    static_hertz_metric(pair) -> IRGHertzMetricMode

Checked `omega=0` compatibility wrapper around [`hertz_metric`](@ref).  It uses the same
fixed-outgoing-`u` Schwarzschild kernel and retains the explicit finite-IVP metadata.
"""
function static_hertz_metric(pair::PhysicalHertzPair{N,T}) where {N,T<:AbstractFloat}
    _validate_physical_hertz_pair(pair)
    pair.mode.omega == zero(Complex{T}) || throw(ArgumentError(
        "static_hertz_metric requires omega = 0; use hertz_metric for finite frequency"))
    _hertz_metric_unchecked(pair)
end

export hertz_metric
