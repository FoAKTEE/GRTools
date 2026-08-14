"An explicit physical `(m, -m)` pairing of solved GHZ middle-corrector modes."
struct PhysicalCorrectorMiddlePair{N,T<:AbstractFloat,
        P<:IRGCorrectorMiddleMode{N,T},Q<:IRGCorrectorMiddleMode{N,T}}
    primary::P
    partner::Q
    pair_fingerprint::UInt128

    function PhysicalCorrectorMiddlePair{N,T,P,Q}(primary::P, partner::Q,
            pair_fingerprint::UInt128, ::Val{:validated}) where {N,T<:AbstractFloat,
            P<:IRGCorrectorMiddleMode{N,T},Q<:IRGCorrectorMiddleMode{N,T}}
        new{N,T,P,Q}(primary, partner, pair_fingerprint)
    end
end

function _validate_middle_source(src::TetradSourceMode{N,T}) where {N,T}
    N isa Integer && !(N isa Bool) && N > 0 ||
        throw(ArgumentError("middle source perturbative order must be positive"))
    _finite(src.M) && src.M > zero(T) ||
        throw(ArgumentError("middle source M must be finite and positive"))
    _require_middle_spin_domain(src.l)
    -src.l <= src.m <= src.l || throw(ArgumentError("invalid middle source mode key"))
    _finite(src.omega) || throw(ArgumentError("middle source omega must be finite"))
    src.source_fingerprint != 0 ||
        throw(ArgumentError("middle source fingerprint must be nonzero"))
    src.normalization isa EinsteinRHS ||
        throw(ArgumentError("middle source normalization is not EinsteinRHS"))
    src.tetrad isa KinnersleyMinusPlus ||
        throw(ArgumentError("middle source tetrad is not Kinnersley (-,+,+,+)"))
    src.slicing isa SchwarzschildOutgoing ||
        throw(ArgumentError("middle source is not on outgoing Schwarzschild slices"))
    src.phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "middle source phase is not exp(-i*omega*u + i*m*phi)"))
    nothing
end

function _validate_corrector_boundary(src::TetradSourceMode{N,T}, bc::InnerRayIVP{N,T},
        label) where {N,T}
    _same_mode(src, bc) || throw(ArgumentError("$label mode key differs from its source"))
    _finite(bc.rb) && bc.rb >= T(2) * src.M ||
        throw(ArgumentError("$label rb must be finite and satisfy rb >= 2M"))
    _finite(bc.x0) && _finite(bc.dx0) ||
        throw(ArgumentError("$label data must be finite"))
    nothing
end

function _validate_corrector_boundary(src::TetradSourceMode{N,T}, bc::MiddleRayIVP{N,T},
        label) where {N,T}
    _require_middle_spin_domain(bc.l)
    _same_mode(src, bc) || throw(ArgumentError("$label mode key differs from its source"))
    _finite(bc.rb) && bc.rb >= T(2) * src.M ||
        throw(ArgumentError("$label rb must be finite and satisfy rb >= 2M"))
    _finite(bc.x0) && _finite(bc.dx0) ||
        throw(ArgumentError("$label data must be finite"))
    nothing
end

function _validate_corrector_middle_mode(mode::IRGCorrectorMiddleMode)
    src = mode.head.source
    _validate_middle_source(src)
    _validate_corrector_boundary(src, mode.head.boundary, "corrector-head boundary")
    _validate_corrector_boundary(src, mode.boundary, "corrector-middle boundary")
    mode.head.boundary.rb == mode.boundary.rb || throw(ArgumentError(
        "corrector head and middle must share the same inner-ray boundary"))
    mode.head.gauge isa GHZIRG || throw(ArgumentError("corrector head is not GHZ IRG"))
    mode.gauge isa GHZIRG || throw(ArgumentError("corrector middle is not GHZ IRG"))
    sll = _checked_complex(typeof(src.M), src.Sll(mode.head.boundary.rb), "S_ll")
    _finite(sll) || throw(ArgumentError("S_ll must be finite at the inner boundary"))
    slm = _checked_complex(typeof(src.M), src.Slm(mode.boundary.rb), "S_lm")
    _finite(slm) || throw(ArgumentError("S_lm must be finite at the inner boundary"))
    nothing
end

_corrector_order_precision(::IRGCorrectorMiddleMode{N,T}) where {N,T} = (N, T)

function PhysicalCorrectorMiddlePair(primary::IRGCorrectorMiddleMode,
        partner::IRGCorrectorMiddleMode; pair_fingerprint=nothing)
    pair_fp = _checked_fingerprint(pair_fingerprint, "pair_fingerprint")
    N, T = _corrector_order_precision(primary)
    partner_N, partner_T = _corrector_order_precision(partner)
    N === partner_N && T === partner_T || throw(ArgumentError(
        "corrector partners have different perturbative orders or precisions"))
    P, Q = typeof(primary), typeof(partner)
    pair = PhysicalCorrectorMiddlePair{N,T,P,Q}(
        primary, partner, pair_fp, Val(:validated))
    _validate_physical_corrector_middle_pair(pair)
    pair
end

PhysicalCorrectorMiddlePair(::IRGCorrectorMiddleMode; pair_fingerprint=nothing) =
    throw(ArgumentError("an explicit (l,-m,-conj(omega)) middle partner is required"))
PhysicalCorrectorMiddlePair(::IRGCorrectorMiddleMode, ::Nothing;
        pair_fingerprint=nothing) =
    throw(ArgumentError("an explicit (l,-m,-conj(omega)) middle partner is required"))

perturbative_order(::PhysicalCorrectorMiddlePair{N}) where {N} = N

function _validate_physical_corrector_middle_pair(pair::PhysicalCorrectorMiddlePair)
    primary, partner = pair.primary, pair.partner
    _validate_corrector_middle_mode(primary)
    _validate_corrector_middle_mode(partner)
    pair.pair_fingerprint != 0 || throw(ArgumentError("pair_fingerprint must be nonzero"))

    psrc, qsrc = primary.head.source, partner.head.source
    psrc.M == qsrc.M || throw(ArgumentError("corrector partners have different M"))
    psrc.l == qsrc.l || throw(ArgumentError("corrector partners have different l"))
    qsrc.m == -psrc.m ||
        throw(ArgumentError("corrector partner must have azimuthal key -m"))
    qsrc.omega == -conj(psrc.omega) ||
        throw(ArgumentError("corrector partner must have frequency key -conj(omega)"))
    primary.boundary.rb == partner.boundary.rb ||
        throw(ArgumentError("corrector partners have different rb"))
    typeof(psrc.normalization) === typeof(qsrc.normalization) ||
        throw(ArgumentError("corrector partner normalizations differ"))
    typeof(psrc.tetrad) === typeof(qsrc.tetrad) ||
        throw(ArgumentError("corrector partner tetrads differ"))
    typeof(psrc.slicing) === typeof(qsrc.slicing) ||
        throw(ArgumentError("corrector partner slicings differ"))
    typeof(psrc.phase) === typeof(qsrc.phase) ||
        throw(ArgumentError("corrector partner phases differ"))
    typeof(primary.head.gauge) === typeof(partner.head.gauge) ||
        throw(ArgumentError("corrector-head partner gauges differ"))
    typeof(primary.gauge) === typeof(partner.gauge) ||
        throw(ArgumentError("corrector-middle partner gauges differ"))
    nothing
end

_barred_middle_parity(::Type{T}, m::Int) where {T} = isodd(m) ? one(T) : -one(T)

function _x_nmb_unchecked(pair::PhysicalCorrectorMiddlePair{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    parity = _barred_middle_parity(T, pair.primary.head.source.m)
    parity * conj(x_nm(pair.partner, r; atol, rtol))
end

function _dr_x_nmb_unchecked(pair::PhysicalCorrectorMiddlePair{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    parity = _barred_middle_parity(T, pair.primary.head.source.m)
    parity * conj(dr_x_nm(pair.partner, r; atol, rtol))
end

"Evaluate barred `x_nm` using only the explicitly supplied partner solution."
function x_nmb(pair::PhysicalCorrectorMiddlePair{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    _validate_physical_corrector_middle_pair(pair)
    _x_nmb_unchecked(pair, r; atol, rtol)
end

"Evaluate the radial derivative of barred `x_nm` using only the explicit partner."
function dr_x_nmb(pair::PhysicalCorrectorMiddlePair{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    _validate_physical_corrector_middle_pair(pair)
    _dr_x_nmb_unchecked(pair, r; atol, rtol)
end
