"""
    TwoRealT0Dagger

Normalization marker for a caller-supplied adjusted source in the project's mostly-plus
`(-,+,+,+)` convention, including the explicit `2 Re(T_0^dagger[...])` normalization.  This
marker records a convention; this bounded kernel does not construct `J - E*x` and does not solve
the separate adjoint Teukolsky equation `O^dagger Phi = eta`.
"""
struct TwoRealT0Dagger end

"""
    ReducedAdjointSourceMode(pair, Smbmb; adjusted_source_fingerprint, ...)

A supplied adjusted `S_mbar_mbar` radial mode keyed to a revalidated, caller-asserted GHZ
reality-partner carrier.  The carrier proves mode-key consistency, not functional reality or
physical boundary data.
The independent adjusted-source fingerprint and a complete snapshot of the pair provenance keep
this carrier distinct from both of its upstream Einstein-RHS modes.  The source is reduced input:
no corrector tail, `J - E*x`, or adjoint-Hertz solve is hidden in this type.
"""
struct ReducedAdjointSourceMode{N,T<:AbstractFloat,
        P<:PhysicalCorrectorMiddlePair{N,T},F}
    Smbmb::F
    pair::P
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    primary_source_fingerprint::UInt128
    partner_source_fingerprint::UInt128
    pair_fingerprint::UInt128
    primary_head_x0::Complex{T}
    primary_head_dx0::Complex{T}
    primary_middle_x0::Complex{T}
    primary_middle_dx0::Complex{T}
    partner_head_x0::Complex{T}
    partner_head_dx0::Complex{T}
    partner_middle_x0::Complex{T}
    partner_middle_dx0::Complex{T}
    pair_rb::T
    adjusted_source_fingerprint::UInt128
    normalization::TwoRealT0Dagger
    gauge::GHZIRG
    tetrad::KinnersleyMinusPlus
    slicing::SchwarzschildOutgoing
    phase::ExpMinusIomegaUPlusImPhi
end

perturbative_order(::ReducedAdjointSourceMode{N}) where {N} = N
_adjoint_source_order_precision(::ReducedAdjointSourceMode{N,T}) where {N,T} = (N, T)

function _require_eta_spin_domain(l::Integer)
    lI = _checked_mode_index(l, "eta mode l")
    lI >= 2 || throw(ArgumentError("reduced adjoint eta modes require l >= 2"))
    lI
end

function _adjoint_pair_snapshot(pair::PhysicalCorrectorMiddlePair)
    primary, partner = pair.primary, pair.partner
    psrc, qsrc = primary.head.source, partner.head.source
    phb, pmb = primary.head.boundary, primary.boundary
    qhb, qmb = partner.head.boundary, partner.boundary
    (psrc.M, psrc.l, psrc.m, psrc.omega, psrc.source_fingerprint,
        qsrc.source_fingerprint, pair.pair_fingerprint,
        phb.x0, phb.dx0, pmb.x0, pmb.dx0,
        qhb.x0, qhb.dx0, qmb.x0, qmb.dx0, pmb.rb)
end

function _adjusted_Smbmb(src::ReducedAdjointSourceMode{N,T}, r::T) where {N,T}
    value = _checked_complex(T, src.Smbmb(r), "adjusted S_mbar_mbar")
    _finite(value) || throw(ArgumentError(
        "adjusted S_mbar_mbar must be finite along the integration ray"))
    value
end

function _validate_reduced_adjoint_source(src::ReducedAdjointSourceMode{N,T}) where {N,T}
    _validate_physical_corrector_middle_pair(src.pair)
    _require_eta_spin_domain(src.l)
    snapshot = _adjoint_pair_snapshot(src.pair)
    carried = (src.M, src.l, src.m, src.omega, src.primary_source_fingerprint,
        src.partner_source_fingerprint, src.pair_fingerprint,
        src.primary_head_x0, src.primary_head_dx0,
        src.primary_middle_x0, src.primary_middle_dx0,
        src.partner_head_x0, src.partner_head_dx0,
        src.partner_middle_x0, src.partner_middle_dx0, src.pair_rb)
    carried == snapshot || throw(ArgumentError(
        "reduced adjoint source provenance differs from its caller-asserted corrector pair"))
    src.adjusted_source_fingerprint != 0 || throw(ArgumentError(
        "adjusted_source_fingerprint must be nonzero"))
    src.adjusted_source_fingerprint ∉ (src.primary_source_fingerprint,
        src.partner_source_fingerprint, src.pair_fingerprint) || throw(ArgumentError(
        "adjusted_source_fingerprint must be independent of upstream fingerprints"))
    src.normalization isa TwoRealT0Dagger || throw(ArgumentError(
        "adjusted source is not in the mostly-plus TwoRealT0Dagger convention"))
    src.gauge isa GHZIRG || throw(ArgumentError("adjusted source is not keyed to GHZ IRG"))
    src.tetrad isa KinnersleyMinusPlus || throw(ArgumentError(
        "adjusted source is not in the Kinnersley (-,+,+,+) tetrad"))
    src.slicing isa SchwarzschildOutgoing || throw(ArgumentError(
        "adjusted source is not on outgoing Schwarzschild slices"))
    src.phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "adjusted source phase is not exp(-i*omega*u + i*m*phi)"))
    _finite(src.pair_rb) && src.pair_rb >= T(2) * src.M || throw(ArgumentError(
        "adjusted source pair boundary must be finite and satisfy rb >= 2M"))
    _adjusted_Smbmb(src, src.pair_rb)
    nothing
end

function ReducedAdjointSourceMode(pair::PhysicalCorrectorMiddlePair{N,T}, Smbmb;
        adjusted_source_fingerprint, normalization=TwoRealT0Dagger(), gauge=GHZIRG(),
        tetrad=KinnersleyMinusPlus(), slicing=SchwarzschildOutgoing(),
        phase=ExpMinusIomegaUPlusImPhi()) where {N,T}
    _validate_physical_corrector_middle_pair(pair)
    _require_eta_spin_domain(pair.primary.head.source.l)
    normalization isa TwoRealT0Dagger || throw(ArgumentError(
        "adjusted source requires the mostly-plus TwoRealT0Dagger convention"))
    gauge isa GHZIRG || throw(ArgumentError("adjusted source requires GHZ IRG"))
    tetrad isa KinnersleyMinusPlus || throw(ArgumentError(
        "adjusted source requires the Kinnersley (-,+,+,+) tetrad"))
    slicing isa SchwarzschildOutgoing || throw(ArgumentError(
        "adjusted source requires outgoing Schwarzschild slicing"))
    phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "adjusted source phase must be exp(-i*omega*u + i*m*phi)"))
    fp = _checked_fingerprint(adjusted_source_fingerprint,
        "adjusted_source_fingerprint")
    snapshot = _adjoint_pair_snapshot(pair)
    P, F = typeof(pair), typeof(Smbmb)
    src = ReducedAdjointSourceMode{N,T,P,F}(Smbmb, pair, snapshot..., fp,
        normalization, gauge, tetrad, slicing, phase)
    _validate_reduced_adjoint_source(src)
    src
end

"Finite-ray free data for the reduced eta equation; neither retarded nor physical."
struct EtaRayIVP{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    primary_source_fingerprint::UInt128
    partner_source_fingerprint::UInt128
    pair_fingerprint::UInt128
    adjusted_source_fingerprint::UInt128
    rb::T
    eta0::Complex{T}
    deta0::Complex{T}
    boundary_fingerprint::UInt128
    policy::FiniteRayFreeIVP
end

perturbative_order(::EtaRayIVP{N}) where {N} = N
_eta_boundary_order_precision(::EtaRayIVP{N,T}) where {N,T} = (N, T)

function EtaRayIVP(src::ReducedAdjointSourceMode{N,T}; rb=src.pair_rb, eta,
        dr_eta, boundary_fingerprint, policy) where {N,T}
    _validate_reduced_adjoint_source(src)
    policy isa FiniteRayFreeIVP || throw(ArgumentError(
        "eta data require the finite-ray free-IVP policy, not a retarded/physical claim"))
    rbT = _checked_real(T, rb, "eta boundary rb")
    eta0 = _checked_complex(T, eta, "eta boundary value")
    deta0 = _checked_complex(T, dr_eta, "dr_eta boundary value")
    _finite(rbT) && rbT >= T(2) * src.M || throw(ArgumentError(
        "eta boundary rb must be finite and satisfy rb >= 2M"))
    rbT == src.pair_rb || throw(ArgumentError(
        "eta and corrector-middle data must share the same finite-ray boundary"))
    _finite(eta0) && _finite(deta0) ||
        throw(ArgumentError("eta boundary data must be finite"))
    boundary_fp = _checked_fingerprint(boundary_fingerprint, "boundary_fingerprint")
    EtaRayIVP{N,T}(src.M, src.l, src.m, src.omega,
        src.primary_source_fingerprint, src.partner_source_fingerprint,
        src.pair_fingerprint, src.adjusted_source_fingerprint, rbT, eta0, deta0,
        boundary_fp, policy)
end

function _validate_eta_boundary(src::ReducedAdjointSourceMode, bc::EtaRayIVP)
    _validate_reduced_adjoint_source(src)
    (src.M, src.l, src.m, src.omega, src.primary_source_fingerprint,
        src.partner_source_fingerprint, src.pair_fingerprint,
        src.adjusted_source_fingerprint) ==
        (bc.M, bc.l, bc.m, bc.omega, bc.primary_source_fingerprint,
         bc.partner_source_fingerprint, bc.pair_fingerprint,
         bc.adjusted_source_fingerprint) || throw(ArgumentError(
            "eta boundary mode or provenance key differs from its adjusted source"))
    _require_eta_spin_domain(bc.l)
    T = typeof(bc.M)
    _finite(bc.rb) && bc.rb >= T(2) * bc.M || throw(ArgumentError(
        "eta boundary rb must be finite and satisfy rb >= 2M"))
    bc.rb == src.pair_rb || throw(ArgumentError(
        "eta and corrector-middle data must share the same finite-ray boundary"))
    _finite(bc.eta0) && _finite(bc.deta0) ||
        throw(ArgumentError("eta boundary data must be finite"))
    bc.boundary_fingerprint != 0 ||
        throw(ArgumentError("boundary_fingerprint must be nonzero"))
    bc.policy isa FiniteRayFreeIVP ||
        throw(ArgumentError("eta boundary is not a finite-ray free IVP"))
    _adjusted_Smbmb(src, bc.rb)
    nothing
end

"A bounded `l >= 2` finite-ray solution of the supplied reduced eta equation."
struct IRGEtaMode{N,T<:AbstractFloat,S<:ReducedAdjointSourceMode{N,T}}
    source::S
    boundary::EtaRayIVP{N,T}
    gauge::GHZIRG
end

perturbative_order(::IRGEtaMode{N}) where {N} = N

function solve_irg_eta(src::ReducedAdjointSourceMode, bc::EtaRayIVP)
    source_N, source_T = _adjoint_source_order_precision(src)
    boundary_N, boundary_T = _eta_boundary_order_precision(bc)
    source_N === boundary_N && source_T === boundary_T || throw(ArgumentError(
        "adjusted source and eta boundary have different orders or precisions"))
    _validate_eta_boundary(src, bc)
    S = typeof(src)
    IRGEtaMode{source_N,source_T,S}(src, bc, GHZIRG())
end

function _validate_irg_eta_mode(mode::IRGEtaMode{N,T}) where {N,T}
    source_N, source_T = _adjoint_source_order_precision(mode.source)
    boundary_N, boundary_T = _eta_boundary_order_precision(mode.boundary)
    N === source_N === boundary_N && T === source_T === boundary_T ||
        throw(ArgumentError("eta mode, source, and boundary orders or precisions differ"))
    _validate_eta_boundary(mode.source, mode.boundary)
    mode.gauge isa GHZIRG || throw(ArgumentError("eta mode is not GHZ IRG"))
    nothing
end

function _eta_radius(mode::IRGEtaMode{N,T}, r) where {N,T}
    rT = _checked_real(T, r, "eta evaluation radius")
    _finite(rT) && rT >= mode.boundary.rb && rT >= T(2) * mode.boundary.M ||
        throw(ArgumentError(
            "eta radius must be finite, outside the horizon, and satisfy r >= rb"))
    rT
end

function _eta_scaled_ray_gap(q::T, r::T) where {T<:AbstractFloat}
    delta = r - q
    if q <= delta
        # Here delta/r >= 1/2, so the division cannot lose a nonzero ratio
        # before it scales the possibly tiny q.
        return q * (delta / r)
    end
    # Otherwise q/r > 1/2, and multiplying that bounded ratio by the shorter
    # ray gap avoids forming a potentially overflowing q*delta numerator.
    (q / r) * delta
end

_eta_ratio_square_scale(ratio::T, value) where {T<:AbstractFloat} =
    ratio * (ratio * value)

function _eta_real_muladd(scale::T, value::Complex{T}, offset::Complex{T}) where {T}
    Complex{T}(fma(scale, real(value), real(offset)),
        fma(scale, imag(value), imag(offset)))
end

function _eta_composition_risk(scale::T, eta0::Complex{T}, deta0::Complex{T},
        sourced::Complex{T}, value::Complex{T}) where {T<:AbstractFloat}
    threshold = sqrt(eps(T))
    for component in (real, imag)
        product = scale * component(deta0)
        isfinite(product) || return true
        magnitude = max(abs(product), abs(component(eta0)), abs(component(sourced)))
        magnitude != zero(T) && abs(component(value)) <= threshold * magnitude && return true
    end
    false
end

_eta_fallback_precision(::Type{Float16}) = 256
_eta_fallback_precision(::Type{Float32}) = 1024
_eta_fallback_precision(::Type{Float64}) = 8192
_eta_fallback_precision(::Type{BigFloat}) = max(256, precision(BigFloat))
_eta_fallback_precision(::Type{T}) where {T} = 256

_eta_has_ieee_subnormal(value::T) where {T<:Union{Float16,Float32,Float64}} =
    !iszero(value) && issubnormal(value)
_eta_has_ieee_subnormal(value::Complex{T}) where
        {T<:Union{Float16,Float32,Float64}} =
    _eta_has_ieee_subnormal(real(value)) || _eta_has_ieee_subnormal(imag(value))
_eta_has_ieee_subnormal(value) = false

function _eta_widened_d2(source::Complex{T}, derivative::Complex{T},
        r::T) where {T<:AbstractFloat}
    bits = max(_eta_fallback_precision(T), precision(BigFloat))
    setprecision(BigFloat, bits) do
        r_big = BigFloat(r)
        real_value = -BigFloat(2) * (BigFloat(real(source)) +
            BigFloat(real(derivative)) / r_big)
        imag_value = -BigFloat(2) * (BigFloat(imag(source)) +
            BigFloat(imag(derivative)) / r_big)
        Complex{T}(T(real_value), T(imag_value))
    end
end

function _eta_widened_residual(r, derivative, d2, source, value::T) where
        {T<:AbstractFloat}
    bits = max(_eta_fallback_precision(T), precision(BigFloat))
    setprecision(BigFloat, bits) do
        result = BigFloat(d2) + BigFloat(2) *
            (BigFloat(source) + BigFloat(derivative) / BigFloat(r))
        T(result)
    end
end

function _eta_widened_residual(r, derivative, d2, source,
        value::Complex{T}) where {T<:AbstractFloat}
    bits = max(_eta_fallback_precision(T), precision(BigFloat))
    setprecision(BigFloat, bits) do
        r_big = BigFloat(r)
        real_value = BigFloat(real(d2)) + BigFloat(2) *
            (BigFloat(real(source)) + BigFloat(real(derivative)) / r_big)
        imag_value = BigFloat(imag(d2)) + BigFloat(2) *
            (BigFloat(imag(source)) + BigFloat(imag(derivative)) / r_big)
        Complex{T}(T(real_value), T(imag_value))
    end
end

function _eta_widened_composition(scale::T, eta0::Complex{T}, deta0::Complex{T},
        sourced::Complex{T}) where {T<:AbstractFloat}
    bits = max(_eta_fallback_precision(T), precision(BigFloat))
    setprecision(BigFloat, bits) do
        scale_big = BigFloat(scale)
        real_value = BigFloat(real(eta0)) + scale_big * BigFloat(real(deta0)) -
            BigFloat(2) * BigFloat(real(sourced))
        imag_value = BigFloat(imag(eta0)) + scale_big * BigFloat(imag(deta0)) -
            BigFloat(2) * BigFloat(imag(sourced))
        Complex{T}(T(real_value), T(imag_value))
    end
end

function _eta_compose(scale::T, eta0::Complex{T}, deta0::Complex{T},
        sourced::Complex{T}) where {T<:AbstractFloat}
    offset = eta0 - sourced
    first = _eta_real_muladd(scale, deta0, offset)
    value = first - sourced
    if !_finite(offset) || !_finite(first) || !_finite(value) ||
            _eta_composition_risk(scale, eta0, deta0, sourced, value)
        return _eta_widened_composition(scale, eta0, deta0, sourced)
    end
    value
end

"Evaluate eta from the one-integral Green formula for `(r*eta)'' = -2r*S_mbar_mbar`."
function eta(mode::IRGEtaMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    _validate_irg_eta_mode(mode)
    rT = _eta_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    bc = mode.boundary
    rT == bc.rb && return bc.eta0
    _adjusted_Smbmb(mode.source, rT)
    # `ratio * (r-rb)` is algebraically `rb * (1-rb/r)`, but preserves an
    # adjacent-float radial increment instead of recovering it by subtracting
    # two nearly equal O(1) numbers.  Form the bounded real coefficient before
    # multiplying the possibly complex boundary derivative.
    homogeneous_step = _eta_scaled_ray_gap(bc.rb, rT)
    sourced = _quad(q -> _eta_scaled_ray_gap(q, rT) *
        _adjusted_Smbmb(mode.source, q), bc.rb, rT; atol=at, rtol=rt)
    # Fuse one source subtraction into the boundary-data product, then retry
    # the full three-term composition at wider precision whenever a finite
    # result is hidden by overflow or severe cancellation.
    value = _eta_compose(homogeneous_step, bc.eta0, bc.deta0, sourced)
    _finite(value) || throw(ArgumentError("eta must be finite"))
    value
end

"Evaluate the first radial derivative of the reduced eta Green solution."
function dr_eta(mode::IRGEtaMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    _validate_irg_eta_mode(mode)
    rT = _eta_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    bc = mode.boundary
    rT == bc.rb && return bc.deta0
    _adjusted_Smbmb(mode.source, rT)
    ratio = bc.rb / rT
    moment = _quad(q -> begin
            ray_ratio = q / rT
            _eta_ratio_square_scale(ray_ratio, _adjusted_Smbmb(mode.source, q))
        end, bc.rb, rT; atol=at, rtol=rt)
    homogeneous = _eta_ratio_square_scale(ratio, bc.deta0)
    # As in `eta`, subtract the two source moments successively so finite
    # cancellation is not destroyed by an eager factor of two.
    value = (homogeneous - moment) - moment
    _finite(value) || throw(ArgumentError("dr_eta must be finite"))
    value
end

"Evaluate the second radial derivative using the reduced eta equation itself."
function d2r_eta(mode::IRGEtaMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    _validate_irg_eta_mode(mode)
    rT = _eta_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    derivative = dr_eta(mode, rT; atol=at, rtol=rt)
    source = _adjusted_Smbmb(mode.source, rT)
    combined = source + derivative / rT
    retry_wide = !_finite(combined)
    if retry_wide
        # For small r, `derivative/r` can overflow before cancelling the
        # supplied source.  Move the division outside that cancellation.
        combined = (source * rT + derivative) / rT
    end
    value = -T(2) * combined
    zero_cancellation = iszero(value) && (!iszero(source) || !iszero(derivative))
    if retry_wide || !_finite(value) || _eta_has_ieee_subnormal(value) ||
            zero_cancellation
        value = _eta_widened_d2(source, derivative, rT)
    end
    _finite(value) || throw(ArgumentError("d2r_eta must be finite"))
    value
end

"Pure local residual `eta'' + 2eta'/r + 2S_mbar_mbar` for independent inputs."
function irg_eta_residual(r, deta, d2eta, Smbmb)
    _require_nonbool_real(r, "eta residual radius")
    _require_nonbool_number(deta, "dr_eta")
    _require_nonbool_number(d2eta, "d2r_eta")
    _require_nonbool_number(Smbmb, "adjusted S_mbar_mbar")
    _finite(r) && r > zero(r) ||
        throw(ArgumentError("eta residual radius must be finite and positive"))
    all(_finite, (deta, d2eta, Smbmb)) ||
        throw(ArgumentError("eta residual inputs must be finite"))
    combined = Smbmb + deta / r
    retry_wide = !_finite(combined)
    if retry_wide
        combined = (Smbmb * r + deta) / r
    end
    value = (d2eta + combined) + combined
    zero_cancellation = iszero(value) &&
        (!iszero(deta) || !iszero(d2eta) || !iszero(Smbmb))
    if retry_wide || !_finite(value) || _eta_has_ieee_subnormal(value) ||
            zero_cancellation
        value = _eta_widened_residual(r, deta, d2eta, Smbmb, value)
    end
    _finite(value) || throw(ArgumentError("eta residual is not finite"))
    value
end
