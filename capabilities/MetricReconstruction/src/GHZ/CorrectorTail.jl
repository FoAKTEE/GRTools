"Finite outgoing-ray free-IVP data; this policy is neither retarded nor physical."
struct FiniteRayFreeIVP end

struct TailRayIVP{N,T<:AbstractFloat}
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
    rb::T
    x0::Complex{T}
    boundary_fingerprint::UInt128
    policy::FiniteRayFreeIVP
end

function TailRayIVP(pair::PhysicalCorrectorMiddlePair{N,T};
        rb=pair.primary.boundary.rb, x_nn, boundary_fingerprint, policy) where {N,T}
    _validate_physical_corrector_middle_pair(pair)
    policy isa FiniteRayFreeIVP || throw(ArgumentError(
        "tail data require the finite-ray free-IVP policy, not a retarded/physical claim"))
    rbT = _checked_real(T, rb, "tail boundary rb")
    rbT == pair.primary.boundary.rb ||
        throw(ArgumentError("tail, primary, and partner boundaries must share the same rb"))
    x0 = _checked_complex(T, x_nn, "x_nn boundary value")
    _finite(x0) || throw(ArgumentError("x_nn boundary value must be finite"))
    boundary_fp = _checked_fingerprint(boundary_fingerprint, "boundary_fingerprint")
    psrc, qsrc = pair.primary.head.source, pair.partner.head.source
    phb, pmb = pair.primary.head.boundary, pair.primary.boundary
    qhb, qmb = pair.partner.head.boundary, pair.partner.boundary
    TailRayIVP{N,T}(psrc.M, psrc.l, psrc.m, psrc.omega,
        psrc.source_fingerprint, qsrc.source_fingerprint, pair.pair_fingerprint,
        phb.x0, phb.dx0, pmb.x0, pmb.dx0, qhb.x0, qhb.dx0, qmb.x0, qmb.dx0,
        rbT, x0, boundary_fp, policy)
end

perturbative_order(::TailRayIVP{N}) where {N} = N
_tail_order_precision(::TailRayIVP{N,T}) where {N,T} = (N, T)

function _validate_tail_boundary(pair::PhysicalCorrectorMiddlePair,
        bc::TailRayIVP)
    _validate_physical_corrector_middle_pair(pair)
    psrc, qsrc = pair.primary.head.source, pair.partner.head.source
    bc.M == psrc.M || throw(ArgumentError("tail boundary has a different M"))
    bc.l == psrc.l || throw(ArgumentError("tail boundary has a different l"))
    bc.m == psrc.m || throw(ArgumentError("tail boundary has a different m"))
    bc.omega == psrc.omega || throw(ArgumentError("tail boundary has a different omega"))
    bc.primary_source_fingerprint == psrc.source_fingerprint || throw(ArgumentError(
        "tail boundary has a different primary source fingerprint"))
    bc.partner_source_fingerprint == qsrc.source_fingerprint || throw(ArgumentError(
        "tail boundary has a different partner source fingerprint"))
    bc.pair_fingerprint == pair.pair_fingerprint ||
        throw(ArgumentError("tail boundary has a different pair fingerprint"))
    phb, pmb = pair.primary.head.boundary, pair.primary.boundary
    qhb, qmb = pair.partner.head.boundary, pair.partner.boundary
    (bc.primary_head_x0, bc.primary_head_dx0) == (phb.x0, phb.dx0) ||
        throw(ArgumentError("tail boundary has different primary head IVP data"))
    (bc.primary_middle_x0, bc.primary_middle_dx0) == (pmb.x0, pmb.dx0) ||
        throw(ArgumentError("tail boundary has different primary middle IVP data"))
    (bc.partner_head_x0, bc.partner_head_dx0) == (qhb.x0, qhb.dx0) ||
        throw(ArgumentError("tail boundary has different partner head IVP data"))
    (bc.partner_middle_x0, bc.partner_middle_dx0) == (qmb.x0, qmb.dx0) ||
        throw(ArgumentError("tail boundary has different partner middle IVP data"))
    bc.rb == pair.primary.boundary.rb || throw(ArgumentError(
        "tail, primary, and partner boundaries must share the same rb"))
    _finite(bc.rb) && bc.rb >= typeof(bc.M)(2) * bc.M ||
        throw(ArgumentError("tail boundary rb must be finite and satisfy rb >= 2M"))
    _finite(bc.x0) || throw(ArgumentError("x_nn boundary value must be finite"))
    bc.boundary_fingerprint != 0 ||
        throw(ArgumentError("boundary_fingerprint must be nonzero"))
    bc.policy isa FiniteRayFreeIVP || throw(ArgumentError(
        "tail boundary is not a finite-ray free IVP"))
    nothing
end

"Bounded `l >= 1` scalar-QuadGK correctness-oracle tail; `l=0` and physical data remain future."
struct IRGCorrectorTailMode{N,T<:AbstractFloat,
        P<:PhysicalCorrectorMiddlePair{N,T}}
    pair::P
    boundary::TailRayIVP{N,T}
    gauge::GHZIRG
end

perturbative_order(::IRGCorrectorTailMode{N}) where {N} = N

function _Sln(src::TetradSourceMode{N,T}, r::T) where {N,T}
    value = _checked_complex(T, src.Sln(r), "S_ln")
    _finite(value) || throw(ArgumentError("S_ln must be finite along the integration ray"))
    value
end

function _tail_angular_factors(::Type{T}, l::Integer) where {T<:AbstractFloat}
    lI = _require_middle_spin_domain(l)
    ell = T(lI)
    _finite(ell) || throw(ArgumentError("tail l is not finite when represented in $T"))
    ellp1 = ell + one(T)
    ellp1 > ell || throw(ArgumentError(
        "precision $T cannot resolve the angular labels l and l+1 for l=$lI"))
    c = sqrt(ell / T(2)) * sqrt(ellp1)
    _finite(c) || throw(ArgumentError(
        "tail angular factors are not finite in $T for l=$lI"))
    ell, ellp1, c
end

function solve_irg_corrector_tail(pair::PhysicalCorrectorMiddlePair,
        bc::TailRayIVP)
    N, T = _corrector_order_precision(pair.primary)
    boundary_N, boundary_T = _tail_order_precision(bc)
    N === boundary_N && T === boundary_T || throw(ArgumentError(
        "corrector pair and tail boundary have different orders or precisions"))
    _validate_tail_boundary(pair, bc)
    ell, ellp1, c = _tail_angular_factors(T, pair.primary.head.source.l)
    P = typeof(pair)
    mode = IRGCorrectorTailMode{N,T,P}(pair, bc, GHZIRG())
    at, rt = zero(T), sqrt(eps(T))
    _tail_Q_unchecked(mode, bc.rb, ell, ellp1, c, at, rt)
    mode
end

function _validate_corrector_tail_mode(mode::IRGCorrectorTailMode{N,T}) where {N,T}
    pair_N, pair_T = _corrector_order_precision(mode.pair.primary)
    N === pair_N && T === pair_T || throw(ArgumentError(
        "tail mode and corrector pair have different orders or precisions"))
    _validate_tail_boundary(mode.pair, mode.boundary)
    mode.gauge isa GHZIRG || throw(ArgumentError("corrector tail is not GHZ IRG"))
    ell, ellp1, c = _tail_angular_factors(T, mode.pair.primary.head.source.l)
    at, rt = zero(T), sqrt(eps(T))
    _tail_Q_unchecked(mode, mode.boundary.rb, ell, ellp1, c, at, rt)
    ell, ellp1, c
end

function _tail_radius(mode::IRGCorrectorTailMode{N,T}, r) where {N,T}
    rT = _checked_real(T, r, "tail evaluation radius")
    rT >= mode.boundary.rb && rT >= T(2) * mode.boundary.M && _finite(rT) ||
        throw(ArgumentError(
            "tail radius must be finite, outside the horizon, and satisfy r >= rb"))
    rT
end

function _tail_Q_unchecked(mode::IRGCorrectorTailMode{N,T}, r::T, ell::T, ellp1::T,
        c::T, atol::T, rtol::T) where {N,T}
    pair = mode.pair
    primary, head = pair.primary, pair.primary.head
    src = head.source
    A = x_mmb(head, r; atol, rtol)
    dA = dr_x_mmb(head, r; atol, rtol)
    B = x_nm(primary, r; atol, rtol)
    dB = dr_x_nm(primary, r; atol, rtol)
    C = _x_nmb_unchecked(pair, r; atol, rtol)
    dC = _dr_x_nmb_unchecked(pair, r; atol, rtol)
    all(_finite, (A, dA, B, dB, C, dC)) ||
        throw(ArgumentError("corrector components must be finite along the tail ray"))

    Sll = _Sll(head, r)
    Sln = _Sln(src, r)
    f = one(T) - T(2) * src.M / r
    two_i_omega = Complex{T}(zero(T), T(2)) * src.omega
    four_i_omega = Complex{T}(zero(T), T(4)) * src.omega
    lambda_over_r2 = (ell / r) * (ellp1 / r)
    _finite(lambda_over_r2) || throw(ArgumentError(
        "tail lambda/r^2 is not finite in $T at r=$r"))
    L0 = -f * Sll +
        (two_i_omega + T(2) / r - T(2) * src.M / r / r) * dA +
        (four_i_omega / r + T(2) / r / r - lambda_over_r2) * A
    D = B - C
    dD = dB - dC
    Q = Sln + L0 / T(2) + c * (dD + T(3) * D / r) / (T(2) * r)
    _finite(Q) || throw(ArgumentError("x_nn tail source Q must be finite"))
    Q
end

function _x_nn_unchecked(mode::IRGCorrectorTailMode{N,T}, r::T, ell::T, ellp1::T,
        c::T, atol::T, rtol::T) where {N,T}
    bc = mode.boundary
    r == bc.rb && return bc.x0
    integral = _quad(q -> q * ((q / r) * _tail_Q_unchecked(
        mode, q, ell, ellp1, c, atol, rtol)),
        bc.rb, r; atol, rtol)
    _finite(integral) || throw(ArgumentError("x_nn tail integral must be finite"))
    value = (bc.rb / r) * bc.x0 + integral
    _finite(value) || throw(ArgumentError("x_nn must be finite"))
    value
end

function x_nn(mode::IRGCorrectorTailMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    ell, ellp1, c = _validate_corrector_tail_mode(mode)
    rT = _tail_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    _x_nn_unchecked(mode, rT, ell, ellp1, c, at, rt)
end

function dr_x_nn(mode::IRGCorrectorTailMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    ell, ellp1, c = _validate_corrector_tail_mode(mode)
    rT = _tail_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    value = _x_nn_unchecked(mode, rT, ell, ellp1, c, at, rt)
    Q = _tail_Q_unchecked(mode, rT, ell, ellp1, c, at, rt)
    derivative = rT * Q - value / rT
    _finite(derivative) || throw(ArgumentError("dr_x_nn must be finite"))
    derivative
end

function irg_tail_residual(r, N, dN, Q)
    _require_nonbool_real(r, "tail residual radius")
    _require_nonbool_number(N, "x_nn")
    _require_nonbool_number(dN, "dr_x_nn")
    _require_nonbool_number(Q, "tail source Q")
    _finite(r) && r > zero(r) ||
        throw(ArgumentError("tail residual radius must be finite and positive"))
    all(_finite, (N, dN, Q)) ||
        throw(ArgumentError("tail residual inputs must be finite"))
    dN + N / r - r * Q
end

struct IRGCorrectorTetradComponents{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    primary_source_fingerprint::UInt128
    partner_source_fingerprint::UInt128
    pair_fingerprint::UInt128
    tail_boundary_fingerprint::UInt128
    rb::T
    r::T
    A::Complex{T}
    B::Complex{T}
    C::Complex{T}
    N::Complex{T}
    normalization::EinsteinRHS
    policy::FiniteRayFreeIVP
    gauge::GHZIRG
    tetrad::KinnersleyMinusPlus
    slicing::SchwarzschildOutgoing
    phase::ExpMinusIomegaUPlusImPhi
end

perturbative_order(::IRGCorrectorTetradComponents{N}) where {N} = N

function corrector_tetrad_components(mode::IRGCorrectorTailMode{N,T}, r;
        atol=zero(T), rtol=sqrt(eps(T))) where {N,T}
    ell, ellp1, c = _validate_corrector_tail_mode(mode)
    rT = _tail_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    pair, primary, head = mode.pair, mode.pair.primary, mode.pair.primary.head
    A = x_mmb(head, rT; atol=at, rtol=rt)
    B = x_nm(primary, rT; atol=at, rtol=rt)
    C = _x_nmb_unchecked(pair, rT; atol=at, rtol=rt)
    Nvalue = _x_nn_unchecked(mode, rT, ell, ellp1, c, at, rt)
    all(_finite, (A, B, C, Nvalue)) ||
        throw(ArgumentError("assembled corrector tetrad components must be finite"))
    psrc, qsrc = head.source, pair.partner.head.source
    IRGCorrectorTetradComponents{N,T}(psrc.M, psrc.l, psrc.m, psrc.omega,
        psrc.source_fingerprint, qsrc.source_fingerprint, pair.pair_fingerprint,
        mode.boundary.boundary_fingerprint, mode.boundary.rb, rT, A, B, C, Nvalue,
        psrc.normalization, mode.boundary.policy, mode.gauge, psrc.tetrad, psrc.slicing,
        psrc.phase)
end
