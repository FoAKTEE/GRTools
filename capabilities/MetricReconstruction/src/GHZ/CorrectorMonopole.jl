"Explicit marker that the `l=0` corrector has no spin-one `B` or `C` coefficient."
struct MonopoleSpinOneComponentsAbsent end

"""
    MonopoleTailRayIVP(head; rb=head.boundary.rb, x_nn,
                       boundary_fingerprint, policy)

Finite outgoing-ray data for the `l=m=0` head-only `x_nn` equation.  The data are a
free finite-ray IVP: they are neither retarded/physical boundary data nor a mass
completion.  No partner or spin-one provenance exists in this sector.
"""
struct MonopoleTailRayIVP{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    source_fingerprint::UInt128
    head_x0::Complex{T}
    head_dx0::Complex{T}
    rb::T
    x0::Complex{T}
    boundary_fingerprint::UInt128
    policy::FiniteRayFreeIVP
end

perturbative_order(::MonopoleTailRayIVP{N}) where {N} = N

function _validate_monopole_head(head::IRGCorrectorHeadMode{N,T}) where {N,T}
    N isa Integer && !(N isa Bool) && N > 0 || throw(ArgumentError(
        "monopole head perturbative order must be a positive non-Boolean integer"))
    src, bc = head.source, head.boundary
    src.l == 0 && src.m == 0 ||
        throw(ArgumentError("GHZ monopole tail requires the l=m=0 source mode"))
    src.Slm isa ZeroRadialSource || throw(ArgumentError(
        "GHZ monopole tail requires S_lm to be the explicit ZeroRadialSource"))
    src.normalization isa EinsteinRHS ||
        throw(ArgumentError("GHZ monopole tail requires an EinsteinRHS source"))
    head.gauge isa GHZIRG ||
        throw(ArgumentError("monopole corrector head is not GHZ IRG"))
    _finite(src.M) && src.M > zero(T) ||
        throw(ArgumentError("monopole source M must be finite and positive"))
    _finite(src.omega) ||
        throw(ArgumentError("monopole source omega must be finite"))
    src.source_fingerprint != 0 ||
        throw(ArgumentError("monopole source fingerprint must be nonzero"))
    _same_mode(src, bc) ||
        throw(ArgumentError("monopole head and inner-ray boundary mode keys differ"))
    _finite(bc.rb) && bc.rb >= T(2) * src.M || throw(ArgumentError(
        "monopole head boundary rb must be finite and satisfy rb >= 2M"))
    all(_finite, (bc.x0, bc.dx0)) ||
        throw(ArgumentError("monopole head boundary data must be finite"))
    Sll = _checked_complex(T, src.Sll(bc.rb), "S_ll")
    _finite(Sll) ||
        throw(ArgumentError("S_ll must be finite at the monopole inner boundary"))
    nothing
end

function MonopoleTailRayIVP(head::IRGCorrectorHeadMode{N,T};
        rb=head.boundary.rb, x_nn, boundary_fingerprint, policy) where {N,T}
    _validate_monopole_head(head)
    policy isa FiniteRayFreeIVP || throw(ArgumentError(
        "monopole tail data require the finite-ray free-IVP policy, not a retarded/physical claim"))
    rbT = _checked_real(T, rb, "monopole tail boundary rb")
    rbT == head.boundary.rb || throw(ArgumentError(
        "monopole head and tail boundaries must share the same rb"))
    x0 = _checked_complex(T, x_nn, "monopole x_nn boundary value")
    _finite(x0) ||
        throw(ArgumentError("monopole x_nn boundary value must be finite"))
    boundary_fp = _checked_fingerprint(boundary_fingerprint, "boundary_fingerprint")
    src, hbc = head.source, head.boundary
    MonopoleTailRayIVP{N,T}(src.M, src.l, src.m, src.omega,
        src.source_fingerprint, hbc.x0, hbc.dx0, rbT, x0, boundary_fp, policy)
end

function _validate_monopole_tail_boundary(head::IRGCorrectorHeadMode,
        bc::MonopoleTailRayIVP)
    _validate_monopole_head(head)
    src, hbc = head.source, head.boundary
    bc.M == src.M ||
        throw(ArgumentError("monopole tail boundary has a different M"))
    bc.l == 0 && bc.l == src.l ||
        throw(ArgumentError("monopole tail boundary must carry l=0"))
    bc.m == 0 && bc.m == src.m ||
        throw(ArgumentError("monopole tail boundary must carry m=0"))
    bc.omega == src.omega ||
        throw(ArgumentError("monopole tail boundary has a different omega"))
    bc.source_fingerprint == src.source_fingerprint || throw(ArgumentError(
        "monopole tail boundary has a different source fingerprint"))
    (bc.head_x0, bc.head_dx0) == (hbc.x0, hbc.dx0) || throw(ArgumentError(
        "monopole tail boundary has different head IVP data"))
    bc.rb == hbc.rb || throw(ArgumentError(
        "monopole head and tail boundaries must share the same rb"))
    _finite(bc.rb) && bc.rb >= typeof(bc.M)(2) * bc.M || throw(ArgumentError(
        "monopole tail boundary rb must be finite and satisfy rb >= 2M"))
    _finite(bc.x0) ||
        throw(ArgumentError("monopole x_nn boundary value must be finite"))
    bc.boundary_fingerprint != 0 ||
        throw(ArgumentError("boundary_fingerprint must be nonzero"))
    bc.policy isa FiniteRayFreeIVP ||
        throw(ArgumentError("monopole tail boundary is not a finite-ray free IVP"))
    nothing
end

"""
Bounded scalar-QuadGK `l=m=0` GHZ corrector tail.  This is a finite-ray
head-only correctness kernel, not physical boundary data and not monopole mass
completion.
"""
struct IRGCorrectorMonopoleTailMode{N,T<:AbstractFloat,
        H<:IRGCorrectorHeadMode{N,T}}
    head::H
    boundary::MonopoleTailRayIVP{N,T}
    gauge::GHZIRG
end

perturbative_order(::IRGCorrectorMonopoleTailMode{N}) where {N} = N

function _monopole_tail_Q_unchecked(mode::IRGCorrectorMonopoleTailMode{N,T},
        r::T, atol::T, rtol::T) where {N,T}
    head, src = mode.head, mode.head.source
    A = x_mmb(head, r; atol, rtol)
    dA = dr_x_mmb(head, r; atol, rtol)
    all(_finite, (A, dA)) || throw(ArgumentError(
        "monopole corrector head must be finite along the tail ray"))
    Sll = _Sll(head, r)
    Sln = _Sln(src, r)
    M_over_r = src.M / r
    f = one(T) - T(2) * M_over_r
    i = Complex{T}(zero(T), one(T))
    two_i = Complex{T}(zero(T), T(2))
    A_over_r = A / r
    gradient_term = (i * src.omega) * dA +
        ((one(T) - M_over_r) * dA) / r
    if !_finite(gradient_term)
        # Retry in a cancellation-preserving form.  This path is needed when
        # the two exact terms are finite only after their large parts cancel.
        gradient_term = ((i * (src.omega * r) + (one(T) - M_over_r)) * dA) / r
    end
    A_term = two_i * (src.omega * A_over_r) + A_over_r / r
    if !_finite(A_term)
        A_coefficient = one(Complex{T}) + two_i * (src.omega * r)
        A_term = ((A_coefficient * A) / r) / r
    end
    # Apply the outer factor 1/2 analytically before forming frequency
    # coefficients.  In particular, never materialize 4im*omega: it can
    # overflow at low precision even when the final Q0 term is finite.
    Q = Sln - (f / T(2)) * Sll + gradient_term + A_term
    _finite(Q) ||
        throw(ArgumentError("monopole x_nn tail source Q0 must be finite"))
    Q
end

function solve_irg_corrector_monopole_tail(
        head::IRGCorrectorHeadMode{HN,HT},
        bc::MonopoleTailRayIVP{BN,BT}) where {HN,HT,BN,BT}
    HN === BN && HT === BT || throw(ArgumentError(
        "monopole head and tail boundary have different orders or precisions"))
    _validate_monopole_tail_boundary(head, bc)
    H = typeof(head)
    mode = IRGCorrectorMonopoleTailMode{HN,HT,H}(head, bc, GHZIRG())
    _monopole_tail_Q_unchecked(mode, bc.rb, zero(HT), sqrt(eps(HT)))
    mode
end

function _validate_corrector_monopole_tail_mode(
        mode::IRGCorrectorMonopoleTailMode{N,T}) where {N,T}
    _validate_monopole_tail_boundary(mode.head, mode.boundary)
    mode.gauge isa GHZIRG ||
        throw(ArgumentError("monopole corrector tail is not GHZ IRG"))
    _monopole_tail_Q_unchecked(mode, mode.boundary.rb, zero(T), sqrt(eps(T)))
    nothing
end

function _monopole_tail_radius(mode::IRGCorrectorMonopoleTailMode{N,T}, r) where {N,T}
    rT = _checked_real(T, r, "monopole tail evaluation radius")
    _finite(rT) && rT >= mode.boundary.rb && rT >= T(2) * mode.boundary.M ||
        throw(ArgumentError(
            "monopole tail radius must be finite, outside the horizon, and satisfy r >= rb"))
    rT
end

function _x_nn_monopole_unchecked(mode::IRGCorrectorMonopoleTailMode{N,T},
        r::T, atol::T, rtol::T) where {N,T}
    bc = mode.boundary
    r == bc.rb && return bc.x0
    integral = _quad(q -> q * ((q / r) * _monopole_tail_Q_unchecked(
        mode, q, atol, rtol)), bc.rb, r; atol, rtol)
    _finite(integral) ||
        throw(ArgumentError("monopole x_nn tail integral must be finite"))
    value = (bc.rb / r) * bc.x0 + integral
    _finite(value) || throw(ArgumentError("monopole x_nn must be finite"))
    value
end

function x_nn(mode::IRGCorrectorMonopoleTailMode{N,T}, r;
        atol=zero(T), rtol=sqrt(eps(T))) where {N,T}
    _validate_corrector_monopole_tail_mode(mode)
    rT = _monopole_tail_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    _x_nn_monopole_unchecked(mode, rT, at, rt)
end

function dr_x_nn(mode::IRGCorrectorMonopoleTailMode{N,T}, r;
        atol=zero(T), rtol=sqrt(eps(T))) where {N,T}
    _validate_corrector_monopole_tail_mode(mode)
    rT = _monopole_tail_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    value = _x_nn_monopole_unchecked(mode, rT, at, rt)
    Q = _monopole_tail_Q_unchecked(mode, rT, at, rt)
    derivative = rT * Q - value / rT
    _finite(derivative) ||
        throw(ArgumentError("monopole dr_x_nn must be finite"))
    derivative
end

"""
`l=m=0` corrector coefficient data.  `A=x_mbar_m` and `N=x_nn` are the only
component coefficients: `spin_one_components` records that `B` and `C` do not
exist here.  This carrier is not a coordinate metric or a mass completion.
"""
struct IRGCorrectorMonopoleCoefficients{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    source_fingerprint::UInt128
    tail_boundary_fingerprint::UInt128
    rb::T
    r::T
    A::Complex{T}
    N::Complex{T}
    spin_one_components::MonopoleSpinOneComponentsAbsent
    normalization::EinsteinRHS
    policy::FiniteRayFreeIVP
    gauge::GHZIRG
    tetrad::KinnersleyMinusPlus
    slicing::SchwarzschildOutgoing
    phase::ExpMinusIomegaUPlusImPhi
end

perturbative_order(::IRGCorrectorMonopoleCoefficients{N}) where {N} = N

function corrector_monopole_coefficients(
        mode::IRGCorrectorMonopoleTailMode{N,T}, r;
        atol=zero(T), rtol=sqrt(eps(T))) where {N,T}
    _validate_corrector_monopole_tail_mode(mode)
    rT = _monopole_tail_radius(mode, r)
    at, rt = _tolerances(T, atol, rtol)
    A = x_mmb(mode.head, rT; atol=at, rtol=rt)
    Nvalue = _x_nn_monopole_unchecked(mode, rT, at, rt)
    all(_finite, (A, Nvalue)) || throw(ArgumentError(
        "assembled monopole corrector coefficients must be finite"))
    src, bc = mode.head.source, mode.boundary
    IRGCorrectorMonopoleCoefficients{N,T}(src.M, src.l, src.m, src.omega,
        src.source_fingerprint, bc.boundary_fingerprint, bc.rb, rT, A, Nvalue,
        MonopoleSpinOneComponentsAbsent(), src.normalization, bc.policy,
        mode.gauge, src.tetrad, src.slicing, src.phase)
end
