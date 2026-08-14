"An owned read-only ray snapshot; construction always breaks mutable aliases."
struct _RaySnapshot{T} <: AbstractVector{T}
    _storage::Vector{T}

    function _RaySnapshot{T}(values::AbstractVector{T}) where {T}
        storage = Vector{T}(undef, length(values))
        copyto!(storage, values)
        new{T}(storage)
    end
end

@inline Base.IndexStyle(::Type{<:_RaySnapshot}) = IndexLinear()
@inline Base.size(values::_RaySnapshot) = size(getfield(values, :_storage))
@inline Base.length(values::_RaySnapshot) = length(getfield(values, :_storage))
@inline Base.getindex(values::_RaySnapshot, index::Int) =
    @inbounds getfield(values, :_storage)[index]
@inline Base.iterate(values::_RaySnapshot, state...) =
    iterate(getfield(values, :_storage), state...)
Base.setindex!(::_RaySnapshot, value, indices...) =
    throw(ArgumentError("ray snapshots are read-only"))
Base.copy(values::_RaySnapshot) = copy(getfield(values, :_storage))
Base.propertynames(::_RaySnapshot, private::Bool=false) = ()
function Base.getproperty(::_RaySnapshot, name::Symbol)
    name === :_storage && throw(ArgumentError("ray snapshots do not expose mutable storage"))
    throw(ArgumentError("ray snapshots have no public property $name"))
end

_ray_snapshot(values::AbstractVector{T}) where {T} = _RaySnapshot{T}(values)

"""
    RayMomentQuadratureDiagnostics

Cumulative sums of the local embedded-rule norm indicators returned by `QuadGK` for each
normalized source-moment vector.  These diagnostics use different endpoint-dependent kernels
on different cells.  They are neither propagated field-error bounds nor quantities that may be
meaningfully summed across source streams.
"""
struct RayMomentQuadratureDiagnostics{T<:AbstractFloat}
    primary_sll::_RaySnapshot{T}
    primary_slm::_RaySnapshot{T}
    partner_sll::_RaySnapshot{T}
    partner_slm::_RaySnapshot{T}
    primary_sln::_RaySnapshot{T}

    function RayMomentQuadratureDiagnostics{T}(primary_sll::AbstractVector{T},
            primary_slm::AbstractVector{T}, partner_sll::AbstractVector{T},
            partner_slm::AbstractVector{T}, primary_sln::AbstractVector{T}) where
            {T<:AbstractFloat}
        snapshots = (_ray_snapshot(primary_sll), _ray_snapshot(primary_slm),
            _ray_snapshot(partner_sll), _ray_snapshot(partner_slm),
            _ray_snapshot(primary_sln))
        count = length(snapshots[1])
        count > 0 || throw(ArgumentError("quadrature diagnostics must not be empty"))
        all(values -> length(values) == count, snapshots) ||
            throw(ArgumentError("quadrature diagnostic streams have different lengths"))
        all(values -> all(value -> _finite(value) && value >= zero(T), values),
            snapshots) || throw(ArgumentError(
            "quadrature diagnostics must be finite and nonnegative"))
        all(values -> values[1] == zero(T), snapshots) ||
            throw(ArgumentError("accumulated quadrature diagnostics must start at zero"))
        all(values -> all(j -> values[j] >= values[j - 1], 2:count), snapshots) ||
            throw(ArgumentError("accumulated quadrature diagnostics must be nondecreasing"))
        new{T}(snapshots[1], snapshots[2], snapshots[3], snapshots[4], snapshots[5])
    end
end

"Complete immutable provenance copied from a finite-ray GHZ corrector-tail mode."
struct RayMomentProvenance{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    primary_source_fingerprint::UInt128
    partner_source_fingerprint::UInt128
    pair_fingerprint::UInt128
    tail_boundary_fingerprint::UInt128
    rb::T
    tail_x0::Complex{T}
    primary_head_x0::Complex{T}
    primary_head_dx0::Complex{T}
    primary_middle_x0::Complex{T}
    primary_middle_dx0::Complex{T}
    partner_head_x0::Complex{T}
    partner_head_dx0::Complex{T}
    partner_middle_x0::Complex{T}
    partner_middle_dx0::Complex{T}
    normalization::EinsteinRHS
    policy::FiniteRayFreeIVP
    gauge::GHZIRG
    tetrad::KinnersleyMinusPlus
    slicing::SchwarzschildOutgoing
    phase::ExpMinusIomegaUPlusImPhi
end

perturbative_order(::RayMomentProvenance{N}) where {N} = N

"Local-embedded-rule diagnostic metadata accumulated through one exact grid point."
struct RayMomentDiagnosticPoint{T<:AbstractFloat}
    primary_sll::T
    primary_slm::T
    partner_sll::T
    partner_slm::T
    primary_sln::T

    function RayMomentDiagnosticPoint{T}(primary_sll::T, primary_slm::T,
            partner_sll::T, partner_slm::T, primary_sln::T) where {T<:AbstractFloat}
        values = (primary_sll, primary_slm, partner_sll, partner_slm, primary_sln)
        all(value -> _finite(value) && value >= zero(T), values) ||
            throw(ArgumentError("quadrature diagnostic point must be finite and nonnegative"))
        new{T}(primary_sll, primary_slm, partner_sll, partner_slm, primary_sln)
    end
end

"All cached corrector components at one exact grid point, with full mode provenance."
struct RayMomentComponents{N,T<:AbstractFloat}
    provenance::RayMomentProvenance{N,T}
    grid_index::Int
    r::T
    A::Complex{T}
    dA::Complex{T}
    B::Complex{T}
    dB::Complex{T}
    C::Complex{T}
    dC::Complex{T}
    N::Complex{T}
    quadrature_diagnostic::RayMomentDiagnosticPoint{T}
    atol::T
    rtol::T
end

perturbative_order(::RayMomentComponents{N}) where {N} = N

"""
    RayMomentCache

Materialized, solver-independent values of a bounded `l >= 1` finite-ray GHZ corrector pair.
The cache supports exact-grid lookup only: it deliberately performs no interpolation and keeps
the scalar `QuadGK` evaluators as independent irregular-point oracles.

Each entry of `quadrature_diagnostics` is a cumulative sum of local `QuadGK` embedded-rule norm
indicators for one normalized source-moment vector.  The indicators are integration diagnostics,
not propagated error bounds for reconstructed fields, and no cross-stream total is defined.
"""
struct RayMomentCache{N,T<:AbstractFloat}
    provenance::RayMomentProvenance{N,T}
    radii::_RaySnapshot{T}
    A::_RaySnapshot{Complex{T}}
    dA::_RaySnapshot{Complex{T}}
    B::_RaySnapshot{Complex{T}}
    dB::_RaySnapshot{Complex{T}}
    C::_RaySnapshot{Complex{T}}
    dC::_RaySnapshot{Complex{T}}
    N::_RaySnapshot{Complex{T}}
    quadrature_diagnostics::RayMomentQuadratureDiagnostics{T}
    atol::T
    rtol::T

    function RayMomentCache{N,T}(provenance::RayMomentProvenance{N,T},
            radii::AbstractVector{T}, A::AbstractVector{Complex{T}},
            dA::AbstractVector{Complex{T}}, B::AbstractVector{Complex{T}},
            dB::AbstractVector{Complex{T}}, C::AbstractVector{Complex{T}},
            dC::AbstractVector{Complex{T}}, Nvalues::AbstractVector{Complex{T}},
            diagnostics::RayMomentQuadratureDiagnostics{T}, atol::T, rtol::T,
            ::Val{:validated}) where {N,T<:AbstractFloat}
        cache = new{N,T}(provenance, _ray_snapshot(radii), _ray_snapshot(A),
            _ray_snapshot(dA), _ray_snapshot(B), _ray_snapshot(dB), _ray_snapshot(C),
            _ray_snapshot(dC), _ray_snapshot(Nvalues), diagnostics, atol, rtol)
        _validate_ray_moment_cache(cache)
        cache
    end
end

perturbative_order(::RayMomentCache{N}) where {N} = N

function _ray_moment_provenance(mode::IRGCorrectorTailMode{N,T}) where {N,T}
    pair, bc = mode.pair, mode.boundary
    primary, partner = pair.primary, pair.partner
    psrc, qsrc = primary.head.source, partner.head.source
    phb, pmb = primary.head.boundary, primary.boundary
    qhb, qmb = partner.head.boundary, partner.boundary
    RayMomentProvenance{N,T}(psrc.M, psrc.l, psrc.m, psrc.omega,
        psrc.source_fingerprint, qsrc.source_fingerprint, pair.pair_fingerprint,
        bc.boundary_fingerprint, bc.rb, bc.x0, phb.x0, phb.dx0, pmb.x0, pmb.dx0,
        qhb.x0, qhb.dx0, qmb.x0, qmb.dx0, psrc.normalization, bc.policy, mode.gauge,
        psrc.tetrad, psrc.slicing, psrc.phase)
end

function _ray_grid(::Type{T}, radii, rb::T) where {T<:AbstractFloat}
    radii isa AbstractVector ||
        throw(ArgumentError("ray radii must be an AbstractVector"))
    isempty(radii) && throw(ArgumentError("ray radii must not be empty"))
    grid = Vector{T}(undef, length(radii))
    for (j, radius) in enumerate(radii)
        converted = try
            _checked_real(T, radius, "ray radius")
        catch error
            error isa ArgumentError && rethrow()
            throw(ArgumentError("ray radius $j is not representable in $T"))
        end
        _finite(converted) || throw(ArgumentError("ray radius $j must be finite"))
        grid[j] = converted
    end
    grid[1] == rb ||
        throw(ArgumentError("the first ray radius must equal the tail boundary rb exactly"))
    for j in 2:length(grid)
        grid[j] > grid[j - 1] ||
            throw(ArgumentError("ray radii must be strictly increasing"))
    end
    grid
end

# Stable positive kernels used when a cell is much shorter than its radius.
function _ray_series_cap(y::T) where {T<:AbstractFloat}
    bits = precision(y)
    bits <= typemax(Int) - 64 ||
        throw(ArgumentError("ray-kernel precision is too large for a safe series cap"))
    max(256, bits + 64)
end

function _ray_series_cutoff(y::T) where {T<:AbstractFloat}
    base_cutoff = inv(T(8))
    precision(y) >= 19 && return base_cutoff
    unit_roundoff = ldexp(one(y), 1 - precision(y))
    conditioned_cutoff = T(2) * unit_roundoff^inv(T(6))
    max(base_cutoff, min(T(5) / T(8), conditioned_cutoff))
end

function _ray_delta_minus_log1p(y::T) where {T<:AbstractFloat}
    y >= zero(T) && _finite(y) || throw(ArgumentError("invalid ray cell ratio"))
    y == zero(T) && return zero(T)
    y > _ray_series_cutoff(y) && return y - log1p(y)
    term = y * y / T(2)
    result = term
    for n in 3:_ray_series_cap(y)
        term *= -y * T(n - 1) / T(n)
        next = result + term
        next == result && return result
        result = next
    end
    throw(ArgumentError("ray delta-log1p series did not converge at $(precision(y))-bit precision"))
end

# F(q,b)/q^2, where F=b^2/4+3q^2/4-q^3/b-(3q^2/2)log(b/q).
function _ray_F_kernel(y::T) where {T<:AbstractFloat}
    y >= zero(T) && _finite(y) || throw(ArgumentError("invalid ray cell ratio"))
    y == zero(T) && return zero(T)
    if y > _ray_series_cutoff(y)
        onepy = one(T) + y
        return onepy^2 / T(4) + T(3) / T(4) - inv(onepy) -
               T(3) * log1p(y) / T(2)
    end
    power = y * y * y
    result = power / T(2)
    for n in 4:_ray_series_cap(y)
        power *= -y
        term = T(2n - 3) * power / T(2n)
        next = result + term
        next == result && return result
        result = next
    end
    throw(ArgumentError("ray F-kernel series did not converge at $(precision(y))-bit precision"))
end

# `(q/b) F(q,b)/q^2`, evaluated without first forming the quadratically growing
# unscaled F kernel.  Multiplication by q gives the normalized F/b kernel used below.
function _ray_tF_kernel(y::T) where {T<:AbstractFloat}
    y >= zero(T) && _finite(y) || throw(ArgumentError("invalid ray cell ratio"))
    y == zero(T) && return zero(T)
    if y <= _ray_series_cutoff(y)
        return _ray_F_kernel(y) / (one(T) + y)
    end
    invy = inv(y)
    t = invy / (one(T) + invy)
    (y / T(4) + one(T) / T(4) + T(3) * t / T(4)) - t * t -
        T(3) * (t * log1p(y)) / T(2)
end

function _ray_checked_values(values::Vector{Complex{T}}, label) where {T}
    all(_finite, values) ||
        throw(ArgumentError("$label moments must remain finite on every ray cell"))
    values
end

function _ray_quad_vector(f, a::T, b::T, atol::T, rtol::T, label) where
        {T<:AbstractFloat}
    values, error = quadgk(f, a, b; atol, rtol)
    values isa Vector{Complex{T}} ||
        throw(ArgumentError("$label quadrature returned an unexpected value type"))
    _ray_checked_values(values, label)
    errorT = try
        T(error)
    catch
        throw(ArgumentError("$label quadrature error is not representable in $T"))
    end
    _finite(errorT) && errorT >= zero(T) ||
        throw(ArgumentError("$label quadrature error must be finite and nonnegative"))
    values, errorT
end

function _ray_sll_moments(head::IRGCorrectorHeadMode{N,T}, a::T, b::T,
        atol::T, rtol::T, label) where {N,T}
    _ray_quad_vector(a, b, atol, rtol, label) do q
        source = _Sll(head, q)
        y = (b - q) / q
        t = q / b
        one_minus_t = (b - q) / b
        gap = b - q
        middle_value_kernel = gap * (one_minus_t * (one(T) + T(2) * t) / T(2))
        values = Complex{T}[
            (q * one_minus_t) * source,
            (t * t) * source,
            middle_value_kernel * source,
            one_minus_t * (T(4) * t^2 + t + one(T)) * source / T(6),
            (q * (t * _ray_delta_minus_log1p(y))) * source,
            (q * _ray_tF_kernel(y)) * source,
        ]
        _ray_checked_values(values, label)
    end
end

function _ray_slm_moments(middle::IRGCorrectorMiddleMode{N,T}, a::T, b::T,
        atol::T, rtol::T, label) where {N,T}
    _ray_quad_vector(a, b, atol, rtol, label) do q
        source = _Slm(middle, q)
        t = q / b
        one_minus_t = (b - q) / b
        gap = b - q
        values = Complex{T}[
            (gap * (one(T) + t + t^2)) * source,
            (one(T) + T(2) * t^3) * source / T(3),
            (gap * (one_minus_t * (one(T) + T(2) * t) / T(2))) * source,
        ]
        _ray_checked_values(values, label)
    end
end

function _ray_sln_I2_over_b(src::TetradSourceMode{N,T}, a::T, b::T, atol::T,
        rtol::T) where {N,T}
    value, error = quadgk(q -> begin
        result = (q * (q / b)) * _Sln(src, q)
        _finite(result) || throw(ArgumentError(
            "normalized S_ln second moment must remain finite on every ray cell"))
        result
    end, a, b; atol, rtol)
    _finite(value) || throw(ArgumentError(
        "normalized S_ln second moment must remain finite on every ray cell"))
    errorT = try
        T(error)
    catch
        throw(ArgumentError("S_ln quadrature error is not representable in $T"))
    end
    _finite(errorT) && errorT >= zero(T) ||
        throw(ArgumentError("S_ln quadrature error must be finite and nonnegative"))
    Complex{T}(value), errorT
end

function _ray_head_step(A::Complex{T}, dA::Complex{T}, a::T, b::T,
        moments::Vector{Complex{T}}) where {T<:AbstractFloat}
    h = b - a
    ratio = a / b
    next_A = A + dA * (h * ratio) - moments[1]
    next_dA = ratio^2 * dA - moments[2]
    all(_finite, (next_A, next_dA)) ||
        throw(ArgumentError("head recurrence produced nonfinite cache data"))
    next_A, next_dA
end

function _ray_middle_step(B::Complex{T}, dB::Complex{T}, dA::Complex{T},
        a::T, b::T, sll::Vector{Complex{T}}, slm::Vector{Complex{T}},
        c::T) where {T<:AbstractFloat}
    ratio = a / b
    h = b - a
    H = dB + T(2) * (B / a)
    homogeneous = ratio^2 * B +
        H * (h * ((one(T) + ratio + ratio^2) / T(3)))
    t = ratio
    k_value_scaled = h * ((h / b) * ((one(T) + T(2) * t) / T(2)))
    k_derivative = (one(T) - t) * (T(4) * t^2 + t + one(T)) / T(6)
    next_B = homogeneous - (T(2) / T(3)) * slm[1] +
        (c / T(3)) * (dA * k_value_scaled - sll[3])
    next_dB = H - T(2) * (homogeneous / b) - T(2) * slm[2] +
        c * (dA * k_derivative - sll[4])
    all(_finite, (next_B, next_dB)) ||
        throw(ArgumentError("middle recurrence produced nonfinite cache data"))
    next_B, next_dB
end

function _ray_integral_A_over_b(A::Complex{T}, dA::Complex{T}, a::T, b::T,
        sll::Vector{Complex{T}}) where {T<:AbstractFloat}
    h = b - a
    delta = h / a
    ratio = a / b
    value = A * (h / b) + dA * ((a * ratio) * _ray_delta_minus_log1p(delta)) -
        sll[5]
    _finite(value) || throw(ArgumentError("normalized cell integral of A is nonfinite"))
    value
end

function _ray_integral_B_over_b(B::Complex{T}, dB::Complex{T}, dA::Complex{T},
        a::T, b::T, sll::Vector{Complex{T}}, slm::Vector{Complex{T}},
        c::T) where {T<:AbstractFloat}
    h = b - a
    ratio = a / b
    h_over_b = h / b
    H = dB + T(2) * (B / a)
    homogeneous = B * (ratio * h_over_b) +
        H * (h * (h_over_b * ((ratio + inv(T(2))) / T(3))))
    boundary_F_over_b = a * _ray_tF_kernel(h / a)
    value = homogeneous - (T(2) / T(3)) * slm[3] +
        (c / T(3)) * (dA * boundary_F_over_b - sll[6])
    _finite(value) || throw(ArgumentError("normalized cell integral of B is nonfinite"))
    value
end

function _ray_KA_over_r(M::T, omega::Complex{T}, r::T, A::Complex{T},
        dA::Complex{T}) where {T<:AbstractFloat}
    two_M = T(2) * M
    f = one(T) - two_M / r
    two_i_omega = Complex{T}(zero(T), T(2)) * omega
    value = (r - two_M) * dA + (T(2) * f) * A + (two_i_omega * A) * r
    _finite(value) || throw(ArgumentError("scaled factored tail K_A/r is nonfinite"))
    value
end

function _ray_add_diagnostic!(destination::Vector{T}, j::Int, increment::T) where {T}
    value = destination[j - 1] + increment
    _finite(value) || throw(ArgumentError("accumulated quadrature diagnostic is nonfinite"))
    destination[j] = value
end

function _validate_ray_materialization_mode(mode::IRGCorrectorTailMode{N,T}) where {N,T}
    pair_N, pair_T = _corrector_order_precision(mode.pair.primary)
    N === pair_N && T === pair_T || throw(ArgumentError(
        "tail mode and corrector pair have different orders or precisions"))
    _validate_tail_boundary(mode.pair, mode.boundary)
    mode.gauge isa GHZIRG || throw(ArgumentError("corrector tail is not GHZ IRG"))
    ell, ellp1, c = _tail_angular_factors(T, mode.pair.primary.head.source.l)

    # Recheck the raw boundary without routing through any scalar corrector evaluator.
    pair, bc = mode.pair, mode.boundary
    primary, head = pair.primary, pair.primary.head
    src = head.source
    r = bc.rb
    A, dA = head.boundary.x0, head.boundary.dx0
    B, dB = primary.boundary.x0, primary.boundary.dx0
    parity = _barred_middle_parity(T, src.m)
    C = parity * conj(pair.partner.boundary.x0)
    dC = parity * conj(pair.partner.boundary.dx0)
    Sll, Sln = _Sll(head, r), _Sln(src, r)
    f = one(T) - T(2) * src.M / r
    two_i_omega = Complex{T}(zero(T), T(2)) * src.omega
    four_i_omega = Complex{T}(zero(T), T(4)) * src.omega
    lambda_over_r2 = (ell / r) * (ellp1 / r)
    L0 = -f * Sll +
        (two_i_omega + T(2) / r - T(2) * src.M / r / r) * dA +
        (four_i_omega / r + T(2) / r / r - lambda_over_r2) * A
    D, dD = B - C, dB - dC
    Q = Sln + L0 / T(2) + c * (dD + T(3) * D / r) / (T(2) * r)
    _finite(Q) || throw(ArgumentError("x_nn tail source Q must be finite at rb"))
    ell, ellp1, c
end

"""
    materialize_grid(mode::IRGCorrectorTailMode, radii; atol=0, rtol=sqrt(eps(T)))

Materialize the head, both middle partners, and the factored tail in one forward pass over a
strictly increasing finite grid whose first point is exactly `mode.boundary.rb`.  Every source is
integrated once per cell as a fixed vector of local moments.  No scalar corrector evaluator and no
prefix integral is called inside the pass.  A source jump is supported only when it lies on an
explicit cell boundary.
"""
function materialize_grid(mode::IRGCorrectorTailMode{N,T}, radii;
        atol=zero(T), rtol=sqrt(eps(T))) where {N,T<:AbstractFloat}
    _, _, c = _validate_ray_materialization_mode(mode)
    at, rt = _tolerances(T, atol, rtol)
    grid = _ray_grid(T, radii, mode.boundary.rb)
    count = length(grid)
    primary, partner = mode.pair.primary, mode.pair.partner
    phead, qhead = primary.head, partner.head
    psrc = phead.source
    parity = _barred_middle_parity(T, psrc.m)

    A = Vector{Complex{T}}(undef, count)
    dA = similar(A)
    B = similar(A)
    dB = similar(A)
    C = similar(A)
    dC = similar(A)
    Nvalues = similar(A)
    A[1], dA[1] = phead.boundary.x0, phead.boundary.dx0
    B[1], dB[1] = primary.boundary.x0, primary.boundary.dx0
    C[1] = parity * conj(partner.boundary.x0)
    dC[1] = parity * conj(partner.boundary.dx0)
    Nvalues[1] = mode.boundary.x0

    qA, qdA = qhead.boundary.x0, qhead.boundary.dx0
    qB, qdB = partner.boundary.x0, partner.boundary.dx0

    diagnostic_vectors = ntuple(_ -> zeros(T, count), 5)
    p_sll_diagnostics, p_slm_diagnostics, q_sll_diagnostics, q_slm_diagnostics,
        p_sln_diagnostics = diagnostic_vectors

    for j in 2:count
        a, b = grid[j - 1], grid[j]
        p_sll, ep_sll = _ray_sll_moments(phead, a, b, at, rt, "primary S_ll")
        p_slm, ep_slm = _ray_slm_moments(primary, a, b, at, rt, "primary S_lm")
        q_sll, eq_sll = _ray_sll_moments(qhead, a, b, at, rt, "partner S_ll")
        q_slm, eq_slm = _ray_slm_moments(partner, a, b, at, rt, "partner S_lm")
        I2_sln_over_b, ep_sln = _ray_sln_I2_over_b(psrc, a, b, at, rt)

        next_A, next_dA = _ray_head_step(A[j - 1], dA[j - 1], a, b, p_sll)
        next_qA, next_qdA = _ray_head_step(qA, qdA, a, b, q_sll)
        next_B, next_dB = _ray_middle_step(B[j - 1], dB[j - 1],
            dA[j - 1], a, b, p_sll, p_slm, c)
        next_qB, next_qdB = _ray_middle_step(qB, qdB, qdA,
            a, b, q_sll, q_slm, c)

        integral_A_over_b = _ray_integral_A_over_b(
            A[j - 1], dA[j - 1], a, b, p_sll)
        integral_B_over_b = _ray_integral_B_over_b(
            B[j - 1], dB[j - 1], dA[j - 1],
            a, b, p_sll, p_slm, c)
        integral_qB_over_b = _ray_integral_B_over_b(
            qB, qdB, qdA, a, b, q_sll, q_slm, c)
        previous_D = B[j - 1] - C[j - 1]
        next_C = parity * conj(next_qB)
        next_dC = parity * conj(next_qdB)
        next_D = next_B - next_C
        integral_D_over_b = integral_B_over_b - parity * conj(integral_qB_over_b)

        ratio = a / b
        KA_over_a = _ray_KA_over_r(psrc.M, psrc.omega, a,
            A[j - 1], dA[j - 1])
        KA_over_b = _ray_KA_over_r(psrc.M, psrc.omega, b, next_A, next_dA)
        lambda_integral_over_b = c * (c * integral_A_over_b)
        half_K_jump = KA_over_b / T(2) - ratio * (KA_over_a / T(2))
        half_c = c / T(2)
        half_c_D_jump = half_c * next_D - half_c * (ratio * previous_D)
        next_N = ratio * Nvalues[j - 1] + I2_sln_over_b + half_K_jump -
            lambda_integral_over_b + half_c_D_jump + c * integral_D_over_b
        _finite(next_N) || throw(ArgumentError("tail recurrence produced nonfinite cache data"))

        A[j], dA[j] = next_A, next_dA
        B[j], dB[j] = next_B, next_dB
        C[j], dC[j] = next_C, next_dC
        Nvalues[j] = next_N
        qA, qdA, qB, qdB = next_qA, next_qdA, next_qB, next_qdB

        # Exactly one local embedded-rule indicator is accumulated for each source stream.
        _ray_add_diagnostic!(p_sll_diagnostics, j, ep_sll)
        _ray_add_diagnostic!(p_slm_diagnostics, j, ep_slm)
        _ray_add_diagnostic!(q_sll_diagnostics, j, eq_sll)
        _ray_add_diagnostic!(q_slm_diagnostics, j, eq_slm)
        _ray_add_diagnostic!(p_sln_diagnostics, j, ep_sln)
    end

    provenance = _ray_moment_provenance(mode)
    diagnostics = RayMomentQuadratureDiagnostics{T}(p_sll_diagnostics,
        p_slm_diagnostics, q_sll_diagnostics, q_slm_diagnostics, p_sln_diagnostics)
    RayMomentCache{N,T}(provenance, grid, A, dA, B, dB, C, dC,
        Nvalues, diagnostics, at, rt, Val(:validated))
end

function _validate_ray_cache_header(cache::RayMomentCache{N,T}) where {N,T}
    p = cache.provenance
    N isa Integer && !(N isa Bool) && N > 0 ||
        throw(ArgumentError("cache perturbative order must be positive"))
    _finite(p.M) && p.M > zero(T) || throw(ArgumentError("cache M must be finite and positive"))
    _require_middle_spin_domain(p.l)
    -p.l <= p.m <= p.l || throw(ArgumentError("cache has an invalid mode key"))
    _finite(p.omega) || throw(ArgumentError("cache omega must be finite"))
    all(!=(UInt128(0)), (p.primary_source_fingerprint, p.partner_source_fingerprint,
        p.pair_fingerprint, p.tail_boundary_fingerprint)) ||
        throw(ArgumentError("cache fingerprints must be nonzero"))
    p.normalization isa EinsteinRHS ||
        throw(ArgumentError("cache normalization is not EinsteinRHS"))
    p.policy isa FiniteRayFreeIVP ||
        throw(ArgumentError("cache policy is not a finite-ray free IVP"))
    p.gauge isa GHZIRG || throw(ArgumentError("cache gauge is not GHZ IRG"))
    p.tetrad isa KinnersleyMinusPlus ||
        throw(ArgumentError("cache tetrad is not Kinnersley (-,+,+,+)"))
    p.slicing isa SchwarzschildOutgoing ||
        throw(ArgumentError("cache is not on outgoing Schwarzschild slices"))
    p.phase isa ExpMinusIomegaUPlusImPhi ||
        throw(ArgumentError("cache has the wrong separated phase"))
    _finite(p.rb) && p.rb >= T(2) * p.M ||
        throw(ArgumentError("cache rb must be finite and outside the horizon"))
    boundary_values = (p.tail_x0, p.primary_head_x0, p.primary_head_dx0,
        p.primary_middle_x0, p.primary_middle_dx0, p.partner_head_x0,
        p.partner_head_dx0, p.partner_middle_x0, p.partner_middle_dx0)
    all(_finite, boundary_values) || throw(ArgumentError("cache boundary data must be finite"))
    _tolerances(T, cache.atol, cache.rtol)
    nothing
end

function _validate_ray_moment_cache(cache::RayMomentCache{N,T}) where {N,T}
    _validate_ray_cache_header(cache)
    p = cache.provenance
    count = length(cache.radii)
    count > 0 || throw(ArgumentError("cache grid must not be empty"))
    cache.radii[1] == p.rb ||
        throw(ArgumentError("cache grid must begin exactly at rb"))
    all(_finite, cache.radii) || throw(ArgumentError("cache grid must be finite"))
    all(j -> cache.radii[j] > cache.radii[j - 1], 2:count) ||
        throw(ArgumentError("cache grid must be strictly increasing"))
    component_arrays = (cache.A, cache.dA, cache.B, cache.dB, cache.C, cache.dC,
        cache.N)
    all(values -> length(values) == count, component_arrays) ||
        throw(ArgumentError("cache component lengths differ from the grid"))
    all(values -> all(_finite, values), component_arrays) ||
        throw(ArgumentError("cache components must be finite"))
    parity = _barred_middle_parity(T, p.m)
    (cache.A[1], cache.dA[1], cache.B[1], cache.dB[1], cache.C[1], cache.dC[1],
        cache.N[1]) == (p.primary_head_x0, p.primary_head_dx0,
        p.primary_middle_x0, p.primary_middle_dx0,
        parity * conj(p.partner_middle_x0), parity * conj(p.partner_middle_dx0),
        p.tail_x0) || throw(ArgumentError("cache values disagree with boundary provenance"))

    diagnostic_arrays = (cache.quadrature_diagnostics.primary_sll,
        cache.quadrature_diagnostics.primary_slm,
        cache.quadrature_diagnostics.partner_sll,
        cache.quadrature_diagnostics.partner_slm,
        cache.quadrature_diagnostics.primary_sln)
    all(values -> length(values) == count, diagnostic_arrays) ||
        throw(ArgumentError("cache diagnostic lengths differ from the grid"))
    all(values -> all(value -> _finite(value) && value >= zero(T), values),
        diagnostic_arrays) ||
        throw(ArgumentError("cache diagnostics must be finite and nonnegative"))
    all(values -> values[1] == zero(T), diagnostic_arrays) ||
        throw(ArgumentError("cache accumulated diagnostics must start at zero"))
    all(values -> all(j -> values[j] >= values[j - 1], 2:count),
        diagnostic_arrays) ||
        throw(ArgumentError("cache accumulated diagnostics must be nondecreasing"))
    nothing
end

function _ray_cache_index(cache::RayMomentCache{N,T}, r) where {N,T}
    _validate_ray_cache_header(cache)
    p = cache.provenance
    count = length(cache.radii)
    count > 0 || throw(ArgumentError("cache grid must not be empty"))
    component_arrays = (cache.A, cache.dA, cache.B, cache.dB, cache.C, cache.dC,
        cache.N)
    diagnostic_arrays = (cache.quadrature_diagnostics.primary_sll,
        cache.quadrature_diagnostics.primary_slm,
        cache.quadrature_diagnostics.partner_sll,
        cache.quadrature_diagnostics.partner_slm,
        cache.quadrature_diagnostics.primary_sln)
    all(values -> length(values) == count, component_arrays) ||
        throw(ArgumentError("cache component lengths differ from the grid"))
    all(values -> length(values) == count, diagnostic_arrays) ||
        throw(ArgumentError("cache diagnostic lengths differ from the grid"))
    cache.radii[1] == p.rb ||
        throw(ArgumentError("cache grid must begin exactly at rb"))
    radius = try
        _checked_real(T, r, "cache lookup radius")
    catch error
        error isa ArgumentError && rethrow()
        throw(ArgumentError("cache lookup radius is not representable in $T"))
    end
    _finite(radius) || throw(ArgumentError("cache lookup radius must be finite"))
    index = searchsortedfirst(cache.radii, radius)
    index <= length(cache.radii) && cache.radii[index] == radius ||
        throw(ArgumentError("cache lookup requires an exact stored grid radius; interpolation is disabled"))
    index > 1 && !(cache.radii[index] > cache.radii[index - 1]) &&
        throw(ArgumentError("cache grid must be strictly increasing"))
    index < count && !(cache.radii[index + 1] > cache.radii[index]) &&
        throw(ArgumentError("cache grid must be strictly increasing"))
    values_at_index = (cache.A[index], cache.dA[index], cache.B[index], cache.dB[index],
        cache.C[index], cache.dC[index], cache.N[index])
    all(_finite, values_at_index) || throw(ArgumentError("cache components must be finite"))
    diagnostics_at_index = (cache.quadrature_diagnostics.primary_sll[index],
        cache.quadrature_diagnostics.primary_slm[index],
        cache.quadrature_diagnostics.partner_sll[index],
        cache.quadrature_diagnostics.partner_slm[index],
        cache.quadrature_diagnostics.primary_sln[index])
    all(value -> _finite(value) && value >= zero(T), diagnostics_at_index) ||
        throw(ArgumentError("cache diagnostics must be finite and nonnegative"))
    index
end

function ray_quadrature_diagnostics(cache::RayMomentCache{N,T}, r) where {N,T}
    index = _ray_cache_index(cache, r)
    diagnostics = cache.quadrature_diagnostics
    RayMomentDiagnosticPoint{T}(diagnostics.primary_sll[index],
        diagnostics.primary_slm[index], diagnostics.partner_sll[index],
        diagnostics.partner_slm[index], diagnostics.primary_sln[index])
end

function ray_moment_components(cache::RayMomentCache{N,T}, r) where {N,T}
    index = _ray_cache_index(cache, r)
    diagnostics = cache.quadrature_diagnostics
    point_diagnostic = RayMomentDiagnosticPoint{T}(diagnostics.primary_sll[index],
        diagnostics.primary_slm[index], diagnostics.partner_sll[index],
        diagnostics.partner_slm[index], diagnostics.primary_sln[index])
    RayMomentComponents{N,T}(cache.provenance, index, cache.radii[index],
        cache.A[index], cache.dA[index], cache.B[index], cache.dB[index],
        cache.C[index], cache.dC[index], cache.N[index], point_diagnostic,
        cache.atol, cache.rtol)
end

x_mmb(cache::RayMomentCache, r) = cache.A[_ray_cache_index(cache, r)]
dr_x_mmb(cache::RayMomentCache, r) = cache.dA[_ray_cache_index(cache, r)]
x_nm(cache::RayMomentCache, r) = cache.B[_ray_cache_index(cache, r)]
dr_x_nm(cache::RayMomentCache, r) = cache.dB[_ray_cache_index(cache, r)]
x_nmb(cache::RayMomentCache, r) = cache.C[_ray_cache_index(cache, r)]
dr_x_nmb(cache::RayMomentCache, r) = cache.dC[_ray_cache_index(cache, r)]
x_nn(cache::RayMomentCache, r) = cache.N[_ray_cache_index(cache, r)]
corrector_tetrad_components(cache::RayMomentCache, r) =
    ray_moment_components(cache, r)
