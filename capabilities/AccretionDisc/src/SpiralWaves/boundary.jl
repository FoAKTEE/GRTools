# Arbitrary-frequency radiation boundaries, G26 (66)-(68).
# Include after solver.jl. The real-frequency, unregularized edge operator and its
# frequency derivative are supplied by the caller; neither m nor a pattern radius
# is needed. Keep the legacy scalar-Doppler methods and circular operation order.

"""
    wkb_split(C, dC_dω::AbstractMatrix; outward, bc=:outgoing)
    outgoing_branch(C, dC_dω::AbstractMatrix; outward, bc=:outgoing)
    wkb_amplitudes(C, dC_dω::AbstractMatrix, y; outward, bc=:outgoing)

Arbitrary-frequency forms of the existing boundary and amplitude functions.
`dC_dω` is the derivative at fixed radius, m and disc, using the same state variables
as `C`. Supply the REAL, unregularized edge operator (complex storage is allowed).
An absorbing/Landau layer is not a real-frequency WKB dispersion relation.

For distinct roots λ of F=λ²-tr(C)λ+det(C)=0,
λ_ω=(tr(C)_ω λ-det(C)_ω)/(2λ-tr(C)). With exp(λr-iωt+imφ),
k_r=Im(λ), so v_g=1/Im(λ_ω). Select n v_g>0, n=+1 at the outer
edge and -1 at the inner edge. A common real trace (amplitude gradient) cannot
turn an oscillatory conjugate pair into an evanescent pair. Real roots instead
require n Re(λ)<0: actual decay into the exterior.

This is a frozen-coefficient WKB closure, not an exact transparent boundary for
a varying disc. Coalesced roots, zero group slope, or the absence of one outgoing
and one incoming/decaying and growing branch are rejected. Move such an edge;
exact turning points require an Airy/matching treatment. Finite separated roots
arbitrarily near a turning point are allowed, with their usual conditioning loss.
`:paper` retains the literal Im(λ)<0 comparison for sign reconciliation only.
"""
function wkb_split(C::AbstractMatrix, dC_dω::AbstractMatrix;
                   outward::Bool, bc::Symbol = :outgoing)
    size(C) == size(dC_dω) == (2, 2) || throw(DimensionMismatch("boundary matrices must be 2×2"))
    all(isfinite, C) && all(isfinite, dC_dω) ||
        throw(ArgumentError("boundary matrices must be finite"))
    all(isreal, C) && all(isreal, dC_dω) ||
        throw(ArgumentError("radiation boundaries need a real, unregularized edge operator"))
    k1, k2, Y1, Y2 = eigen2(complex.(C))
    k1 != k2 || throw(DomainError(k1, "coalesced boundary roots: move the edge away from the turning point"))
    if iszero(imag(k1)) && iszero(imag(k2))
        directions = (-real(k1), -real(k2))
    else
        trω = dC_dω[1, 1] + dC_dω[2, 2]
        detω = dC_dω[1, 1] * C[2, 2] + C[1, 1] * dC_dω[2, 2] -
               dC_dω[1, 2] * C[2, 1] - C[1, 2] * dC_dω[2, 1]
        slope1 = imag((trω * k1 - detω) / (k1 - k2))
        slope2 = imag((trω * k2 - detω) / (k2 - k1))
        # v_g and dk_r/dω have the same sign. Avoid reciprocal overflow near
        # zero frequency; only the orientation is needed for branch selection.
        directions = (slope1, slope2)
    end
    first_away = _away_is_first(k1, k2, directions; outward = outward, bc = bc)
    first_away ? (k1, Y1, k2, Y2) : (k2, Y2, k1, Y1)
end

function _away_is_first(k1, k2, directions::Tuple; outward::Bool, bc::Symbol)
    bc === :paper && return imag(k1) ≤ imag(k2)
    bc === :outgoing || throw(ArgumentError("boundary bc must be :outgoing or :paper"))
    d1, d2 = directions
    all(isfinite, directions) && d1 * sign(d2) < zero(d1) ||
        throw(DomainError(directions, "edge needs one outgoing/decaying and one incoming/growing branch"))
    outward ? d1 > zero(d1) : d1 < zero(d1)
end

"""
    radiation_impedance_row(ka, Ya, kt, Yt, dYa) -> row

First adiabatic correction for a HOMOGENEOUS outgoing wave in a slowly varying
eigenbasis. For P=[Ya Yt], N=P⁻¹P′, y=a Ya+b Yt, the Riccati equation for b/a
gives b/a=N₂₁/(kt-ka)+O(ε²), ε measuring basis variation relative to the eigenvalue
gap. Thus row=(−Ya₂,Ya₁)+β(−Yt₂,Yt₁), β=N₂₁/(kt-ka).
The row is covariant under any differentiable nonzero scaling of either vector.
Supply the radial derivative of Ya in the same gauge as Ya. Nonzero edge forcing
requires subtracting a particular solution before applying this homogeneous row.

The legacy `_impedance_row` / `bc=:wkb1` instead freezes b′=0 and is retained for
circular compatibility. Its kt denominator is not the oscillatory WKB gap kt-ka.
Neither expansion is uniform at a turning point or corotation.
"""
function radiation_impedance_row(ka, Ya, kt, Yt, dYa)
    D = Ya[1] * Yt[2] - Ya[2] * Yt[1]
    gap = kt - ka
    if iszero(D) || iszero(gap)
        throw(DomainError(gap, "impedance needs a separated eigenbasis"))
    end
    Wa = (-Ya[2], Ya[1]); Wt = (-Yt[2], Yt[1])
    β = (Wa[1] * dYa[1] + Wa[2] * dYa[2]) / (D * gap)
    (Wa[1] + β * Wt[1], Wa[2] + β * Wt[2])
end
