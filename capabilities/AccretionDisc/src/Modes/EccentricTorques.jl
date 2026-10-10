"""
    pl00_tauL_inv(e, h_p, b)

PL00 angular-momentum damping rate `τ_L⁻¹/τ₀⁻¹` as printed in FR24 v1
(arXiv:2407.20398), equation (52). Inputs are physical eccentricity `e ≥ 0`,
aspect ratio `h_p > 0`, and dimensionless softening `b = ε/H > 0`.
The scaled eccentricity is `e/h_p`; the rate reverses sign at `e/h_p = 1.1`.
Multiply by `τ₀⁻¹` to restore a dimensional rate. This is not `τ_a⁻¹`.
The printed coefficient 7.33 is retained despite the unresolved sub-percent
normalisation offset from PL00 (31). Float64 is the default; BigFloat inputs
carry their working precision through the fit.
"""
function pl00_tauL_inv(e::Real, h_p::Real, b::Real)
    T = float(promote_type(typeof(e), typeof(h_p), typeof(b)))
    et = T(e) / T(h_p)
    T(733//100) * (T(b)/T(1//2))^T(-7//4) * T(h_p) *
        (1 - (et/T(11//10))^4) / (1 + (et/T(13//10))^5)
end

"""
    pl00_taue_inv(e, h_p, b)

PL00 eccentricity damping rate `τ_e⁻¹/τ₀⁻¹`, FR24 v1 equation (53), with
the printed coefficient 4.26. Inputs and precision follow [`pl00_tauL_inv`](@ref).
The `e = 0` value is the continuous limit of the fractional damping rate.
The unresolved PL00 (32) normalisation offset is not absorbed into this fit.
"""
function pl00_taue_inv(e::Real, h_p::Real, b::Real)
    T = float(promote_type(typeof(e), typeof(h_p), typeof(b)))
    et = T(e) / T(h_p)
    T(213//50) * (T(b)/T(1//2))^T(-5//2) / (1 + et^3/4) / T(h_p)
end

"""
    cn06_tauL_inv(e, h_p, p)

CN06 angular-momentum damping rate `τ_L⁻¹/τ₀⁻¹` as printed in FR24 v1
equation (54), with surface density `Σ ∝ r^(-p)`. Inputs `e ≥ 0` and `h_p > 0`
are physical eccentricity and aspect ratio; the fit uses `e/h_p` internally.
There is no softening argument in this printed prescription. Precision follows
the input types, as in [`pl00_tauL_inv`](@ref). This transcribes FR24's CN06
prescription; it does not independently establish the original CN06 provenance.
"""
function cn06_tauL_inv(e::Real, h_p::Real, p::Real)
    T = float(promote_type(typeof(e), typeof(h_p), typeof(p)))
    et = T(e) / T(h_p)
    (T(27//10) + T(11//10)*T(p))/2 * T(h_p) *
        (1 - (et/T(11//10))^4) / (1 + (et/T(13//10))^5)
end

"""
    cn06_taue_inv(e, h_p, Q_e)

CN06 eccentricity damping rate `τ_e⁻¹/τ₀⁻¹` as printed in FR24 v1 equation
(55). `e ≥ 0`, `h_p > 0`, and `Q_e > 0`; Fig. 11 uses `Q_e = 0.1, 0.2, 1`.
This is the uncapped printed fit, including its continuous `e = 0` limit.
Precision and the provenance scope follow [`cn06_tauL_inv`](@ref).
"""
function cn06_taue_inv(e::Real, h_p::Real, Q_e::Real)
    T = float(promote_type(typeof(e), typeof(h_p), typeof(Q_e)))
    et = T(e) / T(h_p)
    T(39//50) / T(Q_e) / (1 + et^3/4) / T(h_p)
end
