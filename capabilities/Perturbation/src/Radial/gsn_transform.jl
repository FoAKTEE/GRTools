# :GSN backend — Sasaki-Nakamura transformation functions α, β (+ analytic r-derivatives) and the
# η function whose 1/r-expansion coefficients c₀..c₋₄ define the SN↔Teukolsky map (GSN paper
# arXiv:2306.16469 §II; Mathematica derivation 10.5281/zenodo.8080241). Bodies below are the
# MECHANICALLY EXTRACTED upstream closed forms (src/Homogeneous/Transformation.jl, v0.8.5
# commit 982b087) with the renames `I→im`, `*→_gsn_*` (scratchpad/gsn-port/extract_gsn.py);
# they are generic rational algebra in (a, ω, λ, r) — BigFloat-clean as-is (audit 2026-07-02).
# The η/η′/η″ evaluators are re-authored over a SINGLE coefficient fetch (upstream recomputed
# all five c's in each of eta/eta_prime/eta_primeprime — GSN_PLAN §3e hot-path fix).

function _gsn_alpha(s::Int, m::Int, a, omega, lambda, r)
    if s == 0
        return 1
    elseif s == +1
        return begin
            (sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*((-im)*a^3*m - im*a*m*r^2 + im*a^4*omega + 
            r^3*(1 + im*r*omega) + a^2*(-2 + r + 2*im*r^2*omega)))/(r^2*sqrt(a^2 + r^2))
        end
    elseif s == -1
        return begin
            (sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*(-r - (im*(a^2 + r^2)*((-a)*m + (a^2 + r^2)*omega))/
            (a^2 + (-2 + r)*r)))/(r^2*sqrt(a^2 + r^2))
        end
    elseif s == +2
        return begin
            (1/(r^2*(a^2 + (-2 + r)*r)))*(4*a^3*m*r*(im + r*omega) + 2*a*m*r^2*(im - 3*im*r + 2*r^2*omega) - 
            2*a^4*(-3 + 2*im*r*omega + r^2*omega^2) + r^3*(-2*lambda + r*(2 + lambda + 10*im*omega) - 2*r^3*omega^2) - 
            a^2*r*(8 + 2*m^2*r - r*lambda + 2*im*r*omega + 4*im*r^2*omega + 4*r^3*omega^2))
        end
    elseif s == -2
        return begin
            (1/(r^2*(a^2 + (-2 + r)*r)))*(4*a^3*m*r*(-im + r*omega) + 2*a*m*r^2*(3*im - im*r + 2*r^2*omega) + 
            a^4*(6 + 4*im*r*omega - 2*r^2*omega^2) + a^2*r*(-24 + r*(12 - 2*m^2 + lambda - 6*im*omega) + 12*im*r^2*omega - 4*r^3*omega^2) + 
            r^2*(24 - 2*r*(12 + lambda) + r^2*(6 + lambda - 18*im*omega) + 8*im*r^3*omega - 2*r^4*omega^2))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_alpha_prime(s::Int, m::Int, a, omega, lambda, r)
    if s == 0
        return 0
    elseif s == +1
        return begin
            -((im*sqrt(a^2 + r^2)*(-2*a^5*m + a^3*m*(5 - 3*r)*r - a*m*(-1 + r)*r^3 + 2*a^6*omega - 
            r^4*(im - 3*r*omega + r^2*omega) - a^2*r*(10*im - 9*im*r + r^2*(im + 2*omega)) + 
            a^4*(4*im + 3*r^2*omega - r*(im + 5*omega))))/(r^3*(a^2 + (-2 + r)*r)^2*
            sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))))
        end
    elseif s == -1
        return begin
            (im*(a^2 + r^2)^(3/2)*(-2*a^5*m + a^3*m*(5 - 3*r)*r - a*m*(-1 + r)*r^3 + 2*a^6*omega + 
            r^3*(-2*im + im*r + 3*r^2*omega - r^3*omega) + a^4*r*(-im + (-5 + 3*r)*omega) + 
            a^2*r^2*(3*im - r*(im + 2*omega))))/(r^3*((a^2 + (-2 + r)*r)*(a^2 + r^2))^(3/2))
        end
    elseif s == +2
        return begin
            -((1/(r^3*(a^2 + (-2 + r)*r)^2))*(2*(2*im*a^5*m*r + a^3*m*r^2*(-8*im + r*(9*im - 4*omega)) + 
            a^6*(6 - 2*im*r*omega) + 2*r^5*(1 + 5*im*omega + (-3 + r)*r^2*omega^2) + 
            2*a^4*r*(-11 + r*(6 + 4*im*omega) + r^3*omega^2 + r^2*omega*(-2*im + omega)) + 
            2*a^2*r^2*(8 + r*(-6 + m^2 + im*omega) - r^2*(1 + m^2 + 6*im*omega) + 2*r^4*omega^2 - 
            r^3*omega*(im + 2*omega)) + a*m*r^3*(-2*im + 2*im*r + r^2*(-3*im + 4*omega)))))
        end
    elseif s == -2
        return begin
            -((1/(r^3*(a^2 + (-2 + r)*r)^2))*(2*(-2*im*a^5*m*r + a^6*(6 + 2*im*r*omega) + 
            2*a^4*r*(-15 + r*(6 - 4*im*omega) + r^2*omega^2 + r^3*omega^2) + 
            2*a^2*r^2*(24 + r*(-18 + m^2 + 3*im*omega) + r^2*(3 - m^2 + 6*im*omega) + 2*r^4*omega^2 - r^3*omega*(3*im + 2*omega)) + 
            2*r^3*(-12 + 12*r + r^2*(-3 - 9*im*omega) + 8*im*r^3*omega + r^5*omega^2 - r^4*omega*(2*im + 3*omega)) + 
            a*m*r^3*(-6*im + 6*im*r + r^2*(-im + 4*omega)) + a^3*m*r^2*(8*im - r*(5*im + 4*omega)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_beta(s::Int, m::Int, a, omega, lambda, r)
    if s == 0
        return 0
    elseif s == +1
        return begin
            ((a^2 + r^2)/(a^2 + (-2 + r)*r))^(3/2)/(r^2*sqrt(a^2 + r^2))
        end
    elseif s == -1
        return begin
            (sqrt(a^2 + r^2)*sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2)))/r^2
        end
    elseif s == +2
        return begin
            (-2*im*a*m*r + a^2*(-4 + 2*im*r*omega) + 2*r*(3 - r + im*r^2*omega))/(r*(a^2 + (-2 + r)*r)^3)
        end
    elseif s == -2
        return begin
            (2*(a^2 + (-2 + r)*r)*(im*a*m*r + a^2*(-2 - im*r*omega) + r*(3 - r - im*r^2*omega)))/r
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_beta_prime(s::Int, m::Int, a, omega, lambda, r)
    if s == 0
        return 0
    elseif s == +1
        return begin
            -((sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*(2*a^4 + 3*(-1 + r)*r^3 + a^2*r*(-7 + 5*r)))/
            (r^3*(a^2 + (-2 + r)*r)^2*sqrt(a^2 + r^2)))
        end
    elseif s == -1
        return begin
            (-2*a^6 - 3*a^4*(-1 + r)*r + 2*a^2*r^3 + (-1 + r)*r^5)/
            (r^3*sqrt(a^2 + r^2)*sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2)))
        end
    elseif s == +2
        return begin
            (1/(r^2*(a^2 + (-2 + r)*r)^4))*(2*(2*a^4 + 6*im*a*m*(-1 + r)*r^2 + 
            a^2*r*(-16 + r*(13 + 6*im*omega) - 4*im*r^2*omega) + 
            r^2*(18 - 22*r + r^2*(5 + 2*im*omega) - 4*im*r^3*omega)))
        end
    elseif s == -2
        return begin
            2*(-6 + 2*im*a*m*(-1 + r) + (2*a^4)/r^2 + 10*r + r^2*(-3 + 6*im*omega) - 4*im*r^3*omega + 
            a^2*(-3 - 2*im*(-1 + 2*r)*omega))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_eta_coefficient(s::Int, m::Int, a, omega, lambda, order)
    if s == 0
        if order == 0
            return 1
        elseif order == -1
            return 0
        elseif order == -2
            return 0
        elseif order == -3
            return 0
        elseif order == -4
            return 0
        else
            return 0
        end
    elseif s == +1
        if order == 0
            return -2 - lambda
        elseif order == -1
            return 2*im*a*m
        elseif order == -2
            return -3*a^2 - 2*a^2*lambda
        elseif order == -3
            return -2*a^2 + 2*im*a^3*m
        elseif order == -4
            return -a^4 - a^4*lambda
        else
            return 0
        end
    elseif s == -1
        if order == 0
            return -lambda
        elseif order == -1
            return -2*im*a*m
        elseif order == -2
            return a^2 - 2*a^2*lambda
        elseif order == -3
            return -2*a^2 - 2*im*a^3*m
        elseif order == -4
            return a^4 - a^4*lambda
        else
            return 0
        end
    elseif s == +2
        if order == 0
            return 24 + 10*lambda + lambda^2 + 12*im*omega + 12*a*m*omega - 12*a^2*omega^2
        elseif order == -1
            return -32*im*a*m - 8*im*a*m*lambda + 8*im*a^2*omega + 8*im*a^2*lambda*omega
        elseif order == -2
            return 12*a^2 - 24*im*a*m - 24*a^2*m^2 + 24*im*a^2*omega + 48*a^3*m*omega - 24*a^4*omega^2
        elseif order == -3
            return -24*a^2 + 24*im*a^3*m - 24*im*a^4*omega
        elseif order == -4
            return 12*a^4
        else
            return 0
        end
    elseif s == -2
        if order == 0
            return 2*lambda + lambda^2 - 12*im*omega + 12*a*m*omega - 12*a^2*omega^2
        elseif order == -1
            return 8*im*a*m*lambda + 24*im*a^2*omega - 8*im*a^2*lambda*omega
        elseif order == -2
            return 12*a^2 + 24*im*a*m - 24*a^2*m^2 - 24*im*a^2*omega + 48*a^3*m*omega - 24*a^4*omega^2
        elseif order == -3
            return -24*a^2 - 24*im*a^3*m + 24*im*a^4*omega
        elseif order == -4
            return 12*a^4
        else
            return 0
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


"""
    _gsn_eta_coeffs(s, m, a, omega, lambda) -> (c0, c1, c2, c3, c4)

The five 1/r-expansion coefficients of η (η = c0 + c1/r + c2/r² + c3/r³ + c4/r⁴), fetched once.
[GSN Transformation.jl eta_coefficient, orders 0,−1,−2,−3,−4]
"""
_gsn_eta_coeffs(s::Int, m::Int, a, omega, lambda) = (
    _gsn_eta_coefficient(s, m, a, omega, lambda, 0),
    _gsn_eta_coefficient(s, m, a, omega, lambda, -1),
    _gsn_eta_coefficient(s, m, a, omega, lambda, -2),
    _gsn_eta_coefficient(s, m, a, omega, lambda, -3),
    _gsn_eta_coefficient(s, m, a, omega, lambda, -4),
)

"η(r) = c0 + c1/r + c2/r² + c3/r³ + c4/r⁴. [GSN Transformation.jl eta]"
function _gsn_eta(s::Int, m::Int, a, omega, lambda, r)
    c0, c1, c2, c3, c4 = _gsn_eta_coeffs(s, m, a, omega, lambda)
    return c0 + c1/r + c2/r^2 + c3/r^3 + c4/r^4
end

"η′(r) = −c1/r² − 2c2/r³ − 3c3/r⁴ − 4c4/r⁵. [GSN Transformation.jl eta_prime]"
function _gsn_eta_prime(s::Int, m::Int, a, omega, lambda, r)
    _, c1, c2, c3, c4 = _gsn_eta_coeffs(s, m, a, omega, lambda)
    return -c1/r^2 - 2*c2/r^3 - 3*c3/r^4 - 4*c4/r^5
end

"η″(r) = 2c1/r³ + 6c2/r⁴ + 12c3/r⁵ + 20c4/r⁶. [GSN Transformation.jl eta_primeprime]"
function _gsn_eta_primeprime(s::Int, m::Int, a, omega, lambda, r)
    _, c1, c2, c3, c4 = _gsn_eta_coeffs(s, m, a, omega, lambda)
    return 2*c1/r^3 + 6*c2/r^4 + 12*c3/r^5 + 20*c4/r^6
end
