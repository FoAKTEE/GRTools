# :GSN backend — Teukolsky↔SN amplitude conversion factors (the ratio Teukolsky-formalism
# coefficient / SN-formalism coefficient, per boundary and spin: Btrans/Binc/Bref for In,
# Ctrans/Cinc/Cref for Up) and the Teukolsky-Starobinsky constants |C|² (s=±2), |B|² (s=±1)
# (GSN paper arXiv:2306.16469 §III; Mathematica derivation 10.5281/zenodo.8080241). Bodies are
# the MECHANICALLY EXTRACTED upstream closed forms (src/Homogeneous/ConversionFactors.jl,
# v0.8.5 commit 982b087) with rename `I→im` and the audited precision fix `sqrt(2) → _sqrt2 =
# sqrt(2*one(real(a)))` in all 10 affected branches (GSN_PLAN §3c-5; upstream evaluated the
# Float64 sqrt(2) into otherwise-BigFloat expressions).

function _gsn_Ctrans(s::Int, m::Int, a, omega, lambda)
    _sqrt2 = sqrt(2*one(real(a)))   # typed: BigFloat-exact (upstream: Float64 sqrt(2))
    if s == 0
        return 1
    elseif s == +1
        inv = 2*im*omega
        return 1/inv
    elseif s == -1
        return (2*im*omega)/lambda
    elseif s == +2
        inv = -(4*omega^2)
        return 1/inv
    elseif s == -2
        return begin
            (4*omega^2)/(-lambda*(2+lambda)+12*im*omega-12*a*m*omega+12*(a*omega)^2)
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_Binc(s::Int, m::Int, a, omega, lambda)
    _sqrt2 = sqrt(2*one(real(a)))   # typed: BigFloat-exact (upstream: Float64 sqrt(2))
    if s == 0
        return 1
    elseif s == +1
        return -((2*im*omega)/(2 + lambda))
    elseif s == -1
        inv = -(2*im*omega)
        return 1/inv
    elseif s == +2
        return begin
            (4*omega^2)/(-24 - 10*lambda - lambda^2 - 12*im*omega - 12*a*m*omega + 12*a^2*omega^2)
        end
    elseif s == -2
        inv = -(4*omega^2)
        return 1/inv
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_Btrans(s::Int, m::Int, a, omega, lambda)
    _sqrt2 = sqrt(2*one(real(a)))   # typed: BigFloat-exact (upstream: Float64 sqrt(2))
    if s == 0
        inv = _sqrt2*sqrt(1 + sqrt(1 - a^2))
        return 1/inv
    elseif s == +1
        return begin
            (_sqrt2*sqrt(1 + sqrt(1 - a^2))*((-a)*m + 2*(1 + sqrt(1 - a^2))*omega))/
            (2*a*m + 2*im*(2 + lambda))
        end
    elseif s == -1
        inv = begin
            (1/(1 + sqrt(1 - a^2))^(5/2))*(4*_sqrt2*(2 - 2*a^2 + 2*sqrt(1 - a^2) - 
            a^2*sqrt(1 - a^2) + 2*im*a*m - im*a^3*m + 2*im*a*sqrt(1 - a^2)*m - 8*im*omega + 
            6*im*a^2*omega - 8*im*sqrt(1 - a^2)*omega + 2*im*a^2*sqrt(1 - a^2)*omega))
        end
        return 1/inv
    elseif s == +2
        return begin
            (2*_sqrt2*(1 + sqrt(1 - a^2))^(3/2)*(2*a*(1 + sqrt(1 - a^2))*m*(1 + 8*im*omega) - 
            8*im*(1 + sqrt(1 - a^2))*omega*(-im + 4*omega) + 
            a^3*m*(-2 - sqrt(1 - a^2) - 4*im*(3 + sqrt(1 - a^2))*omega) + 
            im*a^4*(m^2 + 2*(im - 2*omega)*omega) + 2*a^2*((-im)*(1 + sqrt(1 - a^2))*m^2 + 
            omega*(5 + 3*sqrt(1 - a^2) + 8*im*(2 + sqrt(1 - a^2))*omega))))/
            (16*a*(1 + sqrt(1 - a^2))*m*(11 + 2*lambda + 6*im*omega) + 8*im*(1 + sqrt(1 - a^2))*
            (24 + lambda*(10 + lambda) + 12*im*omega) - 36*im*a^5*m*omega + 12*im*a^6*omega^2 - 
            8*a^3*m*(18 + 7*sqrt(1 - a^2) + (3 + sqrt(1 - a^2))*lambda - 6*im*sqrt(1 - a^2)*omega) + 
            a^4*(24*im*m^2 + im*(4 + lambda)*(6 + lambda) + 32*sqrt(1 - a^2)*omega + 8*sqrt(1 - a^2)*lambda*omega + 
            12*(5 + 2*lambda)*omega + 48*im*omega^2) + 4*a^2*(-24*im*(2 + sqrt(1 - a^2)) - 
            12*im*(1 + sqrt(1 - a^2))*m^2 - 10*im*(2 + sqrt(1 - a^2))*lambda - 
            im*(2 + sqrt(1 - a^2))*lambda^2 - 8*(1 + sqrt(1 - a^2))*lambda*omega + 
            4*omega*(1 - 2*sqrt(1 - a^2) - 6*im*(1 + sqrt(1 - a^2))*omega)))
        end
    elseif s == -2
        inv = begin
            -((1/(1 + sqrt(1 - a^2))^(3/2))*(4*(_sqrt2*(-4 + 6*a^2 - 2*a^4 - 4*sqrt(1 - a^2) + 
            4*a^2*sqrt(1 - a^2) - 6*im*a*m + 6*im*a^3*m - 6*im*a*sqrt(1 - a^2)*m + 
            3*im*a^3*sqrt(1 - a^2)*m + 2*a^2*m^2 - a^4*m^2 + 2*a^2*sqrt(1 - a^2)*m^2 + 
            24*im*omega - 30*im*a^2*omega + 6*im*a^4*omega + 24*im*sqrt(1 - a^2)*omega - 
            18*im*a^2*sqrt(1 - a^2)*omega - 16*a*m*omega + 12*a^3*m*omega - 16*a*sqrt(1 - a^2)*m*omega + 
            4*a^3*sqrt(1 - a^2)*m*omega + 32*omega^2 - 32*a^2*omega^2 + 4*a^4*omega^2 + 
            32*sqrt(1 - a^2)*omega^2 - 16*a^2*sqrt(1 - a^2)*omega^2))))
        end
        return 1/inv
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_Cinc(s::Int, m::Int, a, omega, lambda)
    _sqrt2 = sqrt(2*one(real(a)))   # typed: BigFloat-exact (upstream: Float64 sqrt(2))
    if s == 0
        return 1/(_sqrt2*sqrt(1 + sqrt(1 - a^2)))
    elseif s == +1
        inv = begin
            (4*_sqrt2*((-im)*a*(1 + sqrt(1 - a^2))*m + a^2*(-1 - 2*im*omega) + 
            (1 + sqrt(1 - a^2))*(1 + 4*im*omega)))/(1 + sqrt(1 - a^2))^(3/2)
        end
        return 1/inv
    elseif s == -1
        return begin
            (-4*a*m + 3*a^3*m - 4*a*sqrt(1 - a^2)*m + a^3*sqrt(1 - a^2)*m + 16*omega - 16*a^2*omega + 
            2*a^4*omega + 16*sqrt(1 - a^2)*omega - 8*a^2*sqrt(1 - a^2)*omega)/
            (_sqrt2*sqrt(1 + sqrt(1 - a^2))*(2 - a^2 + 2*sqrt(1 - a^2))*(a*m - im*lambda))
        end
    elseif s == +2
        inv = begin
            -((1/(1 + sqrt(1 - a^2))^(3/2))*(4*_sqrt2*(-2*a*(1 + sqrt(1 - a^2))*m*(-3*im + 8*omega) + 
            a^3*m*(-3*im*(2 + sqrt(1 - a^2)) + 4*(3 + sqrt(1 - a^2))*omega) - 
            a^4*(2 + m^2 + 6*im*omega - 4*omega^2) + 4*(1 + sqrt(1 - a^2))*(-1 - 6*im*omega + 8*omega^2) + 
            2*a^2*(3 + 2*sqrt(1 - a^2) + (1 + sqrt(1 - a^2))*m^2 + 
            omega*(3*im*(5 + 3*sqrt(1 - a^2)) - 8*(2 + sqrt(1 - a^2))*omega)))))
        end
        return 1/inv
    elseif s == -2
        return begin
            -((2*_sqrt2*(1 + sqrt(1 - a^2))^(3/2)*(-2*im*a*m + 2*im*a^3*m - 2*im*a*sqrt(1 - a^2)*m + 
            im*a^3*sqrt(1 - a^2)*m + 2*a^2*m^2 - a^4*m^2 + 2*a^2*sqrt(1 - a^2)*m^2 + 8*im*omega - 
            10*im*a^2*omega + 2*im*a^4*omega + 8*im*sqrt(1 - a^2)*omega - 6*im*a^2*sqrt(1 - a^2)*omega - 
            16*a*m*omega + 12*a^3*m*omega - 16*a*sqrt(1 - a^2)*m*omega + 4*a^3*sqrt(1 - a^2)*m*omega + 
            32*omega^2 - 32*a^2*omega^2 + 4*a^4*omega^2 + 32*sqrt(1 - a^2)*omega^2 - 
            16*a^2*sqrt(1 - a^2)*omega^2))/((1 + sqrt(1 - a^2))^4*(lambda*(2 + lambda) - 12*im*omega) + 
            24*a^3*(1 + sqrt(1 - a^2))*m*(-im + 2*(1 + sqrt(1 - a^2))*omega) + 
            4*a*(1 + sqrt(1 - a^2))^2*m*(6*im + 2*im*(1 + sqrt(1 - a^2))*lambda + 
            3*(1 + sqrt(1 - a^2))^2*omega) - 4*a^2*(1 + sqrt(1 - a^2))*
            (6 + (1 + sqrt(1 - a^2))*(-3 + 6*m^2) + 2*im*(1 + sqrt(1 - a^2))*
            (3 + (1 + sqrt(1 - a^2))*(-3 + lambda))*omega + 3*(1 + sqrt(1 - a^2))^3*omega^2) + 
            12*a^4*(1 - 2*(1 + sqrt(1 - a^2))*omega*(-im + omega + sqrt(1 - a^2)*omega))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_TS_abs_Csq(s::Int, m::Int, a, omega, lambda)
    _sqrt2 = sqrt(2*one(real(a)))   # typed: BigFloat-exact (upstream: Float64 sqrt(2))
    if abs(s) != 2
        throw(AssertionError("|C|^2 is only defined for s = +/-2"))
    end

    # Here we assume that lambda is computed with s = -2
    if s == 2
        # −sAlm = sAlm + 2s
        # Compute what lambda should be for s = -2
        lambda_minus2 = lambda + 4
    else
        lambda_minus2 = lambda
    end

    E = lambda_minus2 + 2*a*m*omega - a^2*omega^2 + 2
    Q = E + a^2*omega^2 - 2*a*m*omega
    return begin
        (Q^2 + 4*a*omega*m - 4*a^2*omega^2)*((Q-2)^2 + 36*a*omega*m - 36*a^2*omega^2) + (2*Q - 1)*(96*a^2*omega^2 - 48a*omega*m) + 144*omega^2*(1-a^2)
    end
end


function _gsn_TS_Bsq(s::Int, m::Int, a, omega, lambda)
    _sqrt2 = sqrt(2*one(real(a)))   # typed: BigFloat-exact (upstream: Float64 sqrt(2))
    if abs(s) != 1
        throw(AssertionError("|B|^2 is only defined for s = +/-1"))
    end

    # Here we assume that lambda is computed with s = -1
    if s == 1
        # −sAlm = sAlm + 2s
        # Compute what lambda should be for s = -1
        lambda_minus1 = lambda + 2
    else
        lambda_minus1 = lambda
    end

    E = lambda_minus1 + 2*a*m*omega - a^2*omega^2
    Q = E + a^2*omega^2 - 2*a*m*omega
    return begin
        Q^2 + 4*a*omega*m - 4*a^2*omega^2
    end
end


"Up-solution reflection conversion factor: Cref ≡ Btrans. [GSN ConversionFactors.jl l.145]"
_gsn_Cref(s::Int, m::Int, a, omega, lambda) = _gsn_Btrans(s, m, a, omega, lambda)
"In-solution reflection conversion factor: Bref ≡ Ctrans. [GSN ConversionFactors.jl l.150]"
_gsn_Bref(s::Int, m::Int, a, omega, lambda) = _gsn_Ctrans(s, m, a, omega, lambda)
