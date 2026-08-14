# :GSN backend — the GSN potentials sF, sU (per spin weight s ∈ {0,±1,±2}), the Teukolsky radial
# potential V_T, and the SN potential V_SN reconstructed from η (GSN paper arXiv:2306.16469 §II,
# Eqs. (13)-(20) family; Mathematica derivation 10.5281/zenodo.8080241). Bodies below are the
# MECHANICALLY EXTRACTED upstream closed forms (src/Homogeneous/Potentials.jl, v0.8.5 commit
# 982b087) with renames `I→im`, `Delta→_gsn_Delta`, `K→_gsn_KK` and the audited fix
# `sF(s=0) = 0.0 → typed zero` (GSN_PLAN §3c-4). The GSN master equation solved downstream is
#     d²X/dr*² − sF(r) dX/dr* − sU(r) X = 0.        [GSN paper Eq. (13); Solutions.jl l.48]

function _gsn_sF(s::Int, m::Int, a, omega, lambda, r)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    if s == 0
        # s = 0
        return zero(im*(a + omega + r) + lambda)   # typed zero (upstream: Float64 0.0)
    elseif s == 1
        # s = +1
        return begin
            -((2*a*(a^2 + (-2 + r)*r)*(-3*im*a^2*m*r - im*m*r^3 + 2*a^3*(1 + lambda) + a*r*(3 + r*(3 + 2*lambda))))/
            (r*(a^2 + r^2)*(-2*im*a^3*m*r - 2*im*a*m*r^3 + a^4*(1 + lambda) + r^4*(2 + lambda) + 
            a^2*r*(2 + r*(3 + 2*lambda)))))
        end
    elseif s == -1
        # s = -1
        return begin
            -((2*a*(a^2 + (-2 + r)*r)*(3*im*a^2*m*r + im*m*r^3 + 2*a^3*(-1 + lambda) + a*r*(3 + r*(-1 + 2*lambda))))/
            (r*(a^2 + r^2)*(2*im*a^3*m*r + 2*im*a*m*r^3 + a^4*(-1 + lambda) + r^4*lambda + a^2*r*(2 + r*(-1 + 2*lambda)))))
        end
    elseif s == 2
        # s = +2
        return begin
            (8*a*(a^2 + (-2 + r)*r)*(im*m*r^2*(6 + r*(4 + lambda)) - 3*a^2*m*r*(3*im + 4*r*omega) + 
            a*r*(9 + r*(-3 + 6*m^2 - 6*im*omega) - im*r^2*(1 + lambda)*omega) + a^3*(-6 + 9*im*r*omega + 6*r^2*omega^2)))/
            (r*(a^2 + r^2)*(r^4*(24 + 10*lambda + lambda^2 + 12*im*omega) + 24*a^3*m*r*(im + 2*r*omega) + 
            4*a*m*r^2*(-6*im - 2*im*r*(4 + lambda) + 3*r^2*omega) - 12*a^4*(-1 + 2*im*r*omega + 2*r^2*omega^2) - 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*im*omega) - 2*im*r^2*(1 + lambda)*omega + 3*r^3*omega^2)))
        end
    elseif s == -2
        # s = -2
        return begin
            (8*a*(a^2 + (-2 + r)*r)*((-im)*m*r^2*(6 + r*lambda) + 3*a^2*m*r*(3*im - 4*r*omega) + 
            a*r*(9 + r*(-3 + 6*m^2 + 6*im*omega) + im*r^2*(-3 + lambda)*omega) + a^3*(-6 - 9*im*r*omega + 6*r^2*omega^2)))/
            (r*(a^2 + r^2)*(r^4*(2*lambda + lambda^2 - 12*im*omega) + 24*a^3*m*r*(-im + 2*r*omega) + 
            4*a*m*r^2*(6*im + 2*im*r*lambda + 3*r^2*omega) + a^4*(12 + 24*im*r*omega - 24*r^2*omega^2) - 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*im*omega) + 2*im*r^2*(-3 + lambda)*omega + 3*r^3*omega^2)))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_sU(s::Int, m::Int, a, omega, lambda, r)
    #=
    We have derived/shown the explicit expression for
    different physically-relevant spin weight (s=0, \pm 1, \pm2)
    =#
    if s == 0
        # s = 0
        return begin
            (1/(a^2 + r^2)^4)*(r^2*(a^2 + (-2 + r)*r)^2 + (a^2 + (-2 + r)*r)*(a^4 - 4*a^2*r - (-4 + r)*r^3) + 
            (a^2 + r^2)^2*((a^2 + (-2 + r)*r)*lambda - (a*m - (a^2 + r^2)*omega)^2))
        end
    elseif s == 1
        # s = +1
        return begin
            (1/(a^2 + r^2)^4)*(2*(a^2 + (-2 + r)*r)*(a^4 - a^2*r - (-3 + r)*r^3) + 
            (r^2*(-3 + 2*r) + a^2*(-1 + 2*r))^2 + (2*a*(a^2 + (-2 + r)*r)*(a^2 + r^2)*
            (r^2*(-3 + 2*r) + a^2*(-1 + 2*r))*(-3*im*a^2*m*r - im*m*r^3 + 2*a^3*(1 + lambda) + 
            a*r*(3 + r*(3 + 2*lambda))))/(r*(-2*im*a^3*m*r - 2*im*a*m*r^3 + a^4*(1 + lambda) + r^4*(2 + lambda) + 
            a^2*r*(2 + r*(3 + 2*lambda)))) - ((a^2 + r^2)^2*(2*a^7*m*r^2*omega*(-lambda - im*r*omega) + 
            a^8*(1 + lambda)*(2 + r^2*omega^2) - 2*im*a^5*m*r^2*(-6 + r*(4 + m^2 - lambda - 4*im*omega) - 
            3*im*r^2*(1 + lambda)*omega + 3*r^3*omega^2) - 2*im*a*m*r^5*(3 + r*(-5 + 2*lambda) + r^2*(2 - lambda + 2*im*omega) - 
            im*r^3*(3 + lambda)*omega + r^4*omega^2) - 2*im*a^3*m*r^3*(13 + r*(-17 + 2*lambda) + 
            r^2*(6 + m^2 - 2*lambda - 2*im*omega) - 3*im*r^3*(2 + lambda)*omega + 3*r^4*omega^2) + 
            r^6*(-3*(2 + lambda) - r^2*lambda*(2 + lambda) + 2*r*(2 + 3*lambda + lambda^2) + r^4*(2 + lambda)*omega^2) + 
            a^6*r*(-16*(1 + lambda) + r*(5 + m^2*(-1 + lambda) + 6*lambda - lambda^2 + 2*im*omega) + r^3*(5 + 4*lambda)*omega^2 + 
            2*r^2*omega*(im + 2*im*m^2 + omega)) + a^4*r^2*(11 + 25*lambda + 
            2*r*(-5 + 3*m^2 - 13*lambda + lambda^2 - 2*im*omega) + r^2*(5 + 4*lambda - 3*lambda^2 + m^2*(3 + 2*lambda) - 
            8*im*omega) + 3*r^4*(3 + 2*lambda)*omega^2 + 4*r^3*omega*(im + 2*im*m^2 + omega)) + 
            a^2*r^3*(30 + 5*r*(-7 + 2*lambda) - 4*r^2*(-2 + m^2 + lambda - lambda^2 - 3*im*omega) + 
            r^3*(2 - 2*lambda - 3*lambda^2 + m^2*(4 + lambda) - 10*im*omega) + r^5*(7 + 4*lambda)*omega^2 + 
            2*r^4*omega*(im + 2*im*m^2 + omega))))/(r^2*(-2*im*a^3*m*r - 2*im*a*m*r^3 + a^4*(1 + lambda) + 
            r^4*(2 + lambda) + a^2*r*(2 + r*(3 + 2*lambda)))))
        end
    elseif s == -1
        # s = -1
        return begin
            (1/(a^2 + r^2)^4)*((a^2 - r^2)^2 + 2*r*(a^2 + (-2 + r)*r)*(-3*a^2 + r^2) + 
            (2*a*(a^2 + (-2 + r)*r)*(a^2 - r^2)*(a^2 + r^2)*(3*im*a^2*m*r + im*m*r^3 + 2*a^3*(-1 + lambda) + 
            a*r*(3 + r*(-1 + 2*lambda))))/(r*(2*im*a^3*m*r + 2*im*a*m*r^3 + a^4*(-1 + lambda) + r^4*lambda + 
            a^2*r*(2 + r*(-1 + 2*lambda)))) - ((a^2 + r^2)^2*(2*a^7*m*r^2*omega*(2 - lambda + im*r*omega) + 
            a^8*(-1 + lambda)*(2 + r^2*omega^2) + 2*im*a^5*m*r^3*(-2 + m^2 - lambda + 4*im*omega - 3*im*r*omega + 
            3*im*r*lambda*omega + 3*r^2*omega^2) + r^6*lambda*(-3 - r^2*lambda + 2*r*(1 + lambda) + r^4*omega^2) + 
            2*im*a*m*r^5*(-1 + r + 2*r*lambda - r^2*(lambda + 2*im*omega) + im*r^3*(1 + lambda)*omega + r^4*omega^2) + 
            2*im*a^3*m*r^3*(1 + r*(3 + 2*lambda) + r^2*(-2 + m^2 - 2*lambda + 2*im*omega) + 3*im*r^3*lambda*omega + 
            3*r^4*omega^2) + a^6*r*(8 - 8*lambda + r*(-1 + m^2*(-3 + lambda) - lambda^2 - 2*im*omega) + 
            r^3*(-3 + 4*lambda)*omega^2 + 2*r^2*omega*(-im - 2*im*m^2 + omega)) + 
            a^4*r^2*(-11 + 9*lambda + 2*r*(1 + 3*m^2 + lambda + lambda^2 + 2*im*omega) + 
            r^2*(1 - 6*lambda - 3*lambda^2 + m^2*(-1 + 2*lambda) + 8*im*omega) + 3*r^4*(-1 + 2*lambda)*omega^2 + 
            4*r^3*omega*(-im - 2*im*m^2 + omega)) + a^2*r^3*(6 - 3*r*(1 + 2*lambda) - 
            4*r^2*(m^2 - 3*lambda - lambda^2 + 3*im*omega) + r^3*(-4*lambda - 3*lambda^2 + m^2*(2 + lambda) + 10*im*omega) + 
            r^5*(-1 + 4*lambda)*omega^2 + 2*r^4*omega*(-im - 2*im*m^2 + omega))))/
            (r^2*(2*im*a^3*m*r + 2*im*a*m*r^3 + a^4*(-1 + lambda) + r^4*lambda + a^2*r*(2 + r*(-1 + 2*lambda)))))
        end
    elseif s == 2
        # s = +2
        return begin
            (1/(a^2 + r^2)^4)*((a^2 + (-2 + r)*r)*(3*a^4 + (8 - 3*r)*r^3) + 
            (r^2*(-4 + 3*r) + a^2*(-2 + 3*r))^2 + (8*a*(a^2 + (-2 + r)*r)*(a^2 + r^2)*
            (r^2*(-4 + 3*r) + a^2*(-2 + 3*r))*(im*m*r^2*(6 + r*(4 + lambda)) - 3*a^2*m*r*(3*im + 4*r*omega) + 
            a*r*(9 + r*(-3 + 6*m^2 - 6*im*omega) - im*r^2*(1 + lambda)*omega) + a^3*(-6 + 9*im*r*omega + 6*r^2*omega^2)))/
            (r*((-r^4)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 24*a^3*m*r*(im + 2*r*omega) + 
            4*a*m*r^2*(6*im + 2*im*r*(4 + lambda) - 3*r^2*omega) + 12*a^4*(-1 + 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*im*omega) - 2*im*r^2*(1 + lambda)*omega + 3*r^3*omega^2))) - 
            ((a^2 + r^2)^2*(24*a^7*m*r^2*omega*(3 + 3*im*r*omega - 4*r^2*omega^2) - r^6*(24 + 10*lambda + lambda^2 + 12*im*omega)*
            (-12 - r^2*lambda + 2*r*(4 + lambda) + r^4*omega^2) - 4*a^5*m*r^2*
            (-54*im - 2*im*r*(-55 + 3*m^2 - 4*lambda - 18*im*omega) + 2*r^2*(4 + 12*m^2 - 8*lambda + 27*im*omega)*omega - 
            6*im*r^3*(8 + lambda)*omega^2 + 45*r^4*omega^3) + 12*a^8*(-2 - 3*r^2*omega^2 - 2*im*r^3*omega^3 + 
            2*r^4*omega^4) + 4*a^6*r*(60 - 3*r*(16 + 3*m^2 + 18*im*omega) + r^3*(1 + 36*m^2 - 8*lambda + 18*im*omega)*
            omega^2 - 2*im*r^4*(7 + lambda)*omega^3 + 15*r^5*omega^4 - 2*r^2*omega*(-40*im + 9*im*m^2 - 4*im*lambda + 9*omega)) + 
            2*a*m*r^4*(144*im + 8*im*r*(-38 + lambda) + r^4*(40 + 20*lambda + lambda^2 + 24*im*omega)*omega + 
            4*im*r^5*(4 + lambda)*omega^2 - 6*r^6*omega^3 + 4*r^2*(16*im + 3*im*lambda + 2*im*lambda^2 + 6*omega) - 
            4*im*r^3*(-12 + lambda + lambda^2 - 14*im*omega - 5*im*lambda*omega)) - 
            2*a^3*m*r^3*(312*im + r^3*(84 + 30*m^2 - 52*lambda - lambda^2 + 36*im*omega)*omega - 
            4*im*r^4*(19 + 4*lambda)*omega^2 + 48*r^5*omega^3 + 12*r*(-63*im + 3*im*m^2 - 3*im*lambda + 32*omega) + 
            4*r^2*(71*im + 9*im*lambda + im*lambda^2 - im*m^2*(10 + lambda) - 104*omega + 16*lambda*omega)) + 
            a^2*r^3*(576 + 96*r*(-12 + 4*m^2 - 3*im*omega) + 8*r^2*(18 - 3*lambda^2 + m^2*(-58 + 8*lambda) + 
            lambda*(-30 - 2*im*omega) + 46*im*omega) - 8*im*r^4*(-lambda^2 + 2*m^2*(7 + lambda) + lambda*(2 + 5*im*omega) + 
            11*im*omega)*omega + 2*r^5*(-34 + 24*m^2 - 20*lambda - lambda^2 - 24*im*omega)*omega^2 - 8*im*r^6*(1 + lambda)*omega^3 + 
            12*r^7*omega^4 + r^3*(lambda^3 + 36*lambda*(4 + im*omega) + 2*lambda^2*(11 - 8*im*omega) - 
            m^2*(-136 + 42*lambda + lambda^2 - 36*im*omega) - 8*(-18 + 19*im*omega + 6*omega^2))) + 
            a^4*r^2*(-672 + r*(960 - 72*m^2 + 624*im*omega) + r^4*(44 + 180*m^2 - 62*lambda - lambda^2 + 36*im*omega)*
            omega^2 - 8*im*r^5*(5 + 2*lambda)*omega^3 + 48*r^6*omega^4 + 8*r^3*omega*(56*im + 6*im*lambda + im*lambda^2 - 
            3*im*m^2*(9 + lambda) - 46*omega + 8*lambda*omega) + 4*r^2*(6*m^4 + m^2*(7 - 8*lambda + 54*im*omega) + 
            2*(-15 + 10*lambda + lambda^2 - 144*im*omega - 9*im*lambda*omega + 48*omega^2)))))/
            (r^2*((-r^4)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 24*a^3*m*r*(im + 2*r*omega) + 
            4*a*m*r^2*(6*im + 2*im*r*(4 + lambda) - 3*r^2*omega) + 12*a^4*(-1 + 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*im*omega) - 2*im*r^2*(1 + lambda)*omega + 3*r^3*omega^2))))
        end
    elseif s == -2
        # s = -2
        return begin
            (1/(a^2 + r^2)^4)*((a^2*(-2 + r) + r^3)^2 + (a^2 + (-2 + r)*r)*(-a^4 - 8*a^2*r + r^4) + 
            (8*a*(a^2 + (-2 + r)*r)*(a^2 + r^2)*(a^2*(-2 + r) + r^3)*((-im)*m*r^2*(6 + r*lambda) + 
            3*a^2*m*r*(3*im - 4*r*omega) + a*r*(9 + r*(-3 + 6*m^2 + 6*im*omega) + im*r^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*im*r*omega + 6*r^2*omega^2)))/(r*(r^4*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*r*(-im + 2*r*omega) + 4*a*m*r^2*(6*im + 2*im*r*lambda + 3*r^2*omega) + 
            a^4*(12 + 24*im*r*omega - 24*r^2*omega^2) - 4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*r^2*(-3 + lambda)*omega + 3*r^3*omega^2))) - 
            ((a^2 + r^2)^2*(-24*a^7*m*r^2*omega*(-3 + 3*im*r*omega + 4*r^2*omega^2) + 
            12*a^8*(-2 - 3*r^2*omega^2 + 2*im*r^3*omega^3 + 2*r^4*omega^4) + 
            4*a^6*r*(12 - 3*r*(-4 + 3*m^2 + 6*im*omega) + r^3*(-39 + 36*m^2 - 8*lambda - 18*im*omega)*omega^2 + 
            2*im*r^4*(3 + lambda)*omega^3 + 15*r^5*omega^4 + 2*r^2*omega*(24*im + 9*im*m^2 - 4*im*lambda + 15*omega)) - 
            4*a^5*m*r^2*(-18*im + 2*r^2*(-36 + 12*m^2 - 8*lambda - 27*im*omega)*omega + 6*im*r^3*(4 + lambda)*omega^2 + 
            45*r^4*omega^3 + 2*r*(9*im + 3*im*m^2 - 4*im*lambda + 30*omega)) - 
            2*a*m*r^4*(-48*im - 24*im*r*(-2 + lambda) + 4*im*r^2*(7*lambda + 2*lambda^2 + 6*im*omega) - 
            r^4*(12*lambda + lambda^2 - 24*im*omega)*omega + 4*im*r^5*lambda*omega^2 + 6*r^6*omega^3 - 
            4*im*r^3*(lambda + lambda^2 + 6*im*omega + 5*im*lambda*omega)) + a^4*r^3*(24*m^4*r + 48*(-4 + 3*im*omega) - 
            r^3*(60 + 54*lambda + lambda^2 + 36*im*omega)*omega^2 + 8*im*r^4*(-3 + 2*lambda)*omega^3 + 48*r^5*omega^4 + 
            8*r*(9 + 2*lambda + lambda^2 - 72*im*omega + 13*im*lambda*omega) + 8*r^2*omega*(24*im - 6*im*lambda - im*lambda^2 + 18*omega + 
            8*lambda*omega) + 4*m^2*(30 - r*(33 + 8*lambda + 54*im*omega) + 6*im*r^2*(5 + lambda)*omega + 45*r^3*omega^2)) + 
            r^6*((-2 + r)*r*lambda^3 + lambda^2*(12 - 12*r + 2*r^2 - r^4*omega^2) + 12*im*omega*(-12 + 8*r + r^4*omega^2) - 
            2*lambda*(-12 + 4*r*(2 - 3*im*omega) + 6*im*r^2*omega + r^4*omega^2)) + 
            a^2*r^4*(96 - 96*im*omega + 2*r^4*(6 + 24*m^2 - 12*lambda - lambda^2 + 24*im*omega)*omega^2 + 
            8*im*r^5*(-3 + lambda)*omega^3 + 12*r^6*omega^4 + 8*r^3*omega*(2*im*lambda - im*lambda^2 + 2*im*m^2*(3 + lambda) + 3*omega + 
            5*lambda*omega) + r^2*(lambda^3 + lambda*(24 - 34*m^2 - 4*im*omega) + lambda^2*(14 - m^2 + 16*im*omega) - 
            12*im*(22 + 3*m^2 - 4*im*omega)*omega) + 8*r*(m^2*(6 + 8*lambda) - 
            3*(2 + 2*lambda + lambda^2 - 22*im*omega + 2*im*lambda*omega))) - 
            2*a^3*m*r^3*(72*im - 4*im*r*(15 + 9*m^2 - 13*lambda) + r^3*(-36 + 30*m^2 - 44*lambda - lambda^2 - 36*im*omega)*
            omega + 4*im*r^4*(3 + 4*lambda)*omega^2 + 48*r^5*omega^3 + 4*r^2*(im*m^2*(6 + lambda) - 
            im*(3 + 9*lambda + lambda^2 + 24*im*omega + 16*im*lambda*omega)))))/
            (r^2*((-r^4)*(2*lambda + lambda^2 - 12*im*omega) + 24*a^3*m*r*(im - 2*r*omega) - 
            4*a*m*r^2*(6*im + 2*im*r*lambda + 3*r^2*omega) + 12*a^4*(-1 - 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*im*omega) + 2*im*r^2*(-3 + lambda)*omega + 3*r^3*omega^2))))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_VT(s::Int, m::Int, a, omega, lambda, r)
    _K = _gsn_KK(m, a, omega, r)
    _Delta = _gsn_Delta(a, r)
    return lambda - 4im*s*omega*r - (_K^2 - 2im*s*(r-1)*_K)/_Delta
end


"""
    _gsn_VSN(s, m, a, omega, lambda, r) -> Complex

Sasaki-Nakamura potential reconstructed from the transformation function η:
V_SN = sU − ( (1/2)·(d²η/dr*²)/η − (3/4)·((dη/dr*)/η)² ), with dr/dr* = Δ/(r²+a²) and
d²r/dr*² = 2(r²−a²)/(r²+a²)². Typed rationals replace the upstream Float64 0.5/0.75
(GSN_PLAN §3c-4). Diagnostic/identity use. [GSN Potentials.jl VSN, v0.8.5]
"""
function _gsn_VSN(s::Int, m::Int, a, omega, lambda, r)
    drdrs = _gsn_Delta(a, r)/(r^2 + a^2)
    d2rdrs2 = (2*(r^2 - a^2))/((r^2 + a^2)^2)
    etav = _gsn_eta(s, m, a, omega, lambda, r)
    etap = _gsn_eta_prime(s, m, a, omega, lambda, r)
    etapp = _gsn_eta_primeprime(s, m, a, omega, lambda, r)
    detadrs = drdrs*etap
    d2etadrs2 = drdrs*(etap*d2rdrs2 + drdrs*etapp)
    return _gsn_sU(s, m, a, omega, lambda, r) - ( (d2etadrs2/etav)/2 - (3//4)*(detadrs/etav)^2 )
end

"""
    _gsn_teukolsky_residual(s, m, a, omega, lambda, r, R, dRdr, d2Rdr2) -> Complex

Residual of the radial Teukolsky equation Δ R″ + 2(s+1)(r−1) R′ − V_T R (M=1) for VALUE inputs
(diagnostic/test gate). [GSN Potentials.jl radial_Teukolsky_equation]
"""
_gsn_teukolsky_residual(s::Int, m::Int, a, omega, lambda, r, R, dRdr, d2Rdr2) =
    _gsn_Delta(a, r)*d2Rdr2 + 2*(s + 1)*(r - 1)*dRdr - _gsn_VT(s, m, a, omega, lambda, r)*R
