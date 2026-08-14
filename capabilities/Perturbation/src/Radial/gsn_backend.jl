# :GSN backend — carrier integration: the Generalized Sasaki-Nakamura backend of the
# `TeukolskyRadial` dispatcher (method=:GSN). Closes the documented radial_sn boundary
# ("s=-2 Up only" / NotYetPorted): s ∈ {0,±1,±2} × {In,Up} homogeneous Teukolsky radial
# solutions at real ω, self-normalized to unit transmission via the GSN asymptotic ansatzes
# (no MST anchor needed — unlike radial_sn's scale-match). Absorption of the vendored
# GeneralizedSasakiNakamura.jl v0.8.5 commit 982b087 (arXiv:2306.16469) per
# progress/emri-disk/GSN_PLAN.md; oracle fixtures: test/fixtures/gsn_oracle.jl.
#
# Teukolsky-frame conversion [upstream main module l.714-749]:
#   raw Teukolsky amplitudes = conversion-factor × SN amplitudes (SN transmission ≡ 1):
#     In: (Trans, Inc, Ref) = (Btrans, Binc_cf·Binc_SN, Bref_cf·Bref_SN)
#     Up: (Trans, Inc, Ref) = (Ctrans, Cinc_cf·Cinc_SN, Cref_cf·Cref_SN)
#   the dispatcher's mkamp then normalizes to unit transmission — identical to upstream's
#   UNIT_TEUKOLSKY_TRANS; carrier norm = the raw transmission (the conversion factor), and
#   R(r) = M(r)·[X(r*); X′(r*)] / norm with M the SN→Teukolsky conversion matrix below.
#
# Upstream latent defects (GSN_PLAN §3d): D5 (ω=0 UndefVarError branch) is obsolete by design —
# ω=0 routes to the :Static backend before any GSN code runs; D6 (unbound rsin/rsout in the
# auto-widening loop) is fixed by construction in `_gsn_homogeneous` (explicit binding); D4
# (broken Phi-Wronskian) dropped in gsn_solve.jl; D1-D3 live in the inhomogeneous subsystem,
# NOT ported here (quarantined to GSN-3); D7 solver-docstring mismatch is moot (solver replaced).

function _gsn_conversion_matrix(s, m, a, omega, lambda, r)
    #=
    Here we use explicit form for the conversion matrix to
    faciliate cancellations
    =#
    _zc = zero(im*(a + omega + r) + lambda)   # typed zero (upstream: Float64 0.0)
    _sqrt2 = sqrt(2*one(real(a)))
    M11 = M12 = M21 = M22 = _zc
    if s == 0
        M11 = 1/sqrt(a^2 + r^2)
        M21 = -(r/(a^2 + r^2)^(3/2))
        M22 = sqrt(a^2 + r^2)/(a^2 + (-2 + r)*r)
    elseif s == +1
        M11 = begin
            -((im*r*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            ((-a^3)*m*r - a*m*r^3 + r^5*omega + 
            a^4*(2*im + r*omega) + 2*a^2*r*(-2*im + im*r + 
            r^2*omega)))/(sqrt(a^2 + r^2)*
            sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*
            (-2*im*a^3*m*r - 2*im*a*m*r^3 + a^4*(1 + lambda) + 
            r^4*(2 + lambda) + a^2*r*(2 + r*(3 + 2*lambda)))))
        end
        M12 = begin
            (r^2*(a^2 + r^2)^(3/2)*
            sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r)))/
            (sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*
            (-2*im*a^3*m*r - 2*im*a*m*r^3 + a^4*(1 + lambda) + 
            r^4*(2 + lambda) + a^2*r*(2 + r*(3 + 2*lambda))))
        end
        M21 = begin
            -((r*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            (-2*a^7*m*(im + r*omega) + a^8*omega*(2*im + r*omega) + 
            a^5*m*r*(2*im - im*r - 6*r^2*omega) + 
            a*m*r^5*(-4*im + 3*im*r - 2*r^2*omega) - 
            2*a^3*m*r^3*(im - 2*im*r + 3*r^2*omega) + 
            r^6*(2*(2 + lambda) - r*(2 + lambda) - im*r^2*omega + 
            r^3*omega^2) + a^6*(-4 + r*(1 + m^2 - lambda - 
            2*im*omega) + 5*im*r^2*omega + 4*r^3*omega^2) + 
            a^2*r^3*(4 + 4*r*lambda + r^2*(-3 + m^2 - 3*lambda - 
            2*im*omega) - im*r^3*omega + 4*r^4*omega^2) + 
            a^4*r*(8 + 2*r*(-4 + lambda) + 
            r^2*(2*m^2 - 3*lambda - 4*im*omega) + 3*im*r^3*omega + 
            6*r^4*omega^2)))/(sqrt(a^2 + r^2)*
            ((a^2 + (-2 + r)*r)*(a^2 + r^2))^(3/2)*
            (-2*im*a^3*m*r - 2*im*a*m*r^3 + a^4*(1 + lambda) + 
            r^4*(2 + lambda) + a^2*r*(2 + r*(3 + 2*lambda)))))
        end
        M22 = begin
            -((im*r^2*((a^2 + r^2)/(a^2 + (-2 + r)*r))^(3/2)*
            ((-a^3)*m - a*m*r^2 + a^4*omega + r^3*(-im + r*omega) + 
            a^2*(2*im - im*r + 2*r^2*omega)))/
            (sqrt(a^2 + r^2)*sqrt((a^2 + (-2 + r)*r)*
            (a^2 + r^2))*(-2*im*a^3*m*r - 2*im*a*m*r^3 + 
            a^4*(1 + lambda) + r^4*(2 + lambda) + 
            a^2*r*(2 + r*(3 + 2*lambda)))))
        end
    elseif s == -1
        M11 = begin
            (im*r*sqrt(a^2 + r^2)*((-a^3)*m*r - a*m*r^3 + 
            r^5*omega + a^4*(-2*im + r*omega) + 
            2*a^2*r*(2*im - im*r + r^2*omega)))/
            (sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*
            (2*im*a^3*m*r + 2*im*a*m*r^3 + a^4*(-1 + lambda) + 
            r^4*lambda + a^2*r*(2 + r*(-1 + 2*lambda))))
        end
        M12 = begin
            (r^2*sqrt(a^2 + r^2)*sqrt((a^2 + r^2)/
            (a^2 + (-2 + r)*r))*sqrt((a^2 + (-2 + r)*r)*
            (a^2 + r^2)))/(2*im*a^3*m*r + 2*im*a*m*r^3 + 
            a^4*(-1 + lambda) + r^4*lambda + 
            a^2*r*(2 + r*(-1 + 2*lambda)))
        end
        M21 = begin
            (r*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            (a^8*omega*(2*im - r*omega) + 2*a^7*m*(-im + r*omega) + 
            a*m*r^5*(-2*im + im*r + 2*r^2*omega) + 
            2*a^3*m*r^3*(im + 3*r^2*omega) + 
            a^5*m*r*(4*im - 3*im*r + 6*r^2*omega) + 
            a^6*r*(1 - m^2 + lambda - 4*im*omega + 7*im*r*omega - 
            4*r^2*omega^2) - a^4*r^2*(2*(2 + lambda) + 
            r*(-2 + 2*m^2 - 3*lambda + 10*im*omega) - 9*im*r^2*omega + 
            6*r^3*omega^2) + a^2*r^3*(4 - 4*r*(1 + lambda) - 
            r^2*(-1 + m^2 - 3*lambda + 8*im*omega) + 5*im*r^3*omega - 
            4*r^4*omega^2) + r^6*((-2 + r)*lambda - 
            r*omega*(2*im - im*r + r^2*omega))))/
            ((a^2 + r^2)^(3/2)*sqrt((a^2 + (-2 + r)*r)*
            (a^2 + r^2))*(2*im*a^3*m*r + 2*im*a*m*r^3 + 
            a^4*(-1 + lambda) + r^4*lambda + 
            a^2*r*(2 + r*(-1 + 2*lambda))))
        end
        M22 = begin
            -((r^2*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r))*
            sqrt((a^2 + (-2 + r)*r)*(a^2 + r^2))*
            (-r - (im*(a^2 + r^2)*((-a)*m + (a^2 + r^2)*omega))/
            (a^2 + (-2 + r)*r)))/(sqrt(a^2 + r^2)*
            (2*im*a^3*m*r + 2*im*a*m*r^3 + a^4*(-1 + lambda) + 
            r^4*lambda + a^2*r*(2 + r*(-1 + 2*lambda)))))
        end
    elseif s == +2
        M11 = begin
            (r^2*(a^2 + (-2 + r)*r)*(-4*a^5*m*r*(im + r*omega) - 
            2*a^3*m*r^2*(-3*im + 2*im*r + 4*r^2*omega) + 
            a*m*(2*im*r^4 - 4*r^6*omega) + 
            2*a^6*(-5 + 2*im*r*omega + r^2*omega^2) + 
            a^4*r*(32 + r*(-24 + 2*m^2 - lambda - 6*im*omega) + 
            10*im*r^2*omega + 6*r^3*omega^2) + 
            r^4*(-12 + 2*r*(9 + lambda) - 
            r^2*(6 + lambda + 6*im*omega) + 2*im*r^3*omega + 
            2*r^4*omega^2) + 2*a^2*r^2*(-12 + r*(23 + lambda) + 
            r^2*(-10 + m^2 - lambda - 6*im*omega) + 4*im*r^3*omega + 
            3*r^4*omega^2)))/
            (((a^2 + (-2 + r)*r)^2*(a^2 + r^2))^(3/2)*
            ((-r^4)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 
            24*a^3*m*r*(im + 2*r*omega) + 4*a*m*r^2*
            (6*im + 2*im*r*(4 + lambda) - 3*r^2*omega) + 
            12*a^4*(-1 + 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*r^2*(1 + lambda)*omega + 3*r^3*omega^2)))
        end
        M12 = begin
            (2*r^3*(a^2 + (-2 + r)*r)*(a^2 + r^2)^2*
            ((-im)*a*m*r + a^2*(-2 + im*r*omega) + 
            r*(3 - r + im*r^2*omega)))/
            (((a^2 + (-2 + r)*r)^2*(a^2 + r^2))^(3/2)*
            ((-r^4)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 
            24*a^3*m*r*(im + 2*r*omega) + 4*a*m*r^2*
            (6*im + 2*im*r*(4 + lambda) - 3*r^2*omega) + 
            12*a^4*(-1 + 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*r^2*(1 + lambda)*omega + 3*r^3*omega^2)))
        end
        M21 = begin
            -((im*r*(-2*a^7*m*r*(-2 + 4*im*r*omega + 3*r^2*omega^2) - 
            2*a^5*m*r^2*(-4 + r*(3 + m^2 - lambda + 2*im*omega) + 
            4*im*r^2*omega + 9*r^3*omega^2) + 
            2*a*m*r^5*(-10 + r*(3 - 2*lambda) + 
            r^2*(2 + lambda + 2*im*omega) + 4*im*r^3*omega - 
            3*r^4*omega^2) - 2*a^3*m*r^3*
            (12 + r^2*(3 + m^2 - 2*lambda) + r*(-17 + 2*lambda) - 
            4*im*r^3*omega + 9*r^4*omega^2) + 
            2*a^8*(-6*im - 2*r*omega + 2*im*r^2*omega^2 + 
            r^3*omega^3) + 2*a^6*r*(16*im + 
            im*r*(-9 + 2*m^2 - 2*lambda + 4*im*omega) + 
            r^2*(-8 + 3*m^2 - lambda + im*omega)*omega + 
            6*im*r^3*omega^2 + 4*r^4*omega^3) + 
            r^6*(-4*im*lambda - 2*r^3*(4 + lambda + 5*im*omega)*omega + 
            2*r^5*omega^3 - 12*r*(im + omega) + 
            r^2*(6*im + im*lambda + 18*omega + 4*lambda*omega)) + 
            a^4*r^2*(-16*im + 2*im*r*(10 + m^2 + 6*lambda - 
            12*im*omega) + 2*r^3*(-14 + 6*m^2 - 3*lambda - 
            3*im*omega)*omega + 12*im*r^4*omega^2 + 12*r^5*omega^3 + 
            r^2*(-2*im - 4*im*m^2 - 7*im*lambda + 22*omega + 
            4*lambda*omega)) + 2*a^2*r^4*(-4*im*lambda + 
            im*r*(-4 + 3*m^2 + 6*lambda + 6*im*omega) + 
            3*r^3*(-4 + m^2 - lambda - 3*im*omega)*omega + 
            2*im*r^4*omega^2 + 4*r^5*omega^3 + 
            r^2*(5*im - 4*im*m^2 - im*lambda + 24*omega + 
            4*lambda*omega))))/
            (((a^2 + (-2 + r)*r)^2*(a^2 + r^2))^(3/2)*
            ((-r^4)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 
            24*a^3*m*r*(im + 2*r*omega) + 4*a*m*r^2*
            (6*im + 2*im*r*(4 + lambda) - 3*r^2*omega) + 
            12*a^4*(-1 + 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*r^2*(1 + lambda)*omega + 3*r^3*omega^2))))
        end
        M22 = begin
            (r^2*(a^2 + r^2)^2*(-4*a^3*m*r*(im + r*omega) - 
            2*a*m*r^2*(im - 3*im*r + 2*r^2*omega) + 
            2*a^4*(-3 + 2*im*r*omega + r^2*omega^2) + 
            r^3*(2*lambda - r*(2 + lambda + 10*im*omega) + 2*r^3*omega^2) + 
            a^2*r*(8 + 2*m^2*r - r*lambda + 2*im*r*omega + 
            4*im*r^2*omega + 4*r^3*omega^2)))/
            (((a^2 + (-2 + r)*r)^2*(a^2 + r^2))^(3/2)*
            ((-r^4)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 
            24*a^3*m*r*(im + 2*r*omega) + 4*a*m*r^2*
            (6*im + 2*im*r*(4 + lambda) - 3*r^2*omega) + 
            12*a^4*(-1 + 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*r^2*(1 + lambda)*omega + 3*r^3*omega^2)))
        end
    elseif s == -2
        M11 = begin
            (r^2*(-4*a^5*m*r*(-im + r*omega) - 
            2*a*m*r^4*(im + 2*r^2*omega) - 2*a^3*m*r^2*
            (3*im - 2*im*r + 4*r^2*omega) + 
            2*a^6*(-5 - 2*im*r*omega + r^2*omega^2) + 
            a^4*r*(32 + r*(-20 + 2*m^2 - lambda + 6*im*omega) - 
            10*im*r^2*omega + 6*r^3*omega^2) + 
            r^4*(-12 + 2*r*(5 + lambda) - 
            r^2*(2 + lambda - 6*im*omega) - 2*im*r^3*omega + 
            2*r^4*omega^2) + 2*a^2*r^2*(-12 + r*(19 + lambda) + 
            r^2*(-6 + m^2 - lambda + 6*im*omega) - 4*im*r^3*omega + 
            3*r^4*omega^2)))/((a^2 + (-2 + r)*r)^3*
            ((a^2 + r^2)/(a^2 + (-2 + r)*r)^2)^(3/2)*
            ((-r^4)*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*r*(im - 2*r*omega) - 4*a*m*r^2*
            (6*im + 2*im*r*lambda + 3*r^2*omega) + 
            12*a^4*(-1 - 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*r^2*(-3 + lambda)*omega + 3*r^3*omega^2)))
        end
        M12 = begin
            (2*r^3*(a^2 + (-2 + r)*r)*
            sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r)^2)*
            (im*a*m*r + a^2*(-2 - im*r*omega) + 
            r*(3 - r - im*r^2*omega)))/
            ((-r^4)*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*r*(im - 2*r*omega) - 
            4*a*m*r^2*(6*im + 2*im*r*lambda + 3*r^2*omega) + 
            12*a^4*(-1 - 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*r^2*(-3 + lambda)*omega + 3*r^3*omega^2))
        end
        M21 = begin
            (im*r*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r)^2)*
            (2*a^7*m*r*(2 + 4*im*r*omega - 3*r^2*omega^2) - 
            2*a^5*m*r^2*(4 + r*(-1 + m^2 - lambda + 6*im*omega) - 
            12*im*r^2*omega + 9*r^3*omega^2) - 
            2*a^3*m*r^4*(-5 + 2*lambda + r*(3 + m^2 - 2*lambda + 
            16*im*omega) - 12*im*r^2*omega + 9*r^3*omega^2) - 
            2*a*m*r^5*(6 + r*(-7 + 2*lambda) + 
            r^2*(2 - lambda + 10*im*omega) - 4*im*r^3*omega + 
            3*r^4*omega^2) + 2*a^8*(6*im - 2*r*omega - 
            2*im*r^2*omega^2 + r^3*omega^3) + 
            2*a^6*r*(-36*im + r^2*(-12 + 3*m^2 - lambda + 3*im*omega)*
            omega - 10*im*r^3*omega^2 + 4*r^4*omega^3 + 
            r*(21*im - 2*im*m^2 + 2*im*lambda + 4*omega)) + 
            r^5*(-48*im + 12*im*r*(6 + lambda) - 
            12*im*r^2*(3 + lambda - 3*im*omega) - 
            2*r^4*(4 + lambda - 9*im*omega)*omega - 8*im*r^5*omega^2 + 
            2*r^6*omega^3 + r^3*(6*im + 3*im*lambda + 34*omega + 
            4*lambda*omega)) + 2*a^2*r^3*(-48*im + 
            4*im*r*(27 + 2*lambda) + im*r^2*(-72 + m^2 - 14*lambda + 
            30*im*omega) + r^4*(-16 + 3*m^2 - 3*lambda + 
            21*im*omega)*omega - 14*im*r^5*omega^2 + 4*r^6*omega^3 + 
            r^3*(15*im + 5*im*lambda + 48*omega + 4*lambda*omega)) + 
            a^4*r^2*(144*im + 2*im*r*(-90 + 3*m^2 - 8*lambda) + 
            2*r^3*(-22 + 6*m^2 - 3*lambda + 15*im*omega)*omega - 
            36*im*r^4*omega^2 + 12*r^5*omega^3 + 
            r^2*(54*im - 4*im*m^2 + 11*im*lambda + 70*omega + 
            4*lambda*omega))))/((a^2 + r^2)^2*
            ((-r^4)*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*r*(im - 2*r*omega) - 4*a*m*r^2*
            (6*im + 2*im*r*lambda + 3*r^2*omega) + 
            12*a^4*(-1 - 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*r^2*(-3 + lambda)*omega + 3*r^3*omega^2)))
        end
        M22 = begin
            (r^2*sqrt((a^2 + r^2)/(a^2 + (-2 + r)*r)^2)*
            (-4*a^3*m*r*(-im + r*omega) - 2*a*m*r^2*
            (3*im - im*r + 2*r^2*omega) + 
            2*a^4*(-3 - 2*im*r*omega + r^2*omega^2) + 
            a^2*r*(24 + r*(-12 + 2*m^2 - lambda + 6*im*omega) - 
            12*im*r^2*omega + 4*r^3*omega^2) + 
            r^2*(-24 + 2*r*(12 + lambda) - 
            r^2*(6 + lambda - 18*im*omega) - 8*im*r^3*omega + 
            2*r^4*omega^2)))/
            ((-r^4)*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*r*(im - 2*r*omega) - 
            4*a*m*r^2*(6*im + 2*im*r*lambda + 3*r^2*omega) + 
            12*a^4*(-1 - 2*im*r*omega + 2*r^2*omega^2) + 
            4*a^2*r*(6 + r*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*r^2*(-3 + lambda)*omega + 3*r^3*omega^2))
        end
    else
        # Throw an error, this spin weight is not supported
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end 
    return (M11, M12, M21, M22)
end


# ------------------------------------------------------------------------------------------------
# Solve driver with the Eq.-(42) auto order-widening loop
# ------------------------------------------------------------------------------------------------

# upstream defaults [GSN main module l.25-28] + widening controls [l.600-603]
const _GSN_DEFAULT_RSIN = -50
const _GSN_DEFAULT_RSOUT = 1000
const _GSN_MAX_ORDER = 50
const _GSN_STEP_ORDER = 5
const _GSN_SANITY_TOL = 1e-6      # upstream check_XinXup_sanity tolerance

"""
    _gsn_homogeneous(s, l, m, a, omega, lambda; rsin, rsout, hor_order, inf_order, tol, form,
                     autowiden=true) -> (Xin::_GSNXSolution, Xup::_GSNXSolution)

Solve both boundary conditions and, when `autowiden`, bump both expansion orders by
+$( _GSN_STEP_ORDER) (up to $( _GSN_MAX_ORDER)) until the Eq.-(42) identity |Binc/Cinc −
p·c₀/(ω·η(r₊))| ≤ $( _GSN_SANITY_TOL) holds — the upstream real-ω auto loop [GSN main module
l.596-624] with the D6 defect (unbound rsin/rsout inside the loop) fixed by explicit binding.
"""
function _gsn_homogeneous(s::Int, l::Int, m::Int, a, omega, lambda;
                          rsin::Real=_GSN_DEFAULT_RSIN, rsout::Real=_GSN_DEFAULT_RSOUT,
                          hor_order::Int=3, inf_order::Int=6,
                          tol::Union{Nothing,Real}=nothing, form::Symbol=:auto,
                          autowiden::Bool=true)
    ho, io = hor_order, inf_order
    Xin = _gsn_solve_X(s, l, m, a, omega, lambda, :In;  rsin=rsin, rsout=rsout,
                       hor_order=ho, inf_order=io, tol=tol, form=form)
    Xup = _gsn_solve_X(s, l, m, a, omega, lambda, :Up; rsin=rsin, rsout=rsout,
                       hor_order=ho, inf_order=io, tol=tol, form=form)
    if autowiden
        while true
            actual, expected = _gsn_sanity_identity(Xin, Xup)
            # RELATIVE criterion: for large-|Binc/Cinc| (potential-barrier / high-l) modes the
            # upstream ABSOLUTE 1e-6 would spuriously widen forever (upstream never exercised
            # this path at high l — its widening loop carried the D6 UndefVarError). Documented
            # semantic fix, not an inherited behavior.
            abs(actual - expected) <= _GSN_SANITY_TOL*max(1, abs(expected)) && break
            (ho >= _GSN_MAX_ORDER || io >= _GSN_MAX_ORDER) &&
                error("GSN: Eq.-(42) sanity identity not met at maximum expansion orders " *
                      "(s=$s,l=$l,m=$m,a=$a,ω=$omega; |Δ|=$(abs(actual-expected))).")
            ho = min(ho + _GSN_STEP_ORDER, _GSN_MAX_ORDER)
            io = min(io + _GSN_STEP_ORDER, _GSN_MAX_ORDER)
            Xin = _gsn_solve_X(s, l, m, a, omega, lambda, :In;  rsin=rsin, rsout=rsout,
                               hor_order=ho, inf_order=io, tol=tol, form=form)
            Xup = _gsn_solve_X(s, l, m, a, omega, lambda, :Up; rsin=rsin, rsout=rsout,
                               hor_order=ho, inf_order=io, tol=tol, form=form)
        end
    end
    return Xin, Xup
end

# ------------------------------------------------------------------------------------------------
# Carrier glue
# ------------------------------------------------------------------------------------------------

const _GSN_DEMAND_LOCK = ReentrantLock()

"""
    _gsn_teukolsky_radial(s, l, m, a, omega, lambda, nu; kwargs...) -> (In=..., Up=...)

Build the unit-transmission `TeukolskyRadialFunction` pair from a GSN solve (the :GSN branch of
`TeukolskyRadial`). ν is carried for cross-backend comparability only — the GSN pipeline never
uses it (OQ-1 ruling: same treatment as the :SasakiNakamura backend, where the dispatcher
resolves ν regardless of backend). Amplitude/solution conversion per the file header.
"""
function _gsn_teukolsky_radial(s::Integer, l::Integer, m::Integer, a::Real, omega::Real, lambda, nu;
                               rsin::Real=_GSN_DEFAULT_RSIN, rsout::Real=_GSN_DEFAULT_RSOUT,
                               hor_order::Int=3, inf_order::Int=6,
                               tol::Union{Nothing,Real}=nothing, form::Symbol=:auto,
                               autowiden::Bool=true)
    lam = real(lambda)
    Xin, Xup = _gsn_homogeneous(Int(s), Int(l), Int(m), a, omega, lam;
                                rsin=rsin, rsout=rsout, hor_order=hor_order, inf_order=inf_order,
                                tol=tol, form=form, autowiden=autowiden)
    Bt = _gsn_Btrans(Int(s), Int(m), a, omega, lam); Bi = _gsn_Binc(Int(s), Int(m), a, omega, lam)
    Br = _gsn_Bref(Int(s), Int(m), a, omega, lam)
    Ct = _gsn_Ctrans(Int(s), Int(m), a, omega, lam); Ci = _gsn_Cinc(Int(s), Int(m), a, omega, lam)
    Cr = _gsn_Cref(Int(s), Int(m), a, omega, lam)
    mk(inc, trans, ref) = (Incidence=inc/trans, Transmission=one(inc/trans), Reflection=ref/trans)
    In = TeukolskyRadialFunction(s, l, m, float(a), float(omega), nu, lambda, :In,
                                 Bt, mk(Bi*Xin.inc_SN, Bt, Br*Xin.ref_SN), :GSN)
    In.cache = Xin
    Up = TeukolskyRadialFunction(s, l, m, float(a), float(omega), nu, lambda, :Up,
                                 Ct, mk(Ci*Xup.inc_SN, Ct, Cr*Xup.ref_SN), :GSN)
    Up.cache = Xup
    return (In=In, Up=Up)
end

"""
    _gsn_value(R::TeukolskyRadialFunction, r, deriv) -> Complex

Unit-transmission Teukolsky radial value (`deriv=0`) / first r-derivative (`deriv=1`) at `r`
for the :GSN backend: R(r), R′(r) = M(r)·[X(r*); X′(r*)]/norm, with the semianalytic GSN
solution from the carrier cache and the SN→Teukolsky conversion matrix M. [upstream
Teukolsky_radial_function_from_Sasaki_Nakamura_function ∘ GSN_solution / conv-factor]
"""
function _gsn_value(R::TeukolskyRadialFunction, r::Real, deriv::Int)
    sol = R.cache
    sol isa _GSNXSolution || error("GSN: carrier cache missing — construct via TeukolskyRadial(...; method=:GSN).")
    rs = _gsn_rstar_from_r(sol.a, oftype(sol.a, r))
    X, dX = lock(_GSN_DEMAND_LOCK) do
        _gsn_X_at(sol, rs)
    end
    M11, M12, M21, M22 = _gsn_conversion_matrix(sol.s, sol.m, sol.a, sol.omega, sol.lambda, oftype(sol.a, r))
    Rv  = (M11*X + M12*dX)/R.norm
    dRv = (M21*X + M22*dX)/R.norm
    return deriv == 0 ? Rv : dRv
end
