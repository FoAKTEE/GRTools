# Perturbed (1,1) matter-ring dynamics (PertGeodesiceqsolution, densitysubstrule, LHsnormaftersubst,
# BLtoNPruleU, Urhoamplitudes, ResultingfunctionsRing; Sourceincoordinates.nb, Polcar & Witzany
# 2025, arXiv:2507.15720). Closed-form rational functions of (MB,R,ω,mi,S1..S4) at equatorial
# Schwarzschild (θ=π/2, M=MB); ported verbatim from the saved notebook output cells.
#
# `S1..S4` are the four NP/GHP source legs (Sl, Sn, Sm, Sm†); the perturbed-orbit four-velocity
# coefficients (U1..U4) on the n/l/m/m† null tetrad solve LHSace{l,n,m,m†}0coord == -S{1,2,3,4}; the
# density response (drho, dtheta, dr), the continuity residual, the NP→Boyer-Lindquist tetrad map,
# and the deformed-ring (t,φ) reconstruction follow by substitution. Type-generic: the closed forms
# carry no element-type pins, so BigFloat / extended-precision inputs flow through unchanged.

# Solve LHSace{l,n,m,m†}0coord == -S{1,2,3,4} for the four-velocity coefficients (U1..U4).
function _pert_geodesic_solution(MB, R, ω, mi, S1, S2, S3, S4)
    Dcommon = 12*MB^2 + MB*(-4 + mi^2)*R - 2*sqrt(MB)*mi*ω*R^(5/2) + ω^2*R^4
    P = sqrt(MB)*mi - ω*R^(3/2)
    rt = sqrt(-3*MB + R)
    u1 = (rt*((8*im)*MB^3*S1 - (2*im)*ω^2*R^5*S2 + im*sqrt(MB)*ω*R^(7/2)*(4*mi*S2 + sqrt(2)*(S3 - S4)) + 2*MB^(3/2)*R^(3/2)*(2*mi*S2 + sqrt(2)*(1 - im*ω*R)*(S3 - S4)) + 4*sqrt(2)*MB^(5/2)*sqrt(R)*(-S3 + S4) - (2*im)*MB^2*R*(4*S1 + 4*S2 + sqrt(2)*mi*(-S3 + S4)) + im*MB*R^2*(2*S1 + (4 - 2*mi^2 + (4*im)*ω*R)*S2 + sqrt(2)*mi*(-S3 + S4)))) / (2*P*Dcommon)
    u2 = (R*rt*((-4*im)*MB^2*S1 - im*ω^2*R^4*S1 + MB*R*((2*im - im*mi^2 + 2*ω*R)*S1 + im*(4*S2 + sqrt(2)*mi*(S3 - S4))) - 2*MB^(3/2)*sqrt(R)*(mi*S1 + sqrt(2)*(-S3 + S4)) + im*sqrt(MB)*ω*R^(5/2)*(2*mi*S1 + sqrt(2)*(-S3 + S4)))) / (P*Dcommon)
    u3 = (rt*(-4*MB^(5/2)*sqrt(R)*S1 + 2*MB^(3/2)*R^(3/2)*(S1 + im*ω*R*S1 + 2*S2) + im*sqrt(2)*ω^2*R^5*S4 - im*sqrt(MB)*ω*R^(7/2)*(S1 - 2*S2 + 2*sqrt(2)*mi*S4) + im*MB*R^2*(mi*(S1 - 2*S2) + sqrt(2)*mi^2*S4 - 2*sqrt(2)*(S3 + S4)) - (2*im)*MB^2*R*(mi*S1 - 2*sqrt(2)*(S3 + 2*S4)))) / (sqrt(2)*P*Dcommon)
    u4 = (rt*(4*sqrt(2)*MB^(5/2)*sqrt(R)*S1 + 2*sqrt(2)*MB^(3/2)*R^(3/2)*((-1 - im*ω*R)*S1 - 2*S2) + (2*im)*ω^2*R^5*S3 + im*sqrt(MB)*ω*R^(7/2)*(sqrt(2)*S1 - 2*sqrt(2)*S2 - 4*mi*S3) + (2*im)*MB^2*R*(sqrt(2)*mi*S1 + 8*S3 + 4*S4) - im*MB*R^2*(sqrt(2)*mi*(S1 - 2*S2) - 2*mi^2*S3 + 4*(S3 + S4)))) / (2*P*Dcommon)
    return (u1, u2, u3, u4)
end

# Density response amplitudes (drho, dtheta, dr) from the substituted four-velocity (md = ring mass).
function _density_subst(MB, R, ω, mi, md, u1, u2, u3, u4)
    drho = (-1/4)*((md*(-3*MB + R)^(3/2)*(2*R*((-4*im)*MB + R*(2*im + ω*R))*u1 - im*(2*MB - R)*((4*MB + R*(-2 - im*ω*R))*u2 + sqrt(2)*mi*R*(-u3 + u4)))) / (pi*R^3*(-2*MB + R)^2*(-(sqrt(MB)*mi) + ω*R^(3/2))))
    dtheta = ((-im)*sqrt(-3*MB + R)*(u3 + u4)) / (sqrt(2)*(sqrt(MB)*mi - ω*R^(3/2)))
    dr = ((-im)*sqrt(-3*MB + R)*(2*R*u1 + (2*MB - R)*u2)) / (2*sqrt(MB)*mi - 2*ω*R^(3/2))
    return (drho, dtheta, dr)
end

# Continuity residual: LHSnormcoord after substituting the geodesic solution (≡ -Snorm if the
# source is exactly conserved). A diagnostic, not zero at machine precision for truncated sources.
function _norm_after_subst(MB, R, ω, mi, S1, S2, S3, S4)
    return ((2*im)*MB*sqrt(R)*S1 - im*R^(3/2)*(S1 + 2*S2) + sqrt(2)*sqrt(MB)*R*(-S3 + S4)) / (-(sqrt(MB)*mi) + ω*R^(3/2))
end

# NP-tetrad four-velocity (U1..U4 on n/l/m/m†) → Boyer-Lindquist components (Ut,Ur,Uθ,Uφ).
function _bl_to_np_U(MB, R, U1, U2, U3, U4)
    Ut = -((R*U1)/(2*MB - R)) + U2/2
    Ur = U1 + (-1/2 + MB/R)*U2
    Uθ = (U3 + U4)/(sqrt(2)*R)
    Uφ = (im*(U3 - U4))/(sqrt(2)*R)
    return (Ut=Ut, Ur=Ur, Uθ=Uθ, Uφ=Uφ)
end

"""
    ring_urho_amplitudes(mode::SourceMode) -> NamedTuple

Perturbed-ring four-velocity and density amplitudes for one source mode, porting `Urhoamplitudes`.
`mode` is a `.dat` source row (M=MB=1 Schwarzschild, equatorial). Returns
`(m, n, k, ω, radius, ringmass, Fourvelocity, density, normalisation)` where
`Fourvelocity = (U1,U2,U3,U4)` are the n/l/m/m† null-tetrad four-velocity coefficients (from
`_pert_geodesic_solution`), `density = (drho, dtheta, dr)` the linearised density response,
and `normalisation = (normSource, normSourcecheck)` pairs the mode's `Snorm` with the continuity
residual `_norm_after_subst`.
"""
function ring_urho_amplitudes(mode::SourceMode)
    R = mode.radius
    MB = one(R)
    ω = mode.ω; mi = mode.m; md = mode.ringmass
    S1 = mode.Sl; S2 = mode.Sn; S3 = mode.Sm; S4 = mode.Smdag
    u1, u2, u3, u4 = _pert_geodesic_solution(MB, R, ω, mi, S1, S2, S3, S4)
    drho, dtheta, dr = _density_subst(MB, R, ω, mi, md, u1, u2, u3, u4)
    normcheck = _norm_after_subst(MB, R, ω, mi, S1, S2, S3, S4)
    return (m=mode.m, n=mode.n, k=mode.k, ω=ω, radius=R, ringmass=md,
            Fourvelocity=(U1=u1, U2=u2, U3=u3, U4=u4),
            density=(drho=drho, dtheta=dtheta, dr=dr),
            normalisation=(normSource=complex(mode.Snorm), normSourcecheck=normcheck))
end

"""
    ring_functions(table) -> NamedTuple

Deformed equatorial-ring observables as functions of `(t, φ)` for a source `table`
(`Vector{SourceMode}`, one ring radius), porting `ResultingfunctionsRing`. For each mode it forms
[`ring_urho_amplitudes`](@ref) and sums the real reconstruction
`Σ_mode A·e^{-i(ωt-mφ)} + conj(A)·e^{+i(ωt-mφ)}` for each observable.

Returns `(radius, ringmass, FourvelocityNP, Fourvelocity, density, normalisation)` where every
field is a callable `(t,φ)`:
  - `FourvelocityNP = (Ul, Un, Um, Umdag)` — the n/l/m/m† null-tetrad four-velocity. As in upstream,
    the `Um`/`Umdag` reconstructions cross-pair under conjugation (the +c.c. term of `Um` uses
    `conj(U4)`, of `Umdag` uses `conj(U3)`), keeping the m/m̄ legs complex-conjugate.
  - `Fourvelocity = (Ut, Ur, Uθ, Uφ)` — Boyer-Lindquist components via `_bl_to_np_U`.
  - `density = (drho, dtheta, dr)` and `normalisation` (the `normSource + normSourcecheck`
    continuity residual reconstruction).
"""
function ring_functions(table::AbstractVector{SourceMode})
    isempty(table) && throw(ArgumentError("ring_functions: empty source table"))
    R = table[1].radius
    ringmass = table[1].ringmass
    all(mode -> mode.radius == R, table) ||
        throw(ArgumentError("ring_functions: all source modes must have the same radius"))
    all(mode -> mode.ringmass == ringmass, table) ||
        throw(ArgumentError("ring_functions: all source modes must have the same ringmass"))
    amps = [ring_urho_amplitudes(mode) for mode in table]
    MB = one(R)

    # Σ get(a)·e^{-iψ} + conj(ccget(a))·e^{+iψ},  ψ = ω t − m φ.
    function recon(get, ccget)
        return (t::Real, φ::Real) -> begin
            s = zero(get(amps[1]) * cis(zero(promote_type(typeof(t), typeof(φ)))))
            for a in amps
                ψ = a.ω * t - a.m * φ
                s += get(a) * cis(-ψ) + conj(ccget(a)) * cis(ψ)
            end
            return s
        end
    end

    U1 = recon(a -> a.Fourvelocity.U1, a -> a.Fourvelocity.U1)
    U2 = recon(a -> a.Fourvelocity.U2, a -> a.Fourvelocity.U2)
    U3 = recon(a -> a.Fourvelocity.U3, a -> a.Fourvelocity.U4)   # +c.c. uses U4
    U4 = recon(a -> a.Fourvelocity.U4, a -> a.Fourvelocity.U3)   # +c.c. uses U3
    drho   = recon(a -> a.density.drho,   a -> a.density.drho)
    dtheta = recon(a -> a.density.dtheta, a -> a.density.dtheta)
    dr     = recon(a -> a.density.dr,     a -> a.density.dr)
    Normdif = recon(a -> a.normalisation.normSource + a.normalisation.normSourcecheck,
                    a -> a.normalisation.normSource + a.normalisation.normSourcecheck)

    bl(t, φ) = _bl_to_np_U(MB, R, U1(t, φ), U2(t, φ), U3(t, φ), U4(t, φ))
    Ut(t, φ) = bl(t, φ).Ut
    Ur(t, φ) = bl(t, φ).Ur
    Uθ(t, φ) = bl(t, φ).Uθ
    Uφ(t, φ) = bl(t, φ).Uφ

    return (radius=R, ringmass=ringmass,
            FourvelocityNP=(Ul=U1, Un=U2, Um=U3, Umdag=U4),
            Fourvelocity=(Ut=Ut, Ur=Ur, Uθ=Uθ, Uφ=Uφ),
            density=(drho=drho, dtheta=dtheta, dr=dr),
            normalisation=Normdif)
end

export ring_urho_amplitudes, ring_functions
