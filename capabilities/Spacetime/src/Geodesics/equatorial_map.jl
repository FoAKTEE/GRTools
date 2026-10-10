"""
    EL_to_pe_jacobian(a, p, e; x=1, max_condition=nothing)

Constrained `∂(p,e)/∂(E,Lz)` at fixed a and equatorial x=±1, in M=1 units.
Contracts the Geodesics owner with dK=2(Lz-aE)(dLz-a dE), so dx=0.
Requires |a|<1, 0<e<1 and p above the separatrix. Circular turning points
are degenerate: this function deliberately does not evaluate the e=0 limit.

Rejects an infinity-norm condition number above `max_condition` (default
`1/sqrt(eps(T))`). This coordinate-dependent bound uses the stated M=1
coordinates, not a relative error bound on a nearly vanishing rate.
"""
function EL_to_pe_jacobian(a::Real, p::Real, e::Real; x::Real=1, max_condition=nothing)
    a,p,e,x = promote(float(a),float(p),float(e),float(x))
    abs(a)<1 && abs(x)==1 && 0<e<1 && p>separatrix(a,e,x) ||
        throw(ArgumentError("EL_to_pe_jacobian requires a stable eccentric equatorial orbit"))
    E,L,_ = constants_of_motion(a,p,e,x)
    z,u = zero(E),one(E)
    chain = [u z; z u; -2a*(L-a*E) 2*(L-a*E)]
    J = (orbital_params_jacobian(a,p,e,x)*chain)[1:2,:]
    limit = max_condition === nothing ? inv(sqrt(eps(typeof(E)))) : max_condition
    κ = opnorm(J,Inf)*opnorm(inv(J),Inf)
    isfinite(κ) && κ <= limit ||
        throw(ArgumentError("EL_to_pe_jacobian condition number $κ exceeds $limit"))
    return J
end

"""
    eccentric_mode_numerator(Ωφ, Ωr, alpha, modes; axisymmetric_power=0, mu=1)
    eccentric_mode_numerator(orbit, alpha, modes; kwargs...)

Stable rewrite of `Edot-alpha*Ldot` for the reaction on the perturber.
Each mode has integer `m>0`, integer `k`, and signed disc `torque`;
`alpha=∂pE/∂pL`. `axisymmetric_power` is the sum of the signed m=0,k>0
disc powers. With Edot=-P/mu and Ldot=-T/mu, returns
`-[Σ(k Ωr/m)T_mk + (Ωφ-alpha)ΣT_mk + ΣP_0k]/mu`.
No division by m=0 or ω occurs. Mode selection and reality-pair weights
belong to the caller. This identity does not resolve the e=0 0/0 in de/e.
"""
function eccentric_mode_numerator(Ωφ, Ωr, alpha, modes; axisymmetric_power=0, mu=1)
    mu > 0 || throw(ArgumentError("mu must be positive"))
    torque = zero(Ωφ*axisymmetric_power)
    sideband = zero(torque)
    for mode in modes
        mode.m isa Integer && mode.m>0 && mode.k isa Integer ||
            throw(ArgumentError("torque modes require integer m>0 and integer k"))
        torque += mode.torque
        sideband += (mode.k*Ωr/mode.m)*mode.torque
    end
    return -(sideband+(Ωφ-alpha)*torque+axisymmetric_power)/mu
end

eccentric_mode_numerator(o::EquatorialOrbit, alpha, modes; kwargs...) =
    eccentric_mode_numerator(orbit_frequency(o,1,0),orbit_frequency(o,0,1),alpha,modes;kwargs...)

"""
    ELdot_to_pedot(orbit, Edot, Ldot; eccentric_numerator=nothing,
                  max_condition=nothing, max_fraction_per_radial_period=0.01)

Map specific constants rates to `(pdot,edot,e2dot,alpha,fraction_per_radial_period)`.
For small e, supply the independently accumulated `eccentric_mode_numerator`
instead of subtracting Edot-alpha*Ldot. The caller must use the same rates,
mass and mode set in both arguments. The orbit has positive Lz (x=1).

Stops when max(|pdot|/p,|edot|/e) times the radial coordinate-time period
exceeds the supplied averaging threshold (default 1%). The condition guard
and the stable-orbit domain check also apply. e=0 is rejected: a circular
continuation needs a separately specified smooth power/torque law.
"""
function ELdot_to_pedot(o::EquatorialOrbit, Edot, Ldot;
                       eccentric_numerator=nothing, max_condition=nothing,
                       max_fraction_per_radial_period=one(o.p)/100)
    J = EL_to_pe_jacobian(o.a,o.p,o.e;max_condition)
    alpha = -J[2,2]/J[2,1]
    numerator = eccentric_numerator === nothing ? Edot-alpha*Ldot : eccentric_numerator
    edot = J[2,1]*numerator
    pdot = J[1,1]*Edot+J[1,2]*Ldot
    fraction = 2oftype(o.p,pi)/o.Ωr*max(abs(pdot)/o.p,abs(edot)/o.e)
    fraction <= max_fraction_per_radial_period ||
        throw(ArgumentError("element change per radial period $fraction exceeds averaging threshold"))
    return (; pdot,edot,e2dot=2o.e*edot,alpha,fraction_per_radial_period=fraction)
end
