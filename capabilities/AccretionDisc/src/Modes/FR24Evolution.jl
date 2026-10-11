# FR24 v1, ref-paper/arXiv-2407.20398/src/Sections/4_eccentric_torques.tex:105-160.
# Shared algebra also owns FD25 2_linear_theory.tex:155-182 (DESIGN v3 §6).

# Kept free of validation/conversion so derivation_fr can certify the runtime expressions.
function _fr24_kepler_rates(Mp, a, Omega_p, e, FJ0, torque, power)
    L = Mp*Omega_p*a^2*sqrt(1-e^2)
    E = -Mp*a^2*Omega_p^2/2
    tau0_inv = FJ0/(Mp*Omega_p*a^2)
    tauL_inv = torque/L
    taua_inv = 2power/(Mp*a^2*Omega_p^2)
    # e=0 is a convention for the circular slice, not a determined limit of -edot/e.
    taue_inv = iszero(e) ? zero(taua_inv) : (1-e^2)/e^2*(taua_inv/2-tauL_inv)
    (; L, E, torque, power, Ldot=-torque, Edot=-power, tau0_inv, tauL_inv,
       taua_inv, taue_inv, adot=-a*taua_inv, edot=-e*taue_inv)
end

"""
    fr24_orbital_rates(diagnostics; Mp=1, a=1, Omega_p=1, e=0)

FR24 (45)-(51): Keplerian rates from a nonempty collection of
`fr24_mode_diagnostics` outputs. Sum their `net` and `power` fields; no torque
or modal power is recomputed. All modes must share the same positive `FJ0` and
units. `Mp`, `a`, `Omega_p` are in those units, with GM=Omega_p²*a³ fixed.
Positive disc torque/power gives negative planet `Ldot`/`Edot`.

Return total planet L,E, disc torque,power, planet Ldot,Edot, inverse timescales
tau0_inv,tauL_inv,taua_inv,taue_inv and element derivatives adot,edot.
Defaults produce Float64; BigFloat inputs retain their type and precision.
At e=0, require power=Omega_p*torque to summation roundoff, and set taue_inv=0
as in the published circular slice. This does not determine the e→0 damping
limit. For small nonzero e, the subtraction in taue_inv can lose accuracy;
this algebra alone is not a small-e convergence certificate.
"""
function fr24_orbital_rates(diagnostics; Mp=1, a=1, Omega_p=1, e=0)
    isempty(diagnostics) && throw(ArgumentError("at least one modal diagnostic required"))
    FJ0 = first(diagnostics).FJ0
    all(d -> d.FJ0 == FJ0, diagnostics) || throw(ArgumentError("mode normalizations differ"))
    Mp,a,Omega_p,e,FJ0,torque,power = promote(float(Mp),float(a),float(Omega_p),
        float(e),float(FJ0),float(sum(d.net for d in diagnostics)),float(sum(d.power for d in diagnostics)))
    Mp>0 && a>0 && Omega_p>0 && FJ0>0 && 0<=e<1 &&
        all(isfinite,(Mp,a,Omega_p,e,FJ0,torque,power)) ||
        throw(ArgumentError("finite positive mass, radius, frequency and FJ0, and 0<=e<1 required"))
    if iszero(e)
        scale = sum(abs(d.power) + abs(Omega_p*d.net) for d in diagnostics)
        abs(power-Omega_p*torque) <= 64eps(typeof(power))*scale ||
            throw(ArgumentError("circular power and torque violate the Keplerian identity"))
    end
    _fr24_kepler_rates(Mp,a,Omega_p,e,FJ0,torque,power)
end
