# Sensitivity-only thermodynamics. Symbolic owner: derivation/src/thermal_balance.jl.
# Explicit import only (no exports): the SpiralWaves export set is gated by test/SpiralWaves/test_scaffold.jl;
# use SpiralWaves.ThermostattedLI, thermostat_state, thermostat_lagrangian and thermostat_bath qualified.

"""
    ThermostattedLI()

Sensitivity tag for conserved rest density rho with p=rho*T and
energy=rho*(1+T/(gamma_gas-1)), in c=1 units. T is a prescribed stationary,
axisymmetric scalar with delta T=0 in the fixed Boyer–Lindquist source gauge.
An isotropic comoving bath has J_bath^mu=Q*u^mu; it changes energy conservation.
This tag is deliberately rejected by primary DiscModel/ModeContext dispatch.
It supplies local thermodynamic diagnostics, not a primary LI mode solver.
"""
struct ThermostattedLI end

_closure_tag(::ThermostattedLI) = throw(ArgumentError(
    "ThermostattedLI is a sensitivity diagnostic; primary LI requires LocallyIsothermalOptionA"))

"""
    thermostat_state(::ThermostattedLI, rho, temperature; gamma_gas=5//3)

Rest density and temperature are positive/nonnegative, with 1<gamma_gas<=2.
Temperature denotes p/rho, not the primary p/energy. The fixed-temperature
compression index (rho/p)*partial(p)/partial(rho) is 1. Its acoustic speed
squared is p/(energy+p), whereas OPTION A at the same p/energy has speed
squared theta. At q=0 the primary compression index is 1+theta; C-LI-0 does
not apply to this variant. Scalars inherit precision; Float64 is the default.
"""
function thermostat_state(::ThermostattedLI, rho::Real, temperature::Real; gamma_gas::Real=5//3)
    rho, temperature, gamma_gas = float.(promote(rho,temperature,gamma_gas))
    isfinite(rho) && rho>0 || throw(DomainError(rho,"rho must be finite and positive"))
    isfinite(temperature) && temperature>=0 || throw(DomainError(temperature,"temperature must be finite and nonnegative"))
    isfinite(gamma_gas) && 1<gamma_gas<=2 || throw(DomainError(gamma_gas,"require 1<gamma_gas<=2"))
    heat_capacity = inv(gamma_gas-one(gamma_gas))
    pressure = rho*temperature
    energy = rho*(one(temperature)+heat_capacity*temperature)
    (;rho,temperature,gamma_gas,heat_capacity,pressure,energy,
      theta=pressure/energy,sound_speed2=pressure/(energy+pressure),gamma1_eff=one(temperature))
end

"""
    thermostat_lagrangian(::ThermostattedLI, state, Delta_rho, xi_r, temperature_r)

Lagrangian variations: Delta T=xi^r*T', Delta p=T*Delta rho+rho*Delta T,
Delta energy=(energy/rho)*Delta rho+rho*c_v*Delta T. Setting xi_r=0 also
returns the Eulerian fixed-temperature variations for Delta_rho=delta rho.
No zero-bath G26 relation (19) is imposed on these variations.
"""
function thermostat_lagrangian(::ThermostattedLI, state, Delta_rho, xi_r, temperature_r)
    Delta_temperature = xi_r*temperature_r
    (;pressure=state.temperature*Delta_rho+state.rho*Delta_temperature,
      energy=state.energy/state.rho*Delta_rho+state.rho*state.heat_capacity*Delta_temperature)
end

"""
    thermostat_bath(::ThermostattedLI, state, uup, ulow;
                    expansion, temperature_advection)

Q=p*div(u)+rho*c_v*u.grad(T); positive Q heats the fluid. For T(r),
temperature_advection=u^r*T'. uup and ulow are normalized contravariant and
covariant four-velocity components in (t,r,phi,z) order, lowered in the same
metric. J_bath=Q*u, comoving heat=-u.J_bath=Q, spatial bath force zero.
Killing energy work is -Q*u_t and angular-momentum work Q*u_phi. These are
reservoir contributions only: gravitational work (-T.Lie_t(g)/2) and torque
(+T.Lie_phi(g)/2) must be added separately. Background circular flow has Q=0.

Assumptions: perfect fluid, conserved baryons, instantaneous prescribed bath,
constant heat capacity, no conduction/viscosity or bath momentum in the rest
frame. This is a finite-GR phenomenological completion, not radiation transport.
"""
function thermostat_bath(::ThermostattedLI, state, uup, ulow; expansion, temperature_advection)
    length(uup)==length(ulow)==4 || throw(DimensionMismatch("four-velocity needs four components"))
    norm = sum(uup[i]*ulow[i] for i in 1:4)
    isapprox(norm,-one(norm);rtol=64eps(typeof(float(norm))),atol=zero(norm)) ||
        throw(ArgumentError("uup and ulow must satisfy u.u=-1"))
    heat = state.pressure*expansion+state.rho*state.heat_capacity*temperature_advection
    force_up = ntuple(i->heat*uup[i],4)
    force_low = ntuple(i->heat*ulow[i],4)
    (;heat,force_up,force_low,killing_energy=-force_low[1],killing_angular_momentum=force_low[3])
end
