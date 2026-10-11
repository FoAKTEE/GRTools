# FR24 v1: ref-paper/arXiv-2407.20398/src/Sections/2_linear_framework.tex:268-317.
# NewtonianDisc units: GM = a = Omega_p = Sigma_p = 1 for the FR24 reproduction.
using ..ForcingModes: mode_power

"""
    fr24_corotation_radius(d, mode)

FR24 (35), the O(h_p²) radius in the disc's length units. `mode.omega` is the
TOTAL frequency; Omega(r_C) = mode.omega / mode.m. This is an approximation,
not the grid's corotation index. Use `bisect(r -> Omega(d,r)-mode.omega/m,...)`
for an independent exact-root diagnostic. No background profile is duplicated.
"""
function fr24_corotation_radius(d::NewtonianDisc, mode::HarmonicMode)
    mode.m > 0 && mode.omega > 0 || throw(ArgumentError("positive m and frequency required"))
    frequency, _ = promote(mode.omega, d.h0)
    r0 = (frequency / mode.m)^(-2//3)
    r0 * (1 - d.h0^2 * (d.p + d.q) * r0^(1-d.q) / 3)
end

"""
    fr24_torque_profile(r, density; a=one(eltype(r)))

FR24 (33): T(r) = integral_a^r density(s) ds. Piecewise-linear density,
trapezoidal cell integrals; the reference a may lie inside a cell. The return
includes the same interpolant's integral at any requested point in the domain.
"""
function fr24_torque_profile(r, density; a=one(eltype(r)))
    length(r) == length(density) >= 2 && all(diff(r) .> 0) ||
        throw(ArgumentError("matching arrays on a strictly increasing grid required"))
    first(r) <= a <= last(r) || throw(ArgumentError("reference radius is outside grid"))
    cumulative = vcat(zero(first(density)*first(r)), cumsum(diff(r) .* (density[1:end-1] .+ density[2:end]) ./ 2))
    integral_at(s) = begin
        first(r) <= s <= last(r) || throw(ArgumentError("integration point is outside grid"))
        j = clamp(searchsortedlast(r,s),1,length(r)-1)
        dx = s-r[j]
        cumulative[j] + dx*density[j] + dx^2*(density[j+1]-density[j])/(2*(r[j+1]-r[j]))
    end
    origin = integral_at(a)
    (; torque=cumulative .- origin, at=s -> integral_at(s)-origin)
end

"""
    fr24_mode_diagnostics(d, mode, r, f, cache; halfwidth=1//200)

Orbit-averaged FR24 (31)-(39), using the SAME enthalpy and PotentialCache as
`solve_lm_mode`. Returns flux, density, torque, net, corotation, lindblad, power,
planet_power, radius and FJ0. Torque and power have the disc sign; reaction
power has the opposite sign. FJ0=h0^(-3) for G M_p=Sigma_p=a=Omega_p=1.
No normalization is applied to the arrays. beta=0 is required.

The printed radius and 0.005a half-window are reproduction conventions.
Corotation uses a piecewise-linear flux interpolant and the exact integral of
piecewise-linear torque density; no extrapolation is permitted. This split is
FR24's finite-window diagnostic, not a proof of pure corotation deposition at
q!=0. The simple balance F_J'=dT/dr applies only at q=0, away from resonances
and the pole displacement. General LI balance includes thermal exchange.
"""
function fr24_mode_diagnostics(d::NewtonianDisc, mode::HarmonicMode, r, f, cache::PotentialCache;
                               halfwidth=1//200)
    iszero(d.beta) || throw(ArgumentError("FR24 requires locally isothermal beta=0"))
    rc = fr24_corotation_radius(d,mode)
    halfwidth > 0 && first(r) <= rc-halfwidth < rc+halfwidth <= last(r) ||
        throw(ArgumentError("the FR24 corotation window must lie inside the grid"))
    ur, up, sig = vars_from_enthalpy(d,f,r,mode,cache)
    phi = first.(potential_at.(Ref(cache),r))
    flux = pi .* r.^2 .* Sigma.(Ref(d),r) .* real.(ur .* conj.(up))
    density = -pi .* r .* mode.m .* phi .* imag.(sig)
    profile = fr24_torque_profile(r,density)
    flux_at(s) = begin
        j = clamp(searchsortedlast(r,s),1,length(r)-1)
        flux[j] + (s-r[j])/(r[j+1]-r[j])*(flux[j+1]-flux[j])
    end
    corotation = profile.at(rc+halfwidth)-profile.at(rc-halfwidth) -
                 (flux_at(rc+halfwidth)-flux_at(rc-halfwidth))
    net = profile.torque[end]-profile.torque[1]
    power = mode_power(mode,net)
    (; flux, density, torque=profile.torque, net, corotation, lindblad=net-corotation,
       power, planet_power=-power, radius=rc, FJ0=inv(d.h0^3))
end

"""
    fr24_instantaneous_torque(r, modes, densities, potentials, t)

FR24 (28), with real Phi amplitudes: cross-l terms survive instantaneously at
fixed m. `densities` and `potentials` are vectors of radial vectors, in mode
order, from the same linear responses and caches. Time is in the mode's units.
This direct pair sum is intended for small diagnostic mode sets.
"""
function fr24_instantaneous_torque(r, modes, densities, potentials, t)
    all(m -> m.m > 0,modes) || throw(ArgumentError("FR24 sums m>=1"))
    out = zeros(promote_type(eltype(r),typeof(real(first(first(densities))))),length(r))
    for (i, mi) in pairs(modes), (j,mj) in pairs(modes)
        mi.m == mj.m || continue
        phase = (mi.omega-mj.omega)*t
        out .-= pi .* r .* mi.m .* potentials[j] .*
            (cos(phase).*imag.(densities[i]) .- sin(phase).*real.(densities[i]))
    end
    out
end
