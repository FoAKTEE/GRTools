using NumericsCore: brent

"Kepler ellipse with semimajor axis `a > 0`, eccentricity `0 ≤ e < 1`, and mean motion `n > 0`."
struct EccentricOrbit{T<:AbstractFloat}
    a::T
    e::T
    n::T
    function EccentricOrbit(a::Real, e::Real, n::Real)
        T = float(promote_type(typeof(a), typeof(e), typeof(n)))
        aa, ee, nn = T(a), T(e), T(n)
        isfinite(aa) && aa > 0 || throw(DomainError(a, "semimajor axis must be finite and positive"))
        _kepler_eccentricity(ee)
        isfinite(nn) && nn > 0 || throw(DomainError(n, "mean motion must be finite and positive"))
        new{T}(aa, ee, nn)
    end
end

function _kepler_eccentricity(e)
    isfinite(e) && 0 <= e < 1 || throw(DomainError(e, "elliptic eccentricity must satisfy 0 ≤ e < 1"))
    e
end

"""
    kepler_E(M_anom, e; maxiter=32)

Solve `E - e*sin(E) = M_anom`, retaining the turns of the mean anomaly.
Newton iteration on the positive principal anomaly preserves odd symmetry;
an escaped or unconverged iterate falls back to NumericsCore's bracketed Brent
solver. `maxiter` is the Newton budget (zero exercises the bracketed path).
All tolerances inherit the promoted floating type and current precision.
Accuracy is absolute: the residual is a few `eps` and the error in `E` is at the
conditioning limit `eps/(1 - e cos E)`, so for `|E| ≪ 1` as `e → 1` it is not relative.
"""
function kepler_E(M_anom::Real, e::Real; maxiter::Integer=32)
    T = float(promote_type(typeof(M_anom), typeof(e)))
    M, ecc = T(M_anom), _kepler_eccentricity(T(e))
    isfinite(M) || throw(DomainError(M_anom, "mean anomaly must be finite"))
    maxiter >= 0 || throw(ArgumentError("Newton iteration budget must be nonnegative"))
    (iszero(ecc) || iszero(M)) && return M
    piT = T(pi)
    principal = rem(M, 2piT, RoundNearest)
    m = abs(principal)
    iszero(m) && return M
    E, steps = _kepler_newton(m, ecc, maxiter)
    if steps === nothing
        E = brent(x -> x - ecc * sin(x) - m, zero(T), piT; atol=eps(T), rtol=eps(T),
                  maxiter=4precision(T))
    end
    (M - principal) + copysign(E, principal)
end

"""
    _kepler_newton(m, e, maxiter) -> (E, steps)

Newton's method for `E - e*sin(E) = m` on the principal interval `0 < m ≤ π`, started
at `m` (`e < 4/5`) or at `π`. Returns the last iterate and the number of steps taken,
or `steps = nothing` when an iterate leaves `[0, π]` or `maxiter` is exhausted; `kepler_E`
then brackets with Brent. The absolute step test ignores the conditioning
`1/(1 - e cos E)`: near pericentre as `e → 1` (e.g. BigFloat-256, `e ≥ 0.9`) the step can
stall at the rounding floor, and the bracket supplies the root.
"""
function _kepler_newton(m::T, ecc::T, maxiter::Integer) where {T<:AbstractFloat}
    piT = T(pi)
    E = ecc < T(4)/5 ? m : piT
    for k in 1:maxiter
        step = (E - ecc * sin(E) - m) / (1 - ecc * cos(E))
        next = E - step
        0 <= next <= piT || return (E, nothing)
        E = next
        abs(step) <= eps(T) * max(one(T), E) && return (E, k)
    end
    (E, nothing)
end

"""
    planet_orbit(M_anom, e; a=1, n=1)
    planet_orbit(M_anom, orbit::EccentricOrbit)

FR22 (24)–(25), with pericentre at `M_anom = 0`. Return `(r_p, phi_p, E,
dr_dt, dphi_dt, dE_dt)`; `E` and `phi_p` are continuous, unwrapped angles.
Time derivatives use `dM_anom/dt = n`. The shared Kepler certificate in
`derivation_fr/src/kepler.jl` specifies the radius, parity and time measure.
"""
function planet_orbit(M_anom::Real, orbit::EccentricOrbit)
    T = float(promote_type(typeof(M_anom), typeof(orbit.a)))
    a, e, n = T(orbit.a), T(orbit.e), T(orbit.n)
    E = kepler_E(T(M_anom), e)
    Er = rem(E, 2T(pi), RoundNearest)
    rho = 1 - e * cos(Er)
    f = 2atan(sqrt(1 + e) * sin(Er/2), sqrt(1 - e) * cos(Er/2))
    phi = iszero(e) ? T(M_anom) : E + (f - Er)
    (r_p=a*rho, phi_p=phi, E=E, dr_dt=a*n*e*sin(Er)/rho,
     dphi_dt=n*sqrt((1-e)*(1+e))/rho^2, dE_dt=n/rho)
end

planet_orbit(M_anom::Real, e::Real; a::Real=1, n::Real=1) =
    planet_orbit(M_anom, EccentricOrbit(a, e, n))
