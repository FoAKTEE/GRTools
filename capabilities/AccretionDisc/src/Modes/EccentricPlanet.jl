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

"""
    lindblad_locations(orbit::EccentricOrbit, m, l) -> (r_minus, r_plus)

Pressureless Newtonian FR22 (26) radii for positive-frequency modes, `m,l ≥ 1`.
The minus branch is inner and the plus branch outer. For `m=1` the inner
resonance is absent (`nothing`), not a radius at zero. This overload does not
include pressure or relativistic frequency shifts. The common mean motion
cancels; radii inherit the orbit's floating type and precision.
"""
function lindblad_locations(orbit::EccentricOrbit{T}, m::Integer, l::Integer) where {T}
    m >= 1 && l >= 1 || throw(ArgumentError("pressureless FR22 radii require m,l ≥ 1"))
    inner = m == 1 ? nothing : orbit.a * cbrt(T(m - 1) / T(l))^2
    outer = orbit.a * cbrt(T(m + 1) / T(l))^2
    inner, outer
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

_zz22_B(x) = 3 - x - 2/sqrt(x)

"""
    zz22_wavelet(t1, t2, orbit, dir; cs)

ZZ22 (4)–(9) for a Keplerian, globally isothermal disc (`alpha=3/2`,
`beta=0`). A wavelet emitted at `t1` retains the planet's instantaneous
angular velocity at emission. `dir=-1/+1` selects inward/outward travel;
`cs` is the constant sound speed in the orbit's length/time units.
Return `(R, theta)` in the inertial frame, or `nothing` when `R ≤ 0`.
No amplitude or disc response is calculated. Angles remain unwrapped.
"""
function zz22_wavelet(t1::Real, t2::Real, orbit::EccentricOrbit, dir::Integer; cs::Real)
    dir in (-1, 1) || throw(ArgumentError("wavelet direction must be -1 or +1"))
    isfinite(cs) && cs > 0 || throw(DomainError(cs, "sound speed must be finite and positive"))
    isfinite(t1) && isfinite(t2) && t2 >= t1 ||
        throw(DomainError((t1, t2), "finite observation time must follow emission"))
    T = float(promote_type(typeof(t1), typeof(t2), typeof(orbit.a), typeof(cs)))
    age = T(t2) - T(t1)
    p = planet_orbit(orbit.n*T(t1), orbit)
    R = p.r_p + dir*T(cs)*age
    R > 0 || return nothing
    Rc = orbit.a * (orbit.n/p.dphi_dt)^(T(2)/3)
    pref = orbit.n*orbit.a * (Rc/orbit.a)^(-one(T)/2) / T(cs)
    theta_c = p.phi_p - sign(p.r_p-Rc)*pref*_zz22_B(p.r_p/Rc)
    (R=R, theta=p.dphi_dt*age + theta_c + sign(R-Rc)*pref*_zz22_B(R/Rc))
end

"""
    zz22_track(orbit, t2, dir; cs, dt, t_start=0)

Wavelets at `t2`, newest first, emitted at `t2-k*dt ≥ t_start` (to endpoint
roundoff). Return `(k, log10(R/a), mod(theta-phi_p(t2),2pi)/pi)` tuples.
Nonpositive radii are culled; no radial plot clipping is applied. Plotters
must break lines at gaps in `k` or jumps in the wrapped ordinate exceeding 1.
The FR22 driver selects the cadence and fifth-orbit snapshots from FR22Setup.
"""
function zz22_track(orbit::EccentricOrbit, t2::Real, dir::Integer;
                    cs::Real, dt::Real, t_start::Real=0)
    isfinite(dt) && dt > 0 || throw(DomainError(dt, "emission interval must be finite and positive"))
    isfinite(t_start) && isfinite(t2) && t2 >= t_start ||
        throw(DomainError((t_start, t2), "invalid emission window"))
    T = float(promote_type(typeof(orbit.a), typeof(t2), typeof(cs), typeof(dt), typeof(t_start)))
    time = T(t2)
    piT = T(pi)
    phi2 = planet_orbit(orbit.n*time, orbit).phi_p
    count = (time-T(t_start))/T(dt)
    kmax = floor(Int, count + 8eps(T)*max(one(T), count))
    out = Tuple{Int,T,T}[]
    for k in 0:kmax
        w = zz22_wavelet(time-k*T(dt), time, orbit, dir; cs=T(cs))
        w === nothing && continue
        push!(out, (k, log10(w.R/orbit.a), mod(w.theta-phi2, 2piT)/piT))
    end
    out
end

# All-mode Newtonian source bank (DESIGN §3). The circular cache/interpolator
# remains owned by NewtonianSolver.jl; no solver or FFT implementation lives here.
using NumericsCore: fft_radix2!
using LinearAlgebra: mul!
using ..ForcingModes: HarmonicMode, SofteningPolicy, FixedSoftening,
    SemimajorScaleHeight, InstantaneousScaleHeight, LocalScaleHeight

"""
Storage returned by [`potential_bank`](@ref). `caches[eindex, modeindex]` owns
ordinary `PotentialCache` vectors; every cache shares the bank's log-radius grid.
`key` records physical inputs, softening policy, grids, precision and backend.
There is no global memoization. Treat the bank and its cache vectors as read-only.
"""
struct NewtonianPotentialBank{T,K}
    caches::Matrix{PotentialCache{T}}
    mode_indices::Dict{HarmonicMode{T},Int}
    key::K
end

"""
    harmonic_cache(bank, eindex, mode::HarmonicMode) -> PotentialCache

Return the bank-owned per-mode cache (the same object on every call). No copy,
quadrature or interpolation is done here. These are zero-copy per-mode handles:
the existing `PotentialCache` has concrete Vector fields, not SubArray fields.
The full frequency participates in lookup, preventing accidental reuse at a
different mean motion. `potential_at(cache,r)` provides the existing interpolation.
"""
function harmonic_cache(bank::NewtonianPotentialBank{T}, eindex::Integer,
                        mode::HarmonicMode) where {T}
    bank.caches[eindex, bank.mode_indices[HarmonicMode{T}(mode.m, mode.k, mode.omega)]]
end

_bank_softvalue(s::FixedSoftening) = s.eps
_bank_softvalue(s::Union{SemimajorScaleHeight,InstantaneousScaleHeight,LocalScaleHeight}) = s.b
_bank_softlength(s::FixedSoftening, h, power, a, rp, r) = s.eps
_bank_softlength(s::SemimajorScaleHeight, h, power, a, rp, r) = s.b*h*a^power
_bank_softlength(s::InstantaneousScaleHeight, h, power, a, rp, r) = s.b*h*rp^power
_bank_softlength(s::LocalScaleHeight, h, power, a, rp, r) = s.b*h*r^power

"""
    potential_bank(d, orbits, rs, modes; softening=FixedSoftening(d.soft),
                   Nt=1024, Nphi=1024, derivatives=2, method=:direct, Gmu=1)

Project the direct Plummer potential and its first two analytic radial derivatives
for every orbit, radius and requested `HarmonicMode`. `modes` uses the shared
Newtonian adapter (`l=m+k`, `omega=l*n`); signed l and m=0 are supported. Orbits
must share the mean motion implied by the requested frequencies. Radii must be
positive and strictly increasing. All inputs are promoted, Float64 by default;
BigFloat inherits ambient precision. `derivatives=2` is the only supported value,
so every returned cache is ready for the existing second-order master equation.

Both methods sample endpoint-exclusive `(M,phi)=2pi*(i/Nt,j/Nphi)`. `:fft` uses
NumericsCore's 2-D transform and returns `2real(F[mod(-l,Nt)+1,m+1])/(Nt*Nphi)`.
`:direct` evaluates those same bins by batched real matrix products in phi and
direct time sums. It is the default for the measured FR22-size masks; choose
`:fft` for dense spectra. FFT sizes must be powers of two; both methods require
strict sub-Nyquist indices (`2abs(l)<Nt`, `2m<Nphi`). Node refinement is the
caller's convergence obligation; this constructor does not assert convergence.

FR22's corrected normalization is `n/(2pi^2)` in dt, or `1/(2pi^2)` in dM
(E-ECC-1/2; derivation_fr/src/fourier.jl and phi_factorisation.jl). The factor 2
also applies to (0,0) as a Fourier-bin convention; a mean term in reconstruction
must use half this value. No sector or paper mask is silently discarded.

`Gmu` is G times perturber mass; lengths use the disc's reference radius (one).
Scale-height policies use `H(x)=d.h0*x^((3-d.q)/2)`. Field-radius softening
includes both radial derivatives of H. FR22's ambiguous policy must be selected
explicitly by its caller. `d.indirect != 0` is rejected: this API implements the
direct potential only, not the extra indirect tables absent from FR24 v1.

Memory is O(Nt*Nphi + Nt*(unique(m)+unique(l)) + length(orbits)*length(rs)*
length(modes)); temporary arrays are reused across radii and orbits. No global
precision or BLAS/thread configuration is changed here.
"""
function potential_bank(d::NewtonianDisc, orbits, rs, modes;
                        softening::SofteningPolicy=FixedSoftening(d.soft),
                        Nt::Integer=1024, Nphi::Integer=1024, derivatives=2,
                        method::Symbol=:direct, Gmu::Real=1)
    obs, radii, requested = collect(orbits), collect(rs), collect(modes)
    isempty(obs) && throw(ArgumentError("at least one orbit is required"))
    isempty(radii) && throw(ArgumentError("at least one radius is required"))
    derivatives == 2 || throw(ArgumentError("PotentialCache requires derivatives=2"))
    method in (:fft, :direct) || throw(ArgumentError("method must be :fft or :direct"))
    Nt > 0 && Nphi > 0 || throw(ArgumentError("node counts must be positive"))
    method == :fft && !(ispow2(Nt) && ispow2(Nphi)) &&
        throw(ArgumentError("FFT node counts must be powers of two"))
    d.indirect == 0 || throw(ArgumentError("the bank implements only the direct potential"))
    T = float(promote_type(typeof(d.h0), typeof(d.q), typeof(Gmu),
        typeof(_bank_softvalue(softening)), mapreduce(typeof, promote_type, radii),
        mapreduce(o -> typeof(o.a), promote_type, obs),
        mapreduce(m -> typeof(m.omega), promote_type, requested; init=Int)))
    rg = T.(radii)
    all(r -> isfinite(r) && r > 0, rg) && all(diff(rg) .> 0) ||
        throw(ArgumentError("radii must be finite, positive and strictly increasing"))
    isfinite(Gmu) && Gmu >= 0 || throw(DomainError(Gmu, "Gmu must be finite and nonnegative"))
    isfinite(d.h0) && d.h0 >= 0 && isfinite(d.q) ||
        throw(ArgumentError("finite q and finite nonnegative h0 are required"))
    os = [EccentricOrbit(T(o.a), T(o.e), T(o.n)) for o in obs]
    ms = [HarmonicMode{T}(m.m, m.k, m.omega) for m in requested]
    ls = [Base.Checked.checked_add(m.m, m.k) for m in ms]
    for (m, l) in zip(ms, ls)
        -Nt/2 < l < Nt/2 && m.m < Nphi/2 ||
            throw(ArgumentError("requested modes must lie strictly below Nyquist"))
        all(o -> isapprox(m.omega, l*o.n; rtol=8eps(T), atol=zero(T)), os) ||
            throw(ArgumentError("mode frequency must equal (m+k)*orbit.n"))
    end
    indices = Dict(m => i for (i,m) in enumerate(ms))
    length(indices) == length(ms) || throw(ArgumentError("duplicate modes"))
    lg = log10.(rg)
    caches = [PotentialCache{T}(lg, zeros(T,length(rg)), zeros(T,length(rg)),
                               zeros(T,length(rg))) for _ in os, _ in ms]
    key = (; scalar_type=T, bits=precision(T), orbits=Tuple((o.a,o.e,o.n) for o in os),
        radii=Tuple(rg), modes=Tuple(ms), softening, h0=T(d.h0), q=T(d.q),
        Gmu=T(Gmu), Nt, Nphi, derivatives, method)
    bank = NewtonianPotentialBank(caches, indices, key)
    isempty(ms) && return bank
    _fill_potential_bank!(bank, d, os, rg, ms, ls, softening, Int(Nt), Int(Nphi), method, T(Gmu))
    bank
end

function _fill_potential_bank!(bank::NewtonianPotentialBank{T}, d, os, rg, ms, ls,
                               softening, Nt, Nphi, method, Gmu) where {T}
    M = [2T(pi)*T(i)/Nt for i in 0:Nt-1]
    phi = [2T(pi)*T(j)/Nphi for j in 0:Nphi-1]
    A = zeros(method == :fft ? Complex{T} : T, Nt, Nphi, 3)
    # Precompute the separable direct weights once for the whole bank.
    um, ul = unique([m.m for m in ms]), unique(ls)
    mi, li = [findfirst(==(m.m), um) for m in ms], [findfirst(==(l), ul) for l in ls]
    W = method == :direct ? hcat([cos(m*p) for p in phi, m in um],
                                 [sin(m*p) for p in phi, m in um]) : zeros(T,0,0)
    C = method == :direct ? [cos(l*t) for t in M, l in ul] : zeros(T,0,0)
    S = method == :direct ? [sin(l*t) for t in M, l in ul] : zeros(T,0,0)
    projected = zeros(T, method == :direct ? Nt : 0, method == :direct ? 2length(um) : 0)
    power = (3-T(d.q))/2
    scale = T(2)/(T(Nt)*T(Nphi))
    for (ei,o) in enumerate(os)
        orbit = planet_orbit.(M, Ref(o))
        # Stable squared separation near conjunction, shared by all radii.
        bend = [2p.r_p*sin((az-p.phi_p)/2)^2 for p in orbit, az in phi]
        for (ri,r) in enumerate(rg)
            for i in 1:Nt
                rp = orbit[i].r_p
                ep = T(_bank_softlength(softening,T(d.h0),power,o.a,rp,r))
                dep = softening isa LocalScaleHeight ? power*ep/r : zero(T)
                d2ep = softening isa LocalScaleHeight ? power*(power-1)*ep/r^2 : zero(T)
                @inbounds for j in 1:Nphi
                    distance2 = (r-rp)^2 + 2r*bend[i,j] + ep^2
                    distance2 > 0 || throw(DomainError(distance2,"singular unsoftened sample"))
                    v = inv(sqrt(distance2))
                    u = (r-rp) + bend[i,j] + ep*dep
                    A[i,j,1] = -Gmu*v
                    A[i,j,2] = Gmu*u*v^3
                    A[i,j,3] = Gmu*((1+dep^2+ep*d2ep)*v^3-3u^2*v^5)
                end
            end
            if method == :fft
                fft_radix2!(A; dims=(1,2))
            end
            for c in 1:3
                method == :direct && mul!(projected, view(A,:,:,c), W)
                for k in eachindex(ms)
                    value = if method == :fft
                        real(A[mod(-ls[k],Nt)+1,ms[k].m+1,c])
                    else
                        acc = zero(T)
                        @inbounds for i in 1:Nt
                            acc += C[i,li[k]]*projected[i,mi[k]] +
                                   S[i,li[k]]*projected[i,mi[k]+length(um)]
                        end
                        acc
                    end
                    cache = bank.caches[ei,k]
                    dest = c == 1 ? cache.P : c == 2 ? cache.dP : cache.d2P
                    dest[ri] = scale*value
                end
            end
        end
    end
    bank
end
