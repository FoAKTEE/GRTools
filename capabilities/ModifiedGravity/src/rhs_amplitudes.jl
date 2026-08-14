# Summed effective-source amplitudes at the ring radius (RHSamplitudescalculation,
# Sourceincoordinates.nb, Polcar & Witzany 2025, arXiv:2507.15720). Completes SourceAmplitudes:
# sums the per-l `source_amplitudes` legs over l = max(|m|,2)..lmax, evaluated at the ring radius,
# into one (m,n,k) NamedTuple — this IS the generator of the shipped `.dat` source-mode tables.

"""
    rhs_amplitudes(m, n, k, a, p, e, x, ringrad, ringmass; lmax=10, Nmax=80) -> NamedTuple

Summed NP/GHP effective-source amplitudes of the `(m,n,k)` mode at the ring radius `ringrad`,
porting `RHSamplitudescalculation`. Sums the per-`l` [`source_amplitudes`](@ref) over
`l = max(|m|,2)..lmax` for the geodesic `(a,p,e,x)` (M=1), each leg evaluated at `ringrad`:

`Sl = Σ_l source_amplitudes(l,m,n,k,…).Sl(ringrad)`  (and likewise Sn, Sm, Sm†, Snorm).

Returns `(m, n, k, ω, radius, ringmass, Sources)` with `Sources = (Sl, Sn, Sm, Smdag, Snorm)` of
summed leg values and `ω` the +mode frequency. This is the exact pipeline that produced the shipped
`*.dat` files (`Snorm ≡ 0`, matching upstream `Sourcesamplitudes`).
"""
function rhs_amplitudes(m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real, ringrad::Real, ringmass::Real;
        lmax::Integer=10, Nmax::Integer=80)
    lmin = max(abs(m), 2)
    lmax >= lmin || throw(ArgumentError(
        "rhs_amplitudes: lmax must be at least max(abs(m), 2) = $lmin (got $lmax)"))
    Nmax >= 1 || throw(ArgumentError("rhs_amplitudes: Nmax must be positive (got $Nmax)"))
    Sl = Sn = Sm = Smdag = Snorm = zero(ComplexF64)
    ω = NaN
    for l in lmin:lmax
        S = source_amplitudes(l, m, n, k, a, p, e, x; Nmax=Nmax)
        ω = S.omega
        Sl    += S.Sl(ringrad)
        Sn    += S.Sn(ringrad)
        Sm    += S.Sm(ringrad)
        Smdag += S.Smdag(ringrad)
        Snorm += S.Snorm(ringrad)
    end
    return (m=Int(m), n=Int(n), k=Int(k), ω=ω,
            radius=Float64(ringrad), ringmass=Float64(ringmass),
            Sources=(Sl=Sl, Sn=Sn, Sm=Sm, Smdag=Smdag, Snorm=Snorm))
end

export rhs_amplitudes
