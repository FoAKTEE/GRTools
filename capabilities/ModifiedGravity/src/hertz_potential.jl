# Radial Hertz potential per (l,m,n,k) mode (RadHzpotmode, Sourceincoordinates.nb,
# Polcar & Witzany 2025): the −2 point-particle ψ4 amplitude times the homogeneous s=+2 "Up"
# solution. Element type is threaded (no float()/ComplexF64 pins) so BigFloat flows through.

using Perturbation: TeukolskyRadial, teukolsky_radial_deriv,
    TeukolskyPointParticleModeCircular, TeukolskyPointParticleModeSpherical

"""
    hertz_potential_mode(l, m, n, k, a, p, e, x; Nmax=80) -> (RadHZ, dRadHZ, ω)

Radial Hertz potential of the `(l,m,n,k)` mode on a Kerr geodesic `(a,p,e,x)` (M=1), porting
`RadHzpotmode` (Sourceincoordinates.nb). `RadHZ`/`dRadHZ` are callables of Boyer-Lindquist radius
`r`; `ω` is the mode frequency.

`RadHZ(r) = (-1)^(l+kloc+m) (2/ω⁴) Z^ℐ_{-2,lmnk} · R^up_{+2,lm}(r)`, with `kloc = 0` for an
equatorial orbit (`x=1`) and `kloc = k` otherwise — a load-bearing convention pinned against the
upstream `Sourcesamplitudes` minus-mode rule. `Z^ℐ` is the ψ4 amplitude at infinity
(`TeukolskyPointParticleMode*` `Zinf`); `R^up_{+2}` is `TeukolskyRadial(2,l,m,a,ω).Up`.

Only circular-equatorial (`e=0,x=1`) and spherical (`e=0,|x|<1`) orbits are ported; eccentric
orbits throw (no migrated grounding yet).
"""
function hertz_potential_mode(l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real; Nmax::Integer=80)
    kloc = x == 1 ? 0 : k
    mode = _orbit_psi4_mode(l, m, n, k, a, p, e, x; Nmax=Nmax)
    ω = mode.omega
    iszero(ω) && throw(DomainError(ω,
        "hertz_potential_mode: static modes require a separate completion/static reconstruction"))
    hzamp = (-1)^(l + kloc + m) * (2 / ω^4) * mode.Zinf
    up = TeukolskyRadial(2, l, m, a, ω; Nmax=Nmax).Up
    RadHZ(r::Real)  = hzamp * up(r)
    dRadHZ(r::Real) = hzamp * teukolsky_radial_deriv(up, r)
    return (RadHZ, dRadHZ, ω)
end

"""
    _orbit_psi4_mode(l, m, n, k, a, p, e, x; Nmax) -> NamedTuple

Dispatch on orbit type to the matching s=-2 point-particle ψ4 mode, returning its `Zinf`/`omega`.
Circular (`e=0,x=1`) and spherical (`e=0,|x|<1`) only.
"""
function _orbit_psi4_mode(l, m, n, k, a, p, e, x; Nmax)
    if e == 0 && x == 1
        return TeukolskyPointParticleModeCircular(-2, l, m, a, p; Nmax=Nmax)
    elseif e == 0 && abs(x) < 1
        return TeukolskyPointParticleModeSpherical(l, m, k, a, p, x; s=-2, Nmax=Nmax)
    else
        throw(ArgumentError("hertz_potential_mode: only e=0 (circular/spherical) orbits are ported (got e=$e, x=$x)."))
    end
end

export hertz_potential_mode
