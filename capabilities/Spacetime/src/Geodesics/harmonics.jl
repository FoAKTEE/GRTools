"""
    EquatorialOrbit(a, p, e)

Bound, stable equatorial Kerr geodesic in units M=1, with positive Lz (x=1).
Negative `a` describes a retrograde orbit relative to the hole. The pericentre
is q=0. Stores E,L, Mino frequencies Υr,Υφ,Γ=⟨dt/dλ⟩, coordinate frequencies
Ωr=Υr/Γ, Ωφ=Υφ/Γ, and the signed radial trajectory. Inputs are promoted to a
common floating type; BigFloat inputs inherit the current precision.
Requires |a|<1, 0≤e<1 and p above the separatrix. At e=0, q is the limiting
radial phase and all nonzero radial harmonics of a phase-independent field vanish.
"""
struct EquatorialOrbit{T<:AbstractFloat,F}
    a::T
    p::T
    e::T
    E::T
    L::T
    Υr::T
    Υφ::T
    Γ::T
    Ωr::T
    Ωφ::T
    trajectory::F
end

function EquatorialOrbit(a::Real, p::Real, e::Real)
    a,p,e = promote(float(a),float(p),float(e))
    abs(a)<1 && 0<=e<1 && p>separatrix(a,e,one(a)) ||
        throw(ArgumentError("EquatorialOrbit requires |a|<1, 0≤e<1 and p above the separatrix"))
    E,L,_ = constants_of_motion(a,p,e,one(a))
    f = kerr_frequencies(a,p,e,one(a); Time="Mino")
    radial = geodesic_trajectory_radial(a,p,e,one(a))
    return EquatorialOrbit(a,p,e,E,L,f.Ur,f.Uph,f.Ut,f.Ur/f.Ut,f.Uph/f.Ut,radial)
end

"""
    worldline_rates(orbit::EquatorialOrbit, q)

Return the radial trajectory `(r,Ur,Δtr,Δφr)` and `dt_dlambda`,
`dphi_dlambda`, `dtau_dlambda`, plus proper-time velocities `ut,ur,uphi`.
The Kerr first integrals give P=E(r²+a²)-aL, Δ=r²-2r+a²,
dt/dλ=(r²+a²)P/Δ+a(L-aE), dφ/dλ=aP/Δ+L-aE, dτ/dλ=r².
`Ur=dr/dλ` is positive outbound and negative inbound; `ur=Ur/r²`.
"""
function worldline_rates(o::EquatorialOrbit, q::Real)
    g = o.trajectory(q)
    r = g.r
    P = o.E*(r^2+o.a^2)-o.a*o.L
    Δ = delta_kerr(o.a,r)
    dt_dlambda = (r^2+o.a^2)*P/Δ+o.a*(o.L-o.a*o.E)
    dphi_dlambda = o.a*P/Δ+o.L-o.a*o.E
    dtau_dlambda = r^2
    return (; g..., dt_dlambda, dphi_dlambda, dtau_dlambda,
            ut=dt_dlambda/r^2, ur=g.Ur/r^2, uphi=dphi_dlambda/r^2)
end

"""
    orbit_frequency(orbit, m::Integer, k::Integer)

ω_mk=mΩφ+kΩr. Since t=Γq/Υr+Δt and φ=Υφq/Υr+Δφ,
ωt-mφ = kq+ωΔt-mΔφ. In the Kepler limit Ωr→Ωφ this has l=m+k.
"""
orbit_frequency(o::EquatorialOrbit, m::Integer, k::Integer) = m*o.Ωφ+k*o.Ωr

"""
    orbit_harmonic(H, orbit, m::Integer, k::Integer; N=64, paired=false)

Coordinate-time Fourier projection of an azimuthal coefficient `H(q, state)`:
`∫₀²π (dt/dλ) H exp(i[kq+ωΔt-mΔφ]) dq / (2πΓ)`.
Here `state=worldline_rates(orbit,q)` and H already contains the azimuthal
projection with convention exp(-imψ). Time and azimuth offsets are removed at
pericentre. The zeroth projection of H=1 is one. N-point Gauss-Legendre nodes
and weights are computed in the orbit's floating type, including BigFloat.

With `paired=true`, integrate I(q)+I(2π-q) on [0,π] using N nodes. Both states
and both callback values are evaluated; no evenness or reality is assumed of H.
Increase N to check convergence, particularly for large |m|, |k| or e.
"""
function orbit_harmonic(H, o::EquatorialOrbit{T}, m::Integer, k::Integer;
                        N::Int=64, paired::Bool=false) where {T}
    nodes, weights = gl_nodes_weights(T,N)
    πT = T(pi)
    ω = orbit_frequency(o,m,k)
    origin = o.trajectory(zero(T))
    function integrand(q)
        s = worldline_rates(o,q)
        phase = k*q+ω*(s.Δtr-origin.Δtr)-m*(s.Δφr-origin.Δφr)
        return s.dt_dlambda*H(q,s)*cis(phase)
    end
    half = paired ? πT/2 : πT
    integral = sum(zip(nodes,weights)) do (z,w)
        q = half*(1+z)
        w*(paired ? integrand(q)+integrand(2πT-q) : integrand(q))
    end
    return half*integral/(2πT*o.Γ)
end
