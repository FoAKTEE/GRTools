"""
Leading G26 singular-field gauge in Boyer–Lindquist (t,r,phi,z) at the worldline event.
The caller supplies the metric and the signed proper-time rates owned by
Spacetime.worldline_rates; no orbit, turning-point or metric solver lives here.
[I] Equal-time slicing, coordinate-component transport, the softening frame,
and far-zone continuation are model choices. No retarded or subleading claim.
"""
module EccentricSource
using LinearAlgebra

"Inverse Jacobian dx/dy at (rp,phip); time is retained in the contraction."
function chart_jacobian(rp, phip)
    s, c = sincos(phip)
    o, z = one(rp), zero(rp)
    [o z z z; z c s z; z -s/rp c/rp z; z z z o]
end

"Full spacetime projector in the Cartesian chart, before any time restriction."
function chart_projector(g, u, rp, phip)
    ulow = g*u
    J = chart_jacobian(rp, phip)
    transpose(J)*(g+ulow*transpose(ulow))*J
end

"Equal-time distance and its first two field-radius derivatives, holding the event fixed."
function distance_jet(g, u, rp, r, psi, softening)
    ulow = g*u
    s, c = sincos(psi)
    dR, dPhi = r*c-rp, r*s/rp
    A = g[2,2]+ulow[2]^2
    B = ulow[2]*ulow[3]
    C = g[3,3]+ulow[3]^2
    s2 = softening^2+A*dR^2+2B*dR*dPhi+C*dPhi^2
    s2_r = 2A*dR*c+2B*(c*dPhi+dR*s/rp)+2C*dPhi*s/rp
    s2_rr = 2A*c^2+4B*c*s/rp+2C*(s/rp)^2
    (; s2, s2_r, s2_rr)
end

"Covariant metric field and analytic radial jets; no reality projection or component zeroing."
function field_jet(g, u, rp, r, psi, softening, mass)
    j = distance_jet(g,u,rp,r,psi,softening)
    ulow = g*u
    K = 2mass*(g+2ulow*transpose(ulow))
    s = sqrt(j.s2)
    h = K ./ s
    h_r = -K*j.s2_r ./ (2s^3)
    h_rr = K*(3j.s2_r^2/(4s^5)-j.s2_rr/(2s^3))
    (; h, h_r, h_rr, j...)
end

"""
For |dr|,|dphi| <= delta <= min(1,rp/2), bound |s_cart^2-s_twist^2|
by bound*delta^3. Taylor inequalities used: |cos y-1|<=y²/2,
|sin y-y|<=|y|³/6, |sin²(z/2)-z²/4|<=z⁴/48.
Use dimensionless BL coordinates with lengths in M=1 units in this bound.
The absolute values make the same proof apply to both signs of u^r.
"""
function chart_remainder_bound(g, u, rp, twist)
    ulow = g*u
    Q = g[2:3,2:3]+ulow[2:3]*transpose(ulow[2:3])
    K = (rp+1)/2+1/rp+one(rp)/6
    maximum(abs,Q)*(4K+K^2)+4Q[2,2]*(1+abs(twist))^4/48
end
end
