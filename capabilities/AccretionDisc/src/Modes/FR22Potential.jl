using NumericsCore: gl_nodes_weights

"""
    fr22_potential_quadrature(d, orbit, r, mode; softening,
                             panels=32, order=32, Nphi=2048, Gmu=1)

Per-mode factorised orbit quadrature returning `(Phi, dPhi_dr, d2Phi_dr2)`.
Production tabulation uses `potential_bank` and `harmonic_cache`; this independent
oracle uses composite Gauss-Legendre in eccentric anomaly and the existing
softened Laplace coefficients. Refine both `panels` and `Nphi` to budget error.

FR22 (22)–(25), source `ref-paper/arXiv-2207.14637/src/Sections/
3_numerical_procedure.tex:26–48`; factorisation certificate:
`derivation_fr/src/phi_factorisation.jl`. The measure is
`dM=(1-e*cos(E))*dE`; the corrected physical-time inverse has `n/(2pi^2)`.
Thus coefficients at fixed `(m,l=m+k)` are independent of mean motion.
The direct potential has `Gmu=G*M_p=1` by default; no wake mass scaling is applied.

`softening` is required: FR22's `epsilon=0.3H` (same source :79) does not
identify the evaluation radius of H (defined in `2_framework.tex:48`). All four
shared policies are supported with the bank's reference-radius units. Local H(r)
includes its first and second radial chain terms; the other policies hold the
softening fixed when differentiating at each orbit point. No policy is inferred
from a fit. Floating types/precision are inherited from the inputs.
"""
function fr22_potential_quadrature(d::NewtonianDisc, orbit::EccentricOrbit,
                                  r::Real, mode::HarmonicMode;
                                  softening::SofteningPolicy,
                                  panels::Integer=32, order::Integer=32,
                                  Nphi::Integer=2048, Gmu::Real=1)
    T = float(promote_type(typeof(d.h0), typeof(d.q), typeof(orbit.a), typeof(r),
                          typeof(_bank_softvalue(softening)), typeof(Gmu)))
    rr, a, e = T(r), T(orbit.a), T(orbit.e)
    l = Base.Checked.checked_add(mode.m, mode.k)
    isfinite(rr) && rr > 0 || throw(DomainError(r, "radius must be finite and positive"))
    panels > 0 && order > 0 && Nphi > 2mode.m ||
        throw(ArgumentError("positive orbit quadrature and sub-Nyquist azimuth required"))
    d.indirect == 0 || throw(ArgumentError("FR22 uses the direct potential only"))
    isfinite(Gmu) && Gmu >= 0 || throw(DomainError(Gmu, "Gmu must be finite and nonnegative"))
    isfinite(d.h0) && d.h0 >= 0 && isfinite(d.q) ||
        throw(ArgumentError("finite q and finite nonnegative h0 are required"))
    isapprox(mode.omega, l*orbit.n; rtol=8eps(T), atol=zero(T)) ||
        throw(ArgumentError("mode frequency must equal (m+k)*orbit.n"))
    nodes, weights = gl_nodes_weights(T, Int(order))
    power, half = (3-T(d.q))/2, one(T)/2
    width = 2T(pi)/panels
    value, deriv, deriv2 = zero(T), zero(T), zero(T)
    for panel in 0:panels-1, j in eachindex(nodes)
        E = (panel+(nodes[j]+1)/2)*width
        sn, cs = sincos(E)
        rho = 1-e*cs
        rp = a*rho
        phi = atan(sqrt((1-e)*(1+e))*sn, cs-e)
        M = E-e*sn
        ep = T(_bank_softlength(softening, T(d.h0), power, a, rp, rr))
        ep > 0 || throw(DomainError(ep, "the softened oracle requires positive softening"))
        alpha, eta = rr/rp, ep/rp
        b = laplace_b_soft(half, mode.m, alpha, eta; nphi=Nphi)
        db = laplace_b_soft_deriv(half, mode.m, alpha, eta; nphi=Nphi)
        d2b = laplace_b_soft_2deriv(half, mode.m, alpha, eta; nphi=Nphi)
        if softening isa LocalScaleHeight
            dep = power*ep/rr
            d2ep = power*(power-1)*ep/rr^2
            b3 = laplace_b_soft(half+1, mode.m, alpha, eta; nphi=Nphi)
            db3 = laplace_b_soft_deriv(half+1, mode.m, alpha, eta; nphi=Nphi)
            b5 = laplace_b_soft(half+2, mode.m, alpha, eta; nphi=Nphi)
            # b_eta=-eta*b_3/2; b_alpha,eta=-eta*b'_3/2;
            # b_eta,eta=-b_3/2+3eta^2*b_5/2.
            db -= dep*eta*b3
            d2b += -2dep*eta*db3 - (dep^2+rp*d2ep*eta)*b3 + 3dep^2*eta^2*b5
        end
        w = -T(Gmu)*weights[j]*rho*cos(mode.m*phi-l*M)/(2panels)
        value += w*b/rp
        deriv += w*db/rp^2
        deriv2 += w*d2b/rp^3
    end
    (value, deriv, deriv2)
end
