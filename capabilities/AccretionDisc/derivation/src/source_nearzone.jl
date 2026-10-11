"""
Independent scalar-distance oracle, transcribed from EffectiveSource text only:
kerr-equatorial.c:98–118 and kerr-equatorial-coeffs.c:207–212.
This is not a full tensor-field certificate. No reference C code is executed.
"""
module SourceNearZone
using ..EccentricSource

"Sheared-phi coefficients in BL coordinates; ur is the signed contravariant rate."
function coefficients(a, M, R, L, ur)
    Delta = a^2+R*(R-2M)
    G = a^2*(2M+R)+R*(L^2+R^2)
    alpha20 = R^2*(a^4*(2M+R)+a^2*R*(L^2-4M^2+2M*R*ur^2+R^2*(2+ur^2))+
        R^2*(L^2*(R-2M)+R^2*(-2M+R+R*ur^2)))/(Delta^2*G)
    alpha02 = R^2
    beta = 4*(L^2+R^2+a^2*(R+2M)/R)
    c = -L*R^3*ur/(Delta*G)
    (; alpha20, alpha02, beta, c)
end

"""
Six dyadic samples with |Delta x| = max(|dr|,|dphi|), in M=1 units.
Require a bounded cubic quotient, a vanishing quadratic quotient and a
log-log slope 3 +/- 0.1. This is a numerical scaling gate, not an exact proof.
"""
function remainder_scaling(g, u, R, coeff, radial, azimuth)
    steps = [one(R)/2^k for k in 7:12]
    remainder = map(steps) do step
        dr, dp = radial*step, azimuth*step
        cart = EccentricSource.distance_jet(g,u,R,R+dr,dp,zero(R)).s2
        ref = coeff.alpha20*dr^2+coeff.beta*sin((dp-coeff.c*dr)/2)^2
        abs(cart-ref)
    end
    cubic, quadratic = remainder ./ steps.^3, remainder ./ steps.^2
    lx, ly = log.(steps), log.(remainder)
    x = lx .- sum(lx)/length(lx)
    slope = sum(x .* (ly .- sum(ly)/length(ly)))/sum(abs2,x)
    bounded = all(isfinite,cubic) && minimum(cubic)>0 && maximum(cubic)/minimum(cubic)<=5//4
    vanishing = all(diff(quadratic) .< 0) && last(quadratic)<=first(quadratic)/16
    pass = bounded && vanishing && abs(slope-3)<=1//10
    (; steps, remainder, cubic, quadratic, slope, bounded, vanishing, pass)
end
end
