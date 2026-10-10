# Geodesics — constants of motion E, Lz, Q for a generic bound Kerr geodesic (a,p,e,x),
# x = cos(inclination), M = 1. The `_*_generic` (e>0) and `_*_circular` (e=0) bodies are the
# Schmidt closed forms — a faithful transcription of the BHPToolkit/KerrGeodesics symbolic
# output (oracle-verified to ~1e-12 vs Mathematica). Arithmetic only, so BigFloat is inherited.
# The public bundle `constants_of_motion` also serves the KerrBackwards delegating-wrapper alias.
# Two endpoints of the (e = 0) circular forms are exact closed forms instead of the Schmidt
# expression, which is a 0/0 there: the polar x = 0 one (BHPToolkit, below) and the equatorial
# |x| = 1 one at p = 2 (Bardeen-Press-Teukolsky, `_circeq_removable`).

@inline _energy_generic(a, p, e, x) = sqrt((2*((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) + (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2)))) - 2*x*sqrt((((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*(((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))^2*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) - (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))^2*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + ((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*(-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))))/x^2))/(4*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + (((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))^2))

@inline _energy_circular(a, p, x) = sqrt(((-3 + p)*(-2 + p)^2*p^5 - 2*a^5*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) + a^4*p^2*(-1 + x^2)*(4 - 5*p*(-1 + x^2) + 3*p^2*(-1 + x^2)) - a^6*(-1 + x^2)^2*(x^2 + p^2*(-1 + x^2) - p*(1 + 2*x^2)) + a^2*p^3*(4 - 4*x^2 + p*(12 - 7*x^2) - 3*p^3*(-1 + x^2) + p^2*(-13 + 10*x^2)) + a*(-2*p^(9/2)*x*sqrt(p^2 + a^2*(-1 + x^2)) + 4*p^3*x*sqrt(p^3 + a^2*p*(-1 + x^2))) + 2*a^3*(2*p*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) - x^3*sqrt(p^7 + a^2*p^5*(-1 + x^2))))/((p^2 - a^2*(-1 + x^2))*((-3 + p)^2*p^4 - 2*a^2*p^2*(3 + 2*p - 3*x^2 + p^2*(-1 + x^2)) + a^4*(-1 + x^2)*(-1 + x^2 + p^2*(-1 + x^2) - 2*p*(1 + x^2)))))

@inline _angular_momentum_generic(a, p, e, x) = ((-2*a*p*sqrt((2*((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) + (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2)))) - 2*x*sqrt((((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*(((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))^2*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) - (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))^2*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + ((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*(-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))))/x^2))/(4*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + (((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))^2)))/(1 - e) + x*sqrt((-((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (((4*a^2*p^2)/(1 - e)^2 + ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))*(2*((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) + (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2)))) - 2*x*sqrt((((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*(((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))^2*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) - (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))^2*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + ((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*(-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))))/x^2)))/(4*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + (((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))^2))/x^2))/((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)

@inline _angular_momentum_circular(a, p, x) = (x^2*(-2*a*p*sqrt(((-3 + p)*(-2 + p)^2*p^5 - 2*a^5*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) + a^4*p^2*(-1 + x^2)*(4 - 5*p*(-1 + x^2) + 3*p^2*(-1 + x^2)) - a^6*(-1 + x^2)^2*(x^2 + p^2*(-1 + x^2) - p*(1 + 2*x^2)) + a^2*p^3*(4 - 4*x^2 + p*(12 - 7*x^2) - 3*p^3*(-1 + x^2) + p^2*(-13 + 10*x^2)) + a*(-2*p^(9/2)*x*sqrt(p^2 + a^2*(-1 + x^2)) + 4*p^3*x*sqrt(p^3 + a^2*p*(-1 + x^2))) + 2*a^3*(2*p*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) - x^3*sqrt(p^7 + a^2*p^5*(-1 + x^2))))/((p^2 - a^2*(-1 + x^2))*((-3 + p)^2*p^4 - 2*a^2*p^2*(3 + 2*p - 3*x^2 + p^2*(-1 + x^2)) + a^4*(-1 + x^2)*(-1 + x^2 + p^2*(-1 + x^2) - 2*p*(1 + x^2))))) + x*sqrt((-(((a^2 + (-2 + p)*p)*((-2 + p)*p - a^2*(-1 + x^2))*(p^2 - a^2*(-1 + x^2)))/x^2) + (((-3 + p)*(-2 + p)^2*p^5 - 2*a^5*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) + a^4*p^2*(-1 + x^2)*(4 - 5*p*(-1 + x^2) + 3*p^2*(-1 + x^2)) - a^6*(-1 + x^2)^2*(x^2 + p^2*(-1 + x^2) - p*(1 + 2*x^2)) + a^2*p^3*(4 - 4*x^2 + p*(12 - 7*x^2) - 3*p^3*(-1 + x^2) + p^2*(-13 + 10*x^2)) + a*(-2*p^(9/2)*x*sqrt(p^2 + a^2*(-1 + x^2)) + 4*p^3*x*sqrt(p^3 + a^2*p*(-1 + x^2))) + 2*a^3*(2*p*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) - x^3*sqrt(p^7 + a^2*p^5*(-1 + x^2))))*(4*a^2*p^2 + (((-2 + p)*p - a^2*(-1 + x^2))*(p^4 + a^2*(p*(2 + p) - (a^2 + (-2 + p)*p)*(-1 + x^2))))/x^2))/((p^2 - a^2*(-1 + x^2))*((-3 + p)^2*p^4 - 2*a^2*p^2*(3 + 2*p - 3*x^2 + p^2*(-1 + x^2)) + a^4*(-1 + x^2)*(-1 + x^2 + p^2*(-1 + x^2) - 2*p*(1 + x^2)))))/x^2)))/((-2 + p)*p - a^2*(-1 + x^2))

@inline _carter_generic(a, p, e, x) = (1 - x^2)*(a^2*(1 - (2*((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) + (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2)))) - 2*x*sqrt((((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*(((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))^2*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) - (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))^2*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + ((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*(-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))))/x^2))/(4*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + (((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))^2)) + ((-2*a*p*sqrt((2*((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) + (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2)))) - 2*x*sqrt((((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*(((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))^2*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) - (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))^2*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + ((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*(-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))))/x^2))/(4*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + (((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))^2)))/(1 - e) + x*sqrt((-((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (((4*a^2*p^2)/(1 - e)^2 + ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))*(2*((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) + (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2)))) - 2*x*sqrt((((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*(((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))^2*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e)) - (-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))^2*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + ((2*a*p*(a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2)))/(1 + e) - (2*a*p*(a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2)))/(1 - e))*(-((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(p^2/(1 + e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(p^2/(1 - e)^2 + a^2*(1 - x^2))*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))*(((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))))/x^2)))/(4*((-2*a*p*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2))/(1 + e) + (2*a*p*((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2))/(1 - e))*((2*a*p*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))))/(1 + e) - (2*a*p*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))/(1 - e)) + (((p*(-2 + p/(1 + e)))/(1 + e) + ((a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))/x^2)*(p^4/(1 - e)^4 + a^2*((p*(2 + p/(1 - e)))/(1 - e) + (a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))) - ((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)*(p^4/(1 + e)^4 + a^2*((p*(2 + p/(1 + e)))/(1 + e) + (a^2 - (2*p)/(1 + e) + p^2/(1 + e)^2)*(1 - x^2))))^2))/x^2))^2/(x^2*((p*(-2 + p/(1 - e)))/(1 - e) + ((a^2 - (2*p)/(1 - e) + p^2/(1 - e)^2)*(1 - x^2))/x^2)^2))

@inline _carter_circular(a, p, x) = (1 - x^2)*(a^2*(1 - ((-3 + p)*(-2 + p)^2*p^5 - 2*a^5*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) + a^4*p^2*(-1 + x^2)*(4 - 5*p*(-1 + x^2) + 3*p^2*(-1 + x^2)) - a^6*(-1 + x^2)^2*(x^2 + p^2*(-1 + x^2) - p*(1 + 2*x^2)) + a^2*p^3*(4 - 4*x^2 + p*(12 - 7*x^2) - 3*p^3*(-1 + x^2) + p^2*(-13 + 10*x^2)) + a*(-2*p^(9/2)*x*sqrt(p^2 + a^2*(-1 + x^2)) + 4*p^3*x*sqrt(p^3 + a^2*p*(-1 + x^2))) + 2*a^3*(2*p*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) - x^3*sqrt(p^7 + a^2*p^5*(-1 + x^2))))/((p^2 - a^2*(-1 + x^2))*((-3 + p)^2*p^4 - 2*a^2*p^2*(3 + 2*p - 3*x^2 + p^2*(-1 + x^2)) + a^4*(-1 + x^2)*(-1 + x^2 + p^2*(-1 + x^2) - 2*p*(1 + x^2))))) + (x^2*(-2*a*p*sqrt(((-3 + p)*(-2 + p)^2*p^5 - 2*a^5*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) + a^4*p^2*(-1 + x^2)*(4 - 5*p*(-1 + x^2) + 3*p^2*(-1 + x^2)) - a^6*(-1 + x^2)^2*(x^2 + p^2*(-1 + x^2) - p*(1 + 2*x^2)) + a^2*p^3*(4 - 4*x^2 + p*(12 - 7*x^2) - 3*p^3*(-1 + x^2) + p^2*(-13 + 10*x^2)) + a*(-2*p^(9/2)*x*sqrt(p^2 + a^2*(-1 + x^2)) + 4*p^3*x*sqrt(p^3 + a^2*p*(-1 + x^2))) + 2*a^3*(2*p*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) - x^3*sqrt(p^7 + a^2*p^5*(-1 + x^2))))/((p^2 - a^2*(-1 + x^2))*((-3 + p)^2*p^4 - 2*a^2*p^2*(3 + 2*p - 3*x^2 + p^2*(-1 + x^2)) + a^4*(-1 + x^2)*(-1 + x^2 + p^2*(-1 + x^2) - 2*p*(1 + x^2))))) + x*sqrt((-(((a^2 + (-2 + p)*p)*((-2 + p)*p - a^2*(-1 + x^2))*(p^2 - a^2*(-1 + x^2)))/x^2) + (((-3 + p)*(-2 + p)^2*p^5 - 2*a^5*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) + a^4*p^2*(-1 + x^2)*(4 - 5*p*(-1 + x^2) + 3*p^2*(-1 + x^2)) - a^6*(-1 + x^2)^2*(x^2 + p^2*(-1 + x^2) - p*(1 + 2*x^2)) + a^2*p^3*(4 - 4*x^2 + p*(12 - 7*x^2) - 3*p^3*(-1 + x^2) + p^2*(-13 + 10*x^2)) + a*(-2*p^(9/2)*x*sqrt(p^2 + a^2*(-1 + x^2)) + 4*p^3*x*sqrt(p^3 + a^2*p*(-1 + x^2))) + 2*a^3*(2*p*x*(-1 + x^2)*sqrt(p^3 + a^2*p*(-1 + x^2)) - x^3*sqrt(p^7 + a^2*p^5*(-1 + x^2))))*(4*a^2*p^2 + (((-2 + p)*p - a^2*(-1 + x^2))*(p^4 + a^2*(p*(2 + p) - (a^2 + (-2 + p)*p)*(-1 + x^2))))/x^2))/((p^2 - a^2*(-1 + x^2))*((-3 + p)^2*p^4 - 2*a^2*p^2*(3 + 2*p - 3*x^2 + p^2*(-1 + x^2)) + a^4*(-1 + x^2)*(-1 + x^2 + p^2*(-1 + x^2) - 2*p*(1 + x^2)))))/x^2))^2)/((-2 + p)*p - a^2*(-1 + x^2))^2)


# Exact-polar (x=0) formulas from the authoritative local BHPToolkit source
# ref-code/bhptoolkit-KerrGeodesics-0.9.0/Kernel/ConstantsOfMotion.m.  The generic Schmidt
# expressions contain explicit 1/x^2 factors and therefore cannot be evaluated at the polar
# endpoint even though their physical limit is finite.
@inline _energy_polar_circular(a, p) =
    sqrt(p * (a^2 - 2p + p^2)^2 /
         ((a^2 + p^2) * (a^2 + a^2*p - 3p^2 + p^3)))

@inline _carter_polar_circular(a, p) =
    p^2 * (a^4 + 2a^2 * (p - 2) * p + p^4) /
    ((a^2 + p^2) * ((p - 3) * p^2 + a^2 * (1 + p)))

@inline _polar_eccentric_denominator(a, p, e) =
    a^4 * (e^2 - 1)^2 * (e^2 - p - 1) + (3 + e^2 - p) * p^4 -
    2a^2 * p^2 * (-1 - e^4 + p + e^2 * (2 + p))

@inline _energy_polar_generic(a, p, e) =
    sqrt(-p * (a^4 * (e^2 - 1)^2 + (-4e^2 + (p - 2)^2) * p^2 +
               2a^2 * p * (p - 2 + e^2 * (2 + p))) /
         _polar_eccentric_denominator(a, p, e))

@inline _carter_polar_generic(a, p, e) =
    -p^2 * (a^4 * (e^2 - 1)^2 + p^4 +
             2a^2 * p * (p - 2 + e^2 * (2 + p))) /
    _polar_eccentric_denominator(a, p, e)

@inline _constants_zero(a, p, e, x) = zero(float(promote_type(typeof(a), typeof(p), typeof(e), typeof(x))))

# --- the equatorial (|x| = 1) removable singularity of the circular closed forms ---------------
#
# `_angular_momentum_circular` carries the factor 1/((p-2)p - a^2(x^2-1)) and `_carter_circular`
# its square. On the equator x^2 = 1 that denominator is exactly p(p-2); at p = 2 the numerator
# vanishes with it, so the ratio is a REMOVABLE 0/0 that Float64 resolves as NaN - and Q, which is
# (1-x^2)*(...)^2, inherits it as 0*NaN. `_energy_circular` carries no such factor and is regular
# there (it agrees with the value below to 1e-77 at 256 bits), so only Lz and Q are repaired.
#
# The removed value is the Bardeen-Press-Teukolsky (1972, ApJ 178 347, Eqs. 2.12-2.13) equatorial
# circular pair, which with u = sqrt(p) and this package's SIGNED spin (a > 0 prograde, a < 0
# retrograde, x = +1 always; arXiv:2606.04998 main.tex:196-198) is
#
#     E = (u^3 - 2u + a)/(u^(3/2) sqrt(u^3 - 3u + 2a)),
#     Lz = (p^2 - 2au + a^2)/(u^(3/2) sqrt(u^3 - 3u + 2a)),      Q = 0,
#
# the x = -1 branch being its mirror E(a,p,-1) = E(-a,p,+1), Lz(a,p,-1) = -Lz(-a,p,+1). That this
# IS the limit and not merely a plausible substitute is measured, not asserted: at 256 bits the
# untouched Schmidt form approaches it LINEARLY from both sides - |dLz|/Lz = 2.6e-21 (a = 0.9) and
# 2.3e-21 (a = 0.998) at |p - 2| = 1e-20, and 1e-31 at 1e-30.
#
# The trigger is the EXACT vanishing of that denominator at |x| = 1, never a tolerance, so it is
# precision-independent and every other radius keeps its pre-fix Schmidt evaluation BITWISE (the
# two forms otherwise differ by up to 10 ulp, which is why the whole branch is not simply
# replaced). Two consequences, both gated in test_geodesics_constants.jl:
#
#  1. r = 2M carries a circular geodesic only for a > 1/sqrt(2) = 0.70710678..., since the shared
#     radicand is u^3 - 3u + 2a = 2a - sqrt(2) at p = 2 - r = 2M IS the prograde photon orbit at
#     a = 1/sqrt(2). Below that spin the pair raises the ordinary `DomainError` of a negative
#     radicand, which is what every neighbouring radius already does; pre-fix, a = 0 silently
#     returned E = 0.0 with Lz = Q = NaN at what is its own horizon.
#  2. This repairs the singularity, NOT the conditioning of the generic form around it, and the
#     distinction matters to any caller that computes p arithmetically rather than writing the
#     literal 2.0. The Schmidt Lz loses about log10|p-2|^-1 digits in a NEIGHBOURHOOD of p = 2,
#     and the loss runs all the way to TOTAL at the immediate float neighbours. Relative error of
#     Lz against a 256-bit evaluation, a = 0.998:
#
#         |p-2| :  1e-2      1e-6      1e-10     1 ulp (nextfloat)   1 ulp (prevfloat)
#         rel   :  1.1e-14   6.2e-11   1.5e-6    1.15e-1 (11 %)      7.7e-1 (77 %)
#
#     i.e. `geodesic_angular_momentum(0.998, nextfloat(2.0), 0, 1)` = 1.4999999999999998 and
#     `prevfloat(2.0)` = 3.0000000000000004 against a true 1.6941511059695233 - SILENTLY FINITE
#     and badly wrong, identically so before and after this repair (they are the generic form's
#     own cancellation, and 2.0 +/- 1 ulp is not the exact root, so no exact-trigger branch can
#     reach them). Q is affected the same way; E is not (its closed form has no such factor -
#     4.3e-16 at both neighbours). Recomputing the whole equatorial branch through the pair above
#     would fix it, at the cost of moving every equatorial value by up to 10 ulp.
#
#  3. The INCLINED branch (|x| < 1) carries the same denominator and the same removable root, at
#     p* = 1 + sqrt(1 + a^2(x^2-1)) < 2, and is left untouched here - but note WHY the trigger
#     below cannot serve it. p* is irrational, so in Float64 the denominator at the computed p* is
#     +2.1e-17 (a = 0.998, x = 0.99), +1.1e-16 (0.998, 0.9), +1.6e-16 (0.9, 0.99), -5.8e-17
#     (0.998, 0.999) - near zero but never exactly zero, so an exact-zero predicate can never fire
#     there. What comes out is the cancellation ratio: Lz = 0, 3.24, 2.788, -0.0 against 256-bit
#     values 1.687, 1.637, 2.174, 1.693 (28-100 % wrong) with Q wrong by 62-262 %, silently finite.
#     This is pre-existing (identical at 5243f42) and out of the p = 2 node's scope. A node that
#     takes it on cannot copy the mechanism below: it must either cancel the common factor
#     algebraically in Schmidt's general formula, or detect the near-zero denominator with a
#     SCALED tolerance and evaluate the analytic limit there.
@inline _circular_lz_denominator(a, p, x) = (p - 2) * p - a^2 * (x^2 - 1)

# Only the p = 2 root of the equatorial denominator p(p-2) is REMOVABLE: at p = 0 - the ring
# singularity, which is no orbit's radius - the numerator does not vanish with it, so that root
# is excluded and keeps the Schmidt form's NaN rather than gaining a spurious finite/Inf value.
@inline _circeq_removable(a, p, x) =
    isone(abs(x)) && !iszero(p) && iszero(_circular_lz_denominator(a, p, x))

# Bardeen-Press-Teukolsky equatorial circular (E, Lz) at SIGNED spin `a` (M = 1); u = sqrt(p).
@inline function _bpt_circeq(a, p)
    u = sqrt(p)
    den = u * sqrt(u) * sqrt(u * p - 3 * u + 2 * a)   # r^(3/4) sqrt(r^(3/2) - 3 r^(1/2) + 2a)
    ((u * p - 2 * u + a) / den, (p^2 - 2 * a * u + a^2) / den)
end

# Circular Lz and Q with the equatorial 0/0 removed; the Schmidt forms verbatim everywhere else.
@inline _angular_momentum_circular_sel(a, p, x) =
    _circeq_removable(a, p, x) ?
        (x > 0 ? last(_bpt_circeq(a, p)) : -last(_bpt_circeq(-a, p))) :
        _angular_momentum_circular(a, p, x)

@inline _carter_circular_sel(a, p, e, x) =
    _circeq_removable(a, p, x) ? _constants_zero(a, p, e, x) : _carter_circular(a, p, x)

"""
    geodesic_energy(a, p, e, x) -> E

Orbital energy `E` of a generic bound Kerr geodesic (M = 1), `x = cos(inclination)`. Selects the
circular (e = 0) or generic (e > 0) closed form.
"""
@inline function geodesic_energy(a, p, e, x)
    if _reflect_eccentric_equator(a,e,x)
        a,x = -a,-x
    end
    return iszero(x) ?
    (iszero(e) ? _energy_polar_circular(a, p) : _energy_polar_generic(a, p, e)) :
    (iszero(e) ? _energy_circular(a, p, x) : _energy_generic(a, p, e, x))
end

"""
    geodesic_angular_momentum(a, p, e, x[, E]) -> Lz

Axial angular momentum `Lz` of a generic bound Kerr geodesic (M = 1). A passed `E` is ignored
(present for the legacy delegating-wrapper signature).

Equatorial circular endpoint: at `e = 0`, `|x| = 1`, `p = 2` the Schmidt closed form is a
removable `0/0` (its `1/((p-2)p - a^2(x^2-1))` factor), and `Lz` is instead the exact
Bardeen-Press-Teukolsky value `(p^2 - 2a sqrt(p) + a^2)/(p^(3/4) sqrt(p^(3/2) - 3 sqrt(p) + 2a))`
- the two-sided limit, and the paper's own `r0 = 2M` configuration. Every other `p` is the
Schmidt evaluation bitwise. `r = 2M` is a circular geodesic only for `a > 1/sqrt(2)`; below that
spin the radicand is negative and the call raises `DomainError`, as it already does at every
neighbouring radius.

!!! warning "Only `p == 2` exactly is repaired"
    The generic form is ill-conditioned in a NEIGHBOURHOOD of `p = 2`, and this endpoint does not
    change that. At `a = 0.998` the relative error of `Lz` against a 256-bit evaluation is 1.5e-6
    at `|p-2| = 1e-10` and reaches 11 % at `nextfloat(2.0)` (1.4999999999999998) and 77 % at
    `prevfloat(2.0)` (3.0000000000000004), against a true 1.6941511059695233 - silently finite and
    wrong, both before and after the repair. Pass the literal `2.0`; a radius arrived at by
    arithmetic (an ISCO scaling, a mesh coordinate, a unit conversion) may land a few ulp away and
    get one of those values with no warning. `E` is unaffected there (4.3e-16).
"""
@inline function geodesic_angular_momentum(a, p, e, x)
    reflected = _reflect_eccentric_equator(a,e,x)
    if reflected
        a,x = -a,-x
    end
    Lz = iszero(x) ? _constants_zero(a, p, e, x) :
    (iszero(e) ? _angular_momentum_circular_sel(a, p, x) : _angular_momentum_generic(a, p, e, x))
    return reflected ? -Lz : Lz
end
@inline geodesic_angular_momentum(a, p, e, x, _E) = geodesic_angular_momentum(a, p, e, x)

"""
    carter_constant(a, p, e, x[, E[, Lz]]) -> Q

Carter constant `Q` of a generic bound Kerr geodesic (M = 1). Passed `E`/`Lz` are ignored
(present for the legacy delegating-wrapper signature). `Q` is exactly zero on the equator, the
`p = 2` endpoint included (there the Schmidt form is `0*NaN`; see `geodesic_angular_momentum`).
"""
@inline function carter_constant(a, p, e, x)
    if _reflect_eccentric_equator(a,e,x)
        a,x = -a,-x
    end
    return iszero(x) ?
    (iszero(e) ? _carter_polar_circular(a, p) : _carter_polar_generic(a, p, e)) :
    (iszero(e) ? _carter_circular_sel(a, p, e, x) : _carter_generic(a, p, e, x))
end
@inline carter_constant(a, p, e, x, _E) = carter_constant(a, p, e, x)
@inline carter_constant(a, p, e, x, _E, _Lz) = carter_constant(a, p, e, x)

"""
    constants_of_motion(a, p, e, x) -> (E, Lz, Q)

The `(E, Lz, Q)` bundle for a generic bound Kerr geodesic (M = 1). Canonical replacement for the
KerrGeodesics `KerrGeoConstantsOfMotion` and the KerrBackwards `constants_of_motion` wrapper.
For negative spin on an eccentric equatorial orbit, evaluate in the reflected
chart `(-a,p,e,-x)` and reverse Lz. The scalar accessors use the same rule.
The circular formulas and all a≥0 evaluations retain their operation order.
"""
@inline constants_of_motion(a, p, e, x) =
    (E = geodesic_energy(a, p, e, x), Lz = geodesic_angular_momentum(a, p, e, x), Q = carter_constant(a, p, e, x))

# The Schmidt eccentric square-root branch assumes nonnegative a. Kerr reflection
# preserves E,Q and reverses Lz. Keep circular and inclined arithmetic unchanged.
# Inclined a<0, e>0 is still unnormalized (R(r_p)≠0): pass a≥0 with signed x there.
@inline _reflect_eccentric_equator(a,e,x) = a < 0 && e > 0 && abs(x) == 1
