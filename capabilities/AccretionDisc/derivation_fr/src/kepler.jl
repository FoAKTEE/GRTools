"""
FR22 (24)–(25), using M=E-e sin(E) (MD99 2.52, invoked in FR22).
Assume a,n>0, 0≤e<1, pericentre M=E=phi=0 and continuous unwrapped
angles. The tangent-half-angle chart omits E=pi modulo 2pi; continuity
extends its rational identities there because 1-e*cos(E)≥1-e>0.
"""
function certify_fr22_kepler()
    @variables t u a n E
    # Rational chart of sqrt(1-e²), on the positive root: 0≤u<1.
    e, beta = 2u/(1+u^2), (1-u^2)/(1+u^2)
    C, S = (1-t^2)/(1+t^2), 2t/(1+t^2)
    cp = (C-e)/(1-e*C)
    sp = beta*S/(1-e*C)
    rp = a*(1-e^2)/(1+e*cp)
    M = E-e*sin(E)
    c = Certificate[]
    push!(c, certify_symbolic("FR22 (24) unit circle", "24", cp^2+sp^2-1))
    push!(c, certify_symbolic("FR22 (25) radius in eccentric anomaly", "25",
                              rp-a*(1-e*C)))
    push!(c, certify_symbolic("Kepler time measure dt/dE", "MD99 2.52",
                              Symbolics.derivative(M,E)/n-(1-e*cos(E))/n))
    Mchart = E-e*S
    push!(c, certify_symbolic("Kepler map is odd", "MD99 2.52",
        Symbolics.substitute(Mchart,Dict(E=>-E,t=>-t))+Mchart))
    # Exact decomposition plus cos(E)≤1 proves dM/dE≥1-e>0. Together
    # with M(E)→±infinity and oddness, uniqueness proves E(-M)=-E(M).
    push!(c, certify_symbolic("Kepler inverse odd: monotonicity decomposition", "MD99 2.52",
        Symbolics.derivative(M,E)-(1-e)-e*(1-cos(E))))
    println("ASSUMPTIONS Kepler inverse: 0<=e<1, cos(E)<=1; dM/dE>=1-e>0; ",
            "bounded e*sin(E) gives surjectivity. Uniqueness + odd map => E(-M)=-E(M).")
    push!(c, certify_symbolic("FR22 time reversal: cos(phi) even", "24",
                              Symbolics.substitute(cp,Dict(t=>-t))-cp))
    push!(c, certify_symbolic("FR22 time reversal: sin(phi) odd", "24",
                              Symbolics.substitute(sp,Dict(t=>-t))+sp))
    push!(c, certify_symbolic("FR22 time reversal: radius even", "25",
                              Symbolics.substitute(rp,Dict(t=>-t))-rp))
    push!(c, certify_symbolic("FR22 circular limit: r=a, phi=M", "24–25",
        [Symbolics.substitute(rp,Dict(u=>0))-a,
         Symbolics.substitute(cp,Dict(u=>0))-C,
         Symbolics.substitute(sp,Dict(u=>0))-S,
         Symbolics.substitute(M,Dict(u=>0))-E]))
    println("ASSUMPTIONS true anomaly is the continuous lift with phi(0)=0; ",
            "its cosine/sine parity proves phi(-M)=-phi(M).")
    c
end
