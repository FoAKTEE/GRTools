"Thermal mass from the Bondi-radius equals scale-height condition, cs=h*a*n."
function certify_thermal_mass()
    @variables h a G Ms
    cs2 = h^2*G*Ms/a
    certify_symbolic("thermal mass M_th=h^3 M_star", "FR22 (8); FR24 section 2.2",
                     G*(h^3*Ms)/cs2-h*a)
end

"""
FR22 (3)–(8),(18). Positive r,a,G,Mstar,h,Sigmap; real p,q.
Power laws use the real logarithm, x=log(r/a). Logarithmic derivatives
establish the radial jets without treating noninteger powers as polynomials.
Omega and kappa denote positive roots where their squared values are positive.
"""
function certify_fr22_background()
    @variables r a h p q Sigmap GM x Sigma c2
    c = Certificate[]
    sigma_x = Sigmap*exp(-p*x)
    cs_x = h*exp(-q*x/2) # constant velocity scale sqrt(GM/a) factors out
    push!(c, certify_symbolic("FR22 (3) Sigma power-law derivative", "3",
                              Symbolics.derivative(sigma_x,x)+p*sigma_x))
    push!(c, certify_symbolic("FR22 (4) sound-speed power-law derivative", "4",
                              Symbolics.derivative(cs_x,x)+q*cs_x/2))
    push!(c, certify_symbolic("FR22 reference-radius normalisations", "3–4",
        [Symbolics.substitute(Sigmap/(r/a)^p,Dict(r=>a))-Sigmap,
         Symbolics.substitute(h/(r/a)^(q/2),Dict(r=>a))-h]))
    pressure_r = (-p*Sigma/r)*c2 + Sigma*(-q*c2/r)
    push!(c, certify_symbolic("FR22 (5) pressure product derivative", "5",
                              pressure_r+(p+q)*Sigma*c2/r))
    # For r/a>0, x^(-q)=1/x^q, x^(1-q)=x/x^q. Write the
    # source power laws in this quotient form before the rational normal form.
    sound2 = h^2*GM/a / (r/a)^q
    omega2 = GM/r^3 + pressure_r/(r*Sigma)
    paper_omega2 = GM/r^3*(1-h^2*(p+q)*(r/a)/(r/a)^q)
    push!(c, certify_symbolic("FR22 (6) radial pressure balance", "6",
        Symbolics.substitute(omega2,Dict(c2=>sound2))-paper_omega2))
    # 2Omega/r d(r²Omega)/dr = 4Omega²+r*d(Omega²)/dr;
    # c2'=-q*c2/r is the already certified power-law jet.
    omega2_r = Symbolics.derivative(omega2,r) +
               Symbolics.derivative(omega2,c2)*(-q*c2/r) +
               Symbolics.derivative(omega2,Sigma)*(-p*Sigma/r)
    kappa2 = 4omega2+r*omega2_r
    paper_kappa2 = GM/r^3*(1-h^2*(p+q)*(2-q)*(r/a)/(r/a)^q)
    push!(c, certify_symbolic("FR22 (7) epicyclic derivative", "7",
        Symbolics.substitute(kappa2,Dict(c2=>sound2))-paper_kappa2))
    push!(c, certify_symbolic("FR22 (18) thermal length", "18", (-q*c2/r)/c2+q/r))
    push!(c, certify_thermal_mass())
    push!(c, certify_symbolic("FR22 p+q=0: Omega=kappa=OmegaK", "6–7; section 3.3",
        [Symbolics.substitute(omega2,Dict(p=>-q))-GM/r^3,
         Symbolics.substitute(kappa2,Dict(p=>-q))-GM/r^3]))

    d = certificate_disc(h,p,q)
    # Call the actual runtime functions with Num arguments. Squaring removes
    # square roots algebraically BEFORE evaluation; no floating sqrt is sampled.
    push!(c, certify_symbolic("runtime Sigma = FR22 (3)", "3; a=Sigmap=1",
                              Modes.Sigma(d,r)-1/r^p))
    push!(c, certify_symbolic("runtime csiso = FR22 (4)", "4; a=n=1",
                              Modes.csiso(d,r)-h/r^(q/2)))
    push!(c, certify_symbolic("runtime Omega squared = FR22 (6)", "6; a=GM=1",
                              Modes.Omega(d,r)^2-(1-h^2*(p+q)*r^(1-q))/r^3))
    push!(c, certify_symbolic("runtime kappa squared = FR22 (7)", "7; a=GM=1",
                              Modes.kappa(d,r)^2-(1-h^2*(p+q)*(2-q)*r^(1-q))/r^3))
    push!(c, certify_symbolic("runtime LTinv = FR22 (18)", "18",
                              Modes.LTinv(d,r)+q/r))
    push!(c, certify_symbolic("runtime kappa derivative identity", "7",
        Modes.kappa(d,r)^2-4Modes.Omega(d,r)^2-2r*Modes.Omega(d,r)*Modes.dOmegadr(d,r)))
    println("ASSUMPTIONS r,a,h,GM,Sigmap>0; real power laws; positive Omega,kappa ",
            "branches with positive radicands. Squared symbolic identities imply root equality.")
    c
end
