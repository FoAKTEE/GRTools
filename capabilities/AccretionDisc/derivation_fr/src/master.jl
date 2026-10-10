"""
FR22 (13)–(21) with arbitrary smooth locally isothermal temperature.
Independent jets are constrained only by w'=-m*Omega', kappa²=4Omega²+2rOmegaOmega'
and D=kappa²-w². D' remains free (the next background jet supplies it).
L means 1/L_T; L*log(r*Sigma*L/D)' is written as L*A+L', regular at L=0.
"""
function fr22_master_residuals()
    @variables r S Sr O Or K D Dr w c2 L Lr h h1 h2 P P1 P2 m
    A = 1/r+Sr/S-Dr/D
    B = Sr/S+Or/O-Dr/D
    C1 = A-L
    C0 = -2m*O/(r*w)*(L+B)-L*A-Lr-m^2/r^2-D/c2
    Psi = -P2-A*P1+(2m*O/(r*w)*B+m^2/r^2)*P
    master = h2+C1*h1+C0*h-Psi
    # Cramer's rule for the two momentum equations: determinant is D.
    radial_rhs = -(h1-L*h+P1)
    azimuthal_rhs = -im*m/r*(h+P)
    determinant = (-im*w)^2+K
    ur_solved = (radial_rhs*(-im*w)+2O*azimuthal_rhs)/D
    up_solved = ((-im*w)*azimuthal_rhs-K/(2O)*radial_rhs)/D
    ur = im/D*(w*(h1+P1)-2m*O/r*(h+P)-w*L*h)
    up = (K/(2O)*(h1+P1)-m*w/r*(h+P)-K*L*h/(2O))/D
    ds = S*h/c2
    # Differentiate the reconstructed radial velocity through its actual jets.
    jets = (r=>1, S=>Sr, O=>Or, D=>Dr, w=>-m*Or, L=>Lr,
            h=>h1, h1=>h2, P=>P1, P1=>P2)
    ur_r = im*sum(Symbolics.derivative(imag(ur),v)*dv for (v,dv) in jets)
    continuity = -im*w*ds+S*ur_r+(S/r+Sr)*ur+im*m*S/r*up
    constraints = Dict(K=>D+w^2, Or=>(D+w^2-4O^2)/(2r*O))
    residuals = [Symbolics.substitute(determinant-D,constraints),
        ds*c2-S*h,
        (S*h1+Sr*h)/S-(Sr*c2+S*c2*L)*ds/S^2-(h1-L*h),
        ur_solved-ur, up_solved-up,
        Symbolics.substitute(-im*w*ur-2O*up-radial_rhs,constraints),
        Symbolics.substitute(-im*w*up+K/(2O)*ur-azimuthal_rhs,constraints),
        Symbolics.substitute(imag(continuity)*D/(S*w)-master,constraints)]
    (; symbols=(;r,S,Sr,O,Or,K,D,Dr,w,c2,L,Lr,h,h1,h2,P,P1,P2,m),
       A,B,C1,C0,Psi,master,ur,up,ds,constraints,residuals)
end

"Canonical arbitrary-q master, reconstruction, source crosswalk and free-wave AMF."
function certify_fr22_master()
    s = fr22_master_residuals()
    (;r,S,Sr,O,Or,K,D,Dr,w,c2,L,Lr,h,h1,h2,P,P1,P2,m) = s.symbols
    c = Certificate[]
    names = ("momentum determinant D", "density reconstruction", "locally isothermal pressure force",
             "radial Cramer reconstruction", "azimuthal Cramer reconstruction",
             "radial Euler reconstruction", "azimuthal Euler reconstruction", "continuity elimination to master")
    labels = ("eq:detuning","eq:dsigma_dh","eq:ur_eqn","eq:dur_dh","eq:duphi_dh",
              "eq:ur_eqn","eq:uphi_eqn","eq:master_eqn")
    for j in eachindex(names)
        push!(c,certify_symbolic("FR22 "*names[j],"2207.14637v2::"*labels[j],s.residuals[j]))
    end
    # MR20 (9),(17)–(21), Appendix B: independent notation and transcription.
    @variables phi phi1 phi2 omega_p q ell n
    wt = m*(omega_p-O)
    mr1 = 1/r+Sr/S-Dr/D-L
    mr0 = -2m*O/(r*wt)*(L+Sr/S+Or/O-Dr/D)-L*(1/r+Sr/S-Dr/D)-Lr-m^2/r^2-D/c2
    mrpsi = -phi2-(1/r+Sr/S-Dr/D)*phi1+(2m*O/(r*wt)*(Sr/S+Or/O-Dr/D)+m^2/r^2)*phi
    mrur = im/D*(wt*(h1+phi1)-2m*O/r*(h+phi)-wt*L*h)
    mrup = (K/(2O)*(h1+phi1)-m*wt/r*(h+phi)-K*L*h/(2O))/D
    crosswalk = Dict(omega_p=>O+w/m,phi=>P,phi1=>P1,phi2=>P2)
    mr_residuals = [Symbolics.substitute(a,crosswalk)-b for (a,b) in
                    zip((mr1,mr0,mrpsi,mrur,mrup),(s.C1,s.C0,s.Psi,s.ur,s.up))]
    for (label,residual) in zip(("eq:c1_iso","eq:c0_iso","eq:psi_iso","eq:ur_iso","eq:uphi_iso"),mr_residuals)
        push!(c,certify_symbolic("MR20 / FR22 "*label,"1911.01428::"*label,residual))
    end
    push!(c,certify_symbolic("MR20 (9) master form","1911.01428::eq:master",
        Symbolics.substitute(h2+mr1*h1+mr0*h-mrpsi,crosswalk)-s.master))
    push!(c,certify_symbolic("MR20 (16) closure","1911.01428::eq:eng_iso",c2*s.ds-S*h))
    push!(c,certify_symbolic("MR20 D with eccentric total frequency","1911.01428::eq:D",
        Symbolics.substitute(K-wt^2,Dict(omega_p=>ell*n/m))-(K-(ell*n-m*O)^2)))
    d = certificate_disc(Num(1),Num(1),q)
    thermal = Modes.LTinv(d,r)
    thermal_r = Modes.dLTinvdr(d,r)
    push!(c,certify_symbolic("FR22 (18) thermal length from runtime power law","2207.14637v2::eq:thermal_length",thermal+q/r))
    push!(c,certify_symbolic("MR20 thermal length from runtime power law","1911.01428::eq:LT",thermal+q/r))
    # All six shared FR24 labels; FD25 prints the master and definitions and
    # refers to FR22 for reconstruction (its absent labels are not invented).
    for paper in ("2407.20398v1","2506.19917v3")
        push!(c,certify_symbolic("$paper shared arbitrary-q master","$paper::eq:master_eqn",s.residuals[8]))
        push!(c,certify_symbolic("$paper resonance offset","$paper::eqn:resonance_offset",s.residuals[1]))
        push!(c,certify_symbolic("$paper thermal scale","$paper::eqn:thermal_scale",thermal+q/r))
    end
    for (label,residual) in zip(("eq:dsigma_dh","eq:dur_dh","eq:duphi_dh"),s.residuals[[2,4,5]])
        push!(c,certify_symbolic("FR24 shared reconstruction "*label,"2407.20398v1::"*label,residual))
    end
    qzero = Dict(L=>Symbolics.substitute(thermal,Dict(q=>0)),
                 Lr=>Symbolics.substitute(thermal_r,Dict(q=>0)))
    global_iso = [s.A,-2m*O/(r*w)*s.B-m^2/r^2-D/c2,
                 im/D*(w*(h1+P1)-2m*O/r*(h+P)),
                 (K/(2O)*(h1+P1)-m*w/r*(h+P))/D]
    push!(c,certify_symbolic("q=0 removes every thermal-length term","2207.14637v2::2_framework.tex:175",
        [Symbolics.substitute(a,qzero)-b for (a,b) in zip((s.C1,s.C0,s.ur,s.up),global_iso)]))
    # Phi=0, real coefficients, no corotation displacement/dissipation.
    @variables hr hi jr ji pi_exact
    # Substitute separately in the real linear operators for complex amplitudes.
    ur_free = im/D*(w*(jr+im*ji)-2m*O/r*(hr+im*hi)-w*L*(hr+im*hi))
    up_free = (K/(2O)*(jr+im*ji)-m*w/r*(hr+im*hi)-K*L*(hr+im*hi)/(2O))/D
    W = hi*jr-hr*ji # Im(h*conj(h'))
    F = pi_exact*m*r*S/D*W
    push!(c,certify_symbolic("MR20 free-wave Reynolds flux","1911.01428::eq:fj_iso",
        Symbolics.substitute(pi_exact*r^2*S*real(ur_free*conj(up_free))-F,s.constraints)))
    hrr_re, hrr_im = -s.C1*jr-s.C0*hr, -s.C1*ji-s.C0*hi
    W_r = hi*hrr_re-hr*hrr_im
    F_r = pi_exact*m*r*S/D*(s.A*W+W_r)
    flux_residuals = [F_r-L*F, F_r/c2-F*(c2*L)/c2^2,
                      Symbolics.substitute(F_r,qzero)]
    for (name,label,residual) in zip(("F_J'=L F_J","(F_J/c_s²)'=0","q=0 gives F_J'=0"),
        ("eq:fj_law_iso","eq:fj_law_iso","eq:fj_law_iso_full"),flux_residuals)
        push!(c,certify_symbolic("MR20 "*name,"1911.01428::"*label,residual))
    end
    # Num passes through the existing generated cooling kernel at instantaneous
    # cooling. Its cs is ADIABATIC, hence cs²=gamma*c2, even when tc=0.
    gamma = 7//5
    runtime = Modes.master_equation_cool(r,gamma,sqrt(gamma*c2),O,Or,S,Sr,L,Lr,
        Num(2),Num(3),w,0,0,D,Dr,Num(5),Num(7),0,P,P1,P2,m)
    runtime_residuals = [Symbolics.substitute(a-b,s.constraints) for (a,b) in
                         zip(runtime,(s.C1,0,s.C0,0,s.Psi,0))]
    push!(c,certify_symbolic("runtime Num instantaneous-cooling master and forcing",
        "1911.01428::eq:c1_iso; eq:c0_iso; eq:psi_iso",runtime_residuals))
    # Supplemental exact arithmetic: symbolic proofs above carry admission.
    # Fixed-denominator uniform integer coordinates, no singular rejection.
    # K and Or are eliminated by the rational constraints before evaluation.
    checks = Any[s.residuals...; mr_residuals...; runtime_residuals...; flux_residuals...]
    rng = MersenneTwister(221617)
    variables = vcat(collect(values(s.symbols)),[hr,hi,jr,ji,pi_exact])
    bad = 0
    for j in 1:32
        pt = Dict(v=>RatB(rand(rng,1:1009),17) for v in variables)
        vals = [exact_eval(res,pt) for res in checks]
        bad += count(!iszero,vals)
        println("EXACT_POINT master ",j," residuals=",vals)
    end
    push!(c,certify_bool("master exact-rational supplemental points","2207.14637v2::eq:master_eqn",
        bad==0,"$(length(checks)) residuals × 32 points; seed=221617; S={1,...,1009}/17; no rejection; symbolic proof primary";
        modality=:exact_rational,npoints=32))
    println("ASSUMPTIONS smooth real background; r,Sigma,Omega,c_s² nonzero; w,D nonzero; ",
        "kappa²=4Omega²+2rOmegaOmega'; w'=-mOmega'; arbitrary q. AMF: source-free real ",
        "operator away from resonances, no damping; F_J/c_s² is constant, F_J only at q=0. ",
        "FR24/FD25 typeset c_{s²} is read as c_s², matching their EoS and source definitions. ",
        "Exact samples supplement structural symbolic zeros; no probabilistic-only identity claim.")
    c
end
