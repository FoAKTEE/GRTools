"FR24 v1 (2)–(8), exactly the FR22 profiles in a=Omega_p=Sigma_p=1 units."
function certify_fr24_background()
    @variables r h p q x S c2
    d = certificate_disc(h,p,q)
    c = Certificate[]
    # Reuse the landed general power-law derivation, including pressure balance
    # and the epicyclic derivative. Here the paper-qualified runtime crosswalk
    # supplies the labels needed by the shared-label owner table.
    append!(c,certify_fr22_background())
    sigma = Modes.Sigma(d,r)
    cs = Modes.csiso(d,r)
    for paper in ("2407.20398v1","2506.19917v3")
        push!(c,certify_symbolic("$paper surface density","$paper::eqn:surface_density",sigma-1/r^p))
        push!(c,certify_symbolic("$paper sound speed","$paper::eqn:sound_speed",cs-h/r^(q/2)))
        # Use x=log(r/a), so arbitrary real exponents have exact derivatives.
        aspect = h*exp(-q*x/2)/(exp(-3x/2)*exp(x))
        push!(c,certify_symbolic("$paper aspect-ratio logarithmic slope","$paper::eqn:aspect_ratio",
            Symbolics.derivative(aspect,x)-(1-q)/2*aspect))
        push!(c,certify_symbolic("$paper aspect-ratio reference value","$paper::eqn:aspect_ratio",
            Modes.csiso(d,Num(1))/(Modes.OmegaK(d,Num(1))*Num(1))-h))
    end
    push!(c,certify_symbolic("FR24 (4) pressure product jet","2407.20398v1::eqn:eos",
        (-p*S/r)*c2+S*(-q*c2/r)+(p+q)*S*c2/r))
    temp = h^2*exp(-q*x)
    push!(c,certify_symbolic("FR24 (5) temperature slope","2407.20398v1::eqn:temperature",
        Symbolics.derivative(temp,x)+q*temp))
    omega2 = (1-h^2*(p+q)*r^(1-q))/r^3
    kappa2 = (1-h^2*(p+q)*(2-q)*r^(1-q))/r^3
    omega_res = Modes.Omega(d,r)^2-omega2
    kappa_res = Modes.kappa(d,r)^2-kappa2
    push!(c,certify_symbolic("FR24 (7) exact runtime Omega squared","2407.20398v1::eqn:fluid_frequency",omega_res))
    push!(c,certify_symbolic("FR24 (8) exact runtime kappa squared","2407.20398v1::unlabelled-eq-8",kappa_res))
    push!(c,certify_symbolic("FR24 arbitrary-q thermal scale","2407.20398v1::eqn:thermal_scale",Modes.LTinv(d,r)+q/r))
    # Exact rational power-law samples use integer p,q; arbitrary real p,q are
    # covered by the symbolic/logarithmic proofs, not by finite sampling.
    rng = MersenneTwister(240708)
    # Remove paired roots before numeric substitution. Otherwise substituting
    # a rational radius evaluates fractional powers with floating arithmetic.
    squared_profiles = Symbolics.simplify.([Modes.Omega(d,r)^2,Modes.kappa(d,r)^2]; expand=true)
    # The simplifier keeps (3//2)*2 as the rational exponent 3//1;
    # Julia's Rational^Rational evaluates in floating arithmetic even here.
    squared_profiles = Symbolics.substitute.(squared_profiles,Ref(Dict(r^(3//1)=>r^3)))
    bad = 0
    for j in 1:32
        pt = Dict{Num,Any}(r=>RatB(rand(rng,1:17),7),h=>RatB(rand(rng,1:7),100),
                  p=>rand(rng,0:3),q=>rand(rng,0:3))
        vals = [exact_eval(lhs,pt)-exact_eval(rhs,pt)
                for (lhs,rhs) in zip(squared_profiles,(omega2,kappa2))]
        bad += count(!iszero,vals)
        println("EXACT_POINT background ",j," p=",pt[p]," q=",pt[q]," residuals=",vals)
    end
    push!(c,certify_bool("FR24 exact-rational squared-frequency samples","2407.20398v1::eqs-7-8",
        bad==0,"2 residuals × 32 points; seed=240708; integer p,q in 0:3; supplementary to symbolic proof";
        modality=:exact_rational,npoints=32))
    println("ASSUMPTIONS FR24 v1, a=Omega_p=Sigma_p=1; positive r,h; arbitrary real p,q; ",
        "positive frequency branches where radicands are positive. Omega²/kappa² are exact, ",
        "not an h expansion. Common FR22 background is executed as a prerequisite.")
    c
end
