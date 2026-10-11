"""
Runtime certificate for MR20 = arXiv:1911.01428, ApJ 892, 65 (2020):
(9), (24)–(27) reduce to (17)–(21) at t_c=t_c'=pole_disp=0.
The reference is the landed FR22 (16) certificate, not another transcription.
`cs` is adiabatic: c_iso²=cs²/gamma. Keep cs and gamma independent and
rational, avoiding a square root in the direct Rational{BigInt} calls.
"""
function mastereq_code_gate_expressions()
    s = fr22_master_residuals()
    (;r,S,Sr,O,Or,D,Dr,w,L,Lr,P,P1,P2,m,c2) = s.symbols
    @variables gamma cs entropy entropy_r buoyancy buoyancy_r
    args = (r,gamma,cs,O,Or,S,Sr,L,Lr,entropy,entropy_r,w,0,0,
            D,Dr,buoyancy,buoyancy_r,0,P,P1,P2,m)
    expected = (s.C1,0,Symbolics.substitute(s.C0,Dict(c2=>cs^2/gamma)),0,s.Psi,0)
    (;s,args,expected,gamma,cs,entropy,entropy_r,buoyancy,buoyancy_r)
end

"Constrained local jets; no Dfreq/omega_pat call and no singular-point rejection."
function mastereq_code_gate_point(e, rng, j)
    v = e.s.symbols
    # Independent uniform draws from S={1,...,1009}/17. The affine maps below
    # give signed slopes and stable kappa²; these points supplement a symbolic
    # proof, so no Schwartz–Zippel probability is asserted for this sampler.
    draw() = RatB(rand(rng,1:1009),17)
    r,S,O,cs,gamma,n = ntuple(_ -> draw(),6)
    p,q = draw()-30,draw()-30
    # Explicit thermal endpoints coexist with noninteger, signed random q.
    j <= 4 && (q = RatB((0,3//4,1,-1//2)[j]))
    j == 1 && (gamma = RatB(1))
    slope = draw()*17/1009
    Or = -slope*O/r
    Orr = (draw()-30)*O/r^2
    base,offset = rand(rng,1:17),rand(rng,1:11)
    m,l = isodd(j) ? (base,base+offset) : (base+offset,base)
    w = l*n-m*O
    K = 4O^2+2r*O*Or
    Kr = 10O*Or+2r*(Or^2+O*Orr)
    D,Dr = K-w^2,Kr+2m*w*Or
    # Failure at a pole is an error, never a discarded/resampled trial.
    @assert all(!iszero,(r,S,O,cs,gamma,w,D,m-m*O,K-(m-m*O)^2))
    disc = certificate_disc(cs,p,q)
    pt = Dict(v.r=>r,v.S=>S,v.Sr=>-p*S/r,v.O=>O,v.Or=>Or,v.K=>K,
              v.D=>D,v.Dr=>Dr,v.w=>w,v.c2=>cs^2/gamma,
              v.L=>Modes.LTinv(disc,r),v.Lr=>Modes.dLTinvdr(disc,r),
              v.P=>draw()-30,v.P1=>draw()-30,v.P2=>draw()-30,v.m=>RatB(m),
              e.cs=>cs,e.gamma=>gamma,e.entropy=>draw()-30,
              e.entropy_r=>draw()-30,e.buoyancy=>draw()-30,e.buoyancy_r=>draw()-30)
    (;pt,r,S,O,Or,Orr,cs,gamma,p,q,m,l,n,w,K,D,Dr)
end

"Negative control: perturb only the real C0 coefficient of the actual kernel."
function mastereq_coefficient_mutant(args...)
    value = Modes.master_equation_cool(args...)
    ntuple(k -> k == 3 ? value[k]+1/args[1]^2 : value[k],6)
end

"Negative control: use total frequency omega_p=m instead of l*n, consistently."
function mastereq_frequency_mutant(args...)
    O,Or,w,D,Dr,m = args[4],args[5],args[12],args[15],args[16],args[23]
    wrong_w = m-m*O
    wrong_D = D+w^2-wrong_w^2
    wrong_Dr = Dr-2m*w*Or+2m*wrong_w*Or
    wrong_args = Base.setindex(Base.setindex(Base.setindex(args,wrong_w,12),wrong_D,15),wrong_Dr,16)
    Modes.master_equation_cool(wrong_args...)
end

"Symbolic identity for all gamma/q plus direct exact-rational runtime and mutants."
function certify_fr22_mastereq_code_gate()
    e = mastereq_code_gate_expressions()
    c = Certificate[]
    names = ("C1R","C1I","C0R","C0I","PsiR","PsiI")
    labels = ("eq:c1_iso","eq:c1_iso","eq:c0_iso","eq:c0_iso","eq:psi_iso","eq:psi_iso")
    runtime = Modes.master_equation_cool(e.args...)
    for k in 1:6
        # Eliminate Or using D+w²=4Omega²+2rOmegaOmega'. This rational
        # parametrization is invertible for r*Omega != 0 and leaves w, D,
        # D', L, L', gamma, cs and all three forcing jets independent.
        residual = Symbolics.substitute(runtime[k]-e.expected[k],e.s.constraints)
        push!(c,certify_symbolic("mastereq-code-gate symbolic "*names[k],
            "1911.01428::"*labels[k],residual))
    end
    rng = MersenneTwister(221620)
    N = 64
    bad = zeros(Int,6)
    exact_types = true
    killed = zeros(Int,2)
    for j in 1:N
        point = mastereq_code_gate_point(e,rng,j)
        args = map(x -> RatB(exact_eval(x,point.pt)),e.args)
        expected = map(x -> RatB(exact_eval(x,point.pt)),e.expected)
        actual = Modes.master_equation_cool(args...)
        residuals = actual .- expected
        bad .+= .!iszero.(residuals)
        exact_types &= all(x -> x isa RatB,actual)
        mutations = (mastereq_coefficient_mutant(args...) .- expected,
                     mastereq_frequency_mutant(args...) .- expected)
        killed .+= map(x -> any(!iszero,x),mutations)
        println("EXACT_POINT mastereq-code-gate ",j," inputs=",Base.structdiff(point,(;pt=point.pt)),
                " residuals=",residuals," coefficient_mutant=",mutations[1],
                " frequency_mutant=",mutations[2])
    end
    for k in 1:6
        push!(c,certify_bool("mastereq-code-gate exact "*names[k],"1911.01428::"*labels[k],
            bad[k]==0,"$N direct Rational{BigInt} calls; residual exactly 0; seed=221620";
            modality=:exact_rational,npoints=N))
    end
    push!(c,certify_bool("mastereq-code-gate exact output types","runtime arithmetic",exact_types,
        "all six runtime outputs are Rational{BigInt} at all $N points";
        modality=:exact_rational,npoints=N))
    for (k,name) in enumerate(("C0 coefficient wrapper","omega_p=m wrapper"))
        push!(c,certify_bool("mastereq-code-gate KILLED "*name,"negative control",killed[k]==N,
            "$(killed[k])/$N points reject the mutant while the original has zero residual";
            modality=:exact_rational,npoints=N))
    end
    println("ACCOUNTING mastereq-code-gate: symbolic equality is primary; six structural zero residuals; ",
        "supplemental N=$N, seed=221620, S={1,...,1009}/17, affine maps and fixed endpoints in point(); ",
        "no rejection; degree bound/SZ bound not required or claimed for a symbolic proof. ",
        "r,Sigma,Omega,cs,gamma,w,D != 0; t_c=t_c'=pole_disp=0; c_iso²=cs²/gamma; ",
        "K=4Omega²+2rOmegaOmega'; w=l*n-m*Omega; D=K-w²; ",
        "D'=10OmegaOmega'+2r(Omega'²+OmegaOmega'')+2m*w*Omega'; ",
        "L=-q/r, L'=q/r²; arbitrary smooth jets, arbitrary gamma and q; l!=m. ",
        "MR20 source=1911.01428 (9),(24)-(27), LI (17)-(21), not 2007.13787.")
    c
end
