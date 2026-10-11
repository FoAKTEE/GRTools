# C-LI-N consumes the certified R3 LI graph and the canonical FR master.
# MR20 reference: ref-paper/arXiv-1911.01428/src/main.tex:170–195.
module LINewtonian
using Symbolics, Random
import Main.R3Derivation as R3
import FRDerivation as FR
const Q = R3.RatB

function expressions()
    e = FR.mastereq_code_gate_expressions()
    f = e.s
    (;r,S,Sr,O,Or,K,D,Dr,w,c2,L,Lr,h,h1,h2,P,P1,P2,m) = f.symbols
    @variables rN ON O1N O2N SN cN pN qN omegaN mN csN sN YN
    @variables nuN psiN muN dragN psi1N mu1N drag1N
    Kval = 4ON^2+2rN*ON*O1N
    wval = omegaN-mN*ON
    Dval = Kval-wval^2
    Drval = 10ON*O1N+2rN*(O1N^2+ON*O2N)+2mN*wval*O1N
    jet = Dict(r=>rN,S=>SN,Sr=>-pN*SN/rN,O=>ON,Or=>O1N,K=>Kval,
        D=>Dval,Dr=>Drval,w=>wval,c2=>cN,L=>-qN/rN,Lr=>qN/rN^2,m=>mN,
        e.cs=>csN,e.gamma=>csN^2/cN,e.entropy=>2,e.entropy_r=>3,
        e.buoyancy=>5,e.buoyancy_r=>7)
    # First-order limit in y0=s²*h, y1=s*Y. Its coefficients will be tested
    # against the R3 graph, not assumed to follow from a gamma substitution.
    a = L+2m*O/(r*w); b = D/w
    c = w/c2-m^2/(r^2*w); d = -1/r-Sr/S-m*K/(2r*w*O)
    f0 = 2m*O*P/(r*w)-P1; f1 = -m^2*P/(r^2*w)
    yp = a*h+b*YN+f0
    zp = c*h+d*YN+f1
    up = m*(h+P)/(r*w)+K*YN/(2O*w)
    rows = (yp,zp,S*h/c2,up/r,r*up)
    # Consume the same on-shell expressions as the certified LI exports.
    raw = (R3.LocallyIsothermal.MASTER_ONSHELL[1]*R3.y0+
           R3.LocallyIsothermal.MASTER_ONSHELL[2]*R3.y1+R3.LocallyIsothermal.MASTER_ONSHELL[5],
           R3.LocallyIsothermal.MASTER_ONSHELL[3]*R3.y0+
           R3.LocallyIsothermal.MASTER_ONSHELL[4]*R3.y1+R3.LocallyIsothermal.MASTER_ONSHELL[6],
           R3.LocallyIsothermal.ELIMINATED_ONSHELL.deps,
           R3.LocallyIsothermal.ELIMINATED_ONSHELL.duphi_up,
           R3.LocallyIsothermal.ELIMINATED_ONSHELL.duphi_low)
    # Rational conic gives v=s*r*Omega+O(s³), W=1+v²/2+O(s⁴) exactly.
    t = sN*rN*ON/2
    vr = 2(1-t^2)/(1+t^2)^2*sN*(ON+rN*O1N)/2
    scaling = Dict(R3.enu=>1+sN^2*nuN,R3.epsi=>rN*(1+sN^2*psiN),
        R3.emu=>1+sN^2*muN,R3.omt=>sN^3*dragN,R3.omt_r=>sN^3*drag1N,
        R3.psi_r=>1/rN+sN^2*psi1N,R3.mu_r=>sN^2*mu1N,
        R3.v=>2t/(1+t^2),R3.v_r=>vr,R3.W=>(1+t^2)/(1-t^2),
        R3.eps=>SN,R3.eps_r=>-pN*SN/rN,R3.p=>sN^2*SN*cN,
        R3.p_r=>-sN^2*(pN+qN)*SN*cN/rN,R3.m=>mN,R3.omega_d=>sN*wval)
    # Extract each amplitude coefficient before scaling; include zeros. The
    # Newtonian scalar source has htt=-2s²Phi, htt_r=-2s²Phi'. The other metric
    # components enter at higher PN orders and are outside this scalar limit.
    cases = []
    for (i,(name,order)) in enumerate(zip(("hprime","Yprime","density","uphi_up","uphi_low"),(2,1,0,1,1)))
        for (x,z,power,gain) in ((R3.y0,h,2,1),(R3.y1,YN,1,1),
                               (R3.htt,P,2,-2),(R3.htt_r,P1,2,-2))
            actual = gain*sN^power*R3.subs(R3.dsym(raw[i],x),scaling)
            expected = Symbolics.substitute(Symbolics.derivative(rows[i],z),jet)*sN^order
            push!(cases,(name=name*"/"*string(z),order,actual,expected))
        end
    end
    # Differentiate the certified first-order limit, then eliminate Y. No new
    # transcription of the FR master/reconstruction is introduced.
    diff(x) = sum(Symbolics.derivative(x,v)*dv for (v,dv) in
        (r=>1,S=>Sr,O=>Or,D=>Dr,w=>-m*Or,L=>Lr,P=>P1,P1=>P2))
    ysolve = (h1-a*h-f0)/b
    second = h2-diff(a)*h-a*h1-diff(b)*ysolve-b*(c*h+d*ysolve+f1)-diff(f0)
    constraints = f.constraints
    residuals = Pair{String,Any}[]
    for z in (h2,h1,h,P2,P1,P)
        push!(residuals,"MR20 coefficient "*string(z)=>Symbolics.substitute(
            Symbolics.derivative(second-f.master,z),constraints))
    end
    for (name,got,want) in (("density",rows[3],f.ds),("radial",im*ysolve,f.ur),
                            ("azimuthal",Symbolics.substitute(up,Dict(YN=>ysolve)),f.up))
        for z in (h,h1,P,P1)
            push!(residuals,"FR reconstruction $name/$(z)"=>R3.subs(R3.dsym(got-want,z),constraints))
        end
    end
    (;e,f,jet,raw,scaling,cases,residuals,independent=[rN,ON,O1N,O2N,SN,cN,pN,qN,omegaN,csN,
        nuN,psiN,muN,dragN,psi1N,mu1N,drag1N],mN,sN,qN)
end

function certify(;outdir,npts=24,seed=191101430,grid=1048573)
    npts>0 || throw(ArgumentError("positive sample count required"))
    mkpath(outdir)
    ex=expressions(); certs=FR.Certificate[]
    for (name,res) in ex.residuals
        push!(certs,FR.certify_symbolic(name,"C-LI-N; MR20/FR22",res))
    end
    flush(stdout)
    degrees = open(joinpath(outdir,"cleared-numerators.txt"),"w") do io
        [begin
            println(io,"CASE ",c.name,"; target coefficient of s^",c.order)
            R3.rational_circuit(io,Meta.parse(string(real(c.actual-c.expected))*" + im*("*string(imag(c.actual-c.expected))*")"))
         end for c in ex.cases]
    end
    rng=MersenneTwister(seed); bad=zeros(Int,length(ex.cases)); runtime_bad=zeros(Int,8)
    # Exact-rational grid, IID coordinates; m is an independent integer draw.
    # No rejected points. q=0,3/4,1 are supplemental endpoints after the IID set.
    open(joinpath(outdir,"exact-points.txt"),"w") do io
        for j in 1:npts+3
            pt=Dict(v=>Q(rand(rng,1:grid),grid) for v in ex.independent)
            pt[ex.mN]=Q(rand(rng,1:grid))
            j>npts && (pt[ex.qN]=Q((0,3//4,1)[j-npts]))
            env=Dict{Any,R3.LSeries}(Symbolics.value(v)=>R3.series_const(x) for (v,x) in pt)
            env[Symbolics.value(ex.sN)]=R3.series_mono(Q(1),1)
            defects=[]
            for (k,c) in enumerate(ex.cases)
                delta=R3.eval_series(c.actual-c.expected,env)
                parts=delta isa R3.LSeries ? (delta,) : delta
                ok=all(a->R3._lim(a)>c.order && (R3.lorder(a)===nothing || R3.lorder(a)>c.order),parts)
                bad[k]+=!ok
                push!(defects,(c.name,[(a.n0, a.n0<=c.order ? a.c : Q[]) for a in parts]))
            end
            # Actual generated cooling kernel, coefficient-by-coefficient in
            # Phi,Phi',Phi''. Response coefficients are independently returned.
            value(x)=Q(FR.exact_eval(Symbolics.substitute(Num(x),ex.jet),pt))
            for (k,z) in enumerate((ex.f.symbols.P,ex.f.symbols.P1,ex.f.symbols.P2))
                basis=Dict(ex.f.symbols.P=>Int(k==1),ex.f.symbols.P1=>Int(k==2),ex.f.symbols.P2=>Int(k==3))
                args=map(x->value(Symbolics.substitute(Num(x),basis)),ex.e.args)
                @assert !iszero(args[12]) && !iszero(args[15])
                got=FR.Modes.master_equation_cool(args...)
                want=map(x->value(Symbolics.substitute(Num(x),basis)),ex.e.expected)
                residual=got.-want
                k==1 && (runtime_bad[1:4].+=.!iszero.(residual[1:4]))
                runtime_bad[4+k]+=!iszero(residual[5]); runtime_bad[8]+=!iszero(residual[6])
                println(io,"point=$j source_basis=$k runtime_residuals=$residual")
            end
            println(io,"point=$j independent=",sort!(collect(pt);by=x->string(first(x)))," scaled_residual_orders=",defects)
            flush(io)
        end
    end
    # Coefficient extraction cannot increase the total degree of a polynomial.
    # The cleared numerator circuit bounds its s-leading coefficient too. A
    # nonzero rational leading defect implies a nonzero numerator coefficient.
    bound=sum(min(Q(1),Q(d[1],grid))^npts for d in degrees)
    for (k,c) in enumerate(ex.cases)
        push!(certs,FR.certify_bool("scaled R3 "*c.name,"C-LI-N",bad[k]==0,
            "N=$npts IID + 3 q endpoints; d<=$(degrees[k][1]); no rejection";
            modality=:exact_rational,npoints=npts))
    end
    for (k,name) in enumerate(("C1R","C1I","C0R","C0I","Psi/Phi","Psi/Phi1","Psi/Phi2","PsiI"))
        push!(certs,FR.certify_bool("runtime beta=0 "*name,"MR20 170–195",runtime_bad[k]==0,
            "independent forcing basis; csadi rational, gamma=csadi²/cT²; tc=tc'=pole=0";
            modality=:exact_rational,npoints=npts))
    end
    open(joinpath(outdir,"sampling.txt"),"w") do io
        println(io,"N=$npts seed=$seed S={1,...,$grid}/$grid; independent m in {1,...,$grid}; no rejection.")
        println(io,"20 scaled residuals; numerator max total degree=$(maximum(first.(degrees))); union bound=$bound (approximately $(Float64(bound))).")
        println(io,"PRNG IID model over Q; bound applies to any nonzero leading numerator. Coefficient in s has degree <= full numerator. Pole hit aborts, never resamples.")
        println(io,"Nonzero r,Sigma,Omega,cT²,csadi,w,D. Other R3 denominator loci are recorded in cleared-numerators.txt. s->0 formal: invert leading denominator coefficients only.")
        println(io,"Scaled enu,emu=1+O(s²), epsi=r[1+O(s²)], dragging=O(s³), v conic=s*r*Omega+O(s³), W exact conic, p=s²Sigma*cT², y0=s²h,y1=sY, htt=-2s²Phi. Scalar Newtonian metric source only.")
        println(io,"Arbitrary rational Omega jets, not irrational power-law runtime profiles. Endpoint q checks do not enter bound. Runtime identities are additionally covered by the separate symbolic code gate.")
    end
    certs
end
end
