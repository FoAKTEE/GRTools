"C-SOURCE: exact local algebra; reference C files are read-only literature inputs."
function source_certificates(; verbose=true, gate=get(ENV,"R3_SOURCE_GATE","all"),
                             source=EccentricSource)
    verbose && (println("[RUN] C-SOURCE chart contraction"); flush(stdout))
    C = Certificate[]
    function check(name, eqref, expr)
        (gate == "all" || occursin(gate,name)) || return
        verbose && (println("[START] ",name); flush(stdout))
        residual = expr()
        entries = residual isa AbstractArray ? residual : [residual]
        # Expand the cleared polynomial first; the generic simplify rewrite
        # engine can stall on the equivalent unexpanded chart contraction.
        # Rational reduction can factor its output. Expand the final numerator
        # as well; compare polynomial coefficients, never invoke trig simplify.
        reduced = [Symbolics.expand(first(_num_den(Symbolics.simplify_fractions(
            Symbolics.expand(e))))) for e in entries]
        c = certify_bool(name,eqref,all(x -> isequal(x,Num(0)),reduced),
            "$(length(entries)) cleared-polynomial components; exact zero = $(all(x -> isequal(x,Num(0)),reduced))";
            verbose,modality=:symbolic)
        verbose && !c.pass && (println("[RESIDUAL] ",reduced); flush(stdout))
        push!(C,c)
        verbose && (println("[END] ",name," ",c.pass ? "PASS" : "FAIL"); flush(stdout))
    end
    @variables R a M L ur dz dt x y r ep mass ut uphi
    @variables gtt gtphi grr gpp
    g = Num[gtt 0 gtphi 0; 0 grr 0 0; gtphi 0 gpp 0; 0 0 0 1]
    u = Num[ut,ur,uphi,0]
    ul = g*u
    dy = Num[dt,x,R*y,dz]
    full = transpose(dy)*source.chart_projector(g,u,R,Num(0))*dy
    # Symbolics retains sin(0),cos(0) as symbolic calls. Apply their exact
    # values before polynomial normalization, avoiding trig rewrite search.
    full = subs(full,Dict(sin(Num(0))=>Num(0),cos(Num(0))=>Num(1)))
    expanded = gtt*dt^2+2gtphi*dt*y+grr*x^2+gpp*y^2+dz^2+
               (ul[1]*dt+ul[2]*x+ul[3]*y)^2
    check("C-SOURCE full chart contraction before equal time", "G26 (52)", () ->full-expanded)
    eqtime = subs(full,Dict(dt=>0,dz=>0))
    check("C-SOURCE equal-time radial cross term", "DESIGN section 2", () ->
        eqtime-((grr+ul[2]^2)*x^2+2ul[2]*ul[3]*x*y+(gpp+ul[3]^2)*y^2))
    # Spatial BL metric is positive outside the horizon. Positive squares are
    # an explicit SOS certificate: eps²+grr*x²+gpp*y²+(u_r*x+u_phi*y)².
    check("C-SOURCE softened positivity SOS", "grr,gpp>0; ep>0", () ->
        ep^2+eqtime-(ep^2+grr*x^2+gpp*y^2+(ul[2]*x+ul[3]*y)^2))
    @variables psi
    verbose && (println("[RUN] C-SOURCE derivative and component gates"); flush(stdout))
    jet = source.distance_jet(g,u,R,r,psi,ep)
    check("C-SOURCE analytic distance derivatives", "fixed event", () ->
        [R^2*(Symbolics.derivative(jet.s2,r)-jet.s2_r),
         R^2*(Symbolics.derivative(jet.s2_r,r)-jet.s2_rr)])
    check("C-SOURCE runtime distance equals chart contraction", "G26 (52); fixed event", () ->
        jet.s2-ep^2-subs(eqtime,Dict(x=>r*cos(psi)-R,y=>r*sin(psi)/R)))
    # Num flows through the real field implementation. Clear positive root
    # powers only after evaluation; the distance-derivative gate above supplies
    # Qr=dQ/dr, Qrr=d²Q/dr² for these chain-rule identities.
    field = source.field_jet(g,u,R,r,psi,ep,mass)
    root = sqrt(field.s2)
    hnum = field.h*root
    check("C-SOURCE analytic field derivatives", "runtime jets; s>0", () ->
        vcat(vec(2root^3*field.h_r+hnum*jet.s2_r),
             vec(4root^5*field.h_rr-hnum*(3jet.s2_r^2-2jet.s2*jet.s2_rr))))
    check("C-SOURCE h_tr contraction and factor four", "G26 (52)", () ->
        hnum[1,2]-4mass*(gtt*ut+gtphi*uphi)*grr*ur)
    check("C-SOURCE h_rphi contraction and factor four", "G26 (52)", () ->
        hnum[2,3]-4mass*grr*ur*(gtphi*ut+gpp*uphi))
    # Contract the runtime tensor, retaining the normalization constraint as
    # an explicit residual: u.u=-1 makes the trace equal to 4mass/s.
    detg = gtt*gpp-gtphi^2
    ginv = Num[gpp/detg 0 -gtphi/detg 0; 0 1/grr 0 0;
               -gtphi/detg 0 gtt/detg 0; 0 0 0 1]
    check("C-SOURCE trace reversal normalization", "runtime contraction; u.u=-1", () ->
        sum(ginv.*hnum)-4mass-4mass*(dot(u,ul)+1))

    # The actual distance Hessian at the event supplies the local quadratic
    # form. This links the independent coefficients below to distance_jet.
    at_event = Dict(r=>R,psi=>0,ep=>0)
    atzero(v) = subs(subs(v,at_event),Dict(sin(Num(0))=>Num(0),cos(Num(0))=>Num(1)))
    check("C-SOURCE runtime local quadratic Hessian", "unsoftened; z=R*dtheta", () ->
        [atzero(jet.s2),atzero(jet.s2_r),atzero(Symbolics.derivative(jet.s2,psi)),
         atzero(jet.s2_rr)/2-(grr+ul[2]^2),
         atzero(Symbolics.derivative(jet.s2_r,psi))-2ul[2]*ul[3],
         atzero(Symbolics.derivative(Symbolics.derivative(jet.s2,psi),psi))/2-(gpp+ul[3]^2)])

    # Independently transcribed EffectiveSource coefficients, C lines 207-212.
    verbose && (println("[RUN] C-SOURCE EffectiveSource coefficients"); flush(stdout))
    Delta = a^2+R*(R-2M)
    gp = R^2+a^2*(R+2M)/R
    ref = SourceNearZone.coefficients(a,M,R,L,ur)
    # ur is a free signed symbol: one identity covers both radial branches.
    # Smoothness outside R*Delta*G=0 then gives an O(|Delta x|^3)
    # remainder by Taylor's theorem; finite-step scaling is checked separately.
    check("C-SOURCE EffectiveSource exact quadratic both ur signs",
        "C:207-212; dphib=dphi-c*dr; R*Delta*G!=0", () ->
        [ref.alpha20+ref.beta*ref.c^2/4-(R^2/Delta+R^4*ur^2/Delta^2),
         -ref.beta*ref.c/2-2R^2*ur*L/Delta,
         ref.beta/4-(gp+L^2), ref.alpha02-R^2])
    # B1 as printed, a_star=A, chi=a/M; use cos(2psi)=2cos²psi-1.
    verbose && (println("[RUN] C-SOURCE circular and Newtonian limits"); flush(stdout))
    @variables A chi Omega cp sp
    denom = (A-2)*A+chi^2
    b = ut^2*r^2*(M*Omega*(A^3+(A+2)*chi^2)-2chi)^2/(2A^4)
    B1 = ep^2+A^4*M^2/denom-2A^3*M*r*cp/denom+b+
        r^2*((A+2)*chi^2/A^3+A^2/denom+1)/2+
        (2cp^2-1)*(r^2*(-(A+2)*chi^2/A^3+A^2/denom-1)/2-b)
    gc = A^2/denom
    gpc = M^2*(A^2+chi^2*(1+2/A))
    gtc = -2M*chi/A
    circg = Num[-(1-2/A) 0 gtc 0; 0 gc 0 0; gtc 0 gpc 0; 0 0 0 1]
    circjet = source.distance_jet(circg,Num[ut,0,ut*Omega,0],A*M,r,psi,ep)
    # Extract the even quadratic coefficient before imposing sp²+cp²=1.
    # Certify the extraction too: no omitted odd or higher powers can pass.
    trig = Dict(sin(psi)=>sp,cos(psi)=>cp)
    circ = Symbolics.expand(subs(circjet.s2,trig))
    circ0 = subs(circ,Dict(sp=>0))
    circ2 = subs(circ,Dict(sp=>1))-circ0
    check("C-SOURCE exact circular B1 distance", "G26 appendix:31-46; polynomial circle", () ->
        [circ-circ0-circ2*sp^2,2A^4*denom*(circ0+circ2*(1-cp^2)-B1)])
    flatg = Matrix(Diagonal(Num[-1,1,R^2,1]))
    flatfield = source.field_jet(flatg,Num[1,0,0,0],R,r,psi,ep,mass)
    flat = Symbolics.expand(subs(flatfield.s2,trig))
    flat0 = subs(flat,Dict(sp=>0))
    flat2 = subs(flat,Dict(sp=>1))-flat0
    check("C-SOURCE Newtonian Plummer distance and h_tt=-2Phi", "weak-field static leading term", () ->
        [flat-flat0-flat2*sp^2,flat0+flat2*(1-cp^2)-(ep^2+r^2+R^2-2r*R*cp),
         flatfield.h[1,1]*sqrt(flatfield.s2)-2mass])
    # Derivation of the cubic constant, after elementary Taylor bounds:
    # |e|_1<=K*delta², |z|_1<=2delta, delta<=1.
    @variables delta k qmax cb beta0
    bound = qmax*(4k+k^2)*delta^3+beta0*(1+cb)^4*delta^3/48
    raw = qmax*(4k*delta^3+k^2*delta^4)+beta0*(1+cb)^4*delta^4/48
    check("C-SOURCE cubic remainder bound nonnegative slack", "0<delta<=min(1,R/2); R>0", () ->
        bound-raw-(1-delta)*delta^3*(qmax*k^2+beta0*(1+cb)^4/48))
    C
end
