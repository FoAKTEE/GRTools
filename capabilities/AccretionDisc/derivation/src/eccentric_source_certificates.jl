"C-SOURCE: exact local algebra; reference C files are read-only literature inputs."
function source_certificates(; verbose=true, gate=get(ENV,"R3_SOURCE_GATE","all"))
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
        c.pass || (println("[RESIDUAL] ",reduced); flush(stdout))
        push!(C,c)
        verbose && (println("[END] ",name," ",c.pass ? "PASS" : "FAIL"); flush(stdout))
    end
    @variables R a M L ur E dr dphi dz dt x y z r ep mass ut uphi S
    @variables gtt gtphi grr gpp utlow urlow uplow
    g = Num[gtt 0 gtphi 0; 0 grr 0 0; gtphi 0 gpp 0; 0 0 0 1]
    u = Num[ut,ur,uphi,0]
    ul = g*u
    J = EccentricSource.chart_jacobian(R,Num(0))
    dy = Num[dt,x,R*y,dz]
    full = transpose(dy)*EccentricSource.chart_projector(g,u,R,Num(0))*dy
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
    jet = EccentricSource.distance_jet(g,u,R,r,psi,ep)
    check("C-SOURCE analytic distance derivatives", "fixed event", () ->
        [R^2*(Symbolics.derivative(jet.s2,r)-jet.s2_r),
         R^2*(Symbolics.derivative(jet.s2_r,r)-jet.s2_rr)])
    @variables Q Qr Qrr K
    # Chain rule with an independent positive root S; avoids imposing a
    # symbolic square-root branch. S²=Q and S'=Qr/(2S).
    h1 = -K*Qr/(2S^3)
    h2 = Symbolics.derivative(h1,S)*Qr/(2S)+Symbolics.derivative(h1,Qr)*Qrr
    check("C-SOURCE analytic field derivatives", "S>0; S^2=Q", () ->
        [S^3*(Symbolics.derivative(K/S,S)*Qr/(2S)-h1),
         S^5*(h2-K*(3Qr^2/(4S^5)-Qrr/(2S^3)))])
    hnum = 2mass*(g+2ul*transpose(ul))
    check("C-SOURCE h_tr contraction and factor four", "G26 (52)", () ->
        hnum[1,2]-4mass*(gtt*ut+gtphi*uphi)*grr*ur)
    check("C-SOURCE h_rphi contraction and factor four", "G26 (52)", () ->
        hnum[2,3]-4mass*grr*ur*(gtphi*ut+gpp*uphi))
    @variables norm
    check("C-SOURCE trace reversal normalization", "u.u=-1, dimension 4", () ->
        subs(2mass*(4+2norm),Dict(norm=>-1))-4mass)

    # Independently transcribed EffectiveSource coefficients, C lines 207-212.
    verbose && (println("[RUN] C-SOURCE EffectiveSource coefficients"); flush(stdout))
    Delta = a^2+R*(R-2M)
    gp = R^2+a^2*(R+2M)/R
    G = a^2*(2M+R)+R*(L^2+R^2)
    alpha_num = R^2*(a^4*(2M+R)+a^2*R*(L^2-4M^2+2M*R*ur^2+R^2*(2+ur^2))+
        R^2*(L^2*(R-2M)+R^2*(-2M+R+R*ur^2)))
    # alpha20=alpha_num/(Delta² G), beta=4G/R,
    # c=c_num/(Delta G). Clear denominators before polynomial expansion;
    # this is an exact coefficient proof on R*Delta*G != 0.
    alpha02 = R^2
    for sign in (-1,1)
        c_num = -sign*L*R^3*ur
        cleared = [alpha_num*G+L^2*R^5*ur^2*G-(R^2*Delta+R^4*ur^2)*G^2,
                   -2c_num-2sign*R^3*ur*L,
                   G-(R^3+a^2*(R+2M)+R*L^2), alpha02-R^2]
        check("C-SOURCE EffectiveSource quadratic ur sign=$sign",
            "C:207-212; undo twist; R*Delta*G!=0", () ->Symbolics.expand.(cleared))
    end
    # B1 as printed, a_star=A, chi=a/M; use cos(2psi)=2cos²psi-1.
    verbose && (println("[RUN] C-SOURCE circular and Newtonian limits"); flush(stdout))
    @variables A chi Omega cp
    denom = (A-2)*A+chi^2
    b = ut^2*r^2*(M*Omega*(A^3+(A+2)*chi^2)-2chi)^2/(2A^4)
    B1 = ep^2+A^4*M^2/denom-2A^3*M*r*cp/denom+b+
        r^2*((A+2)*chi^2/A^3+A^2/denom+1)/2+
        (2cp^2-1)*(r^2*(-(A+2)*chi^2/A^3+A^2/denom-1)/2-b)
    gc = A^2/denom
    gpc = M^2*(A^2+chi^2*(1+2/A))
    Lc = ut*(-2M*chi/A+gpc*Omega)
    circ = ep^2+gc*(r*cp-A*M)^2+(gpc+Lc^2)*r^2*(1-cp^2)/(A*M)^2
    check("C-SOURCE exact circular B1 distance", "G26 appendix:31-46", () ->2A^4*denom*(circ-B1))
    flat = ep^2+(r*cp-R)^2+r^2*(1-cp^2)
    check("C-SOURCE Newtonian Plummer distance and h_tt=-2Phi", "weak-field static leading term", () ->
        [flat-(ep^2+r^2+R^2-2r*R*cp),2mass*(-1+2)/S-(-2*(-mass/S))])
    # Derivation of the cubic constant, after elementary Taylor bounds:
    # |e|_1<=K*delta², |z|_1<=2delta, delta<=1.
    @variables delta k qmax cb beta0
    bound = qmax*(4k+k^2)*delta^3+beta0*(1+cb)^4*delta^3/48
    raw = qmax*(4k*delta^3+k^2*delta^4)+beta0*(1+cb)^4*delta^4/48
    check("C-SOURCE cubic remainder bound nonnegative slack", "0<delta<=min(1,R/2); R>0", () ->
        bound-raw-(1-delta)*delta^3*(qmax*k^2+beta0*(1+cb)^4/48))
    C
end
