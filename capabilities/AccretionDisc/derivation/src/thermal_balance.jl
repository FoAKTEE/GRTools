"∇_μ T^{μν}, evaluated directly from a contravariant stress tensor."
function stress_divergence_up(T, Γ, ∂)
    n = size(T, 1)
    [sum(∂(T[μ, ν], μ) for μ in 1:n) +
     sum(Γ[μ, μ, σ] * T[σ, ν] + Γ[ν, μ, σ] * T[μ, σ]
         for μ in 1:n, σ in 1:n) for ν in 1:n]
end

"Conserved baryons: ∇_μ(ρu^μ) = u^μ∂_μρ + ρ∇_μu^μ."
baryon_conservation(Γ, u, ρ, ∂) =
    sum(u[μ] * ∂(ρ, μ) for μ in eachindex(u)) + ρ * divergence_u(Γ, u, ∂)

"""
Off-shell stress and coordinate-current identities, and the OPTION A solution.
The stress is the existing TUP_JET, without any adiabatic U/V/W tensor.
The coordinate currents obey
  ∇_μ T^μ_a = (∇_μ T^{μν})g_{νa} + T^{μν}∂_a g_{μν}/2,
where the left side is the divergence of the vector with fixed coordinate label a.
For energy use -T^μ_t; for angular momentum use T^μ_phi. The metric work terms
remain for a nonstationary/nonaxisymmetric perturbation; neither is bath work.

Relation (19) selects zero advected Lagrangian integration constants, as in G26.
The mode construction requires u^t*omega_d != 0. It does not settle stationary
corotation modes, a physical p/rho cooling law, or quadratic wave-flux transport.
"""
function li_balance_expressions()
    LI = LocallyIsothermal
    D = stress_divergence_up(TUP_JET, GAMMAJET, ∂jet)
    split = [D[a] - EULER_JET[a] - UJET[a] * ENERGY_JET for a in 1:4]
    currents = [sum(TUP_JET[μ, ν] * GJET[ν, a] for ν in 1:4)
                for μ in 1:4, a in (IT, IPHI)]
    current_div = [divergence_u(GAMMAJET, currents[:, j], ∂jet) for j in 1:2]
    metric_work = [sum(TUP_JET[μ, ν] * ∂jet(GJET[μ, ν], a)
                       for μ in 1:4, ν in 1:4) / 2 for a in (IT, IPHI)]
    Dlow = lower(GJET, D)
    current_identity = [current_div[j] - metric_work[j] - Dlow[a]
                        for (j, a) in enumerate((IT, IPHI))]

    @variables rho_b rho_b_r delta_rho delta_rho_r
    rhojet = Jet(rho_b, delta_rho)
    partial = (x, a) -> isequal(x, rhojet) ?
        Jet(a == IR ? rho_b_r : Num(0),
            a == IR ? delta_rho_r : ∂amp(delta_rho, a)) : ∂jet(x, a)
    mass = pert(baryon_conservation(GAMMAJET, UJET, rhojet, partial))
    # No radial derivative of delta_rho is needed: u^r_background=0.
    advected = delta_rho/rho_b - deps/(eps+p) +
               XI_R*(rho_b_r/rho_b - eps_r/(eps+p))
    mass_energy = mass/rho_b - DENERGY/(eps+p) + im*UT*(omega-m*Ω_expr)*advected

    # Reconstruct the LI displacement using the existing (16) operator/solver.
    # Fg is a temporary azimuthal unknown; Fgc is the arbitrary parallel part.
    du = [LI.lin_subs(LI.eliminate(a), LI.MASTER_SUB) for a in DUUP]
    xir = im*du[IR]/(UT*(omega-m*Ω_expr))
    trial = [-ULOW[IPHI]*Fg/ULOW[IT], xir, Fg, Num(0)]
    azimuth = solve_linear(rhs16(trial)[IPHI]-du[IPHI], Fg)
    xi = [lin_subs(a, Dict(Fg=>azimuth)) + Fgc*UUP[j] for (j,a) in enumerate(trial)]
    xilow = lower(METRIC, xi)
    nabla = [∂amp(xilow[b], a)-sum(GAMMA[s,a,b]*xilow[s] for s in 1:4)
             for a in 1:4, b in 1:4]
    volume = sum(QUP[a,b]*(HMODE[a,b]+nabla[a,b]+nabla[b,a])
                 for a in 1:4, b in 1:4)/2
    delta_eps = LI.lin_subs(LI.DEPS_SOL, LI.MASTER_SUB)
    relation19 = LI.backsub((delta_eps + xir*eps_r)/(eps+p) + volume)
    delta_rho19 = rho_b*(delta_eps+xir*eps_r)/(eps+p) - xir*rho_b_r
    mass19 = LI.backsub(lin_subs(mass, Dict(delta_rho=>delta_rho19)))
    # Arbitrary nonzero density/gradient for rational evaluation; the preceding
    # symbolic identity certifies the unrestricted baryon/energy relation.
    mass19 = subs(mass19, Dict(rho_b=>5//3, rho_b_r=>-7//5))
    pressure19 = LI.backsub((eps+p)*y0 + xir*p_r +
                            (p/eps)*(eps+p)*volume - eps*LI.THETA_R*xir)
    on_solution = vcat([LI.backsub(DENERGY), mass19, relation19, pressure19],
        [LI.backsub(pert(d)) for d in D],
        [LI.backsub(pert(current_div[j]-metric_work[j])) for j in 1:2])
    (; split, current_identity, mass_energy, relation19, pressure19, on_solution)
end

"""
C-LI-HEAT: symbolic tensor identities plus exact constrained-jet regressions.
The local rest-frame proof uses arbitrary first derivatives at one event,
u=(1,0,0,0), ∂_a u^t=0, Γ=0. Tensoriality makes the projections general for
normalized u and a metric-compatible torsion-free connection. No field equation
or equation of state is imposed in that proof.
"""
function li_energy_balance_certificates(; npts::Int=12, seed::Int=191101429, verbose::Bool=true)
    npts > 0 || throw(ArgumentError("npts must be positive"))
    C = Certificate[]
    @variables x[1:4] eg pg rg ge[1:4] gp[1:4] gr[1:4] B[1:3,1:4]
    origin = Dict(x[a]=>0 for a in 1:4)
    at0(a) = subs(a, origin)
    ∂(a, j) = Symbolics.derivative(a, x[j])
    g = Num[-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    Γ = fill(Num(0), 4, 4, 4)
    u = vcat(Num(1), [sum(B[i,a]*x[a] for a in 1:4) for i in 1:3])
    ee = eg + sum(ge[a]*x[a] for a in 1:4)
    pp = pg + sum(gp[a]*x[a] for a in 1:4)
    rr = rg + sum(gr[a]*x[a] for a in 1:4)
    T = [(ee+pp)*u[a]*u[b]+pp*g[a,b] for a in 1:4,b in 1:4]
    D = at0.(stress_divergence_up(T, Γ, ∂))
    E = at0(energy_eq(Γ,u,pp,ee,∂))
    F = at0.(euler_up(g,Γ,u,pp,ee,∂))
    push!(C, certify_symbolic("off-shell timelike projection u.divT=-E", "C-LI-HEAT; G26 (3)", -D[1]+E; verbose))
    push!(C, certify_symbolic("off-shell spatial projection q.divT=Euler", "C-LI-HEAT; G26 (3)",
        vcat(F[1], [D[a]-F[a] for a in 2:4]); verbose))
    push!(C, certify_symbolic("off-shell baryon continuity", "C-LI-HEAT; G26 (3a)",
        at0(baryon_conservation(Γ,u,rr,∂))-gr[1]-rg*sum(B[a,a+1] for a in 1:3); verbose))
    @variables stress[1:4,1:4] dK[1:4,1:4]
    symmetric_stress = [stress[min(a,b),max(a,b)] for a in 1:4,b in 1:4]
    # Product rule: div(T.K)=K.divT+T^{ab}nabla_a K_b.
    # Symmetry converts the last term to T^{ab}(Lie_K g)_{ab}/2.
    push!(C, certify_symbolic("off-shell current work is half T.Lie_K_g", "C-LI-HEAT Killing balance",
        sum(symmetric_stress[a,b]*(dK[a,b]-(dK[a,b]+dK[b,a])/2)
            for a in 1:4,b in 1:4); verbose))

    # Derive the infinitesimal comoving volume from the spatial determinant.
    @variables lam s11 s22 s33 s12 s13 s23 dvol dr de volume0
    s = Num[s11 s12 s13; s12 s22 s23; s13 s23 s33]
    spatial = Matrix{Num}(I,3,3) + lam*s
    volume_variation = subs(Symbolics.derivative(det(spatial),lam),Dict(lam=>0))/2
    push!(C, certify_symbolic("comoving volume variation = q.Delta_g/2", "G26 (19)",
        volume_variation-(s11+s22+s33)/2; verbose))
    massvar = dr*volume0 + rg*dvol
    energyvar = de*volume0 + (eg+pg)*dvol
    push!(C, certify_symbolic("mass and energy work give relation (19)", "G26 (19)",
        [solve_linear(massvar,dr)/rg+dvol/volume0,
         solve_linear(energyvar,de)/(eg+pg)+dvol/volume0]; verbose))

    ex = li_balance_expressions()
    push!(C, certify_exact("chart jets divT=Euler+u*E off shell", "G26 (3)",
        vcat(bg.(ex.split),pert.(ex.split)); npts, seed, verbose))
    push!(C, certify_exact("off-shell Killing current metric-work identities", "C-LI-HEAT", 
        vcat(bg.(ex.current_identity),pert.(ex.current_identity)); npts, seed, verbose))
    push!(C, certify_symbolic("mass-energy difference is advected G26 (19)", "G26 (19)", ex.mass_energy; verbose))

    @variables zeta_r xi_r
    theta = p/eps
    residual = p*LocallyIsothermal.CLOSURE_RESIDUAL
    gauge = subs(residual, Dict(dp=>dp-zeta_r*p_r, deps=>deps-zeta_r*eps_r))-residual
    push!(C, certify_symbolic("gauge shift = -zeta^r eps Theta'", "C-LI-HEAT gauge x->x+zeta",
        gauge+zeta_r*eps*LocallyIsothermal.THETA_R; verbose))
    lagrangian = dp+xi_r*p_r-theta*(deps+xi_r*eps_r)
    push!(C, certify_symbolic("Delta p and Delta eps invariant under xi->xi+zeta", "C-LI-HEAT gauge",
        subs(lagrangian,Dict(dp=>dp-zeta_r*p_r,deps=>deps-zeta_r*eps_r,xi_r=>xi_r+zeta_r))-lagrangian; verbose))
    push!(C, certify_symbolic("LI displacement term matches selected residual", "C-LI-HEAT gauge",
        lagrangian-eps*LocallyIsothermal.THETA_R*xi_r-residual; verbose))
    @variables delta_theta
    full_eos = residual-eps*delta_theta
    push!(C, certify_symbolic("transforming delta Theta restores scalar EOS covariance", "C-LI-HEAT gauge",
        subs(full_eos,Dict(dp=>dp-zeta_r*p_r,deps=>deps-zeta_r*eps_r,
            delta_theta=>delta_theta-zeta_r*LocallyIsothermal.THETA_R))-full_eos; verbose))

    # Explicit q!=0 jets, independently fixed from Theta'/Theta=-q/r.
    # Solve the background hydrostatic identity for nu', not for p': preserve q.
    rng = MersenneTwister(seed)
    for q in (0//1, 3//4, 1//1)
        ok = true
        sensitive = true
        for j in 1:npts
            pt = rational_point(rng; complex_wd=true)
            radius = RatB(j+7,3)
            pt[p_r] = pt[p]*pt[eps_r]/pt[eps] - q*pt[p]/radius
            pt[nu_r] = exact_eval(first(values(ELIM_NU_R)), pt)
            ok &= all(iszero(exact_eval(a,pt)) for a in ex.on_solution)
            grad = exact_eval(LocallyIsothermal.THETA_R,pt)
            sensitive &= grad == -q*pt[p]/(pt[eps]*radius)
            q != 0 && (sensitive &= !iszero(grad))
            gpnt = merge(pt, Dict(zeta_r=>RatB(2,3),xi_r=>RatB(5,7)))
            expected = RatB(2,3)*q*pt[p]/radius
            sensitive &= exact_eval(gauge,gpnt) == expected
            # Dropping eps*Theta'*xi from the Lagrangian law leaves this
            # nonzero defect; q=0 is deliberately the blind control.
            defect = exact_eval(lagrangian-residual,gpnt)
            sensitive &= defect == -RatB(5,7)*q*pt[p]/radius
        end
        push!(C, certify_bool("OPTION A zero bath and G26 (19), q=$q", "C-LI-HEAT",
            ok, "energy, baryons, Lagrangian pressure, four divT components, E/J metric-work balances";
            verbose, modality=:exact_rational, npoints=npts))
        push!(C, certify_bool("Theta' and gauge/displacement negative controls, q=$q", "C-LI-HEAT R3/R9",
            sensitive, "specified gradient; exact nonzero defects for q=3/4,1; q=0 blind control";
            verbose, modality=:exact_rational, npoints=npts))
    end
    C
end
