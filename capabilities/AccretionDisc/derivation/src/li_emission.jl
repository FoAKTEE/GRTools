"OPTION A C/S and reconstruction, consuming the existing LI reduction."
li_response_exports() = LocallyIsothermal.closure_exports(LocallyIsothermal)

"""
OPTION A stress jets from the perfect-fluid TUP_JET, and linear coordinate
currents. These are inputs to diagnostics, not a quadratic wave-AMF formula.
No adiabatic U/V/W tensor or adiabatic AMF is used in this construction.
Compression energy remains in eps; this is the prescribed energy-density
pressure law, not a finite-GR thermostat. delta Theta=0 is gauge-fixed in BL.
"""
function li_stress_exports()
    LI = LocallyIsothermal
    reduced(e) = LI.onshell_nu(to_omega_d(LI.lin_subs(LI.eliminate(e),LI.MASTER_SUB)))
    ex = Pair{String,Any}[]
    labels = ("t","r","phi","z")
    for a in 1:4, b in a:4
        n = labels[a]*labels[b]
        t0 = bg(TUP_JET[a,b])
        dt = reduced(pert(TUP_JET[a,b]))
        iszero_entry(t0) || push!(ex,"T0_$n"=>t0)
        if !iszero_entry(dt)
            push!(ex,"dT_$n"=>dt)
            push!(ex,"dTbar_$n"=>cbar(dt))
        end
    end
    for (axis,label) in ((IT,"t"),(IPHI,"phi")), mu in 1:4
        current = sum(TUP_JET[mu,nu]*GJET[nu,axis] for nu in 1:4)
        dt = reduced(pert(current))
        iszero_entry(dt) || push!(ex,"dJ$(label)_$(labels[mu])"=>dt)
    end
    ex
end

"""
Execute LI generated solutions against the original fluid equations and an
independent product-rule stress expansion. Local jets have Theta'/Theta=-q/r
at q=0,3/4,1; these are constrained local jets, not radial-profile solves.
"""
function li_emitted_physics(path,stresspath; npts=12,seed=260818,accounting_dir,verbose=true)
    C = Certificate[]
    mod = Module(gensym(:LIEmission))
    Base.include(mod,path); Base.include(mod,stresspath)
    equations = emission_tree.([DEULER_UP[IR],DENERGY,DEULER_UP[IPHI],
        DNORM,LocallyIsothermal.CLOSURE_RESIDUAL,DEULER_UP[IT]])
    rng = MersenneTwister(seed)
    labels = ("t","r","phi","z")
    stressnames = Set(first.(li_stress_exports()))
    mkpath(accounting_dir)
    open(joinpath(accounting_dir,"li-qtemp-residuals.txt"),"w") do io
        println(io,"Local jets Theta'/Theta=-q/r; q=0,3/4,1; N=$npts per q; seed=$seed; exact complex rationals.")
        println(io,"Sampler is nonuniform; no statistical identity probability is claimed. No points rejected.")
        for q in (0//1,3//4,1//1)
            solution,stress,negative = true,true,true
            for j in 1:npts
                pt = complex_emission_point(rng)
                r = RatB(j+7,3)
                pt[p_r] = pt[p]*(pt[eps_r]/pt[eps]-q/r)
                pt[nu_r] = exact_eval(ELIM_NU_R[nu_r],pt)
                inp,tbl = emission_inputs(pt)
                vals = Dict(n=>emitted_value(mod,"li",n,inp) for n in
                    (MASTER_NAMES...,"dut","duphi_up","deps","duphi_low"))
                yp0 = vals["C00"]*pt[y0]+vals["C01"]*pt[y1]+vals["S0"]
                yp1 = vals["C10"]*pt[y0]+vals["C11"]*pt[y1]+vals["S1"]
                tbl[:dp]=(pt[eps]+pt[p])*pt[y0]; tbl[:dur]=im*pt[y1]
                tbl[:dp_r]=(pt[eps_r]+pt[p_r])*pt[y0]+(pt[eps]+pt[p])*yp0
                tbl[:dur_r]=im*yp1
                tbl[:dut]=vals["dut"]; tbl[:duphi]=vals["duphi_up"]; tbl[:deps]=vals["deps"]
                residuals = [eval_tree(e,tbl) for e in equations]
                xi = im*tbl[:dur]/(pt[W]/pt[enu]*pt[omega_d])
                theta = pt[p]/pt[eps]
                blind = tbl[:dp]+xi*pt[p_r]-theta*(tbl[:deps]+xi*pt[eps_r])
                lagrangian = blind+pt[eps]*q*theta*xi/r
                solution &= all(iszero,residuals) && iszero(lagrangian)
                negative &= q==0 ? iszero(blind) : !iszero(blind)
                println(io,"q=$q point=$j fluid=$residuals pressure=$lagrangian omitted-gradient-defect=$blind")

                # Direct product rule, independent of the serialized stress expressions.
                g = map(e->eval_tree(emission_tree(e),tbl),METRIC)
                gi = map(e->eval_tree(emission_tree(e),tbl),METRIC_INV)
                h = map(e->eval_tree(emission_tree(e),tbl),HMODE)
                u = map(e->eval_tree(emission_tree(e),tbl),UUP)
                du = [tbl[:dut],tbl[:dur],tbl[:duphi],RatB(0)]
                t0 = (pt[eps]+pt[p])*u*transpose(u)+pt[p]*gi
                dt = (tbl[:deps]+tbl[:dp])*u*transpose(u)+
                    (pt[eps]+pt[p])*(du*transpose(u)+u*transpose(du))+
                    tbl[:dp]*gi-pt[p]*gi*h*gi
                dj = dt*g+t0*h
                for a in 1:4, b in a:4
                    for (prefix,want) in (("T0",t0[a,b]),("dT",dt[a,b]))
                        name = "$(prefix)_$(labels[a])$(labels[b])"
                        if name in stressnames
                            defect = emitted_value(mod,"li_stress",name,inp)-want
                            stress &= iszero(defect)
                            println(io,"q=$q point=$j $name residual=$defect")
                        else
                            stress &= iszero(want)
                        end
                    end
                end
                for (axis,label) in ((IT,"t"),(IPHI,"phi")), mu in 1:4
                    name = "dJ$(label)_$(labels[mu])"
                    if name in stressnames
                        defect = emitted_value(mod,"li_stress",name,inp)-dj[mu,axis]
                        stress &= iszero(defect)
                        println(io,"q=$q point=$j $name residual=$defect")
                    else
                        stress &= iszero(dj[mu,axis])
                    end
                end
            end
            for (name,ok,note) in (("fluid solution",solution,"Euler r/phi/t, energy, normalization, pressure and Lagrangian temperature term"),
                ("stress jets",stress,"independent perfect-fluid product rule and mixed coordinate currents"),
                ("temperature-gradient negative control",negative,"omitted eps*Theta'*xi defect nonzero for q=3/4,1; q=0 blind control"))
                push!(C,certify_bool("LI emitted $name q=$q","OPTION A",ok,note;
                    verbose,modality=:exact_rational,npoints=npts))
                verbose && flush(stdout)
            end
        end
    end
    C
end
