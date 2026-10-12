# Reduced static/retrograde sector measurement; no mirrored code is executed.
# Run from the repository root under the AccretionDisc project:
#   julia -t1 --project=capabilities/AccretionDisc reproduce/arXiv-2407.20398/excluded_sector.jl out.toml [mmax]
module NewtonianExcludedSector
using AccretionDisc, Printf, TOML
using AccretionDisc.ForcingModes
const MO=AccretionDisc.Modes

function diagnostic_profiles(d,s,x,f,cache;regular)
    ur,up,sigma=MO.vars_from_enthalpy(d,f,x,s,cache;lindblad_regular=regular)
    pc=regular ? MO.RegularPotential(cache) : cache
    phi=[first(MO.potential_at(pc,r)) for r in x]
    # The high-level FR24 diagnostic requires corotation. Reuse its shared
    # single-mode torque contraction and the closure-independent Reynolds flux.
    F=[AccretionDisc.SpiralWaves.reynolds_flux(s;measure=r,enthalpy=MO.Sigma(d,r),
        dur=ur[i]/2,du_phi=r*up[i]/2) for (i,r) in enumerate(x)]
    density=MO.fr24_instantaneous_torque(x,[s],[sigma],[phi],zero(s.omega))
    c2=MO.csiso.(Ref(d),x).^2
    # From RegularRHS (Modes/LindbladCrossing.jl:21-39), for real w:
    # (F_J/c_s^2)' = [T' - L_T^-1 F_Phi]/c_s^2,
    # F_Phi = pi*m*r*Sigma*Re(u_r*conj(Phi))/w.
    # The forced term cannot be dropped where Phi is nonzero. Reuse the flux
    # contraction with the explicit potential part of the azimuthal velocity.
    Fphi=[AccretionDisc.SpiralWaves.reynolds_flux(s;measure=r,enthalpy=MO.Sigma(d,r),
        dur=ur[i]/2,du_phi=s.m*phi[i]/(s.omega-s.m*MO.Omega(d,r)+im*d.pole_disp)/2)
        for (i,r) in enumerate(x)]
    weighted_density=(density.-MO.LTinv.(Ref(d),x).*Fphi)./c2
    (;F,density,c2,weighted_density)
end

function diagnostics(d,s,x,f,cache;regular)
    b=diagnostic_profiles(d,s,x,f,cache;regular)
    F=b.F; c2=b.c2
    T=MO.fr24_torque_profile(x,b.density).torque
    weighted=MO.fr24_torque_profile(x,b.weighted_density).torque
    uncorrected=MO.fr24_torque_profile(x,b.density./c2).torque
    weighted_balance=abs((F[end]/c2[end]-F[1]/c2[1]-weighted[end]+weighted[1])/(weighted[end]-weighted[1]))
    uncorrected_weighted_balance=abs((F[end]/c2[end]-F[1]/c2[1]-uncorrected[end]+uncorrected[1])/(uncorrected[end]-uncorrected[1]))
    (;torque=last(T)-first(T),power=mode_power(s,last(T)-first(T)),fin=first(F),fout=last(F),weighted_balance,uncorrected_weighted_balance)
end

function probe(;n=4000,nt=1024)
    rows=Dict{String,Any}[]
    for (label,h,p,q,rin,e,soft,expected) in
        (("FR22",0.1,0.0,0.0,0.1,0.5,FixedSoftening(0.03),-0.004007638037317219),
         ("FR24",0.06,0.5,0.75,0.05,0.12,InstantaneousScaleHeight(0.3),-0.002025447934808053),
         ("FR24",0.06,0.5,0.75,0.05,0.5,InstantaneousScaleHeight(0.3),-0.0005367354759293588))
        d=NewtonianDisc(h0=h,p=p,q=q,beta=0.0,soft=0.3h,indirect=0)
        x=loggrid(d,rin,5.0,n); modes=[NewtonianMode(2,l,1.0) for l in (0,-1)]
        family=[NewtonianMode(2,l,1.0) for l in 1:(label=="FR22" ? 10 : 42)]
        bank=MO.potential_bank(d,[MO.EccentricOrbit(1.0,e,1.0)],loggrid(d,rin,5.0,n-3),vcat(modes,family);
                                softening=soft,Nt=nt,Nphi=1024)
        family_torque=zero(h)
        for s in family
            cache=MO.harmonic_cache(bank,1,s)
            f=MO.solve_lm_mode(d,s,x,cache;MO.FR22_SOLVER_SETTINGS...)
            family_torque+=diagnostics(d,s,x,f,cache;regular=true).torque*h^3
        end
        for s in modes, regular in (false,true)
            cache=MO.harmonic_cache(bank,1,s)
            f=MO.solve_lm_mode(d,s,x,cache;group_velocity_bc=true,lindblad_regular=regular,
                               phgrad_steps=0,reint_steps=2,rtol=1e-10,atol=1e-12,dense=false)
            b=diagnostics(d,s,x,f,cache;regular)
            # Hold the ODE solution fixed and coarsen only the torque quadrature.
            # This distinguishes a trapezoid error from a radiation imbalance.
            ids=unique!(vcat(collect(1:2:length(x)),length(x)))
            coarse=diagnostics(d,s,x[ids],f[ids,:],cache;regular)
            row=Dict{String,Any}("disc"=>label,"e"=>e,"m"=>s.m,"l"=>s.m+s.k,"regular"=>regular,
                "torque_FJ0"=>b.torque*h^3,"fin_FJ0"=>b.fin*h^3,"fout_FJ0"=>b.fout*h^3,
                "balance_relative"=>abs((b.fout-b.fin-b.torque)/b.torque),"weighted_balance"=>b.weighted_balance,
                "uncorrected_weighted_balance"=>b.uncorrected_weighted_balance,
                "weighted_balance_half_grid"=>coarse.weighted_balance,
                "paper_family_torque_FJ0"=>family_torque,"ratio_to_paper_family"=>b.torque*h^3/family_torque)
            s.m+s.k==0 && (row["reference_relative_error"]=abs((row["torque_FJ0"]-expected)/expected))
            push!(rows,row); println("PROBE ",row);flush(stdout)
        end
    end
    rows
end

function campaign(outpath;mmax=4,n=4000,nt=1024)
    rows=Dict{String,Any}[]; totals=Dict{String,Any}[]
    for policy in (:fr24,:fd25), e in (0.06,0.12,0.5)
        d=NewtonianDisc(h0=0.06,p=0.5,q=policy===:fr24 ? 0.75 : 0.0,beta=0.0,soft=0.018,indirect=0)
        soft=policy===:fr24 ? InstantaneousScaleHeight(0.3) : FixedSoftening(0.018)
        orb=MO.EccentricOrbit(1.0,e,1.0)
        modes=[NewtonianMode(m,l,1.0) for m in 1:mmax for l in m-40:m+40]
        # mmax=4 keeps every prograde corotation and OLR inside [0.05,5].
        x=loggrid(d,0.05,5.0,n)
        bank=MO.potential_bank(d,[orb],loggrid(d,0.05,5.0,n-3),modes;softening=soft,Nt=nt,Nphi=1024)
        cell=Dict{String,Any}[]
        for s in modes
            l=s.m+s.k; sector=s.m==1 && l==0 ? "secular" : l<=0 ? "excluded" : "paper"
            cache=MO.harmonic_cache(bank,1,s)
            f=MO.solve_lm_mode(d,s,x,cache;MO.FR22_SOLVER_SETTINGS...,group_velocity_bc=l<=0)
            b=diagnostics(d,s,x,f,cache;regular=true)
            row=Dict{String,Any}("policy"=>String(policy),"e"=>e,"m"=>s.m,"l"=>l,"sector"=>sector,
                "torque_FJ0"=>b.torque*d.h0^3,"power_FJ0_n"=>b.power*d.h0^3,
                "fin_FJ0"=>b.fin*d.h0^3,"fout_FJ0"=>b.fout*d.h0^3,"weighted_balance"=>b.weighted_balance)
            push!(cell,row);push!(rows,row)
            l==s.m+40 && (println("CELL ",policy," e=",e," m=",s.m," done");flush(stdout))
        end
        total=Dict{String,Any}("policy"=>String(policy),"e"=>e,"mmax"=>mmax,"kmax"=>40,"n"=>n,"Nt"=>nt)
        for key in ("torque_FJ0","power_FJ0_n")
            for sector in ("paper","excluded","secular")
                total[sector*"_"*key]=sum(r[key] for r in cell if r["sector"]==sector)
            end
            net=sum(r[key] for r in cell)
            total["all_"*key]=net
            nonsecular=total["paper_"*key]+total["excluded_"*key]
            total["nonsecular_total_"*key]=nonsecular
            total["excluded_share_nonsecular_"*key]=total["excluded_"*key]/nonsecular
            total["excluded_share_"*key]=total["excluded_"*key]/net
            total["excluded_plus_secular_share_"*key]=(total["excluded_"*key]+total["secular_"*key])/net
        end
        push!(totals,total);println("TOTAL ",total);flush(stdout)
        open(outpath,"w") do io
            TOML.print(io,Dict("cells"=>totals,"modes"=>rows,
                "scope"=>"Reduced m<=$mmax, |k|<=40; secular (1,0) reported separately and boundary conditional; no m=0 power."))
        end
    end
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    if length(ARGS)>1 && ARGS[2]=="probe"
        open(ARGS[1],"w") do io
            NewtonianExcludedSector.TOML.print(io,Dict("probes"=>NewtonianExcludedSector.probe()))
        end
    else
        NewtonianExcludedSector.campaign(ARGS[1];mmax=length(ARGS)>1 ? parse(Int,ARGS[2]) : 4,
            n=length(ARGS)>2 ? parse(Int,ARGS[3]) : 4000,nt=length(ARGS)>3 ? parse(Int,ARGS[4]) : 1024)
    end
end
