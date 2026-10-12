# Preflight the full paper/physical seed before any expensive source projection.
using AccretionDisc, Spacetime, TOML, LinearAlgebra
using AccretionDisc.ForcingModes
const SW=AccretionDisc.SpiralWaves
rows=Dict{String,Any}[]
for (p,outer) in ((1000.0,5.0),(10000.0,5.0),(1000.0,20.0))
    orbit=EquatorialOrbit(0.0,p,0.12)
    aref=p/(1-orbit.e^2)
    d=SW.DiscModel(;χ=0.0,r_ref=aref,h_ref=.06,pSigma=.5,qTemp=.75,
                   closure=LocallyIsothermalOptionA,evolution=:fixed,εp=1e-4)
    mask=SW.eccentric_mode_mask(orbit;policy=:fr24,mask=:physical)
    failures=Dict{String,Any}[]
    for s in mask.physical
        c=SW.ModeContext(d,s;closure=d.closure,domain=(.05aref,outer*aref))
        for (r,outward) in ((first(c.domain),false),(last(c.domain),true))
            try
                b=SW.li_boundary(c,r;outward)
                all(isfinite,b.row) || error("nonfinite boundary row")
            catch err
                C=SW.master_operator(c,r;landau=false)
                λ1,λ2,_,_=SW.eigen2(C)
                push!(failures,Dict("m"=>s.m,"k"=>s.k,"l"=>s.m+s.k,
                    "omega"=>s.omega,"omega_d"=>SW.omega_d(c,r),"r"=>r,
                    "outer"=>outward,"lambda_real"=>real.([λ1,λ2]),
                    "lambda_imag"=>imag.([λ1,λ2]),"error"=>sprint(showerror,err)))
            end
        end
    end
    row=Dict("p"=>p,"outer_over_aref"=>outer,"modes"=>length(mask.physical),
             "edges"=>2length(mask.physical),"failures"=>failures)
    push!(rows,row)
    println("BOUNDARY_PREFLIGHT p=",p," outer=",outer," failures=",length(failures))
    foreach(println,failures)
    flush(stdout)
    open(io->TOML.print(io,Dict("cases"=>rows)),ARGS[1],"w")
end
