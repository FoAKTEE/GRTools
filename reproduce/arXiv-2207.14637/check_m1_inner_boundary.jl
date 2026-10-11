# Characterise the strict inner BC residual without changing a gate or solver.
include("linear_solver_audit.jl")
g=GROUPS[2];d=disc(g);mode=FM.NewtonianMode(1,3,1.)
b=deserialize(joinpath(OUT,"bank-2.jls"));c=subcache(MO.harmonic_cache(b,1,mode),2)
x=loggrid(d,.1,5.,4000)
_,off,_=measured(d,x,mode,c;rtol=1e-13,dense=false,FR_OPTIONS...,phgrad_steps=0)
cases=Tuple{String,typeof(x),NamedTuple}[("phase_off",x,off)]
for (name,file) in (("reference","mode-2-1-3-1.0e-13.jls"),
                    ("grid_double","refine-2-1-3-7999-2.jls"),
                    ("cache_double","refine-2-1-3-4000-1.jls"))
    z=deserialize(joinpath(OUT,file));push!(cases,(name,z.x,z.s))
end
rows=Dict{String,Any}[]
for (name,x,s) in cases
    i=s.indices[1]
    residual=abs(s.mode[i,2]-s.C_after[1]*s.mode[i,1])/maximum(abs,s.mode[:,2])
    push!(rows,Dict("case"=>name,"scaled_inner_bc_residual"=>residual,
        "passes_1e8_gate"=>residual<1e-8,"inner_radius"=>x[i]))
end
save("m1-inner-boundary",Dict("rows"=>rows))
display(rows)
