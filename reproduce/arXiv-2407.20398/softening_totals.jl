# Derived observables only; reuse the landed data-only cube reader and keep banks external.
include("../arXiv-2207.14637/softening_bank.jl")
include("fr24_cube.jl")
using DelimitedFiles
root,out,cube_root=ARGS
mkpath(out)
cube=FR24Cube.read_cube(cube_root)
index=Tuple(findfirst(x -> isapprox(x,t;rtol=0,atol=1e-14),g)
            for (g,t) in zip(cube.grids["eccentric"],(0.,0.,.06,.12)))
oracle=cube.tables["eccentric/t_net"].value[index...]
ladder=[]
for rule in ("fixed","moving")
    s=setup("fr24",rule,.12,40)
    rows=Any[]
    for mode in s.modes
        dir=joinpath(root,"fr24-$rule-0.12-s$(mod(mode.m,2))")
        path=joinpath(dir,"response-$(mode.m)-$(mode.m+mode.k).jls")
        isfile(path) || error("missing required checkpoint: $path")
        y=deserialize(path)
        length(y.r)==4000 || error("incomplete published radial grid")
        push!(rows,[mode.m,mode.m+mode.k,y.torque,y.cpu_s,y.rout,y.grid_nodes,y.cache_nodes])
    end
    open(joinpath(out,"$rule-modes.csv"),"w") do io
        println(io,"m,l,torque_FJ0,cpu_s,rout,grid_nodes,cache_nodes")
        writedlm(io,permutedims(hcat(rows...)),',')
    end
    last=nothing
    for mmax in (10,20,30,40)
        chosen=filter(row -> row[1]<=mmax,rows)
        torque=sum(row[3] for row in chosen)
        push!(ladder,(rule,mmax,40,length(chosen),torque,torque-oracle,
            (torque-oracle)/abs(oracle),last===nothing ? "" : torque-last,
            sum(abs(row[3]) for row in chosen),sum(row[4] for row in chosen)/3600))
        last=torque
    end
    println("FR24 $rule: $(length(rows)) modes, $(sum(r[4] for r in rows)/3600) solve core-hours")
end
open(joinpath(out,"torque_truncation.csv"),"w") do io
    println(io,"rule,mmax,offset,modes,computed_torque_FJ0,residual_to_cube,relative_residual_to_cube,last_shell_increment,sum_abs_modal,solve_core_h")
    for row in ladder;println(io,join(row,','));end
end
open(joinpath(out,"cube-sources.toml"),"w") do io
    TOML.print(io,Dict("root"=>cube_root,"sha256"=>cube.hashes,
        "reader"=>"reproduce/arXiv-2407.20398/fr24_cube.jl:13"))
end
