# Selected-family radial refinement of the mixed-grid full weak-field sum.
# This is not a full-sum radial or exterior-boundary convergence certificate.
include("excluded_sector.jl")
using TOML, Serialization
const GE=GRExcludedSector
baseline=TOML.parsefile(ARGS[1])
bank=deserialize(ARGS[2])
selected=(44,100,175)
reused=filter(r->!(r["m"] in selected),baseline["modes"])
refined=GE.solve_inventory(0.0,10000.0,0.12,bank.modes;Nq=256,n=16001,q=.75,soft=.018,
    source_bank=bank,grid_strategy=:core_log,reuse_rows=reused,core_nodes=513)
old=Dict((r["m"],r["k"])=>r for r in baseline["modes"])
rows=Dict{String,Any}[]
for m in selected
    family=filter(r->r["m"]==m,refined)
    length(family)==81 || error("incomplete refinement family")
    row=Dict{String,Any}("m"=>m,"modes"=>length(family))
    for quantity in ("torque","power")
        before=sum(old[(r["m"],r["k"])][quantity] for r in family)
        after=sum(r[quantity] for r in family)
        absolute_change=sum(abs(r[quantity]-old[(r["m"],r["k"])][quantity]) for r in family)
        row[quantity*"_coarse"]=before
        row[quantity*"_refined"]=after
        row[quantity*"_family_relative_change"]=abs((after-before)/after)
        row[quantity*"_sum_absolute_modal_change_over_full_nonsecular_total"]=
            absolute_change/abs(baseline["totals"]["nonsecular_total_"*quantity])
    end
    push!(rows,row)
    println("TAIL_REFINEMENT ",row)
end
open(io->TOML.print(io,Dict("families"=>rows,"base_n"=>8001,"refined_n"=>16001,
    "base_core_nodes"=>257,"refined_core_nodes"=>513,
    "scope"=>"Only three selected families are refined; this does not bound the full tail error.")),ARGS[3],"w")
