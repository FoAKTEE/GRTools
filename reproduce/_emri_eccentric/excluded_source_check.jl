# Pointwise angular-quadrature audit of a retained radial bank. This is not a
# full torque quadrature-convergence certificate.
using AccretionDisc, Serialization, TOML
using AccretionDisc.ForcingModes
const SW=AccretionDisc.SpiralWaves
bank=deserialize(ARGS[1])
p=bank.perturber; o=p.orbit; a=o.p/(1-o.e^2)
d=SW.DiscModel(;χ=o.a,r_ref=a,h_ref=0.06,pSigma=0.5,qTemp=0.75,
               closure=LocallyIsothermalOptionA,evolution=:fixed)
r=[a*x for x in (0.5,1-o.e,1.0,1+o.e,2.0)]
ref=SW.source_modes(SW.source_azimuth_bank(d,p,r;Nq=512,mmax=maximum(s.m for s in bank.modes)),bank.modes)
old=SW.source_modes(bank,bank.modes,r)
rows=Dict{String,Any}[]
for jet in keys(ref[1])
    scale=maximum(abs(x[jet]) for x in ref)
    change=maximum(abs(x[jet]-y[jet]) for (x,y) in zip(ref,old))
    push!(rows,Dict("jet"=>String(jet),"scale"=>scale,"absolute_change"=>change,
                    "relative_to_component_peak"=>scale==0 ? change : change/scale))
end
open(io->TOML.print(io,Dict("rows"=>rows,"Nq_reference"=>512,"radii"=>r)),ARGS[2],"w")
println("SOURCE_CHECK ",rows)
maximum(row["relative_to_component_peak"] for row in rows)<1e-6 || error("pointwise angular audit failed")
println("PASS pointwise Nq=256/512 audit, all retained modes and source jets")
