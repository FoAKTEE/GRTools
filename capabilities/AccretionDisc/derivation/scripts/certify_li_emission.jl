println("[START] OPTION A LI emission"); flush(stdout)
include(joinpath(@__DIR__,"..","src","R3Derivation.jl"))
const R3 = R3Derivation
const HERE = normpath(joinpath(@__DIR__,".."))
const OUT = isempty(ARGS) ? joinpath(HERE,"out","li-emission-check") : ARGS[1]
const NPTS = length(ARGS)>1 ? parse(Int,ARGS[2]) : 12
certs = R3.Certificate[]
paths = String[]
for (tag,ex) in (("li",R3.li_response_exports()),("li_stress",R3.li_stress_exports()))
    art = joinpath(HERE,"out","ECC_$(tag).jls.txt")
    path = joinpath(HERE,"..","src","SpiralWaves","generated","eccentric_$(tag).jl")
    R3.write_exports(art,ex;header="# ECC closure=$tag; general complex amplitudes")
    R3.emit_eccentric(path,art;tag)
    append!(certs,R3.emission_certificates(path,art,ex;tag,npts=NPTS,seed=260818,accounting_dir=OUT))
    used = union((R3._names(e) for (_,e) in ex)...)
    push!(certs,R3.certify_bool("$tag closure and complex source content","OPTION A",
        !("Gamma1" in used) && all(x in used for x in ("htr","hrphi")),
        "no Gamma1; htr/hrphi retained; graph consumes OPTION A reduction"))
    push!(paths,path)
end
append!(certs,R3.li_emitted_physics(paths...;npts=NPTS,accounting_dir=OUT))
all(c.pass for c in certs) || error("LI emission certificate failed")
println("[PASS] OPTION A LI emission: $(length(certs))/$(length(certs))")
