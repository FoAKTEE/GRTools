println("[START] general adiabatic emission"); flush(stdout)
include(joinpath(@__DIR__,"..","src","R3Derivation.jl"))
const R3 = R3Derivation
const HERE = normpath(joinpath(@__DIR__,".."))
const OUT = isempty(ARGS) ? joinpath(HERE,"out","general-emission-check") : ARGS[1]
const NPTS = length(ARGS)>1 ? parse(Int,ARGS[2]) : 12
const ART = joinpath(HERE,"out","ECC_adiabatic.jls.txt")
const TARGET = joinpath(HERE,"..","src","SpiralWaves","generated","eccentric_adiabatic.jl")
ex = R3.general_response_exports()
R3.write_exports(ART,ex;header="# ECC closure=adiabatic; general complex amplitudes")
R3.emit_eccentric(TARGET,ART;tag="adiabatic")
certs = R3.emission_certificates(TARGET,ART,ex;tag="adiabatic",npts=NPTS,accounting_dir=OUT)
used = union((R3._names(e) for (_,e) in ex)...)
push!(certs,R3.certify_bool("adiabatic odd components and independent conjugates","general metric",
    all(x in used for x in ("htr","hrphi","htrc","hrphic")),"htr, hrphi, htrc, hrphic present in the emitted graph"))
all(c.pass for c in certs) || error("general emission certificate failed")
println("[PASS] general adiabatic emission: $(length(certs))/$(length(certs))")
flush(stdout)
