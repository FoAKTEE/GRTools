using Test, FRDerivation, Symbolics
println("[START] C-LI-N; Julia=",VERSION," threads=",Threads.nthreads()); flush(stdout)
include(joinpath(@__DIR__,"../../derivation/src/R3Derivation.jl"))
include(joinpath(@__DIR__,"../src/li_newtonian.jl"))
const OUT=isempty(ARGS) ? error("supply evidence output directory") : ARGS[1]
const CERTS=LINewtonian.certify(outdir=OUT)
include(joinpath(@__DIR__,"../test/li_newtonian.jl"))
all(c.pass for c in CERTS) || error("C-LI-N failed")
println("[PASS] C-LI-N ",count(c->c.pass,CERTS),"/",length(CERTS)," certificates")
