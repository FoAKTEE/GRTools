# Usage: julia -t1 check_cube_algebra.jl DATA_ROOT OUTPUT.toml
using TOML
include("fr24_cube.jl")
length(ARGS)==2 || error("usage: check_cube_algebra.jl DATA_ROOT OUTPUT.toml")
result = FR24Cube.check_cube(FR24Cube.read_cube(ARGS[1]))
open(ARGS[2],"w") do io
    TOML.print(io,result;sorted=true)
end
for (name,check) in sort(collect(result["checks"]);by=first)
    println(name, ": ",check)
end
println("circular cross-run distributions: ",result["cross_run"])
result["passed"] || error("FR24 v1 cube algebra failed; see output")
