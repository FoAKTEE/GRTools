# Run with the Spacetime package project; --write replaces only TRAJ_EQ_REF.
using Spacetime
include(joinpath(@__DIR__, "..", "..", "equatorial_quadrature_oracle.jl"))
using .EquatorialQuadratureOracle
const CASES = (("0.5","10","0.2"),("0.9","8","0.4"),("0","12","0.3"),("-0.7","11","0.25"))
const PHASES = ("0","0.7","1.6","3","4.5","6")
function bake(bits,n)
    setprecision(BigFloat,bits) do
        [[quadrature_orbit(parse.(BigFloat,c)...;n).at_phase(parse(BigFloat,q)) for q in PHASES] for c in CASES]
    end
end
values=bake(256,128)
refined=bake(384,160)
worst=maximum(abs(x-y)/(1+abs(y)) for (case,refcase) in zip(values,refined)
              for (v,r) in zip(case,refcase) for (x,y) in zip(v,r))
@assert worst < big"1e-65"
header="""# TRAJ_EQ_REF: independent 256-bit Darwin-angle quadrature of the Kerr first integrals.
# Generator: test/fixtures/capture/regenerate_equatorial.jl in the Spacetime package.
# Reproduction command and method: spacetime_trajectory_PROVENANCE.md.
# Remaining TRAJ_{RAD,POL,FULL,PT}_REF: original baked KerrGeodesics fixtures, unchanged.
const TRAJ_EQ_REF = [
"""
io=IOBuffer();print(io,header)
for case in values
    println(io,"    [")
    for v in case
        println(io,"        ",map(Float64,v),",")
    end
    println(io,"    ],")
end
println(io,"]")
block=String(take!(io))
fixture=joinpath(@__DIR__,"..","spacetime_trajectory_ref.jl")
old=read(fixture,String)
rest=old[findfirst("const TRAJ_RAD_REF",old).start:end]
expected=block*rest
if "--write" in ARGS
    write(fixture,expected)
else
    @assert old==expected "Fixture is stale; run with --write"
end
println("Fixture quadrature convergence (256-bit N=128 versus 384-bit N=160): ",worst)
println("24 states / 120 scalars checked; equatorial fixture ", "--write" in ARGS ? "written" : "matches")
