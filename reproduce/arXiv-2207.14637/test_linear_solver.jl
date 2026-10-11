# Float64 lane checks. The full package suite remains the lander's registration gate.
const TESTDIR=normpath(joinpath(@__DIR__,"..","..","capabilities","AccretionDisc","test","Modes"))
include(joinpath(TESTDIR,"test_fr22_phase_gradient.jl"))
include(joinpath(TESTDIR,"test_fr22_linear_solver.jl"))
include(joinpath(TESTDIR,"test_newtonian_planet.jl"))
