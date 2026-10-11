using Test
println("[START] thermostatted sensitivity completion"); flush(stdout)
include(joinpath(@__DIR__,"../src/R3Derivation.jl"))
const R3=R3Derivation
include(joinpath(@__DIR__,"../test/test_li_thermostat.jl"))
println("[PASS] thermostat symbolic completion; 15 symbolic certificates; no primary closure/release edge")
