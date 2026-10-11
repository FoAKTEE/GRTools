# Run with --project=capabilities/AccretionDisc.
using TOML  # CHANDRA_OUTPUT capture; not an AccretionDisc test dependency (lander, lane w17-g43)
include(joinpath(@__DIR__,"../../capabilities/AccretionDisc/test/Modes/test_fr24_mode_set.jl"))
