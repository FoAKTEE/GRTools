# oracle_module.jl — the external-CAS ORACLES of the arXiv:2608.17003 master system, as one module.
#
# `master_generated.jl` is the transcription of the authors' supplementary Mathematica notebook
# (route R1) and `r2_generated.jl` the independent Mathematica derivation (route R2).  Neither is in
# the build path any more (EGD::jd-swap): the build source is the Julia derivation of route R3.  They
# are loaded here, into their own module, so that a test can compare the two answers side by side
# without a name clash with the live `AccretionDisc.SpiralWaves._master_entries`.
#
# Include it defensively — several test files use it and they must not redefine the module:
#     isdefined(@__MODULE__, :SpiralWavesOracle) ||
#         include(joinpath(@__DIR__, "fixtures", "oracle", "oracle_module.jl"))
module SpiralWavesOracle
    using AccretionDisc.SpiralWaves: MasterInputs
    include(joinpath(@__DIR__, "master_generated.jl"))
    include(joinpath(@__DIR__, "r2_generated.jl"))
end
