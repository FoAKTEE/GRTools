# Narrow FR24 anchor checks; no production cube or figure run.
# julia -t1 --project=capabilities/AccretionDisc reproduce/arXiv-2407.20398/check_anchor_interfaces.jl
root=normpath(joinpath(@__DIR__,"..",".."))
for file in ("test_fr24_corotation_radius.jl","test_fr24_torque_diagnostics.jl")
    include(joinpath(root,"capabilities","AccretionDisc","test","Modes",file))
end
