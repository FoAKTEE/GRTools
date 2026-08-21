# oracle_module.jl — the external-CAS ORACLE of the arXiv:2608.17003 master system, as one module.
#
# Route R2 = OUR OWN independent Mathematica derivation (EGD::deriv-master-system-R2). It is NOT in
# the build path (EGD::jd-swap): the build source is the Julia derivation of route R3,
# `src/SpiralWaves/{master,r2}_derived.jl`. The two files here are loaded into their own module so a
# test can compare the two routes side by side without a name clash with the live
# `AccretionDisc.SpiralWaves` kernels:
#     `r2_master_generated.jl`  → `_r2_master_entries(inp, omt, vp)`   — C, S of eq (47)
#     `r2_generated.jl`         → `_r2_delta_u_phi`, `_r2_amf`, `_r2_torque_density` — (72a), (72b)
#
# The notebook-derived oracle of route R1 (`master_generated.jl`) was REMOVED from this repository on
# 2026-08-22 by owner decision (EGD::jd-notebook-purge). It is not needed to keep the gate strong:
# EGD::R1-vs-R2-gate certified R1 == R2 on-shell (`Simplify[R1|(11) − R2] = 0`, and all 300 entries of
# the 50-point fixture byte-identical), and the purge re-sourced every R1 comparison onto R2 at
# unchanged tolerance and unchanged measured margin. The R1 lineage and the full R1-vs-R2-vs-R3
# equivalence evidence are archived in the mission repository (emri-disk results/emri-gasdisk
# paper_EMRI-GASDISK, evidence/{R1-vs-R2-gate,jd-equivalence-gate,jd-notebook-purge}).
#
# Include it defensively — a test file must not redefine the module:
#     isdefined(@__MODULE__, :SpiralWavesOracle) ||
#         include(joinpath(@__DIR__, "fixtures", "oracle", "oracle_module.jl"))
module SpiralWavesOracle
    using AccretionDisc.SpiralWaves: MasterInputs
    include(joinpath(@__DIR__, "r2_master_generated.jl"))
    include(joinpath(@__DIR__, "r2_generated.jl"))
end
