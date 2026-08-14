using Test
using SelfForce

@testset "SelfForce" begin
    @testset "smoke" begin
        # Package + submodule namespaces load (LorenzField.*/SecondOrder.*/Puncture.*).
        @test SelfForce isa Module
        @test SelfForce.LorenzField isa Module
        @test SelfForce.SecondOrder isa Module
        @test SelfForce.Puncture isa Module
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module grounded test files (tests-first, characterization vs captured legacy fixtures):
    include("second_order_coupling.jl") # coupling front door === Wigner.coupling + 3j selection zeros
    include("second_order_modes.jl")   # a_il re-export + selection_modes vs Ricci.cc L86-99
    include("second_order_readcorpus.jl") # read_h1 corpus reader + ZeroPadField vs fixtures + analytic f,fp
    include("second_order_h1schema_bridge.jl") # production 4-row H1Schema -> typed h/dh/ddh radial jet
    include("second_order_sparse_modes.jl") # absent whole (l,m) parents are exact-zero R2 contributions
    include("second_order_assembly.jl") # ricci_source vs captured g++/GSL src_self.h5 oracle + write_src round-trip
    include("second_order_plan_schema.jl") # Q0 typed modes/bandwidth/corpus+plan fingerprint invariants
    include("second_order_analytic.jl") # h1_mass/h1_spin/d2hdr2 vs baked g++ oracle fixture + field-eq self-consistency
    include("lorenz_orbit.jl")          # orbit constants (Spacetime) + R1..R10 BLS sources (sources.c)
    include("lorenz_static.jl")         # construct_static_monopole vs baked legacy H1Lorenz oracle
    include("lorenz_boundary.jl")       # in/out BC series vs baked legacy + ODE-residual vs LorenzRadial RHS
    include("lorenz_coupledsets.jl")    # CoupledSet/CSETS wiring + select_coupled_set + γ↔(l,m) enumeration
    include("lorenz_solver.jl")         # VoP/shoot integrator vs genuine legacy h1 fixtures + round-trip + static identity
    include("puncture_static.jl")       # h0_static vs baked Punctures CSV + elliptic collapse + christoffel identity
    include("puncture_alpha_integration.jl") # exactint/modeint/alphaint(+v1-4,t1-2) vs captured alphaint_reference.csv + f21_terminating collapse
    include("puncture_regsum.jl")       # reg_sum/noreg_sum vs captured regsum_reference.csv + digamma_big collapse + analytic tail
    include("puncture_fieldeq.jl")      # ddhret vs captured ddhret_reference.csv (10 components) + Complex{BigFloat} path + even-static predicate
    include("puncture_tensor_basis.jl") # tensor_harmonic/coord_component/coord_deriv vs captured tensor_basis+coord refs + ForwardDiff θ-deriv FD cross-check + reality/parity
    include("puncture_modes.jl")        # hP/dhP vs captured legacy {hP,dhP}_reference.csv (QuadGK + AD) + central-FD/ForwardDiff dhP identity
    include("puncture_mode_mixing.jl")  # h_mode vs captured h_reference.csv + Wigner.wignerD_fixed collapse vs wignerD_reference.csv + reality identity
    include("puncture_singular_field.jl") # hS_mode === h_mode mapped + zero-mode short-circuit + worldline_index/split identity
    include("puncture_subtraction.jl")  # hret_block enum (Data) + read_hret_mode round-trip vs Data.write_mode_h5 + hR=hret−hS identity
    include("puncture_worldline.jl")    # covariant_first/second vs captured {DhR1,D2hR1}_ref + assemble vs {hR1Data,dhR1Data}_ref + ∇g=0 + a=0 background
end
