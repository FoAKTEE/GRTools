using Test
using MetricReconstruction

@testset "MetricReconstruction" begin
    @testset "smoke" begin
        # Package + submodule namespaces load (Radiation.* / Lorenz.* disambiguation present).
        @test MetricReconstruction isa Module
        @test MetricReconstruction.Radiation isa Module
        @test MetricReconstruction.Lorenz isa Module
        @test MetricReconstruction.Lorenz.LorenzRadial isa Module
        @test MetricReconstruction.GHZ isa Module
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module grounded test files (tests-first, characterization vs captured legacy fixtures):
    include("lorenz_operators.jl")     # geometry/tetrad/jet algebra/CG mixing, vs Spacetime/Wigner + KLC
    include("lorenz_background.jl")     # circular-equatorial orbit scalars, vs Bardeen closed forms + KLC
    include("lorenz_jumps.jl")          # App-H jumps_pm/hj1, vs published App-H tables + Harmonics SWSH
    include("lorenz_closure.jl")        # closure_kerr/closure_a0, vs published App-H LOWER tables + regression
    include("lorenz_kerr_radial.jl")   # a≠0 radial reconstruction (reconstruct_kerr), analytic self-consistency G-K10
    include("lorenz_schwarzschild_radial.jl")  # a=0 radial reconstruction + a→0 limit vs reconstruct_kerr(0)
    include("lorenz_schwarzschild_components.jl")  # a=0 tetrad modes → BLS h̄^i, vs legacy KLC to_BLS fixture
    include("lorenz_kerr_smatrices.jl")  # a≠0 closure S-matrices + lplp/mpmp components, vs legacy KLC + a→0/FD
    include("lorenz_cache.jl")            # bounded identity/precision-aware band-cache invariants
    include("lorenz_kerr_components.jl")  # a≠0 tetrad_modes_kerr/to_BLS_kerr, vs legacy KLC fixture + a→0 regression
    include("lorenz_weyl.jl")          # gauge-invariant ψ4 (weyl_psi4_lorenz/_amp), θ-factorization + Teukolsky identity
    include("lorenz_radial.jl")        # even/odd/monopole/odd_dipole _ode! RHS, vs H1Lorenz/DiscModes
    include("radiation_operators.jl")  # SepField/RTJet jet algebra + Dn/Ln ladders (analytic)
    include("radiation_modes.jl")      # teuk_lambda/radial_jet/angular_jet vs GrTool Harmonics/Perturbation
    include("radiation_tetrad.jl")     # six Hertz/CCK H-components, scalar+jet, vs legacy MR fixture
    include("radiation_weyl.jl")       # calD/calC/Dhat + CCK presets vs ExampleUsage cell-21 oracle
    include("radiation_plus_minus.jl") # ±m tetrad-component value assembly, vs legacy MR component_value fixture
    include("radiation_projection.jl") # tetrad→coordinate projection h_coordinate, vs legacy MR + Kinnersley duality
    include("radiation_einstein.jl")   # linearized-Einstein vacuum residual on reconstructed h (G-M4) + delegation
    include("radiation_weyl_scalars.jl") # ψ0/ψ4 regenerated from reconstructed h (H10 gate) vs Teukolsky + legacy MR
    include("radiation_assembly.jl")   # CCK mode-sum h_coord_modesum + einstein_grid vs legacy MR + vacuum gate
    include("radiation_api.jl")        # user-facing metric_tetrad/h_tetrad/h_coord + reconstruct/psi vs legacy MR

    # Transport submodule (Li–Yunes traceful ORG, 2605.11080): the DECISIVE Move-E gates.
    include("transport_shell_example.jl")   # Schwarzschild static thin shell (G-E1..G-E6), F64 + BigFloat-256
    include("transport_vacuum_crossgate.jl")# vacuum reduction vs RLGP-F/D4 operator algebra (G-E-vac)
    include("transport_smooth_source.jl")   # F4: generic SMOOTH static source via DP5 (G-F4a..G-F4e)
    include("transport_dipole.jl")           # F6-B: static ℓ=1 axial angular-momentum completion (G-F6B-*)

    # F7-3: convention-independent symbolic gate on the Ma–Yang ++ 2nd-order source coefficients
    # (quotient R''=FR'+GR double-entry, r-form==Γ-form named entries, full-A_2-row reduction).
    include("f7_3_source_coeff.jl")          # F7-3: G-F7-3a..d (Dn ladders + RTJet + NumericsCore jets)

    include("ghz_corrector_head.jl")         # sourced GHZ IRG Schwarzschild corrector head
    include("ghz_corrector_middle.jl")       # second equation in the sourced GHZ IRG hierarchy
    include("ghz_corrector_pair.jl")         # explicit physical partner for barred middle data
    include("ghz_corrector_tail.jl")         # bounded finite-ray x_nn tail and tetrad data
    include("ghz_corrector_monopole.jl")     # bounded finite-ray l=0 head-only x_nn tail
    include("ghz_ray_moment_cache.jl")       # forward O(n) bounded GHZ ray cache
    include("ghz_adjoint_source.jl")         # supplied adjusted source to bounded eta rays
    include("ghz_hertz_inversion.jl")        # positive-quarter-thorn4 IRG Hertz radial inversion
    include("ghz_hertz_metric.jl")           # finite-frequency outgoing-u GHZ IRG Hertz metric injection
end
