using Test
using Perturbation

@testset "Perturbation" begin
    @testset "smoke" begin
        @test isdefined(Main, :Perturbation)
        @test Perturbation isa Module
    end

    @testset "firewall closure-load gate" begin
        # Design §4.3: Perturbation is the LIGHT-DEP FIREWALL — its runtime closure must carry
        # NONE of the poison axes. This asserts the package was NOT loaded with any of
        # HDF5/JLD2/Symbolics/OrdinaryDiffEq/CUDA/AMDGPU/KernelAbstractions in its dependency
        # closure (they would appear as loaded modules only if a transitive runtime dep pulled
        # them). External-ODE-free / HDF5-free is a structural invariant of this capability.
        poison = ("HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq", "OrdinaryDiffEqVerner",
                  "CUDA", "AMDGPU", "KernelAbstractions")
        loaded = Set(string(m) for m in values(Base.loaded_modules))
        for p in poison
            @test !(p in loaded)
        end
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (TOML-static dep budget + closure-load),
                            # on top of the (untouched) firewall testset above; before heavy tests

    # EXPLICIT test-file list (was `readdir` discovery: an untracked file in a dirty working
    # tree silently entered the suite, and a renamed/deleted file silently left it). Keep the
    # list ALPHABETICAL — that is the order `sort(readdir(...))` used, so the sequence in which
    # the files run is unchanged. Adding a test file means adding it here.
    TEST_FILES_PERT = [
        "test_amplitudes.jl", "test_cache_threadsafety.jl",
        "test_fluxes.jl",
        "test_negative_spin_constants.jl",
        "test_pointparticle.jl",
        "test_published_tables_2303_16798.jl",   # arXiv:2303.16798 Tables I/II/III (slow: ~2.5 min at -t4)
        "test_qnm_fast.jl",
        "test_qnm_frequency.jl", "test_qnm_hyperboloidal.jl",
        "test_qqnm_angular.jl", "test_qqnm_solve.jl",
        "test_qqnm_source.jl", "test_qqnm_source_tables.jl",
        "test_radial_carrier.jl", "test_radial_gsn.jl",
        "test_radial_heunc.jl", "test_radial_hypergeom.jl",
        "test_radial_jet.jl", "test_radial_mst.jl",
        "test_radial_radial.jl", "test_radial_renorm_angmom.jl",
        "test_radial_sn.jl", "test_radial_static.jl",
        "test_radial_static_precision.jl", "test_radial_teuk_deriv.jl",
        "test_scattering.jl", "test_spinfluxes_circular.jl",
        "test_spinfluxes_eccentric.jl", "test_spinfluxes_modes.jl",
        "test_td_hphc.jl", "test_td_hphc_qnm.jl",
        "test_td_hphc_recon.jl", "test_td_hphc_scd.jl",
        "test_td_hphc_scd_solve.jl", "test_timedomain.jl",
        "test_timedomain_spin.jl", "test_timedomain_tailfits.jl",
        "test_timedomain_spinsource.jl",
        "test_timedomain_spinsource_kerr.jl",
    ]

    @testset "test-file registry is complete and exact" begin
        listed = Set(TEST_FILES_PERT)
        ondisk = Set(filter(n -> startswith(n, "test_") && endswith(n, ".jl"), readdir(@__DIR__)))
        # every listed file exists ...
        @test isempty(setdiff(listed, ondisk))
        # ... and no test_*.jl in this directory is unregistered (a new test file must be added
        # to the list above; a stray scratch file must not be left in test/).
        @test isempty(setdiff(ondisk, listed))
    end

    for f in TEST_FILES_PERT
        include(joinpath(@__DIR__, f))
    end
end
