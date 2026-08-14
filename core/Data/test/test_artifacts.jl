# Module: artifacts — override -> declared+provisioned artifact -> loud-fallback resolution chain
# (design §8.2 / §10). GROUNDING: Pkg.Artifacts closed-form round-trip — bind a freshly-`create_artifact`ed
# tree to a TEMP Artifacts.toml, point the resolver at it, and assert tree-hash identity / precedence /
# loud fallback. Pure path + Artifacts-stdlib logic, no external download, no heavy dep.
#   * provisioned-artifact identity: artifact_path(name) == Artifacts.artifact_path(tree_hash).
#   * override precedence: set_datadir!(tmp) beats the bound+provisioned artifact.
#   * lazy-declared-but-unprovisioned + no override: loud ErrorException (datapath fallback).
#   * set_datadir! override of an UNDECLARED name resolves under the override.
#   * stage_datadir!/ensure_local offline pre-staging hook (§8.2).

using Test
using Data
import Pkg
import Artifacts

@testset "artifacts" begin
    # Save/restore the process-global env + datadir + artifacts-toml so this testset is hermetic.
    _saved_env = get(ENV, "GRTOOLS_DATADIR", nothing)
    try
        @testset "Pkg.Artifacts tree-hash resolution (override->artifact->error chain)" begin
            # Hermetic: clear datadir override + env so the artifact branch is reachable.
            Data.set_datadir!(nothing)
            haskey(ENV, "GRTOOLS_DATADIR") && delete!(ENV, "GRTOOLS_DATADIR")

            # (1) Provisioned-artifact identity. Create a content tree, bind it to a TEMP Artifacts.toml,
            #     point the resolver there ⇒ artifact_path resolves via the Artifacts-stdlib tree-hash.
            toml = joinpath(mktempdir(), "Artifacts.toml")
            h = Pkg.Artifacts.create_artifact() do dir
                write(joinpath(dir, "flux.bin"), "kerr-ecc-eq-flux")
            end
            Pkg.Artifacts.bind_artifact!(toml, "kerr_ecc_eq_flux", h; lazy = true)
            Data.set_artifacts_toml!(toml)
            @test Data._artifact_declared("kerr_ecc_eq_flux")
            @test Data.artifact_path("kerr_ecc_eq_flux") == Artifacts.artifact_path(h)
            @test isfile(joinpath(Data.artifact_path("kerr_ecc_eq_flux"), "flux.bin"))

            # (2) Override precedence: set_datadir!(tmp) beats the bound+provisioned artifact.
            ovr = mktempdir()
            Data.set_datadir!(ovr)
            @test Data.artifact_path("kerr_ecc_eq_flux") ==
                  joinpath(abspath(ovr), "kerr_ecc_eq_flux")
            Data.set_datadir!(nothing)

            # (3) Lazy-declared-but-unprovisioned + no override ⇒ loud ErrorException (datapath fallback).
            toml2 = joinpath(mktempdir(), "Artifacts.toml")
            fake = Base.SHA1("0"^40)         # a tree-hash that is not provisioned in the depot
            Pkg.Artifacts.bind_artifact!(toml2, "absent_corpus", fake; lazy = true)
            Data.set_artifacts_toml!(toml2)
            @test Data._artifact_declared("absent_corpus")
            @test !Artifacts.artifact_exists(fake)
            @test_throws ErrorException Data.artifact_path("absent_corpus")

            Data.set_artifacts_toml!(nothing)
        end

        @testset "shipped Artifacts.toml has a real bound kerr_ecc_eq_flux entry" begin
            # The package ships a REAL bound entry (provisioned locally from the cached .fewdata corpus).
            Data.set_datadir!(nothing)
            haskey(ENV, "GRTOOLS_DATADIR") && delete!(ENV, "GRTOOLS_DATADIR")
            @test Data._artifact_declared("kerr_ecc_eq_flux")
            h = Artifacts.artifact_hash("kerr_ecc_eq_flux", Data.artifacts_toml())
            @test h !== nothing
            # When provisioned in this depot it resolves to the tree-hash root (runtime @artifact_str path).
            if Artifacts.artifact_exists(h)
                @test Data.artifact_path("kerr_ecc_eq_flux") == Artifacts.artifact_path(h)
            end
        end

        @testset "set_datadir! override + artifact_path fallback (undeclared)" begin
            tmp = mktempdir()
            Data.set_datadir!(tmp)
            @test Data.datadir() == abspath(tmp)
            # No Artifacts.toml entry for this name ⇒ the override root wins (precedence).
            @test Data.artifact_path("foo.h5") == joinpath(abspath(tmp), "foo.h5")
            Data.set_datadir!(nothing)
        end

        @testset "set_datadir!(nothing) clears to default_datadir (loud, no ref-code escape)" begin
            Data.set_datadir!(nothing)
            haskey(ENV, "GRTOOLS_DATADIR") && delete!(ENV, "GRTOOLS_DATADIR")
            # No override, no env ⇒ fails loudly (deliberately NO implicit ref-code fallback).
            @test_throws ErrorException Data.datadir()
            # The env var is the configured default root.
            tmp = mktempdir()
            ENV["GRTOOLS_DATADIR"] = tmp
            @test Data.datadir() == abspath(tmp)
            @test Data.artifact_path("bar.dat") == joinpath(abspath(tmp), "bar.dat")
            delete!(ENV, "GRTOOLS_DATADIR")
        end

        @testset "stage_datadir! + ensure_local pre-staging (§8.2 offline/HPC hook)" begin
            # A corpus source directory holding a fixture (stands in for the 4.74GB ZNAmps grid).
            src = mktempdir()
            write(joinpath(src, "corpus.bin"), "payload")
            root = mktempdir()
            out = Data.stage_datadir!(src, root)
            @test out == abspath(root)
            @test Data.datadir() == abspath(root)        # staging set the active root
            # ensure_local resolves the staged fixture under the local root.
            p = Data.ensure_local("corpus.bin")
            @test isfile(p)
            @test read(p, String) == "payload"
            # A corpus absent from the local root throws loudly (no silent escape / download).
            @test_throws ErrorException Data.ensure_local("absent.bin")
            Data.set_datadir!(nothing)
        end

        @testset "stage_datadir! single-file source" begin
            f = joinpath(mktempdir(), "single.bin")
            write(f, "x")
            root = mktempdir()
            Data.stage_datadir!(f, root)
            @test isfile(Data.ensure_local("single.bin"))
            Data.set_datadir!(nothing)
        end
    finally
        Data.set_datadir!(nothing)
        Data.set_artifacts_toml!(nothing)
        if _saved_env === nothing
            haskey(ENV, "GRTOOLS_DATADIR") && delete!(ENV, "GRTOOLS_DATADIR")
        else
            ENV["GRTOOLS_DATADIR"] = _saved_env
        end
    end
end
