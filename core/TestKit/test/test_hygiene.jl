# Unit tests for the §7.7 hygiene checkers (src/hygiene.jl) — a temp-dir fake Project.toml
# exercises the allow / forbid / must_declare / weakdeps_allow / legacy-[sources] PASS and FAIL
# paths, plus assert_not_loaded pass+fail.
#
# Fails-as-expected pattern: `@test`-based checkers cannot be probed with `@test_throws` (a red
# `@test` RECORDS a Fail, it does not throw), so each expected-red probe runs inside a DETACHED,
# never-finished `Test.DefaultTestSet` and we count the recorded Pass/Fail/Error results. The
# probe's failures therefore never redden this suite, while a checker that WRONGLY passes (or
# errors) is caught.

@testset "hygiene checkers (§7.7 gates 1/2/5)" begin

    # ---- probe harness -----------------------------------------------------------------------
    function _probe_counts_walk(ts::Test.DefaultTestSet)
        np, nf, ne = ts.n_passed, 0, 0
        for r in ts.results
            if r isa Test.DefaultTestSet
                p, f, e = _probe_counts_walk(r)
                np += p; nf += f; ne += e
            elseif r isa Test.Fail
                nf += 1
            elseif r isa Test.Error
                ne += 1
            end
        end
        return np, nf, ne
    end

    function _probe(f)
        ts = Test.DefaultTestSet("detached hygiene probe")
        Test.push_testset(ts)
        print_enable = Test.TESTSET_PRINT_ENABLE[]
        Test.TESTSET_PRINT_ENABLE[] = false   # expected-red probes: don't spam "Test Failed"
        try
            f()
        finally
            Test.TESTSET_PRINT_ENABLE[] = print_enable
            Test.pop_testset()
        end
        return _probe_counts_walk(ts)   # ts deliberately NOT finished — nothing propagates
    end

    _write_project(dir, body) = (p = joinpath(dir, "Project.toml"); write(p, body); p)

    # a representative fake package manifest surface
    fake_ok = """
    name = "FakePkg"
    uuid = "11111111-2222-3333-4444-555555555555"

    [deps]
    NumericsCore = "54f76139-0f7a-42ab-9084-b9cc78b29132"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

    [weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

    [sources]
    NumericsCore = {path = "../NumericsCore"}
    """

    fake_poisoned = """
    name = "PoisonPkg"
    uuid = "11111111-2222-3333-4444-666666666666"

    [deps]
    NumericsCore = "54f76139-0f7a-42ab-9084-b9cc78b29132"
    OrdinaryDiffEqVerner = "79d7bb75-1356-48c1-b8c0-6832512096c2"

    [weakdeps]
    HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"

    [sources]
    NumericsCore = {path = "../../packages/Teukolsky"}
    OrdinaryDiffEqVerner = {path = "../../ref-code/WaSABI"}
    """

    # ODE-poisoned deps but CLEAN [sources] — isolates the forbid-glob fail from the
    # always-on legacy-sources fail that fake_poisoned ALSO trips.
    fake_ode = """
    name = "OdePkg"
    uuid = "11111111-2222-3333-4444-777777777777"

    [deps]
    NumericsCore = "54f76139-0f7a-42ab-9084-b9cc78b29132"
    OrdinaryDiffEqVerner = "79d7bb75-1356-48c1-b8c0-6832512096c2"

    [sources]
    NumericsCore = {path = "../NumericsCore"}
    """

    @testset "_name_matches: exact + trailing-* glob" begin
        @test TestKit._name_matches("HDF5", "HDF5")
        @test !TestKit._name_matches("HDF5", "HDF5Ext")           # exact means exact
        @test TestKit._name_matches("OrdinaryDiffEq*", "OrdinaryDiffEqVerner")
        @test TestKit._name_matches("OrdinaryDiffEq*", "OrdinaryDiffEq")
        @test !TestKit._name_matches("OrdinaryDiffEq*", "DiffEqBase")
    end

    mktempdir() do dir
        ok_toml  = _write_project(mkpath(joinpath(dir, "ok")),  fake_ok)
        bad_toml = _write_project(mkpath(joinpath(dir, "bad")), fake_poisoned)
        ode_toml = _write_project(mkpath(joinpath(dir, "ode")), fake_ode)

        @testset "assert_dep_budget: full PASS path" begin
            np, nf, ne = _probe() do
                TestKit.assert_dep_budget(ok_toml;
                    allow=["NumericsCore", "LinearAlgebra"],
                    forbid=["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "WignerSymbols"],
                    must_declare=["NumericsCore"],
                    weakdeps_allow=["CUDA"])
            end
            @test nf == 0
            @test ne == 0
            @test np == 5   # allow + forbid + must_declare + weakdeps_allow + sources checks
        end

        @testset "assert_dep_budget: allow violation FAILS" begin
            np, nf, ne = _probe() do
                TestKit.assert_dep_budget(ok_toml; allow=["NumericsCore"])  # LinearAlgebra unsanctioned
            end
            @test nf == 1
            @test ne == 0
        end

        @testset "assert_dep_budget: forbid hits [deps] exactly" begin
            np, nf, ne = _probe() do
                TestKit.assert_dep_budget(ok_toml; forbid=["LinearAlgebra"])
            end
            @test nf == 1
        end

        @testset "assert_dep_budget: forbid glob hits [deps] AND forbid hits [weakdeps]" begin
            np, nf, ne = _probe() do
                # OrdinaryDiffEq* must catch OrdinaryDiffEqVerner in [deps]; CUDA sits in ok's
                # [weakdeps] — both poison axes land in the ONE forbidden-hits assert per call.
                TestKit.assert_dep_budget(ode_toml; forbid=["OrdinaryDiffEq*"])
                TestKit.assert_dep_budget(ok_toml;  forbid=["CUDA"])
            end
            @test nf == 2
        end

        @testset "assert_dep_budget: must_declare (gate 5) missing trigger FAILS" begin
            np, nf, ne = _probe() do
                TestKit.assert_dep_budget(ok_toml; must_declare=["HDF5"])
            end
            @test nf == 1
        end

        @testset "assert_dep_budget: weakdeps_allow violation FAILS" begin
            np, nf, ne = _probe() do
                TestKit.assert_dep_budget(ok_toml; weakdeps_allow=["AMDGPU"])  # CUDA unsanctioned
            end
            @test nf == 1
        end

        @testset "assert_dep_budget: legacy packages/ + ref-code/ [sources] paths FAIL" begin
            np, nf, ne = _probe() do
                TestKit.assert_dep_budget(bad_toml)   # both legacy paths land in the one sources assert
            end
            @test nf == 1
            @test ne == 0
        end
    end

    @testset "assert_not_loaded: pass + fail + glob" begin
        np, nf, ne = _probe() do
            TestKit.assert_not_loaded(["DefinitelyNotALoadedModule_9f2c1b"])
        end
        @test (np, nf, ne) == (1, 0, 0)

        np, nf, ne = _probe() do
            TestKit.assert_not_loaded(["Test"])       # Test IS loaded here — must fail
            TestKit.assert_not_loaded(["TestKi*"])    # glob catches loaded TestKit — must fail
        end
        @test nf == 2
        @test ne == 0
    end
end
