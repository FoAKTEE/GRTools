# Characterization tests for PostNewtonian.Catalog (PN .m discovery / browse layer).
#
# GROUNDING = the filesystem oracle + the legacy wolframscript selector oracle (baked literals from
# packages/PostNewtonianSelfForce/test/catalog_selector_test.jl, where the upstream FileNames /
# StringDelete / Fold / SelectFromPNSeriesList were run via wolframscript against this same bundled
# SeriesData tree, and `find ... -name '*.m' | wc -l == 91` independently confirmed). Pure stdlib —
# no symbolic oracle.

using Test
using PostNewtonian.Catalog

@testset "Catalog: discovery" begin
    cat = pn_series_catalog()
    # Oracle: 91 *.m files (find census + wolframscript ALLCOUNT=91).
    @test length(cat) == 91
    @test allunique(cat)
    @test "/Schwarzschild/Circular/Flux/Infinity-l2m2" in cat
    @test "/Kerr/Spherical/Orbit/Energy" in cat
    # Key shape: leading '/', no '.m', forward slashes.
    @test all(startswith(k, "/") for k in cat)
    @test all(!endswith(k, ".m") for k in cat)
    @test !any(occursin('\\', k) for k in cat)
end

@testset "Catalog: select_pn_series single filter" begin
    cat = pn_series_catalog()
    @test select_pn_series(cat, "") == cat                         # empty pattern is identity
    @test "/Schwarzschild/Circular/Local/Redshift" in select_pn_series(cat, "redshift")  # case-insensitive
    @test all(!occursin("Spin", s) for s in select_pn_series(cat, "!Spin"))              # negation
end

@testset "Catalog: pn_browse vs wolframscript oracle" begin
    # Oracle (c): Sort[Fold[SelectFromPNSeriesList, allSeries, query]].
    q_local = [
        "/Schwarzschild/Circular/Local/Redshift",
        "/Schwarzschild/Circular/Local/Spin",
        "/Schwarzschild/Circular/Local/Spin/dEdt",
        "/Schwarzschild/Circular/Local/Tidal-lambdaB",
        "/Schwarzschild/Circular/Local/Tidal-lambdaE1",
        "/Schwarzschild/Circular/Local/Tidal-lambdaE2",
        "/Schwarzschild/Circular/Local/Tidal-lambdaE3",
    ]
    @test sort(pn_browse(["Schwarzschild", "Circular", "Local"])) == q_local

    q_local_nospin = [
        "/Schwarzschild/Circular/Local/Redshift",
        "/Schwarzschild/Circular/Local/Tidal-lambdaB",
        "/Schwarzschild/Circular/Local/Tidal-lambdaE1",
        "/Schwarzschild/Circular/Local/Tidal-lambdaE2",
        "/Schwarzschild/Circular/Local/Tidal-lambdaE3",
    ]
    @test sort(pn_browse(["Schwarzschild", "Circular", "Local", "!Spin"])) == q_local_nospin

    q_l2m2 = [
        "/Schwarzschild/Circular/Flux/Infinity-l2m2",
        "/Schwarzschild/Circular/Flux/Spin/Infinity-l2m2",
    ]
    @test sort(pn_browse(["Schwarzschild", "Circular", "Flux", "Infinity-l2m2"])) == q_l2m2

    @test length(pn_browse("")) == 91                              # empty query = whole catalog
    @test sort(pn_browse("Local")) == sort(pn_browse(["Local"]))   # single string == 1-vector
end
