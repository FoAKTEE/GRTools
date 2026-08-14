# Grounded tests for SelfForce.LorenzField.construct_static_monopole (analytic l=m=0 static monopole).
#
# GROUNDING — one-shot characterization vs the legacy H1Lorenz oracle + analytic structure:
#
# (1) CHARACTERIZATION: the EXPECTED inhom[k,part,grididx] values below were captured ONE-SHOT from
#     the legacy `H1Lorenz.construct_static_monopole` (a HEAVY HDF5 package — never a live test-dep)
#     on the orbit r0=10 with grid [4,6,8,10,12,14,16] (r0 at 0-based index 3), via
#     scratchpad/capture_mono.jl (Pkg.activate packages/H1Lorenz). The closed forms are exact
#     algebraic/log expressions, so the port must reproduce them to ~1e-13. The h1-l0m0.h5 corpus
#     fixture is SYNTHETIC seed data (read_h1 structure test) and is deliberately NOT used here.
# (2) STRUCTURE: nf+gauge = 4 fields × 4 parts × (gridsize+1); the gauge slot (k=4) and all imaginary
#     parts (part 2,4) are identically zero; the worldline right-limit lives at the last index.
# (3) ANALYTIC: the static monopole is real and asymptotically flat; h1 In/Out agree at r0 (field
#     continuity) while the radial derivative h1' jumps across the worldline (source term).

using Test
using SelfForce.LorenzField: construct_static_monopole, OrbitalParams

const _GRID_M = [4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]
const _ORBIT_M = OrbitalParams(10.0, _GRID_M, 1, 3, 5)

# Baked legacy H1Lorenz.construct_static_monopole(orbit) nonzero entries (k,part,grididx) => value.
const _MONO_EXPECT = Dict(
    (1,1,1)=>1.1783038054635417, (1,1,2)=>3.724021903687242, (1,1,3)=>6.6279589057324211,
    (1,1,4)=>9.6526647743573317, (1,1,5)=>9.9573642944782961, (1,1,6)=>10.195072325496207,
    (1,1,7)=>10.383490316714965, (1,1,8)=>9.6526647743573335,
    (1,3,1)=>1.080111821674913, (1,3,2)=>1.3965082138827154, (1,3,3)=>1.4912907537897946,
    (1,3,4)=>1.5283385892732442, (1,3,5)=>0.13393574427365573, (1,3,6)=>0.10525871404849818,
    (1,3,7)=>0.084230289459470642, (1,3,8)=>0.17250655642384194,
    (2,1,1)=>4.218342020079378, (2,1,2)=>7.985866534104793, (2,1,3)=>11.186059586240354,
    (2,1,4)=>14.222162923244692, (2,1,5)=>13.877511870914597, (2,1,6)=>13.613476709204377,
    (2,1,7)=>13.407057083504574, (2,1,8)=>14.222162923244696,
    (2,3,1)=>2.2328893104833867, (2,3,2)=>1.6801048091548103, (2,3,3)=>1.5455454239629871,
    (2,3,4)=>1.4976277358741363, (2,3,5)=>-0.14994133407555044, (2,3,6)=>-0.11601112807024383,
    (2,3,7)=>-0.091728414714037809, (2,3,8)=>-0.19716230518761657,
    (3,1,1)=>-1.0761985388432678, (3,1,2)=>0.044055695720823893, (3,1,3)=>0.59697846839506175,
    (3,1,4)=>0.98581152593808075, (3,1,5)=>0.89964800445330517, (3,1,6)=>0.81207789649170181,
    (3,1,7)=>0.73382731771230447, (3,1,8)=>0.98581152593808041,
    (3,3,1)=>0.90925417075272474, (3,3,2)=>0.35646966942414937, (3,3,3)=>0.22191028423232542,
    (3,3,4)=>0.17399259614347465, (3,3,5)=>-0.045076960870765681, (3,3,6)=>-0.041769233597031191,
    (3,3,7)=>-0.036410254545178787, (3,3,8)=>-0.037856158989244326,
)

@testset "construct_static_monopole vs baked legacy H1Lorenz oracle (rtol 1e-13)" begin
    inhom = construct_static_monopole(_ORBIT_M)
    @test size(inhom) == (4, 4, length(_GRID_M) + 1)
    for k in 1:4, p in 1:4, j in 1:size(inhom, 3)
        want = get(_MONO_EXPECT, (k, p, j), 0.0)
        if want == 0.0
            @test inhom[k, p, j] == 0.0       # gauge slot k=4, imaginary parts (p=2,4) ≡ 0
        else
            @test inhom[k, p, j] ≈ want rtol = 1e-13
        end
    end
end

@testset "construct_static_monopole structure (real, gauge=0, r0 right-limit)" begin
    inhom = construct_static_monopole(_ORBIT_M)
    @test all(inhom[4, :, :] .== 0)            # gauge field R2 ≡ 0
    @test all(inhom[:, 2, :] .== 0) && all(inhom[:, 4, :] .== 0)   # purely real
    # field continuity at r0: In value at r0 (grid idx 4) ≈ Out right-limit (last idx)
    nlast = size(inhom, 3)
    @test inhom[1, 1, 4] ≈ inhom[1, 1, nlast] rtol = 1e-13
    # derivative jumps across the worldline (source): In vs Out r0' differ
    @test !isapprox(inhom[1, 3, 4], inhom[1, 3, nlast]; rtol = 1e-6)
end

@testset "construct_static_monopole BigFloat precision path (type-generic {T})" begin
    gridb = BigFloat.(_GRID_M)
    orbitb = OrbitalParams(big(10.0), gridb, 1, 3, 5)
    inhom = construct_static_monopole(orbitb)
    @test eltype(inhom) == BigFloat
    @test Float64(inhom[1, 1, 1]) ≈ _MONO_EXPECT[(1, 1, 1)] rtol = 1e-13
    @test Float64(inhom[2, 1, 3]) ≈ _MONO_EXPECT[(2, 1, 3)] rtol = 1e-13
end
