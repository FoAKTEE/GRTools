using Test
include(joinpath(@__DIR__,"../../../../reproduce/arXiv-2407.20398/fr24_cube.jl"))

# The authors' v1 cube (14 .out files, no licence) is never vendored. It is read from
# FR24_CUBE_ROOT, else from the arXiv ancillary mirror of an emri-disk checkout
# (EMRI_DISK_ROOT, else the sibling of this GRTools tree). Absent data is a recorded skip.
cube_root() = get(ENV,"FR24_CUBE_ROOT") do
    disk = get(ENV,"EMRI_DISK_ROOT",normpath(joinpath(@__DIR__,"..","..","..","..","..","emri-disk")))
    joinpath(disk,"ref-paper","arXiv-2407.20398","src","anc")
end

@testset "FR24 v1 data cube algebra" begin
    root = cube_root()
    if !isfile(joinpath(root,"eccentric","t_net.out"))
        @info "FR24 v1 cube not found; skipped" root
        @test_skip "FR24 v1 cube mirror present"
    else
        cube = FR24Cube.read_cube(root)
        @test size(cube.tables["eccentric/t_net"].value) == (4,4,3,13)
        @test length(cube.hashes) == 14
        # Anchor the axis order against the serialized text independently of reshape.
        raw = parse.(Float64,split(read(joinpath(root,"eccentric/t_net.out"),String)))
        @test cube.tables["eccentric/t_net"].value[2,3,2,5] == raw[((2-1)*4*3+(3-1)*3+1)*13+5]
        result = FR24Cube.check_cube(cube)
        # HOLE H-FR24-CUBE-PRECISION: two identities fail at their unchanged preregistered
        # bounds (evidence/fr24-cube-algebra-gate/PREREG.json); the failure counts are
        # measurements of the fixed v1 data, not acceptance thresholds.
        hole = ("angular_momentum","timescale_relation")
        for (name,check) in result["checks"]
            name in hole ? (@test_broken check["passed"]) : (@test check["passed"])
        end
        @test !result["passed"]
        @test result["checks"]["angular_momentum"]["failure_count"] == 158
        @test result["checks"]["timescale_relation"]["failure_count"] == 24
        @test result["checks"]["timescale_relation"]["count"] == 624
        @test result["cross_run"]["net"]["relative"]["count"] == 48
        @test result["cross_run"]["net"]["worst_relative_qph"] == [0.0,0.0,0.06]
        # A real sign/normalization error must be rejected by the gate.
        cube.tables["eccentric/tau_L_inv"].value[1,1,1,2] += 1e-3
        bad = FR24Cube.check_cube(cube)
        @test !bad["passed"]
        for name in hole
            @test bad["checks"][name]["failure_count"] == result["checks"][name]["failure_count"]+1
        end
    end
end
