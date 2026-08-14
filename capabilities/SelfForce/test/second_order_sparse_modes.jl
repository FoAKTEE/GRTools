# Sparse-corpus safety at the SelfForce.SecondOrder assembly boundary.
#
# An angularly allowed parent mode that is wholly absent from the first-order corpus is the zero
# field, so every bilinear R2 contribution involving it must be an exact radial zero.  This gate
# exercises both supported H1 layouts and checks the assembled low target against the only two
# physically present parent pairs.

using Test
using SelfForce

const _SPARSE_SO = SelfForce.SecondOrder
const _SPARSE_LAYOUTS = (
    (name = "production four-row", dir = joinpath(@__DIR__, "fixtures", "h1lorenz")),
    (name = "legacy six-row", dir = joinpath(@__DIR__, "fixtures", "h1")),
)

function _with_sparse_l2m2(f, source_dir)
    mktempdir() do dir
        cp(joinpath(source_dir, "h1-l2m2.h5"), joinpath(dir, "h1-l2m2.h5"))
        return f(_SPARSE_SO.read_h1(dir; l_max = 2))
    end
end

@testset "SecondOrder sparse whole-mode zero contract" begin
    @testset "$(layout.name)" for layout in _SPARSE_LAYOUTS
        _with_sparse_l2m2(layout.dir) do fields
            assembled = _SPARSE_SO.second_order_ricci(
                fields, 1, 0, 0; h1lmax = 2, dlmax = 2)

            # For target m=0, the sparse corpus contains only the reflected
            # (2,-2)x(2,2) and (2,2)x(2,-2) parent pairs.  Every other
            # angularly allowed pair is an exact absent-mode zero.
            explicit = zeros(ComplexF64, length(fields.r))
            for m1 in (-2, 2)
                explicit .+= _SPARSE_SO.r2_term(
                    1, 1.0, Float64(fields.r0), Float64.(fields.r),
                    Float64.(fields.f), Float64.(fields.fp),
                    fields.h, fields.dh, fields.ddh,
                    fields.h, fields.dh, fields.ddh,
                    0, 0, 2, m1, 2, -m1)
            end

            @test all(isfinite, assembled)
            @test assembled == explicit
            @test maximum(abs, assembled) > 0

            control = _SPARSE_SO.second_order_ricci(
                fields, 1, 4, 4; h1lmax = 2, dlmax = 2)
            @test all(isfinite, control)

            empty_mode = empty(fields.h)
            @test_throws ErrorException _SPARSE_SO.r2_term(
                0, 1.0, Float64(fields.r0), Float64.(fields.r),
                Float64.(fields.f), Float64.(fields.fp),
                empty_mode, empty_mode, empty_mode,
                fields.h, fields.dh, fields.ddh,
                0, 0, 2, 0, 2, 0)
        end
    end
end
