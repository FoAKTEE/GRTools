# 4th-order uniform-grid finite differences + 6th-order Kreiss-Oliger dissipation, in-place over a
# 2D (nx, ny) column-swept field. Port of ref-code/TeukEvolution.jl/src/Radial.jl:27-162.
# Interior stencils == the fd1c/fd2c closures (identity asserted in test/test_fd_uniform.jl);
# exact-rational literals so {T}/BigFloat is inherited. Characterized vs the original in
# test/fixtures/numericscore_fd_uniform_ref.jl.
#
# Exact polynomial degrees (verified in test_fd_uniform.jl, matching upstream 4th-order claims):
#   set_d1! interior & all four edge rows : exact for degree ≤ 4  (5-point stencils)
#   set_d2! interior & all four edge rows : exact for degree ≤ 5  (interior 5-point via symmetry;
#                                                                   edges 6-point one-sided)
#   ko_filter! dissipation term           : annihilates degree ≤ 5 interior (6th central difference)

"""
    set_d1!(dv, v, dr)

Radial first derivative of the `(nx, ny)` field `v` with a 4th-order uniform-grid finite-difference
stencil, written in place into `dv` (swept column-by-column over the `ny` angular columns). Interior
rows (`i = 3:nx-2`) use the 5-point central stencil whose weights are IDENTICAL to [`fd1c`](@ref)
(`(-1/12, 2/3, 0, -2/3, 1/12)/dr`); the four boundary rows use fixed 4th-order one-sided closures
(no ghost zones). Type-generic: element type flows from `v`/`dr`, so `Complex{BigFloat}` fields give
`Complex{BigFloat}` derivatives. Port of `TeukEvolution.jl/src/Radial.jl:27-71` (exact rationals).
"""
function set_d1!(dv::AbstractArray{T}, v::AbstractArray{T}, dr::Real) where {T<:Number}
    inv_dr = 1 / dr
    nx, ny = size(dv)
    for j in 1:ny
        for i in 3:(nx - 2)
            dv[i, j] =
                (
                    -(1 // 12) * v[i+2, j] + (2 // 3) * v[i+1, j] -
                    (2 // 3) * v[i-1, j] + (1 // 12) * v[i-2, j]
                ) * inv_dr
        end
        dv[1, j] =
            (
                -3 * v[5, j] + 16 * v[4, j] - 36 * v[3, j] + 48 * v[2, j] -
                25 * v[1, j]
            ) * (1 // 12) * inv_dr
        dv[2, j] =
            (
                +1 * v[5, j] - 6 * v[4, j] + 18 * v[3, j] - 10 * v[2, j] -
                3 * v[1, j]
            ) * (1 // 12) * inv_dr
        dv[nx, j] =
            -(
                -3 * v[nx-4, j] + 16 * v[nx-3, j] - 36 * v[nx-2, j] +
                48 * v[nx-1, j] - 25 * v[nx, j]
            ) * (1 // 12) * inv_dr
        dv[nx-1, j] =
            -(
                +1 * v[nx-4, j] - 6 * v[nx-3, j] + 18 * v[nx-2, j] -
                10 * v[nx-1, j] - 3 * v[nx, j]
            ) * (1 // 12) * inv_dr
    end
    return nothing
end

"""
    set_d2!(dv, v, dr)

Radial second derivative of the `(nx, ny)` field `v` with a 4th-order uniform-grid finite-difference
stencil, in place into `dv` (column-swept). Interior rows (`i = 3:nx-2`) use the 5-point central
stencil whose weights are IDENTICAL to [`fd2c`](@ref) (`(-1/12, 4/3, -5/2, 4/3, -1/12)/dr^2`); the
four boundary rows use fixed 4th-order 6-point one-sided closures. Type-generic (BigFloat-clean).
Port of `TeukEvolution.jl/src/Radial.jl:82-127` (exact rationals).
"""
function set_d2!(dv::AbstractArray{T}, v::AbstractArray{T}, dr::Real) where {T<:Number}
    inv_dr2 = 1 / dr^2
    nx, ny = size(dv)
    for j in 1:ny
        for i in 3:(nx - 2)
            dv[i, j] =
                (
                    -1 * v[i+2, j] + 16 * v[i+1, j] - 30 * v[i, j] +
                    16 * v[i-1, j] - 1 * v[i-2, j]
                ) * (1 // 12) * inv_dr2
        end
        dv[1, j] =
            (
                -10 * v[6, j] + 61 * v[5, j] - 156 * v[4, j] + 214 * v[3, j] -
                154 * v[2, j] + 45 * v[1, j]
            ) * (1 // 12) * inv_dr2
        dv[2, j] =
            (
                +1 * v[6, j] - 6 * v[5, j] + 14 * v[4, j] - 4 * v[3, j] -
                15 * v[2, j] + 10 * v[1, j]
            ) * (1 // 12) * inv_dr2
        dv[nx, j] =
            (
                -10 * v[nx-5, j] + 61 * v[nx-4, j] - 156 * v[nx-3, j] +
                214 * v[nx-2, j] - 154 * v[nx-1, j] + 45 * v[nx, j]
            ) * (1 // 12) * inv_dr2
        dv[nx-1, j] =
            (
                +1 * v[nx-5, j] - 6 * v[nx-4, j] + 14 * v[nx-3, j] -
                4 * v[nx-2, j] - 15 * v[nx-1, j] + 10 * v[nx, j]
            ) * (1 // 12) * inv_dr2
    end
    return nothing
end

"""
    KO6_KERNEL

The 6th-order Kreiss-Oliger dissipation kernel `(1, -6, 15, -20, 15, -6, 1)` — the 6th central
difference `Σ (-1)^k binom(6,k)`. Sums to zero (annihilates constants); applied as
`v ← v + (ε/64) · KO6[v]` in [`ko_filter!`](@ref), scaled by `1/64 = 2^-6`.
"""
const KO6_KERNEL = (1, -6, 15, -20, 15, -6, 1)

"""
    ko_filter!(v, tmp, eps_KO)

6th-order Kreiss-Oliger dissipation applied in place to the `(nx, ny)` field `v`, using `tmp`
(same shape) as scratch for the pre-filter copy. Interior rows `i = 4:nx-3` are updated to
`v[i] + (eps_KO/64) · (v[i-3] - 6v[i-2] + 15v[i-1] - 20v[i] + 15v[i+1] - 6v[i+2] + v[i+3])`; the
three boundary rows on each end are left untouched. The [`KO6_KERNEL`](@ref) dissipation term is a
6th central difference, so it annihilates polynomials of degree ≤ 5 (genuine high-frequency
dissipation). Type-generic (BigFloat-clean). Port of `TeukEvolution.jl/src/Radial.jl:138-162`.
"""
function ko_filter!(v::AbstractArray{T}, tmp::AbstractArray{T}, eps_KO::Real) where {T<:Number}
    nx, ny = size(v)
    for j in 1:ny
        for i in 1:nx
            tmp[i, j] = v[i, j]
        end
    end
    for j in 1:ny
        for i in 4:(nx - 3)
            v[i, j] = (
                (eps_KO / 64) * tmp[i-3, j] +
                (-6 * eps_KO / 64) * tmp[i-2, j] +
                (15 * eps_KO / 64) * tmp[i-1, j] +
                (1 - (20 * eps_KO / 64)) * tmp[i, j] +
                (15 * eps_KO / 64) * tmp[i+1, j] +
                (-6 * eps_KO / 64) * tmp[i+2, j] +
                (eps_KO / 64) * tmp[i+3, j]
            )
        end
    end
    return nothing
end
