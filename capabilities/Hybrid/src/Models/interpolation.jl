# WL-faithful interpolation leaf — the WaSABI self-force models build their on-disk-table
# interpolants over the ONE consolidated NumericsCore.interp suite. Phase7-D DRY collapse: the
# moving-window order-8 `Interp1WL`, the tensor-product global-Lagrange `Interp2Hermite` /
# `Interp3Lagrange`, their analytic derivatives / mixed `partial`s all live in NumericsCore and are
# imported here (see `Interpolation.jl`). This leaf keeps only (a) the historical Hybrid-facing
# names/adapters — `Interp1`/`interp1_order8`/`derivative`/`derivative2` — over the NumericsCore
# `deriv_order` API, and (b) the two composites NumericsCore does not own: `interp_logexp` and the
# Kerr `kerr_grid_coords` abscissa remap (reuses the canonical `Spacetime.kerr_isco`).
#
# Type-generic `{T}` throughout: BigFloat tables and queries are inherited, never hard-Float64'd.

# --- 1D order-8 moving-window Lagrange (NumericsCore.Interp1WL) -------------

"""
    Interp1{TX,TY}

Historical Hybrid alias of the canonical [`NumericsCore.Interp1WL`](@ref) degree-`order`
piecewise-Lagrange interpolant (matches Wolfram `Interpolation[data, InterpolationOrder->order]`).
`x` must be sorted ascending; `y` may be real or complex. Call `itp(x)` for the value and
[`derivative`](@ref)`(itp, x)` for the first derivative.
"""
const Interp1 = Interp1WL

"""
    interp1_order8(x, y; order=8) -> Interp1
    interp1_order8(pairs;  order=8) -> Interp1

Build an [`Interp1`](@ref) via the canonical [`NumericsCore.interp1_movingwindow`](@ref). `pairs` is
an iterable of 2-element `(x, y)` tuples/vectors (the on-disk WL table layout); `order` defaults to
8 (`InterpolationOrder->8`).
"""
interp1_order8(x::AbstractVector{<:Real}, y::AbstractVector{<:Number}; order::Int = 8) =
    interp1_movingwindow(x, y; order = order)

function interp1_order8(pairs; order::Int = 8)
    x = [float(p[1]) for p in pairs]
    y = [p[2] for p in pairs]
    return interp1_movingwindow(x, y; order = order)
end

"""
    derivative(itp::Interp1, x) -> y'

First derivative of the order-`order` interpolant at `x` (matches `Interpolation'[...][x]`); the
`deriv_order=1` channel of the canonical `Interp1WL`.
"""
derivative(itp::Interp1WL, xq::Real) = itp(xq; deriv_order = 1)

"""
    derivative2(itp::Interp1, x) -> y''

Second derivative of the order-`order` interpolant at `x` (matches `Interpolation''[...][x]`); the
`deriv_order=2` channel of the canonical `Interp1WL`.
"""
derivative2(itp::Interp1WL, xq::Real) = itp(xq; deriv_order = 2)

# --- log-exp composite -----------------------------------------------------

"""
    interp_logexp(itp::Interp1) -> Function

Return `r0 -> real(exp(itp(log(r0))))`, the `Re[Exp[interp[Log[r0]]]]` composite used for the
2SF flux table, where `itp` interpolates a complex `{Log[r0], Log[flux]}` table.
"""
interp_logexp(itp::Interp1WL) = r0 -> real(exp(itp(log(r0))))

# --- Kerr abscissa remap ---------------------------------------------------

"""
    kerr_grid_coords(a, r) -> (log(1-a), sqrt(rISCO[a]/r))

The 2D grid abscissa remap applied before [`interp2_hermite`](@ref). Reuses the canonical
`Spacetime.kerr_isco` for `rISCO[a]`.
"""
kerr_grid_coords(a::Real, r::Real) = (log(1 - a), sqrt(kerr_isco(a) / r))
