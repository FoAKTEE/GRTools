# Interpolation library — the analytically-groundable core of the GrTools interp suite, consolidated
# from the four legacy copies into one canonical type-generic `{T}` set with a coefficient/eval split
# (coefficients built once, evaluated batched — the GPU-eval seam). Everything is exact on the
# polynomial degree it interpolates, which is the grounding oracle (no heavy legacy dep needed).
#
# Sources (rewritten, not copied): CubicSplineInterpolant ← KerrBackwards/interpolant.jl (== FEW
# cubicspline.jl, scipy not-a-knot); BicubicNotAKnot ← FEW/ampinterp2d.jl (== KerrBackwards
# NotAKnotBicubic, one canonical primitive); TricubicNotAKnot ← FEW/kerrflux.jl; GslBicubic +
# natural_spline_node_derivs + gsl_bsearch ← KerrBackwards/gslbicubic.jl; cubic_hermite ←
# Teukolsky/radial_numint.jl `_hermite`; fornberg_weights/highorder_deriv/cheb_*/bary_*/local_interp
# ← DiscEMRIForce/numerics.jl + KerrLorenzCirc/radial_kerr.jl; fd1c/fd2c ← KerrLorenzCirc/operators.jl.

# =====================================================================================================
# Not-a-knot cubic spline (scipy-faithful) — the shared 1D primitive behind the bicubic/tricubic tensors
# =====================================================================================================
"""
    CubicSplineInterpolant(t, y)

Not-a-knot cubic spline through strictly-increasing nodes `t` (length `n ≥ 3`). `y` is a length-`n`
vector (single curve) or an `ncurves × n` matrix (curves share the abscissa `t`). On segment
`[t_i, t_{i+1}]`, `S(t) = y_i + c1_i·x + c2_i·x² + c3_i·x³`, `x = t − t_i`. Call `(spline)(tnew[, k];
deriv_order)` (`deriv_order ∈ 0:3`). Matches `scipy.interpolate.CubicSpline` (not-a-knot), exact on
cubic data. Type-generic (Float64/BigFloat). Coefficients `c1/c2/c3` are precomputed (the eval split).
"""
struct CubicSplineInterpolant{T<:Real}
    t::Vector{T}
    y::Matrix{T}
    c1::Matrix{T}
    c2::Matrix{T}
    c3::Matrix{T}
end

function CubicSplineInterpolant(t::AbstractVector, y::AbstractVecOrMat)
    T = float(promote_type(eltype(t), eltype(y)))
    T <: Real || throw(ArgumentError("CubicSplineInterpolant: t and y must have real element types"))
    tt = T.(t)
    n = length(tt)
    n ≥ 3 || throw(ArgumentError("CubicSplineInterpolant: need ≥3 nodes for not-a-knot, got $n"))
    all(diff(tt) .> 0) || throw(ArgumentError("CubicSplineInterpolant: nodes `t` must be strictly increasing"))

    Y = y isa AbstractVector ? reshape(T.(y), 1, n) : Matrix{T}(y)
    size(Y, 2) == n || throw(DimensionMismatch("CubicSplineInterpolant: y has $(size(Y,2)) columns but t has $n nodes"))
    ncurves = size(Y, 1)
    c1 = Matrix{T}(undef, ncurves, n)
    c2 = Matrix{T}(undef, ncurves, n)
    c3 = Matrix{T}(undef, ncurves, n)
    if n == 3
        # The two not-a-knot boundary equations coincide for three nodes. As in
        # SciPy, use the unique interpolating parabola instead of solving the
        # resulting singular tridiagonal system.
        h1 = tt[2] - tt[1]
        h2 = tt[3] - tt[2]
        @inbounds for k in 1:ncurves
            s1 = (Y[k, 2] - Y[k, 1]) / h1
            s2 = (Y[k, 3] - Y[k, 2]) / h2
            q = (s2 - s1) / (h1 + h2)
            c1[k, 1] = s1 - q * h1
            c1[k, 2] = s2 - q * h2
            c2[k, 1] = c2[k, 2] = q
            c3[k, 1] = c3[k, 2] = zero(T)
            c1[k, 3] = c1[k, 2]
            c2[k, 3] = c2[k, 2]
            c3[k, 3] = c3[k, 2]
        end
        return CubicSplineInterpolant{T}(tt, Matrix(Y), c1, c2, c3)
    end

    diag = Vector{T}(undef, n)
    upper = zeros(T, n)
    lower = zeros(T, n)
    dx2 = tt[3] - tt[2]
    diag[1] = dx2; upper[1] = tt[3] - tt[1]
    for i in 2:n-1
        a = tt[i] - tt[i-1]; b = tt[i+1] - tt[i]
        diag[i] = 2 * (a + b); upper[i] = a; lower[i] = b
    end
    diag[n] = tt[n-1] - tt[n-2]; lower[n] = tt[n] - tt[n-2]
    M = Tridiagonal(lower[2:n], diag, upper[1:n-1])

    b = Vector{T}(undef, n)
    for k in 1:ncurves
        yv = @view Y[k, :]
        d1 = tt[2] - tt[1]; d2 = tt[3] - tt[2]; dd = tt[3] - tt[1]
        s1 = (yv[2] - yv[1]) / d1; s2 = (yv[3] - yv[2]) / d2
        b[1] = ((d1 + 2 * dd) * d2 * s1 + d1^2 * s2) / dd
        for i in 2:n-1
            a = tt[i] - tt[i-1]; c = tt[i+1] - tt[i]
            sa = (yv[i] - yv[i-1]) / a; sc = (yv[i+1] - yv[i]) / c
            b[i] = 3 * (c * sa + a * sc)
        end
        d1 = tt[n-1] - tt[n-2]; d2 = tt[n] - tt[n-1]; dd = tt[n] - tt[n-2]
        s1 = (yv[n-1] - yv[n-2]) / d1; s2 = (yv[n] - yv[n-1]) / d2
        b[n] = (d2^2 * s1 + (2 * dd + d2) * d1 * s2) / dd

        s = M \ b
        for i in 1:n-1
            dx = tt[i+1] - tt[i]
            slope = (yv[i+1] - yv[i]) / dx
            tcoef = (s[i] + s[i+1] - 2 * slope) / dx
            c1[k, i] = s[i]
            c2[k, i] = (slope - s[i]) / dx - tcoef
            c3[k, i] = tcoef / dx
        end
        c1[k, n] = c1[k, n-1]; c2[k, n] = c2[k, n-1]; c3[k, n] = c3[k, n-1]
    end
    return CubicSplineInterpolant{T}(tt, Matrix(Y), c1, c2, c3)
end

@inline function _eval_spline(s::CubicSplineInterpolant, k::Integer, tnew::Real, deriv_order::Integer)
    n = length(s.t)
    i = clamp(searchsortedlast(s.t, tnew), 1, n - 1)
    x = tnew - s.t[i]
    y = s.y[k, i]; c1 = s.c1[k, i]; c2 = s.c2[k, i]; c3 = s.c3[k, i]
    if deriv_order == 0
        return y + c1 * x + c2 * x^2 + c3 * x^3
    elseif deriv_order == 1
        return c1 + 2 * c2 * x + 3 * c3 * x^2
    elseif deriv_order == 2
        return 2 * c2 + 6 * c3 * x
    elseif deriv_order == 3
        return 6 * c3
    else
        throw(ArgumentError("CubicSplineInterpolant: deriv_order must be in 0:3, got $deriv_order"))
    end
end

(s::CubicSplineInterpolant)(tnew::Real, k::Integer; deriv_order::Integer=0) =
    _eval_spline(s, k, tnew, deriv_order)
(s::CubicSplineInterpolant)(tnew::Real; deriv_order::Integer=0) =
    [_eval_spline(s, k, tnew, deriv_order) for k in 1:size(s.y, 1)]

# =====================================================================================================
# Tensor-product not-a-knot bicubic / tricubic (the single canonical surface primitive — AHA)
# =====================================================================================================
"""
    BicubicNotAKnot(w, u, vals)

Tensor-product bicubic not-a-knot spline over the rectangular grid `w × u`, `vals[i, j]` the sample
at `(w[i], u[j])`. Equals `scipy RectBivariateSpline(w, u, vals, kx=3, ky=3, s=0)` and the
KerrBackwards `NotAKnotBicubic` (one canonical primitive). Exact on bicubic data.
"""
struct BicubicNotAKnot{T<:Real}
    w::Vector{T}
    u::Vector{T}
    rows::Vector{CubicSplineInterpolant{T}}
end

function BicubicNotAKnot(w::AbstractVector, u::AbstractVector, vals::AbstractMatrix)
    nw, nu = size(vals)
    nw == length(w) || throw(DimensionMismatch("size(vals,1) must equal length(w)"))
    nu == length(u) || throw(DimensionMismatch("size(vals,2) must equal length(u)"))
    T = float(promote_type(eltype(w), eltype(u), real(eltype(vals))))
    uv = collect(T, u)
    rows = [CubicSplineInterpolant(uv, T[vals[i, j] for j in 1:nu]) for i in 1:nw]
    return BicubicNotAKnot{T}(collect(T, w), uv, rows)
end

function (bs::BicubicNotAKnot)(wq::Real, uq::Real)
    col = [row(uq, 1) for row in bs.rows]
    return CubicSplineInterpolant(bs.w, col)(wq, 1)
end

"""
    TricubicNotAKnot(u, w, z, V)

Tensor-product tricubic not-a-knot interpolant, `V[i,j,k]` the value at `(u[i], w[j], z[k])`. Built
by nested not-a-knot evaluation (along `u`, then `w`, then `z`). Exact on tricubic data / at nodes.
"""
struct TricubicNotAKnot{T<:Real}
    u::Vector{T}
    w::Vector{T}
    z::Vector{T}
    V::Array{T,3}
end

function TricubicNotAKnot(u::AbstractVector, w::AbstractVector, z::AbstractVector,
                         V::AbstractArray{<:Real,3})
    T = float(promote_type(eltype(u), eltype(w), eltype(z), eltype(V)))
    return TricubicNotAKnot{T}(collect(T, u), collect(T, w), collect(T, z), Array{T}(V))
end

function (s::TricubicNotAKnot{T})(uq::Real, wq::Real, zq::Real) where {T}
    nu, nw, nz = size(s.V)
    bz = Matrix{T}(undef, nw, nz)
    for k in 1:nz, j in 1:nw
        bz[j, k] = CubicSplineInterpolant(s.u, collect(@view s.V[:, j, k]))(uq, 1)
    end
    cz = Vector{T}(undef, nz)
    for k in 1:nz
        cz[k] = CubicSplineInterpolant(s.w, collect(@view bz[:, k]))(wq, 1)
    end
    return CubicSplineInterpolant(s.z, cz)(zq, 1)
end

# =====================================================================================================
# GSL bicubic (natural-cubic-spline-based Hermite patch) — exact on bilinear data / at nodes
# =====================================================================================================
"""
    natural_spline_node_derivs(t, v) -> Vector

First derivatives `dv/dt` at every knot of the NATURAL cubic spline through `(t, v)` (second
derivative zero at both ends — `gsl_interp_cspline`). Exact slope for linear data.
"""
function natural_spline_node_derivs(t::AbstractVector{T}, v::AbstractVector{T}) where {T<:Real}
    n = length(t)
    n == length(v) || throw(DimensionMismatch("t and v length mismatch"))
    d = Vector{T}(undef, n)
    if n == 2
        d[1] = d[2] = (v[2] - v[1]) / (t[2] - t[1])
        return d
    end
    h = T[t[k+1] - t[k] for k in 1:n-1]
    M = zeros(T, n)
    ni = n - 2
    if ni >= 1
        sub = Vector{T}(undef, ni); dia = Vector{T}(undef, ni)
        sup = Vector{T}(undef, ni); rhs = Vector{T}(undef, ni)
        for k in 1:ni
            i = k + 1
            sub[k] = h[i-1]; dia[k] = 2 * (h[i-1] + h[i]); sup[k] = h[i]
            rhs[k] = 6 * ((v[i+1] - v[i]) / h[i] - (v[i] - v[i-1]) / h[i-1])
        end
        for k in 2:ni
            w = sub[k] / dia[k-1]
            dia[k] -= w * sup[k-1]; rhs[k] -= w * rhs[k-1]
        end
        Mint = Vector{T}(undef, ni)
        Mint[ni] = rhs[ni] / dia[ni]
        for k in ni-1:-1:1
            Mint[k] = (rhs[k] - sup[k] * Mint[k+1]) / dia[k]
        end
        for k in 1:ni
            M[k+1] = Mint[k]
        end
    end
    for i in 1:n-1
        d[i] = (v[i+1] - v[i]) / h[i] - h[i] * (2 * M[i] + M[i+1]) / 6
    end
    d[n] = (v[n] - v[n-1]) / h[n-1] + h[n-1] * (M[n-1] + 2 * M[n]) / 6
    return d
end

"""`gsl_bsearch(arr, x)` — largest index `i` with `arr[i] ≤ x`, clamped to `[1, n-1]` (GSL `gsl_interp_bsearch`)."""
function gsl_bsearch(arr::AbstractVector, x::Real)
    n = length(arr)
    x <= arr[1] && return 1
    x >= arr[n] && return n - 1
    lo = 1; hi = n - 1
    while hi - lo > 0
        mid = (lo + hi + 1) >>> 1
        if arr[mid] <= x
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

"""
    GslBicubic(x, y, Z)

GSL `gsl_interp2d_bicubic` bicubic interpolant over `x × y`, `Z[i,j]` the sample at `(x[i], y[j])`.
Corner `z_x`/`z_y`/`z_xy` come from natural cubic splines (`z_xy = ∂_x z_y`), composed into a bicubic
Hermite patch. Distinct from `BicubicNotAKnot` (natural, not not-a-knot BC). Exact on bilinear data.
"""
struct GslBicubic{T<:Real}
    x::Vector{T}
    y::Vector{T}
    Z::Matrix{T}
    zx::Matrix{T}
    zy::Matrix{T}
    zxy::Matrix{T}
end

function GslBicubic(x::AbstractVector, y::AbstractVector, Z::AbstractMatrix)
    nx = length(x); ny = length(y)
    size(Z) == (nx, ny) || throw(DimensionMismatch("size(Z) must be (length(x), length(y))"))
    T = float(promote_type(eltype(x), eltype(y), eltype(Z)))
    xv = collect(T, x); yv = collect(T, y); Zm = Matrix{T}(Z)
    zx = Matrix{T}(undef, nx, ny); zy = Matrix{T}(undef, nx, ny); zxy = Matrix{T}(undef, nx, ny)
    for j in 1:ny
        @views zx[:, j] .= natural_spline_node_derivs(xv, Zm[:, j])
    end
    for i in 1:nx
        @views zy[i, :] .= natural_spline_node_derivs(yv, Zm[i, :])
    end
    for j in 1:ny
        @views zxy[:, j] .= natural_spline_node_derivs(xv, zy[:, j])
    end
    return GslBicubic{T}(xv, yv, Zm, zx, zy, zxy)
end

function (gb::GslBicubic)(xq::Real, yq::Real)
    xi = gsl_bsearch(gb.x, xq); yi = gsl_bsearch(gb.y, yq)
    xmin = gb.x[xi]; xmax = gb.x[xi+1]; ymin = gb.y[yi]; ymax = gb.y[yi+1]
    zminmin = gb.Z[xi, yi]; zminmax = gb.Z[xi, yi+1]
    zmaxmin = gb.Z[xi+1, yi]; zmaxmax = gb.Z[xi+1, yi+1]
    dx = xmax - xmin; dy = ymax - ymin
    t = (xq - xmin) / dx; u = (yq - ymin) / dy
    dt = 1 / dx; du = 1 / dy
    zxminmin = gb.zx[xi, yi] / dt; zxminmax = gb.zx[xi, yi+1] / dt
    zxmaxmin = gb.zx[xi+1, yi] / dt; zxmaxmax = gb.zx[xi+1, yi+1] / dt
    zyminmin = gb.zy[xi, yi] / du; zyminmax = gb.zy[xi, yi+1] / du
    zymaxmin = gb.zy[xi+1, yi] / du; zymaxmax = gb.zy[xi+1, yi+1] / du
    zxyminmin = gb.zxy[xi, yi] / (dt*du); zxyminmax = gb.zxy[xi, yi+1] / (dt*du)
    zxymaxmin = gb.zxy[xi+1, yi] / (dt*du); zxymaxmax = gb.zxy[xi+1, yi+1] / (dt*du)
    t0 = one(t); t1 = t; t2 = t*t; t3 = t*t2
    u0 = one(u); u1 = u; u2 = u*u; u3 = u*u2
    z = zero(t)
    z += zminmin * t0*u0
    z += zyminmin * t0*u1
    z += (-3*zminmin + 3*zminmax - 2*zyminmin - zyminmax) * t0*u2
    z += (2*zminmin - 2*zminmax + zyminmin + zyminmax) * t0*u3
    z += zxminmin * t1*u0
    z += zxyminmin * t1*u1
    z += (-3*zxminmin + 3*zxminmax - 2*zxyminmin - zxyminmax) * t1*u2
    z += (2*zxminmin - 2*zxminmax + zxyminmin + zxyminmax) * t1*u3
    z += (-3*zminmin + 3*zmaxmin - 2*zxminmin - zxmaxmin) * t2*u0
    z += (-3*zyminmin + 3*zymaxmin - 2*zxyminmin - zxymaxmin) * t2*u1
    z += (9*zminmin - 9*zmaxmin + 9*zmaxmax - 9*zminmax + 6*zxminmin + 3*zxmaxmin - 3*zxmaxmax - 6*zxminmax + 6*zyminmin - 6*zymaxmin - 3*zymaxmax + 3*zyminmax + 4*zxyminmin + 2*zxymaxmin + zxymaxmax + 2*zxyminmax) * t2*u2
    z += (-6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 4*zxminmin - 2*zxmaxmin + 2*zxmaxmax + 4*zxminmax - 3*zyminmin + 3*zymaxmin + 3*zymaxmax - 3*zyminmax - 2*zxyminmin - zxymaxmin - zxymaxmax - 2*zxyminmax) * t2*u3
    z += (2*zminmin - 2*zmaxmin + zxminmin + zxmaxmin) * t3*u0
    z += (2*zyminmin - 2*zymaxmin + zxyminmin + zxymaxmin) * t3*u1
    z += (-6*zminmin + 6*zmaxmin - 6*zmaxmax + 6*zminmax - 3*zxminmin - 3*zxmaxmin + 3*zxmaxmax + 3*zxminmax - 4*zyminmin + 4*zymaxmin + 2*zymaxmax - 2*zyminmax - 2*zxyminmin - 2*zxymaxmin - zxymaxmax - zxyminmax) * t3*u2
    z += (4*zminmin - 4*zmaxmin + 4*zmaxmax - 4*zminmax + 2*zxminmin + 2*zxmaxmin - 2*zxmaxmax - 2*zxminmax + 2*zyminmin - 2*zymaxmin - 2*zymaxmax + 2*zyminmax + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax) * t3*u3
    return z
end

# =====================================================================================================
# Cubic Hermite (single segment + mesh) — exact on cubics
# =====================================================================================================
"""
    cubic_hermite(x0, x1, f0, f1, d0, d1, xq) -> value

Cubic-Hermite interpolation on `[x0, x1]` matching values `f0,f1` and derivatives `d0,d1` at the
endpoints (the unique interpolating cubic ⇒ exact on cubic data). Teukolsky `_hermite` basis.
"""
function cubic_hermite(x0::Real, x1::Real, f0, f1, d0, d1, xq::Real)
    hh = x1 - x0
    t = (xq - x0) / hh
    h00 = (1 + 2t) * (1 - t)^2
    h10 = t * (1 - t)^2
    h01 = t^2 * (3 - 2t)
    h11 = t^2 * (t - 1)
    return h00 * f0 + h10 * hh * d0 + h01 * f1 + h11 * hh * d1
end

"""
    cubic_hermite(x, f, d, xq) -> value

Mesh cubic-Hermite: locate the bracketing interval of `xq` in the sorted nodes `x` (values `f`,
slopes `d`) and Hermite-interpolate there.
"""
function cubic_hermite(x::AbstractVector{<:Real}, f::AbstractVector, d::AbstractVector, xq::Real)
    n = length(x)
    i = clamp(searchsortedlast(x, xq), 1, n - 1)
    return cubic_hermite(x[i], x[i+1], f[i], f[i+1], d[i], d[i+1], xq)
end

# =====================================================================================================
# Fornberg finite-difference weights + high-order grid derivative
# =====================================================================================================
"""
    fornberg_weights(z, xnodes, m) -> Matrix (n × (m+1))

Finite-difference weights for derivatives `0..m` at `z` on arbitrarily-spaced `xnodes`: `C[i, d+1]`
is the weight of `f(xnodes[i])` in `f^(d)(z) ≈ Σ_i C[i,d+1] f(xnodes[i])`. Exact for polynomials of
degree ≤ `n−1`. (Fornberg 1988.)
"""
function fornberg_weights(z::T, xnodes::AbstractVector{T}, m::Integer) where {T<:Real}
    n = length(xnodes)
    n >= m + 1 || throw(ArgumentError("fornberg_weights: need ≥ m+1 nodes for the m-th derivative"))
    C = zeros(T, n, m + 1)
    c1 = one(T); c4 = xnodes[1] - z
    C[1, 1] = one(T)
    @inbounds for i in 2:n
        mn = min(i - 1, m); c2 = one(T); c5 = c4; c4 = xnodes[i] - z
        for j in 1:(i-1)
            c3 = xnodes[i] - xnodes[j]; c2 *= c3
            if j == i - 1
                for k in mn:-1:1
                    C[i, k+1] = c1 * (k * C[i-1, k] - c5 * C[i-1, k+1]) / c2
                end
                C[i, 1] = -c1 * c5 * C[i-1, 1] / c2
            end
            for k in mn:-1:1
                C[j, k+1] = (c4 * C[j, k+1] - k * C[j, k]) / c3
            end
            C[j, 1] = c4 * C[j, 1] / c3
        end
        c1 = c2
    end
    return C
end

"""
    highorder_deriv(f, rgrid, k; npts=7, order=1)

`order`-th derivative of `f` (sampled on non-uniform `rgrid`) at index `k` via Fornberg weights over
a centred `npts`-point stencil (edge-clamped). `npts=7` ⇒ O(h⁶) truncation. `f` real or complex.
"""
function highorder_deriv(f::AbstractVector, rgrid::AbstractVector{T}, k::Integer;
                        npts::Integer=7, order::Integer=1) where {T<:Real}
    N = length(rgrid)
    npts = min(npts, N)
    half = (npts - 1) ÷ 2
    lo = clamp(k - half, 1, N - npts + 1)
    hi = lo + npts - 1
    C = fornberg_weights(rgrid[k], collect(@view rgrid[lo:hi]), order)
    acc = zero(eltype(f))
    @inbounds for (j, i) in enumerate(lo:hi)
        acc += C[j, order+1] * f[i]
    end
    return acc
end

# =====================================================================================================
# Barycentric Lagrange interpolation (general / local) + Chebyshev spectral differentiation
# =====================================================================================================
"""
    bary_weights(x) -> Vector

Barycentric weights `w_i = 1/∏_{j≠i}(x_i − x_j)` for the node set `x` (any distribution).
"""
function bary_weights(x::AbstractVector{T}) where {T<:Real}
    N = length(x); w = ones(T, N)
    @inbounds for j in 1:N, k in 1:N
        k == j && continue
        w[j] /= (x[j] - x[k])
    end
    return w
end

"""
    bary_eval(x, w, f, xq; nodetol=0) -> value

Second-form barycentric Lagrange evaluation of samples `f` at nodes `x` (weights `w`) at the scalar
query `xq`. `f` may be real or complex. Exact for polynomials of degree ≤ `length(x)−1`.

`nodetol` is the node-collision radius: a query within `abs(xq − x_j) ≤ nodetol` returns the nodal
sample `f[j]` verbatim instead of forming the 0/0 second-form quotient. The default `0` is the
exact test `xq == x_j`; a DG probe row whose `ξ` is *computed* (so a node coincidence is only exact
to roundoff) passes a few `eps` — this is the one knob the three former private copies of this
kernel disagreed on.
"""
@inline function bary_eval(x::AbstractVector{<:Real}, w::AbstractVector, f::AbstractVector, xq::Real;
                           nodetol::Real = zero(real(eltype(w))))
    num = zero(eltype(f)); den = zero(eltype(w))
    @inbounds for j in eachindex(x)
        d = xq - x[j]
        abs(d) <= nodetol && return f[j]
        t = w[j] / d
        num += t * f[j]; den += t
    end
    return num / den
end

"""
    bary_row(x, w, xq; nodetol=0) -> Vector

Second-form barycentric Lagrange CARDINAL ROW `[ℓ_1(xq) … ℓ_N(xq)]` for the nodes `x` with weights
`w` (`bary_weights(x)`) — the matrixified companion of [`bary_eval`](@ref), so that
`dot(bary_row(x, w, xq), f) ≈ bary_eval(x, w, f, xq)` while the row itself can be stored as one line
of an interpolation/probe matrix. A node collision (`abs(xq − x_j) ≤ nodetol`, see
[`bary_eval`](@ref)) returns the exact unit row `e_j`, which is what keeps the row finite when the
query sits on (or within roundoff of) a node.

Single home of the row kernel: `EllipticCore` reference-element face/probe rows and the
`Perturbation.TimeDomain` DG probe rows both call it.
"""
function bary_row(x::AbstractVector{<:Real}, w::AbstractVector, xq::Real;
                  nodetol::Real = zero(real(eltype(w))))
    RT = eltype(w)
    N = length(x)
    row = Vector{RT}(undef, N)
    s = zero(RT)
    @inbounds for j in 1:N
        d = xq - x[j]
        if abs(d) <= nodetol
            fill!(row, zero(RT))
            row[j] = one(RT)
            return row
        end
        t = w[j] / d
        row[j] = t
        s += t
    end
    row ./= s
    return row
end

"""
    bary_interp(xq, xnodes, fvals) -> value(s)

Convenience barycentric interpolation (builds the weights internally). `xq` scalar or vector.
"""
function bary_interp(xq, xnodes::AbstractVector{T}, fvals::AbstractVector) where {T<:Real}
    w = bary_weights(xnodes)
    return xq isa AbstractVector ? [bary_eval(xnodes, w, fvals, x) for x in xq] :
           bary_eval(xnodes, w, fvals, xq)
end

"""
    local_interp(xq, xnodes, fvals; npts=8) -> value(s)

Stable interpolation of `f` on sorted ~equispaced `xnodes` via a LOCAL Lagrange polynomial over the
`npts` nearest nodes (Fornberg 0-th-derivative weights) — avoids the Runge instability of global
barycentric interpolation on equispaced grids.
"""
function local_interp(xq, xnodes::AbstractVector{T}, fvals::AbstractVector; npts::Integer=8) where {T<:Real}
    n = length(xnodes)
    npts = min(npts, n)
    _eval(x) = begin
        k = clamp(searchsortedfirst(xnodes, T(x)), 1, n)
        half = (npts - 1) ÷ 2
        lo = clamp(k - half, 1, n - npts + 1); hi = lo + npts - 1
        C = fornberg_weights(T(x), collect(@view xnodes[lo:hi]), 0)
        acc = zero(eltype(fvals))
        @inbounds for (j, i) in enumerate(lo:hi)
            acc += C[j, 1] * fvals[i]
        end
        acc
    end
    return xq isa AbstractVector ? [_eval(x) for x in xq] : _eval(xq)
end

"""
    cheb_nodes(T, N, a, b) -> Vector

Chebyshev–Gauss–Lobatto nodes `x_j = ½(a+b) + ½(b−a)cos(jπ/N)`, `j=0..N`, on `[a,b]` (descending,
endpoints included). `N+1` nodes.
"""
cheb_nodes(::Type{T}, N::Integer, a::T, b::T) where {T<:Real} =
    [(a + b) / 2 + (b - a) / 2 * cos(T(j) * T(pi) / T(N)) for j in 0:N]

"""
    cheb_diffmatrix(T, N, a, b) -> (D, nodes)

Chebyshev (Trefethen) spectral first-derivative matrix `D` on `N+1` Gauss–Lobatto nodes over `[a,b]`
(plus the nodes). `D*fvals` is the spectral derivative; spectrally exact for polynomials of degree ≤ N.
"""
function cheb_diffmatrix(::Type{T}, N::Integer, a::T, b::T) where {T<:Real}
    x = cheb_nodes(T, N, a, b)
    c = T[(j == 0 || j == N) ? T(2) : one(T) for j in 0:N]
    for j in 0:N
        c[j+1] *= isodd(j) ? -one(T) : one(T)
    end
    D = zeros(T, N + 1, N + 1)
    @inbounds for i in 0:N, j in 0:N
        i != j && (D[i+1, j+1] = (c[i+1] / c[j+1]) / (x[i+1] - x[j+1]))
    end
    for i in 0:N
        D[i+1, i+1] = -sum(@view D[i+1, :])
    end
    return D, x
end

"""
    cheb_deriv(f, rgrid; N, resample_npts=14) -> Vector

Spectral first derivative of a SMOOTH field `f` sampled on `rgrid`: stable-local-resample onto `N+1`
Chebyshev nodes, apply the differentiation matrix, map back. Exponentially accurate for smooth `f`.
"""
function cheb_deriv(f::AbstractVector, rgrid::AbstractVector{T};
                   N::Integer=min(length(rgrid) - 1, 48), resample_npts::Integer=14) where {T<:Real}
    a = rgrid[1]; b = rgrid[end]
    D, xc = cheb_diffmatrix(T, N, a, b)
    fc = local_interp(xc, rgrid, f; npts=resample_npts)
    dfc = D * fc
    return bary_interp(collect(rgrid), xc, dfc)
end

# Chebyshev T_0..T_deg basis matrix on x∈[a,b] (ξ = 2(x−a)/(b−a)−1).
function _cheb_vandermonde(x::AbstractVector{T}, deg::Integer, a::T, b::T) where {T<:Real}
    n = length(x)
    ξ = [2 * (x[i] - a) / (b - a) - 1 for i in 1:n]
    V = zeros(T, n, deg + 1)
    @inbounds for i in 1:n
        V[i, 1] = one(T)
        deg >= 1 && (V[i, 2] = ξ[i])
        for j in 2:deg
            V[i, j+1] = 2 * ξ[i] * V[i, j] - V[i, j-1]
        end
    end
    return V
end

# d/dx of T_0..T_deg at x (chain rule dξ/dx; T'_k = k U_{k-1}).
function _cheb_vandermonde_deriv(x::AbstractVector{T}, deg::Integer, a::T, b::T) where {T<:Real}
    n = length(x); sc = T(2) / (b - a)
    ξ = [2 * (x[i] - a) / (b - a) - 1 for i in 1:n]
    Vd = zeros(T, n, deg + 1)
    @inbounds for i in 1:n
        Um1 = zero(T); U0 = one(T)
        Vd[i, 1] = zero(T)
        for k in 1:deg
            Vd[i, k+1] = k * U0 * sc
            U1 = 2 * ξ[i] * U0 - Um1
            Um1 = U0; U0 = U1
        end
    end
    return Vd
end

"""
    cheb_smooth_deriv(f, rgrid; deg=8) -> Vector

Least-squares Chebyshev-smoothing first derivative: fit a degree-`deg` Chebyshev polynomial to `f`
on `rgrid` (deg ≪ length ⇒ averages point noise), then evaluate the analytic derivative of the fit.
Exact derivative when `f` is itself a polynomial of degree ≤ `deg`. `f` may be complex (re/im fit).
"""
function cheb_smooth_deriv(f::AbstractVector, rgrid::AbstractVector{T}; deg::Integer=8) where {T<:Real}
    a = rgrid[1]; b = rgrid[end]
    deg = min(deg, length(rgrid) - 1)
    V = _cheb_vandermonde(rgrid, deg, a, b)
    Vd = _cheb_vandermonde_deriv(rgrid, deg, a, b)
    if eltype(f) <: Complex
        cr = V \ real.(f); ci = V \ imag.(f)
        return complex.(Vd * cr, Vd * ci)
    else
        c = V \ collect(f)
        return Vd * c
    end
end

# =====================================================================================================
# Compact 4th-order central finite differences (closure form) — exact on cubics
# =====================================================================================================
"""`fd1c(f, r, h)` — 5-point 4th-order central first derivative `f'(r)` (KerrLorenzCirc `fd1c`)."""
fd1c(f, r, h) = (-f(r + 2h) + 8f(r + h) - 8f(r - h) + f(r - 2h)) / (12h)

"""`fd2c(f, r, h)` — 5-point 4th-order central second derivative `f''(r)` (KerrLorenzCirc `fd2c`)."""
fd2c(f, r, h) = (-f(r + 2h) + 16f(r + h) - 30f(r) + 16f(r - h) - f(r - 2h)) / (12h^2)

# =====================================================================================================
# Wolfram-Language-faithful moving-window / global-tensor Lagrange interpolants (WaSABI Interp1/2/3)
# =====================================================================================================
"""
    Interp1WL{TX,TY}

Degree-`order` piecewise-Lagrange ("moving window") interpolant reproducing Wolfram
`Interpolation[data, InterpolationOrder->order]` to machine precision. For a query in the bracketing
interval `[x_i, x_{i+1}]` the degree-`order` polynomial through the `order+1` points whose x-span is
*centered on that interval* is used (the interval-centered stencil, sliding to stay in range at the
boundaries). `x` sorted ascending; `y` real or complex. Call `itp(xq; deriv_order=0)` (`0:2` — value,
`Interpolation'`, `Interpolation''`). Exact on polynomials of degree ≤ `order`. Type-generic `{T}`.
"""
struct Interp1WL{TX<:Real,TY<:Number}
    x::Vector{TX}
    y::Vector{TY}
    order::Int
end

"""
    interp1_movingwindow(x, y; order=8) -> Interp1WL

Build an [`Interp1WL`](@ref) over sorted-ascending `x` and samples `y` (`length ≥ order+1`).
`order` defaults to 8 (`InterpolationOrder->8`).
"""
function interp1_movingwindow(x::AbstractVector{<:Real}, y::AbstractVector{<:Number}; order::Int=8)
    length(x) == length(y) || throw(ArgumentError("x and y must have equal length"))
    length(x) ≥ order + 1 || throw(ArgumentError("interp1_movingwindow: need ≥ order+1 = $(order+1) points"))
    issorted(x) || throw(ArgumentError("interp1_movingwindow: x must be sorted ascending"))
    return Interp1WL(collect(float.(x)), collect(y), order)
end

# Interval-centered window of `order+1` consecutive indices (Wolfram's reverse-engineered stencil):
# among the contiguous windows containing `[x_i,x_{i+1}]`, the one whose span midpoint is closest to
# the interval midpoint (ties toward the larger `lo`). Reproduces WL on uniform AND steep log-spaced
# abscissae. The `1e-10`/`1e-15` tie tolerances are dimensionless comparisons (grid-independent).
@inline function _wl_window(itp::Interp1WL, xq::Real)
    x = itp.x; N = length(x); k = itp.order
    i = clamp(searchsortedlast(x, xq), 1, N - 1)
    Xm = (x[i] + x[i+1]) / 2
    lomin = max(1, i + 1 - k); lomax = min(N - k, i)
    best = lomin; bd = abs(Xm - (x[lomin] + x[lomin+k]) / 2)
    @inbounds for lo in (lomin+1):lomax
        d = abs(Xm - (x[lo] + x[lo+k]) / 2)
        if d <= bd * (1 + 1e-10) + 1e-15
            bd = min(bd, d); best = lo
        end
    end
    return best:(best + k)
end

# Value / 1st / 2nd derivative of the local Lagrange polynomial on the window (node-safe product form).
function _wl_eval(itp::Interp1WL, xq::Real, deriv_order::Integer)
    w = _wl_window(itp, xq); x = itp.x; y = itp.y
    acc = zero(eltype(y)) * one(float(xq))
    @inbounds for i in w
        if deriv_order == 0
            L = one(float(xq))
            for j in w
                j == i && continue
                L *= (xq - x[j]) / (x[i] - x[j])
            end
            acc += y[i] * L
        elseif deriv_order == 1
            dLi = zero(float(xq))
            for k in w
                k == i && continue
                term = one(float(xq)) / (x[i] - x[k])
                for j in w
                    (j == i || j == k) && continue
                    term *= (xq - x[j]) / (x[i] - x[j])
                end
                dLi += term
            end
            acc += y[i] * dLi
        elseif deriv_order == 2
            d2Li = zero(float(xq))
            for k in w
                k == i && continue
                for l in w
                    (l == i || l == k) && continue
                    term = one(float(xq)) / ((x[i] - x[k]) * (x[i] - x[l]))
                    for j in w
                        (j == i || j == k || j == l) && continue
                        term *= (xq - x[j]) / (x[i] - x[j])
                    end
                    d2Li += term
                end
            end
            acc += y[i] * d2Li
        else
            throw(ArgumentError("Interp1WL: deriv_order must be in 0:2, got $deriv_order"))
        end
    end
    return acc
end

(itp::Interp1WL)(xq::Real; deriv_order::Integer=0) = _wl_eval(itp, xq, deriv_order)

# First-form Lagrange basis over ALL nodes: L_i(xq) (order 0), L_i'(xq) (1), L_i''(xq) (2). Node-safe.
function _lagrange_all(nodes::AbstractVector, xq::Real, deriv_order::Integer)
    n = length(nodes); B = Vector{float(typeof(xq))}(undef, n)
    @inbounds for i in 1:n
        if deriv_order == 0
            p = one(float(xq))
            for j in 1:n
                j == i && continue
                p *= (xq - nodes[j]) / (nodes[i] - nodes[j])
            end
            B[i] = p
        elseif deriv_order == 1
            dLi = zero(float(xq))
            for k in 1:n
                k == i && continue
                term = one(float(xq)) / (nodes[i] - nodes[k])
                for j in 1:n
                    (j == i || j == k) && continue
                    term *= (xq - nodes[j]) / (nodes[i] - nodes[j])
                end
                dLi += term
            end
            B[i] = dLi
        elseif deriv_order == 2
            d2Li = zero(float(xq))
            for k in 1:n
                k == i && continue
                for l in 1:n
                    (l == i || l == k) && continue
                    term = one(float(xq)) / ((nodes[i] - nodes[k]) * (nodes[i] - nodes[l]))
                    for j in 1:n
                        (j == i || j == k || j == l) && continue
                        term *= (xq - nodes[j]) / (nodes[i] - nodes[j])
                    end
                    d2Li += term
                end
            end
            B[i] = d2Li
        else
            throw(ArgumentError("_lagrange_all: deriv_order must be in 0:2, got $deriv_order"))
        end
    end
    return B
end

"""
    Interp2Hermite{T}

Tensor-product *global* Lagrange interpolant over a full grid `x × y` (`z[i,j] = f(x[i], y[j])`),
matching Wolfram `Interpolation[grid, Method->"Hermite", InterpolationOrder->{All,All}]` (one degree-
`n−1` polynomial per axis through all nodes). Exact on its tensor bidegree. Call `itp(xq, yq)`; mixed
partials via [`partial`](@ref). Well-conditioned on the cosine-spaced Kerr grids it was built for.
"""
struct Interp2Hermite{T<:Real}
    x::Vector{T}
    y::Vector{T}
    z::Matrix{T}
end

"""
    interp2_hermite(grid) -> Interp2Hermite

Build an [`Interp2Hermite`](@ref) from a WL grid table `{{x, y}, z}` (full tensor product).
"""
function interp2_hermite(grid)
    xs = sort(unique(float(p[1][1]) for p in grid))
    ys = sort(unique(float(p[1][2]) for p in grid))
    nx, ny = length(xs), length(ys)
    nx * ny == length(grid) || throw(ArgumentError("interp2_hermite: grid is not a full tensor product"))
    xi = Dict(v => k for (k, v) in enumerate(xs)); yi = Dict(v => k for (k, v) in enumerate(ys))
    T = float(promote_type(eltype(xs), eltype(ys)))
    z = Matrix{T}(undef, nx, ny)
    @inbounds for p in grid
        z[xi[float(p[1][1])], yi[float(p[1][2])]] = p[2]
    end
    return Interp2Hermite{T}(collect(T, xs), collect(T, ys), z)
end

function (itp::Interp2Hermite)(xq::Real, yq::Real)
    Lx = _lagrange_all(itp.x, xq, 0); Ly = _lagrange_all(itp.y, yq, 0); z = itp.z
    acc = zero(eltype(z)) * one(float(xq))
    @inbounds for i in eachindex(Lx), j in eachindex(Ly)
        acc += Lx[i] * z[i, j] * Ly[j]
    end
    return acc
end

"""
    partial(itp::Interp2Hermite, xq, yq, ax, ay) -> value

Mixed partial `∂ₓ^ax ∂_y^ay` of the tensor-product interpolant (`ax, ay ∈ 0:2`), matching Wolfram
`Derivative[ax, ay][itp][x, y]`.
"""
function partial(itp::Interp2Hermite, xq::Real, yq::Real, ax::Integer, ay::Integer)
    Bx = _lagrange_all(itp.x, xq, ax); By = _lagrange_all(itp.y, yq, ay); z = itp.z
    acc = zero(eltype(z)) * one(float(xq))
    @inbounds for i in eachindex(Bx), j in eachindex(By)
        acc += Bx[i] * z[i, j] * By[j]
    end
    return acc
end

"""
    Interp3Lagrange{T}

Tensor-product *global* Lagrange interpolant over a full 3D grid `x × y × z` (`f[i,j,k]`), matching
Wolfram `Interpolation[grid, InterpolationOrder->All]` — the 3D analogue of [`Interp2Hermite`](@ref).
Exact on its tensor tridegree. Call `itp(xq, yq, zq)`.
"""
struct Interp3Lagrange{T<:Real}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    f::Array{T,3}
end

"""
    interp3_lagrange(grid) -> Interp3Lagrange

Build an [`Interp3Lagrange`](@ref) from a WL grid table `{{x, y, z}, f}` (full tensor product).
"""
function interp3_lagrange(grid)
    xs = sort(unique(float(p[1][1]) for p in grid))
    ys = sort(unique(float(p[1][2]) for p in grid))
    zs = sort(unique(float(p[1][3]) for p in grid))
    nx, ny, nz = length(xs), length(ys), length(zs)
    nx * ny * nz == length(grid) || throw(ArgumentError("interp3_lagrange: grid is not a full tensor product"))
    xi = Dict(v => k for (k, v) in enumerate(xs)); yi = Dict(v => k for (k, v) in enumerate(ys))
    zi = Dict(v => k for (k, v) in enumerate(zs))
    T = float(promote_type(eltype(xs), eltype(ys), eltype(zs)))
    f = Array{T,3}(undef, nx, ny, nz)
    @inbounds for p in grid
        f[xi[float(p[1][1])], yi[float(p[1][2])], zi[float(p[1][3])]] = p[2]
    end
    return Interp3Lagrange{T}(collect(T, xs), collect(T, ys), collect(T, zs), f)
end

function (itp::Interp3Lagrange)(xq::Real, yq::Real, zq::Real)
    Lx = _lagrange_all(itp.x, xq, 0); Ly = _lagrange_all(itp.y, yq, 0); Lz = _lagrange_all(itp.z, zq, 0)
    f = itp.f
    acc = zero(eltype(f)) * one(float(xq))
    @inbounds for k in eachindex(Lz), j in eachindex(Ly), i in eachindex(Lx)
        acc += Lx[i] * Ly[j] * Lz[k] * f[i, j, k]
    end
    return acc
end

# =====================================================================================================
# Log-domain Chebyshev-barycentric interpolants (KerrLorenzCirc LogChebInterp / BigLogRInterp)
# =====================================================================================================
# A single Chebyshev panel on a radial Teukolsky function R(r) fails when |R| spans many decades. The
# fix interpolates the COMPLEX LOG L=log R (phase-unwrapped along the ascending nodes), whose range is
# tiny, then reconstructs analytically: R=exp(L), R'=R·L', R''=R·(L''+L'²). L', L'' come from the
# Chebyshev spectral differentiation matrix (reused), so both interpolants are EXACT when log-magnitude
# is a polynomial (degree ≤ N in the Cheb variable) and phase is affine. Type-generic `{T}`.

# Unwrap the imaginary (phase) part of a complex log-sample vector along its (monotone) node order.
function _unwrap_phase!(L::AbstractVector{Complex{T}}) where {T}
    twoπ = 2 * T(π)
    @inbounds for k in 2:length(L)
        δ = imag(L[k]) - imag(L[k-1])
        m = round(δ / twoπ)
        L[k] -= Complex(zero(T), twoπ * m)
    end
    return L
end

"""
    LogChebInterp{T}

Log-domain Chebyshev-barycentric interpolant of a SMOOTH, NONZERO complex function `R(r)` over the
Cheb variable `r ∈ [lo,hi]`. Stores `L=log R`, `L'`, `L''` at the `N+1` Gauss–Lobatto nodes and
reconstructs `R=exp(L)`, `R'=R·L'`, `R''=R·(L''+L'²)`. Call `itp(r; deriv_order=0)` (`0:2`). Exact
when `log|R|` is a polynomial of degree ≤ `N` and the phase is affine. Build with [`log_cheb_interp`](@ref).
"""
struct LogChebInterp{T<:Real}
    lo::T
    hi::T
    x::Vector{T}
    w::Vector{T}
    L::Vector{Complex{T}}
    Lp::Vector{Complex{T}}
    Lpp::Vector{Complex{T}}
end

"""
    log_cheb_interp(f, lo, hi; N=30) -> LogChebInterp

Sample the smooth nonzero function `f` at `N+1` Gauss–Lobatto nodes on `[lo,hi]`, take the
phase-unwrapped complex log, and build the value/1st/2nd-derivative log-domain interpolant. `L'`,
`L''` are obtained from the Chebyshev differentiation matrix ([`cheb_diffmatrix`](@ref)).
"""
function log_cheb_interp(f, lo::T, hi::T; N::Integer=30) where {T<:Real}
    D, x = cheb_diffmatrix(T, N, lo, hi)
    R = Complex{T}[Complex{T}(f(xi)) for xi in x]
    L = log.(R)
    _unwrap_phase!(L)
    Lp = D * L
    Lpp = D * Lp
    return LogChebInterp{T}(lo, hi, x, bary_weights(x), L, Lp, Lpp)
end

function (ci::LogChebInterp{T})(r::Real; deriv_order::Integer=0) where {T}
    rr = T(r)
    Lr = bary_eval(ci.x, ci.w, ci.L, rr)
    Rv = exp(Lr)
    deriv_order == 0 && return Rv
    Lpr = bary_eval(ci.x, ci.w, ci.Lp, rr)
    deriv_order == 1 && return Rv * Lpr
    Lppr = bary_eval(ci.x, ci.w, ci.Lpp, rr)
    deriv_order == 2 && return Rv * (Lppr + Lpr * Lpr)
    throw(ArgumentError("LogChebInterp: deriv_order must be in 0:2, got $deriv_order"))
end

"""
    BigLogRInterp{T}

Log-domain Chebyshev-barycentric interpolant of a smooth nonzero `R(r)` sampled in LOG-`r` (Cheb
variable `u=log r`, nodes over `[log lo, log hi]`). Stores `L=log R(e^u)` and its `u`-derivatives;
reconstructs `R=exp(L)`, and via `r-`chain rule `R'=R·L_u/r`, `R''=R·(L_uu+L_u²−L_u)/r²`. Call
`itp(r; deriv_order=0)` (`0:2`). Exact when `log|R|` is a polynomial in `u` and phase affine in `u`.
The log-`r` abscissa suits functions varying over many `r`-decades. Build with [`biglogr_interp`](@ref).
"""
struct BigLogRInterp{T<:Real}
    lo::T
    hi::T
    u::Vector{T}
    w::Vector{T}
    L::Vector{Complex{T}}
    Lu::Vector{Complex{T}}
    Luu::Vector{Complex{T}}
end

"""
    biglogr_interp(f, lo, hi; N=30) -> BigLogRInterp

Sample `f` at `N+1` Gauss–Lobatto nodes in `u=log r` over `[log lo, log hi]`, take the phase-unwrapped
complex log, and build the log-`r` value/1st/2nd-derivative interpolant (BigFloat-inherited via `T`).
"""
function biglogr_interp(f, lo::T, hi::T; N::Integer=30) where {T<:Real}
    D, u = cheb_diffmatrix(T, N, log(lo), log(hi))
    R = Complex{T}[Complex{T}(f(exp(ui))) for ui in u]
    L = log.(R)
    _unwrap_phase!(L)
    Lu = D * L
    Luu = D * Lu
    return BigLogRInterp{T}(lo, hi, u, bary_weights(u), L, Lu, Luu)
end

function (ci::BigLogRInterp{T})(r::Real; deriv_order::Integer=0) where {T}
    rr = T(r); uq = log(rr)
    Lr = bary_eval(ci.u, ci.w, ci.L, uq)
    Rv = exp(Lr)
    deriv_order == 0 && return Rv
    Lur = bary_eval(ci.u, ci.w, ci.Lu, uq)
    deriv_order == 1 && return Rv * Lur / rr
    Luur = bary_eval(ci.u, ci.w, ci.Luu, uq)
    deriv_order == 2 && return Rv * (Luur + Lur * Lur - Lur) / (rr * rr)
    throw(ArgumentError("BigLogRInterp: deriv_order must be in 0:2, got $deriv_order"))
end
