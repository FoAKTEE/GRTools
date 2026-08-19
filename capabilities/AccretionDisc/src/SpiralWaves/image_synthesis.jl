# SpiralWaves/image_synthesis.jl — the mode-sum reconstruction of the spiral wake of arXiv:2608.17003
# (Figs. 1-2, main.tex:853-... ): with the Fourier convention (45a)
#
#     f(t, r, φ) = Σ_{m=-∞}^{∞} f̂_m(r) e^{-iω t + i m φ},   ω = mΩ•,
#
# and a REAL field (f̂_{−m} = f̂_m^*, ω_{−m} = −ω_m), every quadratic-free reconstruction collapses onto
# the co-rotating angle Ψ = φ − Ω• t:
#
#     f(r, Ψ) = f̂_0(r) + 2 Σ_{m≥1} Re[ f̂_m(r) e^{imΨ} ] ,
#
# i.e. the SAME ×2 convention (m = 0 unweighted, m ≥ 1 doubled) that the quadratic identity (73)
# carries. `solve_mode` produces m ≥ 1 only (the forced problem has no m = 0 response), so the m = 0
# term is present in the API but empty in practice.
#
# Different m are solved on DIFFERENT adaptive grids (the domain (C1)-(C3) is per-mode), so the modes
# are put on a common radial grid by the NumericsCore cubic spline (per mode, real and imaginary parts
# separately) and set to ZERO outside their own domain. For the wake IMAGES use grids built with
# `extent = :fiducial` (r_ISCO … 5a•, App. C (C1)) so that every mode shares one domain — that is what
# Figs. 1-2 show, with the spiral truncated at the ISCO. Type-generic `{T}`; BigFloat inherits.
# Additive (EGD::sw-image-synthesis).

using NumericsCore: CubicSplineInterpolant

"""
    ModeInterp{T}

Cubic-spline interpolant of ONE mode amplitude f̂_m(r) (real and imaginary parts) together with its
radial domain `[r_lo, r_hi]` and its azimuthal number `m`. Built by [`mode_interp`](@ref); evaluate
with `mi(r)` (zero outside the domain).
"""
struct ModeInterp{T<:Real}
    m::Int
    r_lo::T
    r_hi::T
    re::CubicSplineInterpolant{T}
    im::CubicSplineInterpolant{T}
end

(mi::ModeInterp{T})(r) where {T} =
    (r < mi.r_lo || r > mi.r_hi) ? zero(Complex{T}) : Complex{T}(mi.re(T(r), 1), mi.im(T(r), 1))

"""
    mode_interp(sol; component = 1) -> ModeInterp

Interpolant of one solved mode's amplitude: `component = 1` gives δÛ_m = δp̂/(ε+p) (48) — the
"enthalpy perturbation" of Figs. 1-2 — and `component = 2` gives the master variable −iδû^r (49).
`sol` may be a [`ModeSolution`](@ref) or any object carrying `m`, `r` and `y` (a 2×N matrix).
Exact (interpolatory) at the mode's own nodes.
"""
function mode_interp(sol; component::Integer = 1)
    component in (1, 2) || throw(ArgumentError("mode_interp: component must be 1 (δÛ) or 2 (−iδû^r)"))
    r = sol.r
    T = eltype(r)
    vals = view(sol.y, component, :)
    ModeInterp{T}(Int(sol.m), first(r), last(r),
                  CubicSplineInterpolant(r, T[real(v) for v in vals]),
                  CubicSplineInterpolant(r, T[imag(v) for v in vals]))
end

"""
    synthesize(sols, rq, Ψq; component = 1) -> Matrix

Reconstruct the real perturbation field on the polar grid `rq` × `Ψq` (Ψ = φ − Ω•t, the co-rotating
azimuth) from a collection of solved modes:

    f(r, Ψ) = Σ_m w_m Re[ f̂_m(r) e^{imΨ} ],    w_0 = 1, w_{m≥1} = 2,

returning a `length(rq) × length(Ψq)` matrix. Each mode contributes only inside its own radial domain.
`sols` may be `ModeSolution`s or anything with `m`, `r`, `y`. Duplicate `m` values are summed (so a
resolution study can be superposed), and the modes are accumulated in the order given, which makes the
result bitwise reproducible.
"""
function synthesize(sols, rq::AbstractVector, Ψq::AbstractVector; component::Integer = 1)
    isempty(sols) && throw(ArgumentError("synthesize: no modes given"))
    interps = [mode_interp(s; component = component) for s in sols]
    T = float(promote_type(eltype(rq), eltype(Ψq), typeof(interps[1].r_lo)))
    out = zeros(T, length(rq), length(Ψq))
    for mi in interps
        w = mi.m == 0 ? one(T) : T(2)
        for (i, r) in enumerate(rq)
            f = mi(r)
            iszero(f) && continue
            for (j, Ψ) in enumerate(Ψq)
                out[i, j] += w * real(f * cis(mi.m * T(Ψ)))
            end
        end
    end
    out
end

"""
    polar_to_cartesian(rq, Ψq) -> (X, Y)

Cartesian node coordinates of the polar grid used by [`synthesize`](@ref): X[i,j] = r_i cos Ψ_j,
Y[i,j] = r_i sin Ψ_j — the (x, y) plane of Figs. 1-2, with the small object on the +x axis (Ψ = 0).
"""
function polar_to_cartesian(rq::AbstractVector, Ψq::AbstractVector)
    T = float(promote_type(eltype(rq), eltype(Ψq)))
    X = Matrix{T}(undef, length(rq), length(Ψq))
    Y = similar(X)
    for (i, r) in enumerate(rq), (j, Ψ) in enumerate(Ψq)
        X[i, j] = T(r) * cos(T(Ψ)); Y[i, j] = T(r) * sin(T(Ψ))
    end
    (X, Y)
end

"""
    image_grid(bg::DiscBackground; nr = 512, nΨ = 1024, r_lo = r_isco(bg), r_hi = 5bg.aEMRI) -> (rq, Ψq)

The fiducial image grid of Figs. 1-2: `nr` radii, uniform on [r_lo, r_hi] (the fiducial extent (C1)
r_ISCO … 5a•), and `nΨ` azimuths uniform on [0, 2π) (periodic, endpoint excluded).
"""
function image_grid(bg::DiscBackground{T}; nr::Integer = 512, nΨ::Integer = 1024,
                    r_lo = r_isco(bg), r_hi = 5 * bg.aEMRI) where {T}
    (nr ≥ 2 && nΨ ≥ 2) || throw(ArgumentError("image_grid: need nr, nΨ ≥ 2"))
    rq = collect(range(T(r_lo), T(r_hi); length = nr))
    Ψq = collect(range(zero(T), 2 * T(π); length = nΨ + 1)[1:nΨ])
    (rq, Ψq)
end

"""
    synthesize_cartesian(sols, xs, ys; component = 1, r_lo = nothing, r_hi = nothing,
                         blank = NaN) -> Matrix

Reconstruct the real perturbation field of [`synthesize`](@ref) directly on the CARTESIAN grid
`xs` × `ys` of Figs. 1-2 — the (r cos φ, r sin φ) plane, with the small object on the +x axis:

    F[i, j] = Σ_m w_m Re[ f̂_m(r_ij) e^{i m Ψ_ij} ],   r_ij = √(x_i² + y_j²), Ψ_ij = atan(y_j, x_i),

so the map is EXACT at every pixel (each mode is evaluated through its own cubic spline at the
pixel's radius — nothing is interpolated a second time in angle or in x, y). Returns a
`length(xs) × length(ys)` matrix, first index = x, matching [`polar_to_cartesian`](@ref).

Pixels outside `[r_lo, r_hi]` — by default the union of the modes' own radial domains, i.e. the
plotted disc r_ISCO … 5a• of (C1) — are set to `blank` (`NaN`, which plotting back-ends leave
transparent: the white region inside the ISCO and outside the disc edge of Figs. 1-2). Modes are
accumulated in the order given, so the result is bitwise reproducible.
"""
function synthesize_cartesian(sols, xs::AbstractVector, ys::AbstractVector; component::Integer = 1,
                              r_lo = nothing, r_hi = nothing, blank = NaN)
    isempty(sols) && throw(ArgumentError("synthesize_cartesian: no modes given"))
    interps = [mode_interp(s; component = component) for s in sols]
    T = float(promote_type(eltype(xs), eltype(ys), typeof(interps[1].r_lo)))
    rlo = r_lo === nothing ? minimum(mi.r_lo for mi in interps) : T(r_lo)
    rhi = r_hi === nothing ? maximum(mi.r_hi for mi in interps) : T(r_hi)
    rlo ≤ rhi || throw(ArgumentError("synthesize_cartesian: need r_lo ≤ r_hi, got ($(rlo), $(rhi))"))
    ws = T[mi.m == 0 ? one(T) : T(2) for mi in interps]
    out = fill(T(blank), length(xs), length(ys))
    for (j, y) in enumerate(ys)
        yT = T(y)
        for (i, x) in enumerate(xs)
            xT = T(x)
            r = hypot(xT, yT)
            (r < rlo || r > rhi) && continue
            Ψ = atan(yT, xT)
            acc = zero(T)
            for (k, mi) in enumerate(interps)
                f = mi(r)
                iszero(f) && continue
                acc += ws[k] * real(f * cis(mi.m * Ψ))
            end
            out[i, j] = acc
        end
    end
    out
end

"""
    arm_phase(F, Ψq; period = 2π) -> Vector

Winding curve of the spiral arm: the azimuth Ψ of the field CREST at each radius of a wake image
`F` (`length(rq) × length(Ψq)`, from [`synthesize`](@ref) on a uniform periodic `Ψq`), refined to
sub-cell accuracy by a parabolic fit through the three samples around the row maximum and then
UNWRAPPED in radius (each value is shifted by the multiple of `period` nearest to its predecessor).

`period` is the azimuthal period of the PATTERN, i.e. the ambiguity of "the crest": `2π` for the
one-armed wake of Figs. 1-2 (the default, where the mode sum leaves a single dominant arm), `2π/m`
for a pure single-`m` field, whose m crests are indistinguishable.

The result is the quantitative form of the caption statements of Fig. 1: the arm winds more
tightly when |dΨ/dr| is larger, and `(Ψ[i] − Ψ[j])/2π` is the number of turns between two radii.
Rows that are entirely non-finite or constant give the plain grid maximum.
"""
function arm_phase(F::AbstractMatrix, Ψq::AbstractVector; period = 2π)
    n, nΨ = size(F)
    nΨ == length(Ψq) || throw(DimensionMismatch("arm_phase: size(F, 2) ≠ length(Ψq)"))
    nΨ ≥ 3 || throw(ArgumentError("arm_phase: need at least 3 azimuths"))
    T = float(promote_type(eltype(F), eltype(Ψq), typeof(period)))
    dΨ = T(Ψq[2]) - T(Ψq[1])
    per = T(period)
    per > 0 || throw(ArgumentError("arm_phase: need period > 0"))
    out = Vector{T}(undef, n)
    prev = T(NaN)
    for i in 1:n
        row = view(F, i, :)
        j = argmax(row)
        f0, fm, fp = T(row[j]), T(row[mod1(j - 1, nΨ)]), T(row[mod1(j + 1, nΨ)])
        den = fm - 2 * f0 + fp
        δ = isfinite(den) && !iszero(den) ? (fm - fp) / (2 * den) : zero(T)
        abs(δ) ≤ 1 || (δ = zero(T))                       # parabola not peaked inside the cell
        ψ = T(Ψq[j]) + δ * dΨ
        isnan(prev) || (ψ += per * round((prev - ψ) / per))
        out[i] = ψ
        prev = ψ
    end
    out
end

"""
    wake_modes(bg::DiscBackground, ms; grid_kwargs = (;), solve_kwargs = (;), threads = true,
               verbose = false) -> Vector{ModeSolution}

Solve every azimuthal mode in `ms` for the wake IMAGES of Figs. 1-2 and return the solutions in
the order of `ms` (the mode sum of [`synthesize`](@ref) then needs no further bookkeeping).

Unlike [`torque_survey`](@ref) — which reduces each mode to four scalars and throws the solution
away — the full `ModeSolution`s are kept, and the grids default to `extent = :fiducial`, the
image domain r_ISCO … 5a• of (C1): every mode then shares one domain, so the synthesized map is
continuous and the spiral is truncated at the ISCO exactly as in the published figures.
`grid_kwargs` (which may override `extent`) is forwarded to [`build_grid`](@ref) and
`solve_kwargs` to [`solve_mode`](@ref).

`threads = true` distributes the modes over the available threads, LONGEST FIRST (the cost of a
mode grows with `m` because its adaptive grid does, so the descending-`m` order keeps the tail of the
sweep from serialising). Each task owns everything it touches and writes one disjoint slot, and the
mode sum of [`synthesize`](@ref) is taken afterwards in the fixed order of `ms`, so the result is
BITWISE identical to `threads = false`. `verbose = true` prints the [`compute_budget`](@ref) note.
"""
function wake_modes(bg::DiscBackground{T}, ms::AbstractVector{<:Integer};
                    grid_kwargs = NamedTuple(), solve_kwargs = NamedTuple(),
                    threads::Bool = true, verbose::Bool = false) where {T}
    isempty(ms) && throw(ArgumentError("wake_modes: `ms` must be non-empty"))
    all(m -> m ≥ 1, ms) || throw(ArgumentError("wake_modes: every m must be ≥ 1"))
    gk = merge((extent = :fiducial,), grid_kwargs)
    verbose && @info compute_budget().note
    n = length(ms)
    out = Vector{ModeSolution{T}}(undef, n)
    if threads
        perm = sortperm(collect(ms); rev = true)          # longest job first (cost grows with m)
        Threads.@threads :dynamic for t in 1:n
            i = perm[t]
            @inbounds out[i] = solve_mode(bg, ms[i]; grid = build_grid(bg, ms[i]; gk...),
                                          threads = false, solve_kwargs...)
        end
    else
        for i in 1:n
            @inbounds out[i] = solve_mode(bg, ms[i]; grid = build_grid(bg, ms[i]; gk...),
                                          threads = false, solve_kwargs...)
        end
    end
    out
end
