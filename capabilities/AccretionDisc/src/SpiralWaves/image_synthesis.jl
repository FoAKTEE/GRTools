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
