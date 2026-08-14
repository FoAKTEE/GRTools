# HPHC — Horizon-Penetrating Hyperboloidal-Compactified spectral-collocation Teukolsky **evolution**
# backend (Move D). Port of the TeukEvolution.jl `linear_field` path: the 2nd-order-in-time
# Teukolsky equation for a spin-`s` field about Kerr, reduced to first order via `(f, p) = (ψ, ∂_T ψ)`
# and evolved by classical RK4 on the compactified radial grid `R = cl²/r ∈ [0, 1/r₊]` × the
# Gauss-Legendre spin-weighted angular collocation `y = −cosϑ` (Ripley, JLRipley314/TeukEvolution.jl;
# same RLGP-N chart as Move B, arXiv:2010.00162 §V — `R=cl²/r`, `y=−cosϑ`, `cl≡L`).
#
# PROVENANCE (all coefficient/step forms cited to the ref-code line):
#   - coefficient closed forms A_pp,A_pq,B_pp,B_pq,B_pf,pre : Evolution.jl:62-114
#   - RHS assembly (set_kp)                                 : Evolution.jl:151-181
#   - RK4 + post-step filters                               : Evolution.jl:186-303
#   - grid + dt formula                                     : TeukEvolution.jl:63-66
#   - Gaussian initial data (incl. ingoing pulse)          : Id.jl:32-100
# The `pre` factor is the PAPER-VERIFIED inverse of the ∂_T²ψ coefficient (RLGP_equations.md P1 =
# RLGP-N Eq. eq:Teukolsky_equation, the `[8M(2M−a²R/L²)(1+2MR/L²) − a²sin²ϑ]` bracket) — its presence
# on the ∂_T p RHS is mandatory for the `p ≡ ∂_T ψ` reduction. All forms are transcribed VERBATIM.
#
# REUSE (no numerics duplicated here):
#   - radial 4th-order FD `set_d1!/set_d2!` + 6th-order KO `ko_filter!`  ← NumericsCore.fd_uniform (D1)
#   - GL angular grids `Y_vals/cos_vals/sin_vals`, the spin-weighted `swal`, and the collocation
#     `swal_laplacian_matrix` / `swal_killtop_matrix` / `angular_matrix_mult!`  ← Harmonics.SphereCollocation (D2)
#   - the RK4 integration loop `rk4!`  ← NumericsCore.dg  (the upstream hand-unrolled RK4 is the
#     SAME classical (1/6,1/3,1/3,1/6) scheme; the fixture gate G-D3a protects the equivalence).
#
# Everything is `{T}`-generic (element type flows from the `T` keyword), so `Complex{BigFloat}` fields
# evolve at inherited precision. Characterized end-to-end against the ORIGINAL TeukEvolution run in
# test/fixtures/td_hphc_upstream.jl (config a=0.7, s=−2, m=2, nx=48, ny=24, nt=600).

# ---------------------------------------------------------------------------------------------------
# (1) coefficient closed forms — Evolution.jl:62-114 (transcribed verbatim; {T}-generic)
# ---------------------------------------------------------------------------------------------------

"""
    hphc_A_pp(R, bhm, bhs, cl)

Radial coefficient of `∂_R p` in the reduced ψ-momentum RHS (Evolution.jl:86-91).
"""
@inline hphc_A_pp(R, bhm, bhs, cl) =
    (2 / (cl^4)) * ((cl^6) + (cl^2) * (bhs^2 - 8 * (bhm^2)) * (R^2) + 4 * (bhs^2) * bhm * (R^3))

"""
    hphc_A_pq(R, bhm, bhs, cl)

Radial coefficient of `∂_R² f` (∝ R²·Δ_Kerr/cl⁴, vanishes at the horizon; Evolution.jl:93-94).
"""
@inline hphc_A_pq(R, bhm, bhs, cl) =
    ((R^2) / (cl^4)) * ((cl^4) - 2 * (cl^2) * bhm * R + (bhs^2) * (R^2))

"""
    hphc_B_pp(R, cy, bhm, bhs, cl, mv, spin)

Coefficient of `p` (carries the `im·a·m` azimuthal and `im·a·s·cosϑ` terms; Evolution.jl:96-101).
`cy = cosϑ`.
"""
@inline hphc_B_pp(R, cy, bhm, bhs, cl, mv, spin) = -(
    2 * im * bhs * mv * ((cl^2) + 4 * bhm * R) / (cl^2) -
    2 * (bhs^2) * (cl^2 + 6 * bhm * R) * R / (cl^4) - 4 * bhm * spin +
    8 * (bhm^2) * (2 + spin) * R / (cl^2) +
    2 * im * bhs * spin * cy
)

"""
    hphc_B_pq(R, bhm, bhs, cl, mv, spin)

Coefficient of `∂_R f` (Evolution.jl:103-107).
"""
@inline hphc_B_pq(R, bhm, bhs, cl, mv, spin) =
    (2 * R / (cl^4)) * (
        2 * (bhs^2) * (R^2) + (1 + spin) * (cl^4) -
        (im * bhs * mv + (3 + spin) * bhm) * ((cl^2) * R)
    )

"""
    hphc_B_pf(R, bhm, bhs, cl, mv, spin)

Coefficient of `f` (Evolution.jl:109-111).
"""
@inline hphc_B_pf(R, bhm, bhs, cl, mv, spin) =
    -(2 * R / (cl^4)) * (im * bhs * (cl^2) * mv - (bhs^2) * R + (cl^2) * bhm * (1 + spin))

"""
    hphc_pre(R, sy, bhm, bhs, cl)

The PAPER-VERIFIED inverse of the `∂_T²ψ` coefficient (Evolution.jl:67-73). `1/pre` equals the
RLGP-N master-equation bracket `[8M(2M−a²R/L²)(1+2MR/L²) − a²sin²ϑ]` (RLGP_equations.md P1); the
`p ≡ ∂_T ψ` reduction REQUIRES this factor on the momentum RHS. `sy = sinϑ`.
"""
@inline hphc_pre(R, sy, bhm, bhs, cl) =
    1 / (8 * bhm * (2 * (cl^2) * bhm - (bhs^2) * R) * ((cl^2) + 2 * bhm * R) / (cl^4) - (bhs * sy)^2)

# ---------------------------------------------------------------------------------------------------
# (2) Gaussian initial data — Id.jl:32-100 (incl. the ingoing-pulse branch :79-87)
# ---------------------------------------------------------------------------------------------------

# In-place Gaussian ID for `(f, p) = (ψ, ∂_Tψ)`: a Gaussian in the areal radius `r = cl²/R` times the
# spin-weighted harmonic `swal(spin, m, l_ang, y)`, pinned to 0 at scri (`i=1`, R=0). Optional ingoing
# pulse `p += −R²/(2+4R)·∂_R f`. Finally rescaled by `amp/max|f|`. Transcribes Id.set_gaussian!.
function _hphc_set_gaussian!(f::AbstractMatrix{Complex{T}}, p::AbstractMatrix{Complex{T}},
        spin::Integer, mv::Integer, l_ang::Integer, ru::T, rl::T, width::T, amp::Complex{T},
        cl::T, Rv::AbstractVector{T}, Yv::AbstractVector{T}, ingoing::Bool, dr::T,
        nx::Integer, ny::Integer) where {T<:AbstractFloat}
    max_val = zero(T)
    @inbounds for j in 1:ny
        f[1, j] = 0
        p[1, j] = 0
        for i in 2:nx
            r = (cl^2) / Rv[i]
            g = exp(-(r - (rl + ru) / 2)^2 / width^2)
            fv = g * swal(spin, mv, l_ang, Yv[j])
            f[i, j] = fv
            p[i, j] = 0
            max_val = max(abs(fv), max_val)
        end
    end
    if ingoing
        f_rd1 = zeros(Complex{T}, nx, ny)
        set_d1!(f_rd1, f, dr)
        @inbounds for j in 1:ny, i in 1:nx
            p[i, j] += -Rv[i] * Rv[i] / (2 + 4 * Rv[i]) * f_rd1[i, j]
        end
    end
    @inbounds for j in 1:ny, i in 1:nx
        f[i, j] *= amp / max_val
        p[i, j] *= amp / max_val
    end
    return nothing
end

# ---------------------------------------------------------------------------------------------------
# (3) problem assembly — grid (TeukEvolution.jl:63-66) + coefficient arrays + collocation matrices + ID
# ---------------------------------------------------------------------------------------------------

"""
    HPHCTeukolskyProblem{T}

Precomputed single-`m` HPHC Teukolsky evolution problem (Move D): the compactified radial grid `Rv`,
the Gauss-Legendre angular grid `Yv`, the `Evo_lin_f` operator coefficient arrays
(`A_pp,A_pq,B_pp,B_pq,B_pf,pre`), the angular Laplacian / killtop-2 collocation matrices, and the
Gaussian initial data `(f0, p0)`. Build with [`HPHCTeukolskyProblem(; …)`](@ref); evolve with
[`solve_hphc`](@ref). Multi-`m` runs are embarrassingly parallel — one problem per `m`
([`solve_hphc_mmodes`](@ref)).
"""
struct HPHCTeukolskyProblem{T<:AbstractFloat}
    nx::Int
    ny::Int
    mv::Int
    spin::Int
    bhm::T
    bhs::T
    cl::T
    dr::T
    dt::T
    Rv::Vector{T}
    Yv::Vector{T}
    A_pp::Matrix{T}
    A_pq::Matrix{T}
    B_pp::Matrix{Complex{T}}
    B_pq::Matrix{Complex{T}}
    B_pf::Matrix{Complex{T}}
    pre::Matrix{T}
    laplM::Matrix{T}
    fltrM::Matrix{T}
    f0::Matrix{Complex{T}}
    p0::Matrix{Complex{T}}
end

"""
    HPHCTeukolskyProblem(; nx, ny, mv, spin=-2, bhm=1, bhs=0.7, cl=1, cfl=0.5,
                           id_l_ang=2, id_ru=4.0, id_rl=2.5, id_width=1.5,
                           id_amp=complex(0.1,0.0), ingoing=true, T=Float64)

Assemble the HPHC problem. The grid follows TeukEvolution EXACTLY (incl. the `cl` conventions):
`r₊ = M(1+√(1−(a/M)²))`, `maxR = 1/r₊`, `dr = maxR/(nx−1)`, `Rv[i] = dr·(i−1)` (so `R=0` is scri⁺
and `R=maxR` the horizon), and `dt = min(cfl·dr·M², 6/ny²)`. The angular grid, `swal` values, and the
Laplacian / killtop-2 matrices come from `Harmonics.SphereCollocation`; the coefficients from the
verbatim closed forms above.
"""
function HPHCTeukolskyProblem(; nx::Integer, ny::Integer, mv::Integer, spin::Integer=-2,
        bhm=1.0, bhs=0.7, cl=1.0, cfl=0.5,
        id_l_ang::Integer=2, id_ru=4.0, id_rl=2.5, id_width=1.5,
        id_amp=complex(0.1, 0.0), ingoing::Bool=true,
        T::Type{<:AbstractFloat}=Float64)
    bhmT, bhsT, clT, cflT = T(bhm), T(bhs), T(cl), T(cfl)
    minr = bhmT * (1 + sqrt(1 - (bhsT / bhmT)^2))   # r₊ (uncompactified horizon), TeukEvolution.jl:63
    maxR = 1 / minr                                  # dt is cl-independent by design, :64
    dr = maxR / (nx - 1)                             # :65
    dt = min(cflT * dr * bhmT^2, 6 / T(ny)^2)        # :66
    Rv = T[dr * (i - 1) for i in 1:nx]               # R=0 scri⁺ … R=maxR horizon
    Yv = Y_vals(T, ny)                               # GL roots, y = −cosϑ
    Cv = cos_vals(T, ny)                             # cosϑ = −y
    Sv = sin_vals(T, ny)                             # sinϑ
    A_pp = Matrix{T}(undef, nx, ny)
    A_pq = Matrix{T}(undef, nx, ny)
    B_pp = Matrix{Complex{T}}(undef, nx, ny)
    B_pq = Matrix{Complex{T}}(undef, nx, ny)
    B_pf = Matrix{Complex{T}}(undef, nx, ny)
    pre = Matrix{T}(undef, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        R, cy, sy = Rv[i], Cv[j], Sv[j]
        A_pp[i, j] = hphc_A_pp(R, bhmT, bhsT, clT)
        A_pq[i, j] = hphc_A_pq(R, bhmT, bhsT, clT)
        B_pp[i, j] = hphc_B_pp(R, cy, bhmT, bhsT, clT, mv, spin)
        B_pq[i, j] = hphc_B_pq(R, bhmT, bhsT, clT, mv, spin)
        B_pf[i, j] = hphc_B_pf(R, bhmT, bhsT, clT, mv, spin)
        pre[i, j] = hphc_pre(R, sy, bhmT, bhsT, clT)
    end
    laplM = swal_laplacian_matrix(T, ny, spin, mv)          # S_lapl, Evolution.jl:77
    fltrM = swal_killtop_matrix(T, ny, spin, mv, 2)         # S_fltr (killtop N=2), :78
    f0 = zeros(Complex{T}, nx, ny)
    p0 = zeros(Complex{T}, nx, ny)
    _hphc_set_gaussian!(f0, p0, spin, mv, id_l_ang, T(id_ru), T(id_rl), T(id_width),
                        Complex{T}(id_amp), clT, Rv, Yv, ingoing, dr, nx, ny)
    return HPHCTeukolskyProblem{T}(nx, ny, mv, spin, bhmT, bhsT, clT, dr, dt,
        Rv, Yv, A_pp, A_pq, B_pp, B_pq, B_pf, pre, laplM, fltrM, f0, p0)
end

# ---------------------------------------------------------------------------------------------------
# (4) evolution — RK4 (NumericsCore.rk4!) with the upstream-ordered post-step filters
# ---------------------------------------------------------------------------------------------------

"""
    solve_hphc(prob; nt, ts) -> (times, f, p)

Evolve the HPHC problem for `nt` RK4 steps of size `prob.dt`, recording snapshots every `ts` steps
(and the initial state). Returns `times` and vectors of `nx×ny` field matrices `f` (= ψ) and `p`
(= ∂_Tψ), one per recorded snapshot (snapshot `k` is at step `(k−1)·ts`).

Design: the state `(f, p)` is packed into one `Vector{Complex{T}}` of length `2·nx·ny`
(`reshape → (nx,ny,2)`), and `NumericsCore.rk4!` drives the classical (1/6,1/3,1/3,1/6) RK4 through
an `rhs!(du,u,t)` closure that writes `∂_T f = p` and `∂_T p = set_kp(…)` (Evolution.jl:169 assembly).
The per-step post-processing (Evolution.jl:291-300) — `ko_filter!(f)`, `ko_filter!(p)` at ε=0.5, then
the killtop angular multiply on both — is applied in the `rk4!` callback AFTER each step (never on the
initial record), reproducing the upstream `Evolve_lin_f!` step exactly. Snapshots are keyed by integer
step count, so the recorded states are the post-filtered fields after 0/ts/2ts/… uniform-`dt` steps.
"""
function solve_hphc(prob::HPHCTeukolskyProblem{T}; nt::Integer, ts::Integer) where {T}
    nx, ny = prob.nx, prob.ny
    dr, dt = prob.dr, prob.dt
    u = Vector{Complex{T}}(undef, 2 * nx * ny)
    U = reshape(u, nx, ny, 2)
    @views U[:, :, 1] .= prob.f0
    @views U[:, :, 2] .= prob.p0
    # rhs! scratch
    f_rd1 = zeros(Complex{T}, nx, ny)
    f_rd2 = zeros(Complex{T}, nx, ny)
    p_rd1 = zeros(Complex{T}, nx, ny)
    f_lap = zeros(Complex{T}, nx, ny)
    # filter scratch
    tmpf = zeros(Complex{T}, nx, ny)
    tmpp = zeros(Complex{T}, nx, ny)
    A_pp, A_pq = prob.A_pp, prob.A_pq
    B_pp, B_pq, B_pf = prob.B_pp, prob.B_pq, prob.B_pf
    pre, laplM, fltrM = prob.pre, prob.laplM, prob.fltrM

    # RHS: ∂_T f = p ; ∂_T p = pre·(A_pp·∂_R p + A_pq·∂_R²f + B_pp·p + B_pq·∂_R f + B_pf·f + Δ_S f)
    function rhs!(du, uu, ::Real)
        UU = reshape(uu, nx, ny, 2)
        dUU = reshape(du, nx, ny, 2)
        ff = @view UU[:, :, 1];  pp = @view UU[:, :, 2]
        df = @view dUU[:, :, 1]; dp = @view dUU[:, :, 2]
        set_d1!(f_rd1, ff, dr)
        set_d1!(p_rd1, pp, dr)
        set_d2!(f_rd2, ff, dr)
        angular_matrix_mult!(f_lap, ff, laplM)
        @inbounds for j in 1:ny, i in 1:nx
            df[i, j] = pp[i, j]
            dp[i, j] = pre[i, j] * (
                A_pp[i, j] * p_rd1[i, j] + A_pq[i, j] * f_rd2[i, j] +
                B_pp[i, j] * pp[i, j] + B_pq[i, j] * f_rd1[i, j] +
                B_pf[i, j] * ff[i, j] + f_lap[i, j])
        end
        return du
    end

    # post-step: KO(ε=0.5) on f then p, then killtop angular multiply on both (Evolution.jl:291-300)
    epsKO = T(1) / 2
    function postfilter!(uu)
        UU = reshape(uu, nx, ny, 2)
        ff = @view UU[:, :, 1]; pp = @view UU[:, :, 2]
        ko_filter!(ff, tmpf, epsKO)
        ko_filter!(pp, tmpp, epsKO)
        tmpf .= ff; angular_matrix_mult!(ff, tmpf, fltrM)
        tmpp .= pp; angular_matrix_mult!(pp, tmpp, fltrM)
        return nothing
    end

    snaps_f = Matrix{Complex{T}}[]
    snaps_p = Matrix{Complex{T}}[]
    snaps_t = T[]
    function cb(uu, t, step)
        step >= 1 && postfilter!(uu)         # never filter the initial record
        if step % ts == 0
            UU = reshape(uu, nx, ny, 2)
            push!(snaps_f, copy(@view UU[:, :, 1]))
            push!(snaps_p, copy(@view UU[:, :, 2]))
            push!(snaps_t, T(t))
        end
        return nothing
    end

    rk4!(rhs!, u, (zero(T), nt * dt), dt; callback=cb)
    return (times=snaps_t, f=snaps_f, p=snaps_p)
end

"""
    solve_hphc_mmodes(m_vals; nt, ts, kwargs...) -> Dict(m => (times,f,p))

Build one [`HPHCTeukolskyProblem`](@ref) per `m ∈ m_vals` (sharing all keyword construction args) and
evolve them in parallel (`Threads.@threads` over the embarrassingly-parallel `m` axis — each `m` is an
independent Teukolsky solve with its own scratch). Returns a `Dict` keyed by `m`.
"""
function solve_hphc_mmodes(m_vals; nt::Integer, ts::Integer, kwargs...)
    ms = collect(m_vals)
    probs = [HPHCTeukolskyProblem(; mv=m, kwargs...) for m in ms]
    results = Vector{Any}(undef, length(ms))
    Threads.@threads for k in eachindex(ms)
        results[k] = solve_hphc(probs[k]; nt=nt, ts=ts)
    end
    return Dict(ms[k] => results[k] for k in eachindex(ms))
end
