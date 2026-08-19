# SpiralWaves/solver.jl — the implicit-midpoint global solve of arXiv:2608.17003 Sec. II.E
# (main.tex:737-782): the master system (47) dy/dr = C·y + S is discretised on the adaptive grid of
# grid.jl with the second-order implicit midpoint (Crank–Nicolson) rule (64)
#
#     (I − Δr_i/2 C_{i+1/2})·y_{i+1} − (I + Δr_i/2 C_{i+1/2})·y_i = Δr_i S_{i+1/2},   Δr_i = r_{i+1} − r_i,
#
# assembled for ALL cells at once into a 2N×2N `SparseMatrixCSC{Complex{T}}` (block bidiagonal, lower
# bandwidth 2 and upper bandwidth 2) and closed by ONE boundary row at each edge (66)-(68): the
# eigen-decomposition C·Y = kY of the edge matrix selects the wave travelling AWAY from the orbit and
# the Neumann-type row y·Y_⊥ = 0 forces y ∥ Y. The paper solves this system globally (M:737-741:
# single and multiple shooting were tested and rejected as less robust — the transfer matrix across
# the corotation block is exponentially ill-conditioned).
#
# `C_{i+1/2}` is evaluated AT THE MIDPOINT RADIUS r_{i+1/2} = (r_i + r_{i+1})/2 (the decomposition
# records this as the only reading of "C_{i+1/2}" consistent with the second-order local error; using
# (C_i + C_{i+1})/2 instead changes the scheme at O(Δr²), i.e. within its own truncation error).
#
# ────────────────────────────────────────────────────────────────────────────────────────────────
# BOUNDARY-CONDITION BRANCH — ERRATUM E-M4 [FINDING, this node]
# ────────────────────────────────────────────────────────────────────────────────────────────────
# M:774-776 says: "choose the eigenvector for which Im[k] < 0. This eigenvector serves as the outgoing
# solution". With the paper's OWN Fourier convention (45a) f = Σ_m f̂_m e^{−iωt+imφ} and y ∝ e^{kr},
# a propagating WKB branch has k = i k_r with the radial wavenumber k_r of the local dispersion
# relation ω_d² = κ² + c_s² k_r² (grid.jl `wkb_wavenumber2`), whose group velocity is
#
#     v_g = ∂ω/∂k_r|_m = c_s² k_r / ω_d .
#
# Outside corotation ω_d > 0 and "away from the orbit" means v_g > 0 ⇒ k_r > 0; inside corotation
# ω_d < 0 and "away" means v_g < 0 ⇒ again k_r > 0. So the branch that carries energy and angular
# momentum away from the orbit has k_r > 0, i.e. Im[k] > 0 at BOTH edges — the OPPOSITE sign to the
# printed criterion. The same conclusion follows from the spiral geometry: k_r > 0 gives
# dφ/dr|_phase = −k_r/m < 0, a TRAILING arm, which is what Fig. 1/2 of the paper show and what
# Goldreich–Tremaine theory requires. The printed Im[k] < 0 is therefore read here as a sign slip (or
# an implicit e^{+iωt−imφ} convention in the authors' code).
#
# `solve_mode(...; bc = :outgoing)` (DEFAULT) selects by the PHYSICAL criterion — the sign of the
# group velocity, evaluated from the actual eigenvalues and the actual ω_d, never hard-coded — and in
# an evanescent edge region (both eigenvalues real, k² < 0) it selects the branch that DECAYS away
# from the orbit instead. `bc = :paper` reproduces the printed Im[k] < 0 literally, for the erratum
# comparison. `bc = :custom` takes caller-supplied rows (used by the manufactured-solution tests).
#
# Linear algebra: for `Complex{Float64}` the sparse system is solved by `\` (UMFPACK, `linsolve =
# :sparse`); for every other element type (BigFloat!) SparseArrays has no generic LU, so the SAME
# assembled matrix is copied into LAPACK-style band storage and factorised by a generic partial-pivot
# banded LU (`linsolve = :banded`, O(N) work and memory). `linsolve = :auto` (default) picks the first
# for `Complex{Float64}` and the second otherwise; the two agree to roundoff (tested).
# Type-generic `{T}`; BigFloat precision is inherited. Additive (EGD::sw-mode-solver).

# ---------------------------------------------------------------------------------------------
# generic banded LU (LAPACK `gbtrf`/`gbtrs` layout, partial pivoting) — the BigFloat path
# ---------------------------------------------------------------------------------------------
# [FUTURE] this is a general numerical primitive; its home is NumericsCore once that package is free
# of concurrent edits. It is kept private here (leading underscore) so nothing outside SpiralWaves
# can depend on it in the meantime.

# band storage: `AB[kl + ku + 1 + i - j, j] == A[i, j]`, size (2kl + ku + 1) × n (the extra `kl` rows
# hold the fill-in that partial pivoting pushes into the upper triangle).
@inline _bidx(kl::Int, ku::Int, i::Int, j::Int) = kl + ku + 1 + i - j

function _band_from_sparse(A::SparseMatrixCSC{TC}, kl::Int, ku::Int) where {TC}
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("_band_from_sparse: matrix must be square"))
    AB = zeros(TC, 2 * kl + ku + 1, n)
    rows = rowvals(A); vals = nonzeros(A)
    @inbounds for j in 1:n, idx in nzrange(A, j)
        i = rows[idx]
        (-ku ≤ i - j ≤ kl) ||
            throw(ArgumentError("_band_from_sparse: entry ($i, $j) lies outside the (kl = $kl, ku = $ku) band"))
        AB[_bidx(kl, ku, i, j), j] = vals[idx]
    end
    AB
end

"""
    _band_lu!(AB, kl, ku) -> ipiv

In-place banded LU with partial pivoting (LAPACK `gbtrf` layout, see [`_band_from_sparse`](@ref)).
Overwrites `AB` with the `L\\U` band factors and returns the pivot vector. Type-generic (BigFloat).
Throws `SingularException` on a zero pivot column.
"""
function _band_lu!(AB::AbstractMatrix{TC}, kl::Int, ku::Int) where {TC}
    n = size(AB, 2)
    ipiv = Vector{Int}(undef, n)
    @inbounds for j in 1:n
        kmax = min(j + kl, n)
        p = j; best = abs(AB[_bidx(kl, ku, j, j), j])
        for i in (j + 1):kmax
            a = abs(AB[_bidx(kl, ku, i, j), j])
            if a > best
                best = a; p = i
            end
        end
        ipiv[j] = p
        iszero(best) && throw(SingularException(j))
        jmax = min(j + kl + ku, n)
        if p != j
            for jj in j:jmax
                t = AB[_bidx(kl, ku, j, jj), jj]
                AB[_bidx(kl, ku, j, jj), jj] = AB[_bidx(kl, ku, p, jj), jj]
                AB[_bidx(kl, ku, p, jj), jj] = t
            end
        end
        piv = AB[_bidx(kl, ku, j, j), j]
        for i in (j + 1):kmax
            l = AB[_bidx(kl, ku, i, j), j] / piv
            AB[_bidx(kl, ku, i, j), j] = l
            iszero(l) && continue
            for jj in (j + 1):jmax
                AB[_bidx(kl, ku, i, jj), jj] -= l * AB[_bidx(kl, ku, j, jj), jj]
            end
        end
    end
    ipiv
end

"""
    _band_solve!(AB, ipiv, kl, ku, b) -> b

Solve `A x = b` in place from the band factors of [`_band_lu!`](@ref) (`b` is overwritten by `x`).
"""
function _band_solve!(AB::AbstractMatrix, ipiv::Vector{Int}, kl::Int, ku::Int, b::AbstractVector)
    n = length(b)
    @inbounds for j in 1:n                                   # forward: P and L
        p = ipiv[j]
        if p != j
            b[j], b[p] = b[p], b[j]
        end
        bj = b[j]
        iszero(bj) && continue
        for i in (j + 1):min(j + kl, n)
            b[i] -= AB[_bidx(kl, ku, i, j), j] * bj
        end
    end
    @inbounds for j in n:-1:1                                # backward: U
        s = b[j]
        for jj in (j + 1):min(j + kl + ku, n)
            s -= AB[_bidx(kl, ku, j, jj), jj] * b[jj]
        end
        b[j] = s / AB[_bidx(kl, ku, j, j), j]
    end
    b
end

"Solve the assembled system; see the file header for the `linsolve` policy."
function _solve_linear(A::SparseMatrixCSC{TC}, b::Vector{TC}, kl::Int, ku::Int, linsolve::Symbol) where {TC}
    mode = linsolve === :auto ? (TC === ComplexF64 ? :sparse : :banded) : linsolve
    if mode === :sparse
        TC === ComplexF64 || throw(ArgumentError("solve_mode: linsolve = :sparse needs ComplexF64 (UMFPACK); got $TC"))
        return A \ b
    elseif mode === :banded
        AB = _band_from_sparse(A, kl, ku)
        ipiv = _band_lu!(AB, kl, ku)
        return _band_solve!(AB, ipiv, kl, ku, copy(b))
    end
    throw(ArgumentError("solve_mode: linsolve must be :auto, :sparse or :banded, got $(repr(linsolve))"))
end

# ---------------------------------------------------------------------------------------------
# the outgoing-wave boundary rows (66)-(68)
# ---------------------------------------------------------------------------------------------

"""
    eigen2(C) -> (k1, k2, Y1, Y2)

Closed-form eigenvalues and eigenvectors of a 2×2 matrix `C` (no `LinearAlgebra.eigen`, which has no
generic method for `Complex{BigFloat}`): k± = [tr ± √(tr² − 4det)]/2 with tr = C₁₁ + C₂₂ and
det = C₁₁C₂₂ − C₁₂C₂₁; the eigenvector of k is (C₁₂, k − C₁₁) or (k − C₂₂, C₂₁), whichever has the
larger norm (both are exact null vectors of C − kI; the larger one is the numerically safe choice),
normalised to unit ∞-norm. For a diagonal `C` the columns of the identity are returned.
"""
function eigen2(C::AbstractMatrix{TC}) where {TC<:Number}
    a, c, b, d = C[1, 1], C[2, 1], C[1, 2], C[2, 2]
    tr = a + d
    disc = sqrt((a - d)^2 + 4 * b * c)
    k1 = (tr + disc) / 2
    k2 = (tr - disc) / 2
    (k1, k2, _eigvec2(a, b, c, d, k1), _eigvec2(a, b, c, d, k2))
end

function _eigvec2(a::TC, b::TC, c::TC, d::TC, k::TC) where {TC<:Number}
    v1 = (b, k - a)                       # from row 1: (a − k) Y₁ + b Y₂ = 0
    v2 = (k - d, c)                       # from row 2: c Y₁ + (d − k) Y₂ = 0
    n1 = max(abs(v1[1]), abs(v1[2]))
    n2 = max(abs(v2[1]), abs(v2[2]))
    if n1 == 0 && n2 == 0
        return (one(TC), zero(TC))        # C = kI: every vector is an eigenvector
    end
    v, n = n1 ≥ n2 ? (v1, n1) : (v2, n2)
    (v[1] / n, v[2] / n)
end

"""
    outgoing_branch(C, ωd; outward, bc = :outgoing) -> (k, Y, Y⊥)

Boundary datum (66)-(68) at one domain edge: the eigenpair of the edge matrix `C` selected as the wave
travelling AWAY from the orbit, and the perpendicular vector Y_⊥ = (−Y₂, Y₁) of the Neumann-type row
y·Y_⊥ = 0 (the unconjugated bilinear form, so the row is exactly "y ∥ Y").

`outward = true` at the outer edge, `false` at the inner edge; `ωd` is the (real part of the) Doppler
frequency (51) there. Selection (see the file header, ERRATUM E-M4):

  * `:outgoing` (default) — PHYSICAL. If the edge is in a wave zone (|Im k| > |Re k| for either root,
    i.e. the eigenvalues are a conjugate pair k = i k_r) pick the root with sign(Im k) = ±sign(ω_d)
    (+ outward, − inward), which is sign(group velocity c_s² k_r/ω_d) = ± as required. If the edge is
    evanescent (both roots essentially real) pick the root that DECAYS away from the orbit
    (Re k < 0 outward, Re k > 0 inward).
  * `:paper` — the literal printed criterion Im[k] < 0 at both edges.
"""
function outgoing_branch(C::AbstractMatrix{TC}, ωd::Real; outward::Bool, bc::Symbol = :outgoing) where {TC<:Number}
    k1, k2, Y1, Y2 = eigen2(C)
    first = if bc === :paper
        imag(k1) ≤ imag(k2)
    elseif bc === :outgoing
        s = outward ? 1 : -1
        propagating = abs(imag(k1)) > abs(real(k1)) || abs(imag(k2)) > abs(real(k2))
        if propagating
            w = s * (signbit(ωd) ? -1 : 1)
            w * imag(k1) ≥ w * imag(k2)
        else
            -s * real(k1) ≥ -s * real(k2)
        end
    else
        throw(ArgumentError("outgoing_branch: bc must be :outgoing or :paper, got $(repr(bc))"))
    end
    k, Y = first ? (k1, Y1) : (k2, Y2)
    (k, Y, (-Y[2], Y[1]))
end

"""
    wkb_wavenumber_gr(bg::DiscBackground, m, r) -> k_r
    wkb_wavenumber2_gr(bg::DiscBackground, m, r) -> k_r²

RELATIVISTIC local WKB radial wavenumber of the master system — the fluid-frame form of the local
dispersion relation ω̃² = κ̃² + c_s² k̃² with the PROPER Doppler frequency ω̃ = u^t ω_d, the proper
epicyclic frequency κ̃ = u^t κ and the PROPER radial wavenumber k̃ = e^{−μ}k_r:

    k_r = e^{μ} u^t √(ω_d² − κ²) / c_s ,        k_r² = e^{2μ} (u^t)² (ω_d² − κ²)/c_s² .

FINDING (this node, measured): this is what the eigenvalue of the relativistic C at the domain edges
actually is. At the reference configuration (χ = 0, a• = 15, h0 = 0.05, 𝔭 = 𝔮 = 0, m = 3) the selected
branch has Im k = 20.4468 at r_in = r_ISCO = 6 and 3.70629 at r_out = 66.264, while this formula gives
20.4579 and 3.70664 — 0.05 % and 0.009 % (the remainder is the sub-leading WKB amplitude term, the
trace of C). The NEWTONIAN [`wkb_wavenumber`](@ref) of grid.jl (App. C) gives 11.8114 and 3.56657, i.e.
it is too small by exactly e^{μ}u^t = 1.7321 and 1.0392. The App.-C grid therefore UNDER-RESOLVES the
inner wave zone by the redshift/curvature factor e^{μ}u^t (→ 1 far from the hole, √3 at the χ = 0
ISCO). [FUTURE] feeding this wavenumber into `build_grid` (a strict refinement of the App.-C recipe)
is left to a separate node — it changes every grid and therefore every landed grid test.
[PRELIMINARY] whether κ² enters as (u^tκ)² or as κ² is below the 1e-3 agreement level here.
"""
function wkb_wavenumber2_gr(bg::DiscBackground, m::Integer, r)
    rp = rotation_profile(bg, r)
    pot = kerr_equatorial_potentials(bg.χ, r)
    ωd = omega_d(bg, m, r)
    (exp(pot.μ) * rp.ut)^2 * (ωd^2 - kappa2(bg, r)) / csound(bg, r)^2
end

function wkb_wavenumber_gr(bg::DiscBackground, m::Integer, r)
    k2 = wkb_wavenumber2_gr(bg, m, r)
    sqrt(max(k2, zero(k2)))
end

# ---------------------------------------------------------------------------------------------
# assembly of (64) + the two boundary rows
# ---------------------------------------------------------------------------------------------

"""
    midpoint_system(r, Cmid, Smid, row_in, rhs_in, row_out, rhs_out) -> (A, b)

Assemble the implicit-midpoint system (64) on the ascending node vector `r` (N nodes) as a
`SparseMatrixCSC{Complex{T}}` of size 2N×2N and its right-hand side, with the unknown ordering
u[2i−1] = y_i[1], u[2i] = y_i[2]:

  * row 1: the inner boundary row  `row_in · y_1 = rhs_in`;
  * rows 2i, 2i+1 (i = 1 … N−1): the two rows of cell i,
    (I − Δr_i/2 C_i)·y_{i+1} − (I + Δr_i/2 C_i)·y_i = Δr_i S_i with `Cmid[i]`, `Smid[i]` the matrix
    and source AT THE CELL MIDPOINT;
  * row 2N: the outer boundary row `row_out · y_N = rhs_out`.

Lower and upper bandwidths are 2 (see [`_band_from_sparse`](@ref)). Used directly by the
manufactured-solution tests, which pass their own boundary rows.
"""
function midpoint_system(r::AbstractVector, Cmid::AbstractVector, Smid::AbstractVector,
                         row_in, rhs_in, row_out, rhs_out)
    N = length(r)
    N ≥ 2 || throw(ArgumentError("midpoint_system: need at least 2 nodes"))
    (length(Cmid) == N - 1 && length(Smid) == N - 1) ||
        throw(DimensionMismatch("midpoint_system: need N−1 = $(N-1) midpoint matrices/sources"))
    TC = promote_type(eltype(first(Cmid)), eltype(first(Smid)), typeof(rhs_in), typeof(rhs_out),
                      eltype(row_in), eltype(row_out))
    n = 2 * N
    nz = 4 + 8 * (N - 1)
    I = Vector{Int}(undef, nz); J = Vector{Int}(undef, nz); V = Vector{TC}(undef, nz)
    b = zeros(TC, n)
    q = 0
    @inbounds begin
        I[q+1] = 1; J[q+1] = 1; V[q+1] = TC(row_in[1])
        I[q+2] = 1; J[q+2] = 2; V[q+2] = TC(row_in[2])
        q += 2
        b[1] = TC(rhs_in)
        for i in 1:(N - 1)
            Δ = r[i+1] - r[i]
            C = Cmid[i]; S = Smid[i]
            h = Δ / 2
            for a in 1:2
                row = 2 * i + a - 1
                for c in 1:2
                    δ = (a == c) ? one(TC) : zero(TC)
                    I[q+1] = row; J[q+1] = 2 * (i - 1) + c; V[q+1] = -(δ + h * C[a, c])   # −(I + Δr/2 C) y_i
                    I[q+2] = row; J[q+2] = 2 * i + c;       V[q+2] =  (δ - h * C[a, c])   # +(I − Δr/2 C) y_{i+1}
                    q += 2
                end
                b[row] = Δ * S[a]
            end
        end
        I[q+1] = n; J[q+1] = n - 1; V[q+1] = TC(row_out[1])
        I[q+2] = n; J[q+2] = n;     V[q+2] = TC(row_out[2])
        q += 2
        b[n] = TC(rhs_out)
    end
    (sparse(view(I, 1:q), view(J, 1:q), view(V, 1:q), n, n), b)
end

"""
    solve_midpoint(r, Cmid, Smid, row_in, rhs_in, row_out, rhs_out; linsolve = :auto) -> y::Matrix

Assemble [`midpoint_system`](@ref) and solve it, returning the 2×N solution matrix
(`y[:, i]` = y_i). The entry point of [`solve_mode`](@ref) and of the manufactured-solution tests
(which supply exact boundary rows instead of the WKB ones).
"""
function solve_midpoint(r::AbstractVector, Cmid::AbstractVector, Smid::AbstractVector,
                        row_in, rhs_in, row_out, rhs_out; linsolve::Symbol = :auto)
    A, b = midpoint_system(r, Cmid, Smid, row_in, rhs_in, row_out, rhs_out)
    u = _solve_linear(A, b, 2, 2, linsolve)
    N = length(r)
    y = Matrix{eltype(u)}(undef, 2, N)
    @inbounds for i in 1:N
        y[1, i] = u[2i-1]; y[2, i] = u[2i]
    end
    y
end

"""
    midpoint_residual(r, Cmid, Smid, y) -> (abs, rel)

Maximum ∞-norm residual of the discrete equations (64) over all cells, absolute and relative (divided
by the largest term entering that cell's row). A correctly solved system leaves only roundoff here;
the diagnostic is kept because it is the cheapest guard against a mis-assembled band.
"""
function midpoint_residual(r::AbstractVector, Cmid::AbstractVector, Smid::AbstractVector,
                           y::AbstractMatrix)
    N = length(r)
    ra = zero(real(eltype(y))); rr = zero(ra)
    @inbounds for i in 1:(N - 1)
        Δ = r[i+1] - r[i]
        C = Cmid[i]; S = Smid[i]
        for a in 1:2
            lhs = (y[a, i+1] - y[a, i]) / Δ
            rhs = (C[a, 1] * (y[1, i] + y[1, i+1]) + C[a, 2] * (y[2, i] + y[2, i+1])) / 2 + S[a]
            res = abs(lhs - rhs)
            scale = max(abs(lhs), abs(rhs), abs(S[a]))
            ra = max(ra, res)
            rr = max(rr, iszero(scale) ? res : res / scale)
        end
    end
    (ra, rr)
end

# ---------------------------------------------------------------------------------------------
# the mode solution
# ---------------------------------------------------------------------------------------------

"""
    ModeSolution{T}

Solution of one azimuthal mode `m` of the master system (47) on the grid `r` (ascending, N nodes):

  * `y::Matrix{Complex{T}}` (2×N) — `y[1, i]` = δÛ(r_i) = δp̂/(ε+p) (48), `y[2, i]` = −iδû^r(r_i) (49);
  * `k_in`, `k_out`, `Y_in`, `Y_out` — the selected boundary eigenvalues/eigenvectors (66)-(68);
  * `residual_abs`, `residual` — [`midpoint_residual`](@ref) of the returned solution;
  * `bc_residual` — |y_1·Y_⊥,in| + |y_N·Y_⊥,out|, the boundary rows evaluated on the solution;
  * `params` — the settings the solution was produced with (`landau`, `bc`, `linsolve`, `mEMRI`,
    `hmodes_N`, `threads`, `rc` = the corotation radius used by the Landau shift), so that node
    quantities (δû_φ, the diagnostics of `diagnostics.jl`) can be recomputed consistently;
  * `seconds` — wall time of the solve (assembly + factorisation), for the survey timings.
"""
struct ModeSolution{T<:Real}
    m::Int
    r::Vector{T}
    y::Matrix{Complex{T}}
    k_in::Complex{T}
    k_out::Complex{T}
    Y_in::NTuple{2,Complex{T}}
    Y_out::NTuple{2,Complex{T}}
    residual_abs::T
    residual::T
    bc_residual::T
    params::NamedTuple{(:landau, :bc, :linsolve, :mEMRI, :hmodes_N, :threads, :rc),
                       Tuple{Bool,Symbol,Symbol,T,Union{Nothing,Int},Bool,T}}
    seconds::Float64
end

Base.length(sol::ModeSolution) = length(sol.r)

"""
    solve_mode(bg::DiscBackground, m; grid = build_grid(bg, m), landau = true, bc = :outgoing,
               linsolve = :auto, mEMRI = 1, hmodes_N = nothing, threads = false) -> ModeSolution

Solve the relativistic master system (47) for azimuthal number `m` on the adaptive grid of
[`build_grid`](@ref) by the implicit-midpoint scheme (64) with the outgoing-wave boundary rows
(66)-(68) — the numerical scheme of Sec. II.E.

The matrices and sources come from [`master_matrices`](@ref) evaluated AT THE CELL MIDPOINTS with the
covariant singular-field modes of [`source_modes`](@ref) (`hmodes_N` overrides the trapezoid node
count, `mEMRI` = m•, `threads` distributes the source projection); the Landau regularisation (62)-(63)
is applied unless `landau = false` (then the system is singular at r_c and only usable away from it).
`bc` selects the boundary branch (`:outgoing`, the physical group-velocity criterion — DEFAULT — or
`:paper`, the literal Im[k] < 0 of M:776; see the file header, ERRATUM E-M4). `linsolve` selects the
linear-algebra backend (`:auto`, `:sparse` = UMFPACK, `:banded` = the generic banded LU).

Deterministic: the result is bitwise independent of the thread count.
"""
function solve_mode(bg::DiscBackground{T}, m::Integer;
                    grid::SpiralGrid{T} = build_grid(bg, m),
                    landau::Bool = true, bc::Symbol = :outgoing, linsolve::Symbol = :auto,
                    mEMRI = 1, hmodes_N::Union{Nothing,Integer} = nothing,
                    threads::Bool = false) where {T}
    m ≥ 1 || throw(ArgumentError("solve_mode: m ≥ 1 required"))
    t0 = time()
    r = grid.r
    N = length(r)
    rc = grid.r_c
    rmid = [(r[i] + r[i+1]) / 2 for i in 1:(N - 1)]
    hmid = source_modes(bg, m, rmid; threads = threads, N = hmodes_N, mEMRI = mEMRI)
    Cmid = Vector{Matrix{Complex{T}}}(undef, N - 1)
    Smid = Vector{Vector{Complex{T}}}(undef, N - 1)
    for i in 1:(N - 1)
        Cmid[i], Smid[i] = master_matrices(bg, m, rmid[i]; hmodes = hmid[i], landau = landau, rc = rc)
    end
    # boundary rows from the edge matrices (the source does not enter (66)-(68))
    Cin, _ = master_matrices(bg, m, r[1]; hmodes = singular_field_mode(bg, m, r[1]; N = hmodes_N, mEMRI = mEMRI),
                             landau = landau, rc = rc)
    Cout, _ = master_matrices(bg, m, r[N]; hmodes = singular_field_mode(bg, m, r[N]; N = hmodes_N, mEMRI = mEMRI),
                              landau = landau, rc = rc)
    k_in, Y_in, W_in = outgoing_branch(Cin, omega_d(bg, m, r[1]); outward = false, bc = bc)
    k_out, Y_out, W_out = outgoing_branch(Cout, omega_d(bg, m, r[N]); outward = true, bc = bc)
    y = solve_midpoint(r, Cmid, Smid, W_in, zero(Complex{T}), W_out, zero(Complex{T}); linsolve = linsolve)
    ra, rr = midpoint_residual(r, Cmid, Smid, y)
    bcres = abs(W_in[1] * y[1, 1] + W_in[2] * y[2, 1]) + abs(W_out[1] * y[1, N] + W_out[2] * y[2, N])
    params = (landau = landau, bc = bc, linsolve = linsolve, mEMRI = convert(T, mEMRI),
              hmodes_N = hmodes_N === nothing ? nothing : Int(hmodes_N), threads = threads, rc = rc)
    ModeSolution{T}(Int(m), collect(r), y, k_in, k_out, Y_in, Y_out, ra, rr, bcres, params, time() - t0)
end

"""
    delta_u_phi(bg::DiscBackground, sol::ModeSolution) -> Vector{Complex{T}}

The covariant δû_φ (route R2, [`delta_u_phi`](@ref)) at every node of a solved mode, with the source
modes re-evaluated at the NODES with the same settings the solve used (`sol.params`).
"""
function delta_u_phi(bg::DiscBackground{T}, sol::ModeSolution{T}) where {T}
    p = sol.params
    hm = source_modes(bg, sol.m, sol.r; threads = p.threads, N = p.hmodes_N, mEMRI = p.mEMRI)
    [delta_u_phi(bg, sol.m, sol.r[i], (sol.y[1, i], sol.y[2, i]);
                 hmodes = hm[i], landau = p.landau, rc = p.rc) for i in eachindex(sol.r)]
end
