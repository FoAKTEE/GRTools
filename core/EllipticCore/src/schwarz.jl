# schwarz.jl — element-centered additive Schwarz smoother (CRP main.tex:676): one
# subdomain per DG element covering the element plus `overlap` points into each face
# neighbor (spectre element-centered subdomain geometry, Schwarz.hpp/Vu_2022; DG
# information only crosses faces, so bands are face-neighbor only). Subdomain matrices
# are built PATCH-LOCALLY by RESTRICTED column probes of the global matrix-free
# LINEARIZED apply (the W5-B seam): unit vectors supported on patch dofs, residuals
# evaluated only on the patch's elements — correct off-diagonal couplings within the
# patch, zero outside. Per-ELEMENT matrices, NOT translation-invariance classes: curved
# maps, variable coefficients, hooks, mortars, and boundary ghosts are in the patch by
# construction (the SpecCore local_assembly class-cache trap — audit-speccore risk 2 —
# cannot arise). Factorizations are complex LU, cached under the explicit key
# (element id, subdomain shape, coefficient epoch, shift) with builds/hits counters so
# they persist across AMR iterations (main.tex:676 "cached across AMR iterations").
# One smoothing pass = solve ALL subdomains + weighted overlap combine; the threaded
# sweep uses a gather-then-ordered-fold reduction that is BITWISE identical to serial.
# A diagonal complex shift −i·s (spectre ExplicitInverse.hpp:249-254) is threaded
# through for the schwarz-complex-shift node.
#
# Combine weights (ported PATTERN: spectre Weighting.hpp quintic smoothstep; SpecCore
# schwarz.jl normalization idea): raw weights are the smoothstep self profile on the
# central element × the subdomain's own weight tail continued across each face into
# the overlap bands, then normalized per dof to an EXACT partition of unity by
# quantizing every non-anchor weight to the 2^-30 grid and giving the largest-raw
# anchor the exact residual 1 − Σ others. All quantities are ≤31-bit dyadics in [0,1],
# so per-dof weight sums are EXACT in every fold order (gated == 1.0), and the 1e-9
# quantization only perturbs the smoothing heuristic.
#
# ERRATUM found during pre-red calibration (evidence/schwarz-smoother/): SpecCore
# schwarz.jl:62 `_wintr(ξ, δ, sign) = _extruding(ξ - sign*2, δ, -sign)` evaluates the
# smoothstep on its identically-zero far side for every band point (δ ≤ 1), so ALL
# SpecCore raw band weights are 0 and its normalization silently degenerates the
# combine to RAS (central-indicator weights). The correct intruding weight is the own
# subdomain's extruding tail at the band point's continued coordinate,
# `_extruding(ξ + 2·sign, δ, sign)` — 0.5 at the face, → 0 one overlap-width into the
# neighbor. Measured on the 4x4/N=5/ov=2 Poisson model: 2-pass white-noise residual
# damping 90.0 (corrected) vs 7.6 (RAS-degenerate) — the >10x damping gate falsifies
# the degenerate weighting.
#
# DAG nodes: schwarz-smoother (this file), schwarz-complex-shift (exercises `shift`).
# Threading: never for non-isbits scalars (BigFloat — binding owner directive).

export SchwarzSmoother, SubdomainCache, smooth!, subdomain_matrix

# ---- quintic smoothstep weighting (spectre Weighting.hpp pattern, RT-generic) ----------
function _sw_smoothstep(x0::RT, x1::RT, x::RT) where {RT<:AbstractFloat}
    t = clamp((x - x0) / (x1 - x0), zero(RT), one(RT))
    return t * t * t * (t * (6t - 15) + 10)
end
# extruding weight toward side s ∈ (-1, +1): 1 deep inside, 1/2 at the ±1 face, → 0 one
# overlap-width δ outside the face
_sw_extruding(ξ::RT, δ::RT, s::Int) where {RT} =
    (1 + s) / RT(2) - s * _sw_smoothstep(RT(s) - δ, RT(s) + δ, ξ)
# two-sided central self weight
_sw_wself(ξ::RT, δ::RT) where {RT} =
    _sw_extruding(ξ, δ, -1) + _sw_extruding(ξ, δ, 1) - 1
# intruding weight of the subdomain on a band point of the face neighbor on side s,
# evaluated at the NEIGHBOR's logical coordinate ξ: the own extruding tail continued
# across the face (neighbor coordinate ξ ↔ continued coordinate ξ + 2s). See the
# header ERRATUM — this is the corrected SpecCore `_wintr`.
_sw_wintr(ξ::RT, δ::RT, s::Int) where {RT} = _sw_extruding(ξ + RT(2s), δ, s)
# logical overlap width: distance from the face to the first node OUTSIDE the overlap
_sw_delta(nodes::Vector{RT}, ov::Int) where {RT} =
    ov >= length(nodes) ? RT(2) : nodes[ov+1] + 1

"""
    SubdomainCache()

Cache of Schwarz subdomain LU factorizations, shared across [`SchwarzSmoother`](@ref)
builds (the CRP main.tex:676 "matrix inverses remain cached across AMR iterations"
mechanism). `entries` maps the explicit key `(element id, subdomain shape, coefficient
epoch, shift)` to a factorization; `builds` counts factorizations constructed, `hits`
counts reuses. The subdomain shape is `(central orders, overlap, per-face neighbor
(ElementId, orders) tuples)`, so any mesh change around an element changes its key.
The coefficient `epoch` is CALLER-managed: bump it whenever the operator's
coefficients, hooks, boundary conditions, or penalty change — the cache cannot see
those. Share a cache only across smoothers of the same operator family (same scalar
type); mismatched entries fail the factorization type assertion at lookup.
"""
mutable struct SubdomainCache
    entries::Dict{Any,Any}
    builds::Int
    hits::Int
end
SubdomainCache() = SubdomainCache(Dict{Any,Any}(), 0, 0)

# One element-centered subdomain: patch dofs (central element block first, in layout
# order, then the ov-deep face-neighbor bands in (d, sidx, mortar) order), the
# per-patch-dof (element, local dof) addresses for the restricted probes, the row
# window (central + band elements), and the exact-PoU combine weights.
struct SchwarzSubdomain{RT<:AbstractFloat}
    gdofs::Vector{Int}
    eldof::Vector{NTuple{2,Int}}
    window::Vector{Int}
    weights::Vector{RT}
end

# geometry of subdomain e: dofs + raw smoothstep weights (normalized later)
function _sw_subdomain_dofs(mesh::Mesh2D{RT}, e::Int, ov::Int) where {RT}
    Nx, Ny = mesh.orders[e]
    ξc = (refelement(RT, Nx).nodes, refelement(RT, Ny).nodes)
    δc = (_sw_delta(ξc[1], ov), _sw_delta(ξc[2], ov))
    npatch = Nx * Ny
    gdofs = Int[]; eldof = NTuple{2,Int}[]; raw = RT[]
    sizehint!(gdofs, 2npatch); sizehint!(eldof, 2npatch); sizehint!(raw, 2npatch)
    off = mesh.offsets[e]
    for j in 1:Ny, i in 1:Nx
        l = tensor_lin(i, j, Nx)
        push!(gdofs, off + l)
        push!(eldof, (e, l))
        push!(raw, _sw_wself(ξc[1][i], δc[1]) * _sw_wself(ξc[2][j], δc[2]))
    end
    window = [e]
    for d in 1:2, sidx in 1:2
        iface = mesh.faces[e][d][sidx]
        iface === nothing && continue
        s = sidx == 1 ? -1 : 1
        for nb in iface.nbs
            Nn = mesh.orders[nb]
            ξn = (refelement(RT, Nn[1]).nodes, refelement(RT, Nn[2]).nodes)
            ovd = min(ov, Nn[d])
            ovd == 0 && continue
            δn = _sw_delta(ξn[d], ov)
            δt = _sw_delta(ξn[3-d], ov)
            nrange = sidx == 2 ? (1:ovd) : (Nn[d]-ovd+1:Nn[d])
            ir = d == 1 ? nrange : (1:Nn[1])
            jr = d == 2 ? nrange : (1:Nn[2])
            offn = mesh.offsets[nb]
            for j in jr, i in ir
                l = tensor_lin(i, j, Nn[1])
                push!(gdofs, offn + l)
                push!(eldof, (nb, l))
                ξd = d == 1 ? ξn[1][i] : ξn[2][j]
                ξt = d == 1 ? ξn[2][j] : ξn[1][i]
                push!(raw, _sw_wintr(ξd, δn, s) * _sw_wself(ξt, δt))
            end
            push!(window, nb)
        end
    end
    return gdofs, eldof, window, raw
end

# all subdomains + the per-dof contribution CSR (subdomain-ascending — BOTH the
# exact-PoU normalization and the threaded ordered fold consume it) + exact weights
function _sw_geometry(mesh::Mesh2D{RT}, ov::Int) where {RT}
    nel = length(mesh.elements)
    gd = Vector{Vector{Int}}(undef, nel)
    ed = Vector{Vector{NTuple{2,Int}}}(undef, nel)
    win = Vector{Vector{Int}}(undef, nel)
    raws = Vector{Vector{RT}}(undef, nel)
    for e in 1:nel
        gd[e], ed[e], win[e], raws[e] = _sw_subdomain_dofs(mesh, e, ov)
    end
    N = mesh.ndofs
    counts = zeros(Int, N)
    for e in 1:nel, g in gd[e]
        counts[g] += 1
    end
    dof_ptr = Vector{Int}(undef, N + 1)
    dof_ptr[1] = 1
    for i in 1:N
        dof_ptr[i+1] = dof_ptr[i] + counts[i]
    end
    dof_s = Vector{Int}(undef, dof_ptr[N+1] - 1)
    dof_k = Vector{Int}(undef, dof_ptr[N+1] - 1)
    cursor = copy(dof_ptr)
    for e in 1:nel                      # e ascending ⇒ every dof's list is s-ordered
        for (k, g) in enumerate(gd[e])
            p = cursor[g]
            dof_s[p] = e
            dof_k[p] = k
            cursor[g] = p + 1
        end
    end
    # exact-PoU weights: quantize non-anchor weights to the 2^-30 grid; the anchor
    # (largest raw, first on tie) takes the exact dyadic residual (header rationale)
    wts = [Vector{RT}(undef, length(g)) for g in gd]
    for g in 1:N
        lo, hi = dof_ptr[g], dof_ptr[g+1] - 1     # never empty: the central block covers
        if lo == hi
            wts[dof_s[lo]][dof_k[lo]] = one(RT)
            continue
        end
        total = zero(RT)
        anchor = lo
        for p in lo:hi
            rw = raws[dof_s[p]][dof_k[p]]
            total += rw
            rw > raws[dof_s[anchor]][dof_k[anchor]] && (anchor = p)
        end
        acc = zero(RT)
        for p in lo:hi
            p == anchor && continue
            w = ldexp(round(ldexp(raws[dof_s[p]][dof_k[p]] / total, 30)), -30)
            wts[dof_s[p]][dof_k[p]] = w
            acc += w                    # exact: sums of ≤31-bit dyadics in [0,1]
        end
        wts[dof_s[anchor]][dof_k[anchor]] = one(RT) - acc      # exact dyadic residual
    end
    subs = [SchwarzSubdomain{RT}(gd[e], ed[e], win[e], wts[e]) for e in 1:nel]
    return subs, dof_ptr, dof_s, dof_k
end

# subdomain shape for the cache key: central orders + overlap + per-face neighbor
# (ElementId, orders) — any surrounding mesh change changes it
_sw_shape(mesh::Mesh2D, e::Int, ov::Int) =
    (mesh.orders[e], ov,
     ntuple(d -> ntuple(sidx -> begin
                            iface = mesh.faces[e][d][sidx]
                            iface === nothing ? () :
                                ntuple(k -> (mesh.elements[iface.nbs[k]],
                                             mesh.orders[iface.nbs[k]]),
                                       length(iface.nbs))
                        end, 2), 2))

# the patch-local subdomain matrix by RESTRICTED column probes: for each patch dof,
# apply the LINEARIZED operator with a unit vector supported there, evaluating states
# only where needed (probe element + zero states over the window's neighbor closure)
# and residuals only on the window elements that carry patch rows. Bitwise identical
# to restricting build_matrix(op) to the patch (gated), at patch-local cost.
function _sw_patch_matrix(op::MeshOperator{S,RT,TC}, sub::SchwarzSubdomain,
                          ::Type{TS}, shift::RT) where {S,RT,TC,TS}
    mesh = op.mesh
    win = sub.window
    w2 = copy(win)                      # state closure: window + its face neighbors
    for w in win, d in 1:2, sidx in 1:2
        iface = mesh.faces[w][d][sidx]
        iface === nothing && continue
        for nb in iface.nbs
            nb in w2 || push!(w2, nb)
        end
    end
    Tp = promote_type(S, TC, RT)
    sts = Vector{ElementState{Tp}}(undef, length(mesh.elements))
    zst = Dict{Int,ElementState{Tp}}()
    for f in w2
        st = element_state(op.ems[f], op.coeffs[f], zeros(S, prod(mesh.orders[f])))
        zst[f] = st
        sts[f] = st
    end
    y = zeros(S, mesh.ndofs)
    n = length(sub.gdofs)
    Ass = Matrix{TS}(undef, n, n)
    for j in 1:n
        ej, lj = sub.eldof[j]
        u = zeros(S, prod(mesh.orders[ej]))
        u[lj] = 1
        sts[ej] = element_state(op.ems[ej], op.coeffs[ej], u)
        for w in win
            _residual!(y, op, sts, w, true)
        end
        sts[ej] = zst[ej]
        for i in 1:n
            Ass[i, j] = y[sub.gdofs[i]]
        end
    end
    if !iszero(shift)
        sh = im * shift                 # spectre ExplicitInverse.hpp:249-254 convention
        for i in 1:n
            Ass[i, i] -= sh
        end
    end
    return Ass
end

"""
    SchwarzSmoother(op::MeshOperator; overlap = 2, shift = 0, epoch = 0,
                    cache = SubdomainCache()) -> SchwarzSmoother

Element-centered additive Schwarz smoother for `op` (CRP main.tex:676): one subdomain
per mesh element, covering the element plus `overlap ≥ 0` points into each face
neighbor (full transverse extent, including both halves across a 2:1 face). Each
subdomain matrix is the restriction of the LINEARIZED global matrix to the patch
(built by restricted column probes — hooks, mortars, ghosts, curved maps, and variable
coefficients included by construction), optionally shifted by `−i·shift` on the
diagonal (the Helmholtz damping hook of the schwarz-complex-shift node, spectre
ExplicitInverse `ComplexShift`), and LU-factored ONCE into `cache` under the key
`(element id, subdomain shape, epoch, shift)` — see [`SubdomainCache`](@ref) for the
epoch discipline. Combine weights are the quintic-smoothstep partition of unity, exact
on every dof (file header).

Public fields: `subdomains` (per element, layout order: `gdofs`, `weights`), `fact`
(cached factorizations), `cache`, `overlap`, `shift`, `epoch`. Apply passes with
[`smooth!`](@ref); inspect patch matrices with [`subdomain_matrix`](@ref).
"""
struct SchwarzSmoother{S<:Number,RT<:AbstractFloat,TS<:Number,FT,
                       OT<:MeshOperator}
    op::OT
    overlap::Int
    shift::RT
    epoch::Int
    cache::SubdomainCache
    subdomains::Vector{SchwarzSubdomain{RT}}
    fact::Vector{FT}
    dof_ptr::Vector{Int}
    dof_s::Vector{Int}
    dof_k::Vector{Int}
end

function SchwarzSmoother(op::MeshOperator{S,RT}; overlap::Integer = 2,
                         shift::Real = 0, epoch::Integer = 0,
                         cache::SubdomainCache = SubdomainCache()) where {S,RT}
    overlap >= 0 ||
        throw(ArgumentError("SchwarzSmoother: overlap must be ≥ 0, got $overlap"))
    ov = Int(overlap)
    mesh = op.mesh
    subs, dof_ptr, dof_s, dof_k = _sw_geometry(mesh, ov)
    sh = convert(RT, shift)
    TS = iszero(sh) ? S : promote_type(S, Complex{RT})
    FT = typeof(lu(ones(TS, 1, 1)))
    fact = Vector{FT}(undef, length(subs))
    for s in eachindex(subs)
        key = (mesh.elements[s], _sw_shape(mesh, s, ov), Int(epoch), sh)
        f = get(cache.entries, key, nothing)
        if f === nothing
            f = lu!(_sw_patch_matrix(op, subs[s], TS, sh))
            cache.entries[key] = f
            cache.builds += 1
        else
            cache.hits += 1
        end
        fact[s] = f::FT
    end
    return SchwarzSmoother{S,RT,TS,FT,typeof(op)}(op, ov, sh, Int(epoch), cache,
                                                  subs, fact, dof_ptr, dof_s, dof_k)
end

"""
    subdomain_matrix(sm::SchwarzSmoother, s::Int) -> Matrix

The subdomain matrix of subdomain `s` (element `sm.op.mesh.elements[s]`), rebuilt by
the restricted column probes (NOT read from the cache — the cache stores only the
factorization): the restriction of `build_matrix(sm.op)` to `sm.subdomains[s].gdofs`
with the diagonal `−i·shift` applied (gated at 1e-13 in test_schwarz.jl).
"""
subdomain_matrix(sm::SchwarzSmoother{S,RT,TS}, s::Int) where {S,RT,TS} =
    _sw_patch_matrix(sm.op, sm.subdomains[s], TS, sm.shift)

# one additive sweep: acc += Σ_s scatter(w_s ⊙ A_ss⁻¹ r|_s). Serial scatters in
# subdomain order; the threaded path solves subdomains in parallel into private
# buffers (phase 1) then folds every dof's contributions in the SAME subdomain order
# over disjoint dofs (phase 2) — bitwise identical to serial for any thread count.
function _sw_sweep!(acc::AbstractVector, sm::SchwarzSmoother{S,RT,TS},
                    r::AbstractVector, threaded::Bool) where {S,RT,TS}
    subs = sm.subdomains
    ns = length(subs)
    if threaded
        ds = Vector{Vector{TS}}(undef, ns)
        Threads.@threads for s in 1:ns
            g = subs[s].gdofs
            rs = Vector{TS}(undef, length(g))
            @inbounds for k in eachindex(g)
                rs[k] = r[g[k]]
            end
            ds[s] = ldiv!(sm.fact[s], rs)
        end
        dof_ptr, dof_s, dof_k = sm.dof_ptr, sm.dof_s, sm.dof_k
        Threads.@threads for g in 1:(length(dof_ptr) - 1)
            @inbounds begin
                a = acc[g]
                for p in dof_ptr[g]:(dof_ptr[g+1] - 1)
                    s = dof_s[p]
                    k = dof_k[p]
                    a += subs[s].weights[k] * ds[s][k]
                end
                acc[g] = a
            end
        end
    else
        for s in 1:ns
            g = subs[s].gdofs
            w = subs[s].weights
            rs = Vector{TS}(undef, length(g))
            @inbounds for k in eachindex(g)
                rs[k] = r[g[k]]
            end
            ldiv!(sm.fact[s], rs)
            @inbounds for k in eachindex(g)
                acc[g[k]] += w[k] * rs[k]
            end
        end
    end
    return acc
end

# BigFloat is never threaded (binding owner directive; matrixfree precedent)
_sw_default_threaded(sm::SchwarzSmoother{S,RT,TS}, u, b) where {S,RT,TS} =
    Threads.nthreads() > 1 && isbitstype(promote_type(TS, eltype(u), eltype(b)))

"""
    smooth!(u, sm::SchwarzSmoother, b; passes = 1, threaded = auto) -> u

Run `passes` Schwarz smoothing passes on `A·u = b` with the LINEARIZED operator
(`b` must already carry any operator-on-zero correction). Each pass recomputes the
residual `r = b − A·u`, solves every subdomain with its cached factorization, and adds
the weighted overlap combination `u += Σ_s R_sᵀ (w_s ⊙ A_ss⁻¹ r|_s)` (main.tex:676 —
the paper runs 2 passes per multigrid level). `threaded` defaults to
`Threads.nthreads() > 1` for bits scalar types (never for BigFloat — owner directive)
and is BITWISE identical to serial (gated).
"""
function smooth!(u::AbstractVector, sm::SchwarzSmoother{S,RT,TS}, b::AbstractVector;
                 passes::Integer = 1,
                 threaded::Bool = _sw_default_threaded(sm, u, b)) where {S,RT,TS}
    n = sm.op.mesh.ndofs
    length(u) == n ||
        throw(ArgumentError("smooth!: u has length $(length(u)), the mesh has $n dofs"))
    length(b) == n ||
        throw(ArgumentError("smooth!: b has length $(length(b)), the mesh has $n dofs"))
    Tw = promote_type(eltype(u), eltype(b), TS)
    y = Vector{Tw}(undef, n)
    r = Vector{Tw}(undef, n)
    δ = Vector{Tw}(undef, n)
    for _ in 1:passes
        apply!(y, sm.op, u; threaded = threaded)
        @. r = b - y
        fill!(δ, zero(Tw))
        _sw_sweep!(δ, sm, r, threaded)
        @. u = u + δ
    end
    return u
end
