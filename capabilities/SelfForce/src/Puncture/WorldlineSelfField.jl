# Worldline regular self-field assembly (the h1_regular/hR1 design surface) + its covariant
# derivatives. Port of Punctures/worldline_self_field.jl (hR.m L154-174, Wardell et al.
# arXiv:2508.00087), REDESIGNED:
#
#   • the Schwarzschild worldline connection Γ, ∂_rΓ, ∂_θΓ is CONSUMED from `Spacetime.christoffel`/
#     `christoffel_dr`/`christoffel_dtheta` (re-exported into this Puncture namespace by Static.jl —
#     AHA: one connection primitive; the a=0 method returns Γ[d,a,b] with the same index convention as
#     the legacy local `_christoffel`). No local Christoffel is re-owned.
#   • the per-l regularized l-sum routing REUSES `reg_sum`/`noreg_sum` (RegSum) and the reconstructed
#     coordinate components REUSE `coord_deriv` (TensorBasis) — both already in the Puncture namespace.
#
# `assemble_hR1Data`/`assemble_dhR1Data`/`assemble_d2hR1Data` build the coordinate-frame regular field
# h^R_ab and its 1st/2nd partials at the worldline; `covariant_first`/`covariant_second` lift the
# partials to covariant derivatives on the a=0 background. The second covariant derivative is NOT
# symmetric in (c,d) — the two derivatives do not commute (Riemann) — matching the upstream array.

# index map 1=t,2=r,3=θ,4=φ.  hComps (hR.m L154).
const _HCOMPS = Dict{Tuple{Int,Int},Symbol}(
    (1, 1) => :htt, (1, 2) => :htr, (1, 3) => :htth, (1, 4) => :htph,
    (2, 2) => :hrr, (2, 3) => :hrth, (2, 4) => :hrph,
    (3, 3) => :hthth, (3, 4) => :hthph, (4, 4) => :hphph,
)
_comp_of(a::Int, b::Int) = _HCOMPS[(min(a, b), max(a, b))]

_counts(idxs) = (count(==(1), idxs), count(==(2), idxs),
                 count(==(3), idxs), count(==(4), idxs))

# RegSum vs NoRegSum selection (hR.m L156-160): RegSum when #(t,φ indices) is even.
function _assemble_entry(base, r0, M, lmax, a, b, deriv_idxs)
    nt, nr, nθ, nφ = _counts(deriv_idxs)
    comp = _comp_of(a, b)
    table = [coord_deriv(comp, l, base, r0, M; nt = nt, nr = nr, nθ = nθ, nφ = nφ)
             for l in 0:lmax]
    all_idx = (a, b, deriv_idxs...)
    cnt14 = count(==(1), all_idx) + count(==(4), all_idx)
    return iseven(cnt14) ? reg_sum(table).value : noreg_sum(table).value
end

"""
    assemble_hR1Data(base, r0, M; lmax=40) -> 4×4 Matrix

Coordinate-frame regular field `h^R_ab` at the worldline (hR.m L156). `base(i,l,m)` supplies the
(hR,dhR,ddhR) Lorenz amplitudes for valid `m ≥ 0`. Symmetric in (a,b).
"""
function assemble_hR1Data(base, r0::Real, M::Real; lmax::Integer = 40)
    H = zeros(float(typeof(r0)), 4, 4)
    for a in 1:4, b in a:4
        v = _assemble_entry(base, r0, M, lmax, a, b, ())
        H[a, b] = v; H[b, a] = v
    end
    return H
end

"""
    assemble_dhR1Data(base, r0, M; lmax=40) -> 4×4×4 Array

Partial derivatives `∂_c h^R_ab` at the worldline (hR.m L158), `c` the 3rd index.
"""
function assemble_dhR1Data(base, r0::Real, M::Real; lmax::Integer = 40)
    D = zeros(float(typeof(r0)), 4, 4, 4)
    for a in 1:4, b in a:4, c in 1:4
        v = _assemble_entry(base, r0, M, lmax, a, b, (c,))
        D[a, b, c] = v; D[b, a, c] = v
    end
    return D
end

"""
    assemble_d2hR1Data(base, r0, M; lmax=40) -> 4×4×4×4 Array

Second partial derivatives `∂_c∂_d h^R_ab` at the worldline (hR.m L160), symmetric in (a,b) and (c,d).
"""
function assemble_d2hR1Data(base, r0::Real, M::Real; lmax::Integer = 40)
    D = zeros(float(typeof(r0)), 4, 4, 4, 4)
    for a in 1:4, b in a:4, c in 1:4, d in c:4
        v = _assemble_entry(base, r0, M, lmax, a, b, (c, d))
        for (p, q) in ((a, b), (b, a)), (s, t) in ((c, d), (d, c))
            D[p, q, s, t] = v
        end
    end
    return D
end

# ---- consumed Schwarzschild connection on the worldline (θ=π/2) --------------------------------
# Spacetime.christoffel(r,θ;M) returns Γ[d,a,b] (a=0 method), re-exported as christoffel here.
_worldline_Γ(r0, M) = christoffel(r0, π / 2; M = M)
_worldline_dΓ(r0, M) = (nothing, christoffel_dr(r0, π / 2; M = M),
                        christoffel_dtheta(r0, π / 2; M = M), nothing)

"""
    covariant_first(hR1, dhR1, r0, M) -> 4×4×4 Array

Covariant first derivative `DhR1Data[a,b,c] = ∇_c h^R_ab = ∂_c h_ab − Γ^e_{ac} h_eb − Γ^e_{bc} h_ae`
on the Schwarzschild worldline (hR.m L162).
"""
function covariant_first(hR1::AbstractMatrix, dhR1::AbstractArray, r0::Real, M::Real)
    Γ = _worldline_Γ(r0, M)
    T = promote_type(eltype(hR1), eltype(dhR1), typeof(float(r0)))
    Dh = zeros(T, 4, 4, 4)
    for a in 1:4, b in 1:4, c in 1:4
        acc = dhR1[a, b, c]
        for e in 1:4
            acc -= Γ[e, a, c] * hR1[e, b] + Γ[e, b, c] * hR1[a, e]
        end
        Dh[a, b, c] = acc
    end
    return Dh
end

"""
    covariant_second(hR1, dhR1, d2hR1, r0, M) -> 4×4×4×4 Array

Covariant second derivative `D2hR1Data[a,b,c,d] = ∇_d ∇_c h^R_ab` on the Schwarzschild worldline
(hR.m L165); the 4th index `d` is the outer derivative. Not symmetric in (c,d): covariant derivatives
do not commute.
"""
function covariant_second(hR1::AbstractMatrix, dhR1::AbstractArray,
                          d2hR1::AbstractArray, r0::Real, M::Real)
    Γ = _worldline_Γ(r0, M)
    dΓ = _worldline_dΓ(r0, M)
    Dh = covariant_first(hR1, dhR1, r0, M)
    T = promote_type(eltype(hR1), eltype(dhR1), eltype(d2hR1), typeof(float(r0)))
    D2 = zeros(T, 4, 4, 4, 4)
    for a in 1:4, b in 1:4, c in 1:4, d in 1:4
        dG = dΓ[d]                              # ∂_d Γ
        pc = d2hR1[a, b, c, d]                  # ∂_d(∇_c h_ab)
        for e in 1:4
            pc -= Γ[e, a, c] * dhR1[e, b, d] + Γ[e, b, c] * dhR1[a, e, d]
            if dG !== nothing
                pc -= dG[e, a, c] * hR1[e, b] + dG[e, b, c] * hR1[a, e]
            end
        end
        acc = pc                               # connection on the (0,3) tensor ∇_c h, direction d
        for e in 1:4
            acc -= Γ[e, a, d] * Dh[e, b, c] + Γ[e, b, d] * Dh[a, e, c] +
                   Γ[e, c, d] * Dh[a, b, e]
        end
        D2[a, b, c, d] = acc
    end
    return D2
end
