# dg_operator.jl — elementwise strong-form first-order-flux DG residual of CRP
# eq:dgres_strong (main.tex:565-574) in its matrix form eq:dgres_strong_matrix
# (main.tex:636-654), assembled over the wave-3 per-element matrix layer (dg_matrices.jl):
#
#   R = − MD_i·F^i + M·(β Ψ + γ^i (DΨ)_i)
#       + MLD^T_i · (𝔄^{ij} n_j (Ψ* − Ψ))  −  ML · ((n_i F^i)* − n_i F^i),
#
# with F^i = 𝔄^{ii} ∂_i Ψ (the flux matrix 𝔄 is DIAGONAL, main.tex:575), so the discrete
# equation is R(Ψ) = M·S (massive formulation). Physics coefficients 𝔄/β/γ are CALLER
# callables of the physical coordinates (they live app-side, DESIGN.md §2.1); complex
# coefficients and fields are first-class. DAG node: dg-strong-residuals (elliptic-ssf
# GLOBAL_DAG row 11, TOPOLOGICAL_SCHEDULE W4-A, CRITICAL PATH).
#
# ================= BINDING sign / normal / exchange conventions (owned here) =============
# Downstream nodes (dg-penalty-flux, dg-mortars, boundary-data-hooks, matrix-free-operator)
# inherit these verbatim:
#
#  N1. Tensor-product block maps with positive Jacobians ⇒ element faces are axis-aligned.
#      The outward unit normal of face (dim d, side s ∈ {-1,+1}) is n = s·ê_d in physical
#      coordinates; hence on that face n_i F^i = s F^d and 𝔄^{ij} n_j = s 𝔄^{dd} δ^i_d.
#  N2. A `FaceState` carries (ψ, φ = n·F) w.r.t. the OWNING element's outward normal, at
#      the face's transverse collocation points (the wave-3 face-matrix column ordering).
#      Volume fluxes are formed at volume nodes (F^d = 𝔄^{dd}(nodes)·∂_dΨ, collocation
#      aliasing as in the CRP/spectre) and INTERPOLATED to the face: φ = s·(Iᵀ F^d).
#      The 𝔄 factor of the MLD^T (symmetrization) term is instead evaluated AT the face
#      points (`ElementCoefficients.Aface`), consistent with the wave-3 MLD^T convention.
#  N3. Exchange is SENDER-convention (spectre's): the exterior state received from a
#      neighbor is the neighbor's own FaceState, whose nflux refers to the NEIGHBOR's
#      outward normal — i.e. minus the receiver's. The numerical flux converts:
#      n_i^int F^i_ext = − φ_ext.
#  N4. `CentralFlux` (the σ→0 limit of eq:numflux, main.tex:584-593):
#          Ψ*    = (ψ_int + ψ_ext)/2,
#          (nF)* = (φ_int − φ_ext)/2      [w.r.t. the RECEIVER's outward normal].
#      On a shared face this gives (nF)*_L = −(nF)*_R exactly (face-normal antisymmetry):
#      lifted interface contributions telescope in the global sum for any data, and both
#      corrections vanish spectrally for smooth fields.
#  N5. Face-sign of the lifted terms in the residual, per face (d, s):
#          R += MLD^T · (s·𝔄^{dd}_face · (Ψ* − ψ_int))  −  ML · ((nF)* − φ_int).
#  N6. Interface plumbing is a PURE FUNCTION FaceState → FaceState applied to the SENT
#      neighbor state before it is used as the exterior state (`exterior_state` = identity
#      default). Worldtube/vtu jump hooks and external BCs later substitute their own
#      functions with the same contract (spectre GetModifyBoundaryData pattern,
#      audit-spectre-ssf §2) — they modify ext DATA, never this assembly.
#  N7. The transmissive/no-op exterior is FaceState(ψ_int, −φ_int): both face corrections
#      vanish identically (used by no-op BCs; gated in test_dg_operator.jl).
# ==========================================================================================
#
# Structural identities inherited from wave 3 and gated through this layer
# (test_dg_operator.jl, constants calibrated pre-red in scratchpad calibrate_dg_operator.jl):
# strong ↔ weak equivalence −MD_d F = MD_dᵀ F − Σ_s ML_s φ_s (2D Jacobian-weighted SBP),
# the discrete divergence theorem Σ 1ᵀMD_d F^d = Σ_faces 1ᵀ ML φ (MD annihilates constants),
# and MLD^T column sums vanishing (Σ_k ℓ′_k ≡ 0) which makes the symmetrization term
# conservation-neutral. MMS: the residual of the exact interpolant of a smooth solution
# converges spectrally in N (measured log10 slope −1.01/order over N = 4..10).

export FaceState, ElementCoefficients, ElementState,
       element_coefficients, element_state,
       NumericalFlux, CentralFlux, numerical_flux, exterior_state,
       strong_residual, strong_residual!

"""
    FaceState{T<:Number}

The DG trace data of ONE element face: `psi` (Ψ interpolated to the face points) and
`nflux` (φ = n·F at the face points), both length `nf` = the face's transverse collocation
count, **w.r.t. the owning element's outward normal** n = s·ê_d (convention N1/N2 of the
dg_operator.jl header). Exchanged between neighbors in SENDER convention (N3): a received
exterior `nflux` refers to the sender's outward normal. Construct manually (BC ghost data)
or via [`element_state`](@ref) (own faces). Vectors must have equal length.
"""
struct FaceState{T<:Number}
    psi::Vector{T}
    nflux::Vector{T}
    function FaceState{T}(psi::AbstractVector, nflux::AbstractVector) where {T<:Number}
        length(psi) == length(nflux) ||
            throw(ArgumentError("FaceState: psi and nflux lengths differ " *
                                "($(length(psi)) vs $(length(nflux)))"))
        return new{T}(convert(Vector{T}, psi), convert(Vector{T}, nflux))
    end
end
function FaceState(psi::AbstractVector, nflux::AbstractVector)
    T = promote_type(eltype(psi), eltype(nflux))
    return FaceState{T}(psi, nflux)
end

"""
    ElementCoefficients{T<:Number}

Per-element nodal values of the first-order-flux-form coefficients of eq:first-order /
eq:dgres_strong (CRP main.tex:369-373, 565-575), evaluated once from caller callables by
[`element_coefficients`](@ref) (matrix-free discipline: coefficients are cached per element,
never assembled globally):

- `A::NTuple{2,Vector{T}}` — diagonal flux-matrix entries 𝔄^{11}, 𝔄^{22} at the volume
  nodes (`F^d = A[d] .* ∂_dΨ`).
- `beta::Vector{T}`, `gamma::NTuple{2,Vector{T}}` — source coefficients at the volume
  nodes (`sources = beta.*Ψ + Σ_d gamma[d].*∂_dΨ`).
- `Aface::NTuple{2,NTuple{2,Vector{T}}}` — `Aface[d][s]` = 𝔄^{dd} at the face points of
  face (d, side s: 1 ↔ −1, 2 ↔ +1), consumed by the MLD^T symmetrization term (which
  evaluates 𝔄 AT the face, convention N2) and by penalty-type numerical fluxes.

`T` promotes the element's geometry real type with every evaluated value — complex
coefficients (the CRP m-mode operator) make `T` complex.
"""
struct ElementCoefficients{T<:Number}
    A::NTuple{2,Vector{T}}
    beta::Vector{T}
    gamma::NTuple{2,Vector{T}}
    Aface::NTuple{2,NTuple{2,Vector{T}}}
end

"""
    ElementState{T<:Number}

Precomputed per-element DG apply state (one derivative evaluation per apply — the natural
two-pass mesh sweep: every element builds its state, THEN residuals consume neighbor face
states through the [`exterior_state`](@ref) plumbing):

- `psi::Vector{T}` — nodal field values (may alias the input vector when no promotion is
  needed).
- `dpsi::NTuple{2,Vector{T}}` — physical derivatives ∂_dΨ at the volume nodes (chain rule
  through the diagonal tensor-product Jacobian, 1D differentiation from the refelement
  cache).
- `F::NTuple{2,Vector{T}}` — volume fluxes `F^d = A[d] .* dpsi[d]`.
- `faces::NTuple{2,NTuple{2,FaceState{T}}}` — own-face trace data `faces[d][s]`
  (s: 1 ↔ side −1, 2 ↔ side +1) with `psi = Iᵀψ`, `nflux = s·(Iᵀ F^d)` (conventions
  N1/N2).

Built by [`element_state`](@ref); consumed by [`strong_residual!`](@ref).
"""
struct ElementState{T<:Number}
    psi::Vector{T}
    dpsi::NTuple{2,Vector{T}}
    F::NTuple{2,Vector{T}}
    faces::NTuple{2,NTuple{2,FaceState{T}}}
end

# physical derivatives of a nodal vector on the element: Dphys_d = diag(1/j_d)(D_d ⊗ δ),
# with the 1D differentiation matrices from the (RT, N, precision)-keyed refelement cache
function _physical_derivatives(em::ElementMatrices{RT}, ψ::Vector{T}) where {RT,T}
    Nx, Ny = em.N
    U = reshape(ψ, Nx, Ny)
    D1 = refelement(RT, Nx).D
    D2 = refelement(RT, Ny).D
    du1 = vec((D1 * U) ./ em.jac[1])
    du2 = vec((U * transpose(D2)) ./ transpose(em.jac[2]))
    return (du1, du2)
end

"""
    element_coefficients(em, blk::Block, id::ElementId, A, β, γ) -> ElementCoefficients

Evaluate the caller's coefficient callables on element `id` of block `blk` (the element
whose matrices are `em`): `A = (A¹¹, A²²)`, `β`, and `γ = (γ¹, γ²)`, each a callable
`(x₁, x₂) -> value` of the PHYSICAL coordinates (values may be complex — the CRP m-mode
coefficients are). Volume values are taken at the tensor-product collocation nodes
(i-fastest ordering); `Aface` values at the face points — the face's normal-dimension
block-map endpoint × the transverse volume nodes (Gauss faces are off the volume grid).
Throws `ArgumentError` if `id` does not belong to `blk`.
"""
function element_coefficients(em::ElementMatrices{RT}, blk::Block, id::ElementId,
                              A::Tuple{Any,Any}, β, γ::Tuple{Any,Any}) where {RT}
    id.block == blk.id ||
        throw(ArgumentError("element $(id) does not belong to block $(blk.id)"))
    x1, x2 = em.x
    evalvol(f) = vec([f(ξ, η) for ξ in x1, η in x2])       # i-fastest, matches tensor_lin
    Av = (evalvol(A[1]), evalvol(A[2]))
    βv = evalvol(β)
    γv = (evalvol(γ[1]), evalvol(γ[2]))
    Afc = ntuple(d -> ntuple(sidx -> begin
                                 lo, hi = logical_interval(RT, id, d)
                                 xf = RT(map_forward(blockmap(blk, d), sidx == 1 ? lo : hi))
                                 xt = em.x[3 - d]
                                 d == 1 ? [A[1](xf, t) for t in xt] :
                                          [A[2](t, xf) for t in xt]
                             end, 2), 2)
    T = promote_type(RT, eltype(Av[1]), eltype(Av[2]), eltype(βv),
                     eltype(γv[1]), eltype(γv[2]),
                     eltype(Afc[1][1]), eltype(Afc[1][2]),
                     eltype(Afc[2][1]), eltype(Afc[2][2]))
    cv(v) = convert(Vector{T}, v)
    return ElementCoefficients{T}((cv(Av[1]), cv(Av[2])), cv(βv), (cv(γv[1]), cv(γv[2])),
                                  ((cv(Afc[1][1]), cv(Afc[1][2])),
                                   (cv(Afc[2][1]), cv(Afc[2][2]))))
end

"""
    element_state(em, coeffs::ElementCoefficients, Ψ) -> ElementState

Build the per-element apply state for nodal field `Ψ` (length `Nx·Ny`, i-fastest;
real or complex): physical derivatives, volume fluxes `F^d = A^d.*∂_dΨ`, and the four
own-face [`FaceState`](@ref)s `(Iᵀψ, s·IᵀF^d)` (conventions N1/N2). The scalar type
promotes `eltype(Ψ)` with the coefficient and geometry types. This is the "send" pass of
the two-pass DG sweep — neighbors consume `state.faces[d][s]` through the
[`exterior_state`](@ref) plumbing before [`strong_residual!`](@ref) runs.
"""
function element_state(em::ElementMatrices{RT}, c::ElementCoefficients{TC},
                       Ψ::AbstractVector) where {RT,TC}
    n = em.N[1] * em.N[2]
    length(Ψ) == n ||
        throw(ArgumentError("element_state: Ψ has length $(length(Ψ)), element has $n nodes"))
    T = promote_type(eltype(Ψ), TC, RT)
    ψ = convert(Vector{T}, Ψ)
    du = _physical_derivatives(em, ψ)
    F = (c.A[1] .* du[1], c.A[2] .* du[2])
    faces = ntuple(d -> ntuple(sidx -> begin
                                   fm = em.faces[d][sidx]
                                   s = sidx == 1 ? -1 : 1
                                   FaceState{T}(transpose(fm.I) * ψ,
                                                s * (transpose(fm.I) * F[d]))
                               end, 2), 2)
    return ElementState{T}(ψ, du, F, faces)
end

"""
    NumericalFlux

Abstract supertype of the pluggable interface numerical flux (Ψ*, (nF)*) of
eq:dgres_strong / eq:numflux (CRP main.tex:570-571, 583-595). Instances implement

    numerical_flux(flux, fm::FaceMatrices, Aface, int::FaceState, ext::FaceState)
        -> (Ψstar, nFstar)

with `ext` in SENDER convention (N3) and `nFstar` returned w.r.t. the RECEIVER's outward
normal. `fm` and `Aface` (𝔄^{dd} at the face points) carry the face geometry that
penalty-type fluxes need (σ from `fm.njac`/`Aface`; the dg-penalty-flux node) — the
central flux ignores them. This node ships [`CentralFlux`](@ref) only.
"""
abstract type NumericalFlux end

"""
    CentralFlux()

The central/average numerical flux — the trivial (σ → 0) instance of the CRP internal
penalty flux (eq:numflux, main.tex:584-593):

    Ψ*    = {{Ψ}}  = (ψ_int + ψ_ext)/2,
    (nF)* = {{nF}} = (φ_int − φ_ext)/2   [receiver normal; the minus converts the
                                          sender-convention φ_ext, convention N3/N4].

Antisymmetric across a shared face: `(nF)*_L == −(nF)*_R` and `Ψ*_L == Ψ*_R`, so lifted
interface contributions telescope in the global sum (gated in test_dg_operator.jl).
"""
struct CentralFlux <: NumericalFlux end

function numerical_flux(::CentralFlux, ::FaceMatrices, ::AbstractVector,
                        int::FaceState, ext::FaceState)
    return (int.psi .+ ext.psi) ./ 2, (int.nflux .- ext.nflux) ./ 2
end

"""
    exterior_state(fs::FaceState) -> FaceState

The interface-data plumbing seam (convention N6): a PURE FUNCTION mapping the neighbor's
sent face state to the exterior state the receiver feeds into the numerical flux. This
default is the identity (conforming interior interface, no jumps). Worldtube/vtu jump
conditions and external boundary conditions are later imposed by substituting functions
with the same contract that modify the received data (spectre GetModifyBoundaryData
pattern) — e.g. Ψ_ext ± Ψ^P worldtube shifts, −imΩ{{Ψ}} vtu flux jumps, BC ghosts such
as the transmissive `FaceState(ψ_int, −φ_int)` (N7).
"""
exterior_state(fs::FaceState) = fs

"""
    strong_residual!(R, em, coeffs, state, ext, flux=CentralFlux()) -> R

Elementwise strong-form DG residual (CRP eq:dgres_strong_matrix, main.tex:637-646):

    R = − Σ_d MD_d·F^d + M·(β ψ + Σ_d γ^d ∂_dψ)
        + Σ_{faces (d,s)} [ MLD^T·(s·𝔄^{dd}_face·(Ψ* − ψ_face)) − ML·((nF)* − φ_face) ],

overwriting `R` (length `Nx·Ny`). `ext[d][s]` (any nested indexable, s: 1 ↔ −1, 2 ↔ +1)
supplies the four exterior [`FaceState`](@ref)s in SENDER convention — produced from
neighbor states via the [`exterior_state`](@ref) plumbing, or synthesized as BC ghost
data. The interface numerical flux is pluggable ([`NumericalFlux`](@ref)); the discrete
equation is `R(Ψ) = M·S` (massive formulation). All conventions N1-N7: file header.
"""
function strong_residual!(R::AbstractVector, em::ElementMatrices,
                          c::ElementCoefficients, st::ElementState, ext,
                          flux::NumericalFlux = CentralFlux())
    n = em.N[1] * em.N[2]
    length(R) == n ||
        throw(ArgumentError("strong_residual!: R has length $(length(R)), element has $n nodes"))
    src = c.beta .* st.psi .+ c.gamma[1] .* st.dpsi[1] .+ c.gamma[2] .* st.dpsi[2]
    R .= em.M * src .- em.MD[1] * st.F[1] .- em.MD[2] * st.F[2]
    for d in 1:2, sidx in 1:2
        s = sidx == 1 ? -1 : 1
        fm = em.faces[d][sidx]
        fs = st.faces[d][sidx]
        Af = c.Aface[d][sidx]
        ψstar, φstar = numerical_flux(flux, fm, Af, fs, ext[d][sidx])
        R .+= fm.MLDT * ((s .* Af) .* (ψstar .- fs.psi))    # + (∂ψ_p, 𝔄 n (Ψ*−Ψ))_∂Ω
        R .-= fm.ML * (φstar .- fs.nflux)                   # − (ψ_p, (nF)*−nF)_∂Ω
    end
    return R
end

"""
    strong_residual(em, coeffs, state, ext, flux=CentralFlux()) -> Vector

Allocating variant of [`strong_residual!`](@ref): the result eltype promotes the state,
coefficient, and exterior-state scalar types.
"""
function strong_residual(em::ElementMatrices, c::ElementCoefficients{TC},
                         st::ElementState{TS}, ext,
                         flux::NumericalFlux = CentralFlux()) where {TC,TS}
    T = promote_type(TS, TC)
    for d in 1:2, sidx in 1:2
        e = ext[d][sidx]
        T = promote_type(T, eltype(e.psi), eltype(e.nflux))
    end
    return strong_residual!(Vector{T}(undef, em.N[1] * em.N[2]), em, c, st, ext, flux)
end
