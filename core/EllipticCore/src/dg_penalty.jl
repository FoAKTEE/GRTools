# dg_penalty.jl — internal penalty numerical flux of CRP eq:numflux (main.tex:583-595)
# over the dg_operator.jl flux interface (conventions N1-N7 inherited verbatim):
#
#   Ψ*    = {{Ψ}}  = (ψ_int + ψ_ext)/2,
#   (nF)* = {{nF}} − σ n_i 𝔄^{ij} n_j (ψ_int − ψ_ext),
#
# with the penalty law of Vu:2024 (ref-paper/arXiv-2405.06120v2/paper.tex:594-619, the
# CRP's cited σ definition):
#
#   σ = C (max(p_int, p_ext) + 1)² / min(h_int, h_ext),      C = 1.5 (CRP main.tex:595,
#                                                             yaml PenaltyParameter),
#
# where p is the polynomial order PERPENDICULAR to the face and h = 2𝒥/𝒥^surf is the
# element size perpendicular to the face. h choice for mapped elements (documented per the
# node contract): Vu:2024 defines h through the volume/surface Jacobian ratio; we evaluate
# it AT the face. For the tensor-product diagonal block maps (mesh2d.jl) the transverse
# Jacobian cancels exactly, 2𝒥/𝒥^surf = 2·j_normal(face) = 2·njac — on affine maps this is
# EXACTLY the element's physical extent along the face-normal dimension, and on curved
# (inverse-radial) maps it is the face-local normal scale. See perpendicular_size.
#
# Also ships the FIRST external-boundary data instance (node deliverable 2): the
# zero-Dirichlet ghost of Vu:2024 eq:dirichlet — u imposed on the AVERAGE, u_ext = 2u_D −
# u_int with u_D = 0, and flux equality F^i_ext = F^i_int (the audit-speccore "ubc ghost"
# convention u_ext = 2u_D − u_int, du_ext = du_int) — expressed in the SENDER convention
# N3, where the ghost's own outward normal is −n_int: zero_dirichlet_ghost(fs) =
# FaceState(−ψ_int, −φ_int). Through the flux this yields Ψ* = 0, {{nF}} = φ_int, and a
# doubled penalty jump 2ψ_int. (General Dirichlet/Neumann/Bayliss-Turkel data are the
# boundary-data-hooks + dg-boundary-conditions nodes.)
#
# SIPG structure (gated in test_dg_penalty.jl, constants calibrated pre-red in scratchpad
# calibrate_dg_penalty{,2}.jl): on affine meshes the global Poisson matrix is symmetric to
# machine precision and positive definite for C = 1.5 (min sym-eig 1.34 on the 2×1 N=3
# probe) while the σ = 0 central flux is indefinite (−20.5); A(C) is affine-linear in C
# with slope exactly the lifted penalty blocks σ(ΔC)·ML·diag(𝔄face)·Iᵀ. On CURVED maps
# Gauss collocation aliasing leaves a small adjoint inconsistency (rel 1.4e-4 on the
# inverse-radial probe — NOT machine-zero; positivity still holds). First end-to-end
# solves (manufactured Poisson + complex Helmholtz through the wave-3 GMRES) are gated in
# the test file; mesh orchestration stays test-side until matrix-free-operator.
# DAG node: dg-penalty-flux (elliptic-ssf GLOBAL_DAG row 12, TOPOLOGICAL_SCHEDULE W4-B,
# CRITICAL PATH).

export PenaltyFlux, penalty_sigma, perpendicular_size, zero_dirichlet_ghost

"""
    penalty_sigma(C, p_int, p_ext, h_int, h_ext) -> σ

The Vu:2024 penalty factor (arXiv:2405.06120 paper.tex:594-619, the σ of CRP eq:numflux):

    σ = C · (max(p_int, p_ext) + 1)² / min(h_int, h_ext),

with `p` the polynomial order perpendicular to the face on each side and `h` the
perpendicular element size ([`perpendicular_size`](@ref)). Symmetric under exchanging the
int/ext pairs — both sides of a face compute the SAME σ, which SIPG symmetry requires.
This is the general pairwise law; [`PenaltyFlux`](@ref) evaluates it with receiver-side
data on conforming interfaces (bitwise-identical across the face), and the h/p-mortar
node feeds it genuinely two-sided values. The result type follows Julia promotion
(`BigFloat` geometry ⇒ `BigFloat` σ).
"""
penalty_sigma(C::Real, p_int::Integer, p_ext::Integer, h_int::Real, h_ext::Real) =
    C * (max(p_int, p_ext) + 1)^2 / min(h_int, h_ext)

"""
    perpendicular_size(fm::FaceMatrices) -> h

The element size perpendicular to face `fm` in the Vu:2024 penalty sense, h = 2𝒥/𝒥^surf,
evaluated AT the face: for the tensor-product diagonal block maps the transverse Jacobian
cancels, leaving `h = 2·fm.njac` (the composed face-normal Jacobian of the wave-3 layer).
On affine maps this is exactly the element's physical extent along the face-normal
dimension; on curved maps it is the face-local normal scale, bitwise-equal between
conforming same-level neighbors (both sides compose the same block map at the same exact
dyadic endpoint — gated in test_dg_penalty.jl).
"""
perpendicular_size(fm::FaceMatrices) = 2 * fm.njac

"""
    PenaltyFlux(C = 1.5)

The internal penalty numerical flux of CRP eq:numflux (main.tex:584-595) — the
[`NumericalFlux`](@ref) instance coupling neighboring elements in the elliptic DG
operator:

    Ψ*    = (ψ_int + ψ_ext)/2,
    (nF)* = (φ_int − φ_ext)/2 − σ 𝔄face (ψ_int − ψ_ext)     [receiver normal; the minus
                                                             converts the sender-convention
                                                             φ_ext, conventions N3/N4],

with σ = [`penalty_sigma`](@ref)`(C, p, p, h, h)` from the receiver's face data: `p` =
order perpendicular to the face (from the face-matrix shape) and `h` =
[`perpendicular_size`](@ref). On conforming interfaces (this node's scope) both sides'
`(p, h)` are bitwise-identical, so the receiver-side σ IS the pairwise max/min law and the
two sides penalize identically; nonconforming pairing happens on mortars (dg-mortars
node). `C = 1.5` is the CRP's penalty parameter (main.tex:595). `PenaltyFlux(0)`
reproduces [`CentralFlux`](@ref) exactly; C > 0 stabilizes SIPG (positive definite on the
affine Poisson probe, gated).
"""
struct PenaltyFlux{RT<:Real} <: NumericalFlux
    C::RT
end
PenaltyFlux() = PenaltyFlux(1.5)

function numerical_flux(pf::PenaltyFlux, fm::FaceMatrices, Aface::AbstractVector,
                        int::FaceState, ext::FaceState)
    nvol, nf = size(fm.I)
    p = nvol ÷ nf - 1                      # order perpendicular to the face (N_normal − 1)
    h = perpendicular_size(fm)
    σ = penalty_sigma(pf.C, p, p, h, h)    # receiver-side == pairwise on conforming faces
    ψstar = (int.psi .+ ext.psi) ./ 2
    φstar = (int.nflux .- ext.nflux) ./ 2 .- σ .* (Aface .* (int.psi .- ext.psi))
    return ψstar, φstar
end

"""
    zero_dirichlet_ghost(fs::FaceState) -> FaceState

Ghost exterior state imposing the homogeneous Dirichlet condition Ψ = 0 through the
numerical flux (Vu:2024 eq:dirichlet with u_D = 0; the spectre/audit-speccore ghost
convention u_ext = 2u_D − u_int, F_ext = F_int): in the SENDER convention N3 the ghost's
outward normal is −n_int, so

    zero_dirichlet_ghost(FaceState(ψ, φ)) = FaceState(−ψ, −φ).

Fed to [`PenaltyFlux`](@ref)/[`CentralFlux`](@ref) this gives Ψ* = 0 exactly, {{nF}} =
φ_int, and the penalty riding on the doubled jump 2ψ_int. Use as the `ext` of a
block-boundary face (the first boundary-condition instance of the exterior_state plumbing
seam N6; richer boundary data arrive with boundary-data-hooks/dg-boundary-conditions).
"""
zero_dirichlet_ghost(fs::FaceState{T}) where {T} = FaceState{T}(-fs.psi, -fs.nflux)
