# boundary_hooks.jl — the generic ModifyBoundaryData seam (spectre GetModifyBoundaryData
# pattern, audit-spectre-ssf §1.2/§2) over the N6 interface plumbing of dg_operator.jl:
# typed hooks post-process the RECEIVED exterior FaceState before the numerical flux, in
# TWO modes with distinct roles in the linear solve,
#
#   (a) LINEARIZED hooks — linear-homogeneous in the trial-field trace data, applied in
#       EVERY operator application (they belong to the GMRES/Schwarz/multigrid matrix;
#       spectre DgOperator.hpp:708-728). Downstream instance: the vtu −imΩ⟨⟨Ψ⟩⟩ flux jump
#       (spectre Equations.cpp:92-116; dg-jump-conditions node).
#   (b) FIXED-DATA hooks — additive shifts of the received data INDEPENDENT of the trial
#       field, dropped from the linearized application and moved to the RHS via
#       OPERATOR-ON-ZERO: applying the FULL operator to zero yields exactly the fixed
#       source contribution, b -= A_hooks(0), keeping the matrix linear-homogeneous
#       (spectre DgOperator.hpp:1081-1130). Downstream instance: the worldtube ±Ψ^P jumps
#       (spectre Equations.cpp:70-91).
#
# The pair of generic functions `modify_boundary_data` / `modify_boundary_data_linearized`
# is the seam contract (spectre's ModifyBoundaryData::apply / ::apply_linearized); every
# hook type implements both. Also owns the external ghost-BC seam generalized from the
# W4-B zero-Dirichlet ghost: `ghost_state(bc, int)` builds the ghost exterior from the
# interior trace + boundary data, `ghost_state_linearized` is its data-zeroed twin so
# inhomogeneous BC data ride the same operator-on-zero path. Conventions N1-N7
# (dg_operator.jl header) inherited verbatim — hooks and ghosts produce SENDER-convention
# FaceStates and modify ext DATA only, never the assembly.
# DAG node: boundary-data-hooks (elliptic-ssf GLOBAL_DAG row 14, TOPOLOGICAL_SCHEDULE
# W5-B). Built TDD — gates calibrated pre-red in scratchpad calibrate_boundary_hooks.jl.

export BoundaryDataHook, NoHook, LinearizedHook, FixedDataHook, ComposedHook,
       modify_boundary_data, modify_boundary_data_linearized,
       GhostBoundary, DirichletGhost, NeumannGhost,
       ghost_state, ghost_state_linearized

"""
    BoundaryDataHook

Abstract supertype of the ModifyBoundaryData seam (spectre GetModifyBoundaryData pattern,
audit-spectre-ssf §1.2/§2): a hook post-processes the RECEIVED exterior [`FaceState`](@ref)
of one interior interface before the numerical flux, through the pair of generic functions

    modify_boundary_data(hook, int, ext)             -> FaceState   (full application)
    modify_boundary_data_linearized(hook, int, ext)  -> FaceState   (matrix application)

with `int` the receiver's own face trace and `ext` the received neighbor state in SENDER
convention (N3). The linearized function keeps exactly the parts linear-homogeneous in the
trial-field data — it runs inside every matrix application; trial-independent (fixed) data
appear only in the full function and reach the RHS via operator-on-zero (`b -= A_full(0)`,
spectre DgOperator.hpp:1081-1130). Instances: [`NoHook`](@ref) (identity default),
[`LinearizedHook`](@ref) (mode a), [`FixedDataHook`](@ref) (mode b),
[`ComposedHook`](@ref) (both simultaneously). Hooks modify ext DATA only (seam N6) — the
residual assembly and the numerical flux never change.
"""
abstract type BoundaryDataHook end

"""
    NoHook()

The identity hook — the default ModifyBoundaryData seam: the received exterior state is
used unchanged (returned as the SAME object), in both the full and the linearized
application. With `NoHook` on every interface the hook-routed sweep IS the wave-4
`exterior_state` behavior, bitwise (gated in test_boundary_hooks.jl).
"""
struct NoHook <: BoundaryDataHook end

"""
    LinearizedHook(f)

Mode (a) of the seam: a hook LINEAR-HOMOGENEOUS in the trial-field trace data, applied
identically in the full AND the linearized application — it is part of the operator
matrix (spectre `ModifyBoundaryData::apply_linearized`, DgOperator.hpp:708-728). `f` is a
callable `(int::FaceState, ext::FaceState) -> FaceState` returning the modified exterior
state in SENDER convention; the CALLER guarantees linearity `f(αu + v) = αf(u) + f(v)`
in the trace data jointly (gated through the operator in test_boundary_hooks.jl at
1e-13). Downstream instance (dg-jump-conditions node): the vtu slicing-junction flux jump
`φ_ext -= imΩ (ψ_int + ψ_ext)/2` (spectre Equations.cpp:92-116; paper eq:jump_vt/
eq:jump_tu junction conditions [∂_{r*}Ψ] = −imΩ Ψ).
"""
struct LinearizedHook{F} <: BoundaryDataHook
    f::F
end

"""
    FixedDataHook(delta::FaceState)
    FixedDataHook(dpsi::AbstractVector, dnflux::AbstractVector)

Mode (b) of the seam: an additive, trial-field-INDEPENDENT shift of the received data
(SENDER convention, both components at the face points),

    full:        ext -> FaceState(ext.psi + delta.psi, ext.nflux + delta.nflux),
    linearized:  ext -> ext                (identity — fixed data never enter the matrix).

The shift reaches the RHS via operator-on-zero: with zero trial field the full operator
residual is exactly the lifted fixed-data contribution, so `b -= A_full(0)` reproduces the
inhomogeneous problem with a linear-homogeneous matrix (spectre
`impose_inhomogeneous_boundary_conditions_on_source`, DgOperator.hpp:1081-1130).
Downstream instance (dg-jump-conditions node): the worldtube Ψ ↔ Ψ^R = Ψ − Ψ^P transform,
`ψ_ext ± Ψ^P`, `φ_ext ∓ (−n·F(Ψ^P))` (spectre Equations.cpp:70-91) — note the ψ shift has
OPPOSITE sign on the two sides of the interface while the φ shift has the SAME sign,
because n·F flips with the sender's normal (conventions N1/N3; the prescribed-jump sign
table is derived and gated in test_boundary_hooks.jl: `[Ψ] = J, [∂Ψ] = K` ⇔
`(Δψ, Δφ) = (−J, +K)` left / `(+J, +K)` right).
"""
struct FixedDataHook{T<:Number} <: BoundaryDataHook
    delta::FaceState{T}
end
FixedDataHook(dpsi::AbstractVector, dnflux::AbstractVector) =
    FixedDataHook(FaceState(dpsi, dnflux))

"""
    ComposedHook(hooks::BoundaryDataHook...)

Sequential composition of hooks — linearized and fixed modifications acting on the SAME
interface simultaneously. Stage `k` receives the ORIGINAL interior trace `int` (hooks
modify received exterior data only, seam N6) and stage `k−1`'s output as `ext`, in the
order given. The linearized application composes the stages' linearized actions (fixed
stages drop out), so the affine decomposition `A_full(u) = A_lin(u) + A_full(0)` holds for
the composition (gated at 1e-13 in test_boundary_hooks.jl). NOTE order matters when a
later stage reads ψ a previous stage shifted; the composed junction conditions gated in
the manufactured test apply the linearized stage to the RAW received data FIRST (the
spectre vtu convention, Equations.cpp:105-107), then the fixed shifts.
"""
struct ComposedHook{H<:Tuple} <: BoundaryDataHook
    hooks::H
end
ComposedHook(hooks::BoundaryDataHook...) = ComposedHook(hooks)

"""
    modify_boundary_data(hook, int::FaceState, ext::FaceState) -> FaceState

FULL application of a [`BoundaryDataHook`](@ref) to the received exterior state `ext`
(SENDER convention N3), given the receiver's own trace `int`: linearized parts AND fixed
data both act (spectre `ModifyBoundaryData::apply`). Used when applying the affine
operator — in particular on ZERO trial data to harvest the fixed sources
(`b -= A_full(0)`).
"""
modify_boundary_data(::NoHook, int::FaceState, ext::FaceState) = ext
modify_boundary_data(h::LinearizedHook, int::FaceState, ext::FaceState) = h.f(int, ext)
modify_boundary_data(h::FixedDataHook, int::FaceState, ext::FaceState) =
    FaceState(ext.psi .+ h.delta.psi, ext.nflux .+ h.delta.nflux)
modify_boundary_data(c::ComposedHook, int::FaceState, ext::FaceState) =
    foldl((e, h) -> modify_boundary_data(h, int, e), c.hooks; init = ext)

"""
    modify_boundary_data_linearized(hook, int::FaceState, ext::FaceState) -> FaceState

LINEARIZED application of a [`BoundaryDataHook`](@ref): only the parts linear-homogeneous
in the trial-field data act; fixed data are dropped (spectre
`ModifyBoundaryData::apply_linearized`). This is the seam every MATRIX application routes
through — GMRES operator, Schwarz subdomain builds, multigrid levels, `build_matrix`
column probes all see the linearized hooks and never the fixed data.
"""
modify_boundary_data_linearized(::NoHook, int::FaceState, ext::FaceState) = ext
modify_boundary_data_linearized(h::LinearizedHook, int::FaceState, ext::FaceState) =
    h.f(int, ext)
modify_boundary_data_linearized(::FixedDataHook, int::FaceState, ext::FaceState) = ext
modify_boundary_data_linearized(c::ComposedHook, int::FaceState, ext::FaceState) =
    foldl((e, h) -> modify_boundary_data_linearized(h, int, e), c.hooks; init = ext)

"""
    GhostBoundary

Abstract supertype of the external ghost-BC seam — the generalization of the W4-B
`zero_dirichlet_ghost`: at a block-boundary face the exterior state is a GHOST built from
the interior trace and boundary data,

    ghost_state(bc, int)             -> FaceState   (full: data included)
    ghost_state_linearized(bc, int)  -> FaceState   (data zeroed — the matrix ghost),

both in SENDER convention (the ghost's own outward normal is −n_int, N3). Inhomogeneous
boundary data are trial-independent, so they ride the same operator-on-zero path as
[`FixedDataHook`](@ref): the matrix uses the linearized ghost, `b -= A_full(0)` carries
the data lift to the RHS. Instances: [`DirichletGhost`](@ref), [`NeumannGhost`](@ref);
richer condition families (Bayliss-Turkel via flux, no-op None) are app-side consumers of
this seam (dg-boundary-conditions node).
"""
abstract type GhostBoundary end

"""
    DirichletGhost(g::AbstractVector)

Dirichlet ghost with boundary data `g` at the face points — Vu:2024 eq:dirichlet, u
imposed on the AVERAGE: u_ext = 2g − u_int with flux equality F_ext = F_int. In the
SENDER convention N3 (the ghost's outward normal is −n_int):

    ghost_state(bc, int)            = FaceState(2g − ψ_int, −φ_int),
    ghost_state_linearized(bc, int) = FaceState(−ψ_int, −φ_int)      (g ≡ 0),

so the linearized ghost IS the W4-B [`zero_dirichlet_ghost`](@ref) exactly, and `g = 0`
reproduces it in the full application too (gated). Through the penalty flux: Ψ* = g,
{{nF}} = φ_int, penalty riding on the jump 2(ψ_int − g).
"""
struct DirichletGhost{T<:Number} <: GhostBoundary
    g::Vector{T}
end

"""
    NeumannGhost(h::AbstractVector)

Neumann ghost with boundary data `h` = the imposed n·F at the face points, w.r.t. the
RECEIVER's outward normal — Vu:2024 eq:neumann, n·F imposed on the average: u_ext = u_int
(zero jump — the penalty is off) and (n·F)_ext = 2h − (n·F)_int w.r.t. n_int. In the
SENDER convention N3 the ghost's own-normal flux flips sign:

    ghost_state(bc, int)            = FaceState(ψ_int, φ_int − 2h),
    ghost_state_linearized(bc, int) = FaceState(ψ_int, φ_int)        (h ≡ 0).

Through the flux: Ψ* = ψ_int exactly and (nF)* = h, so the residual face correction is
the classical weak Neumann term −ML·(h − φ_int).
"""
struct NeumannGhost{T<:Number} <: GhostBoundary
    h::Vector{T}
end

"""
    ghost_state(bc::GhostBoundary, int::FaceState) -> FaceState

FULL ghost exterior of a block-boundary face: boundary data included. Feed to the
numerical flux as the face's `ext`; on zero trial data the residual then yields the
boundary-data lift that `b -= A_full(0)` moves to the RHS (Vu:2024 §inhomogeneous /
spectre DgOperator.hpp:1081-1130).
"""
ghost_state(bc::DirichletGhost, int::FaceState) =
    FaceState(2 .* bc.g .- int.psi, -int.nflux)
ghost_state(bc::NeumannGhost, int::FaceState) =
    FaceState(int.psi, int.nflux .- 2 .* bc.h)

"""
    ghost_state_linearized(bc::GhostBoundary, int::FaceState) -> FaceState

LINEARIZED ghost exterior: the boundary DATA are zeroed, keeping only the (linear)
dependence on the interior trace — the ghost every matrix application uses. For
[`DirichletGhost`](@ref) this is exactly the W4-B [`zero_dirichlet_ghost`](@ref).
"""
ghost_state_linearized(::DirichletGhost, int::FaceState) =
    FaceState(-int.psi, -int.nflux)
ghost_state_linearized(::NeumannGhost, int::FaceState) =
    FaceState(int.psi, int.nflux)
