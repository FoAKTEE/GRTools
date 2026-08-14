# mortars.jl — nonconforming (h/p) interface coupling of CRP main.tex:678-686: where
# elements meet two-to-one (h) or with different resolution along the interface (p), a
# "mortar" covers the full face with the LARGER of the two abutting resolutions; BOTH
# neighbors project their face data to the mortar (Lagrange interpolation P), compute the
# numerical flux there, and restrict the result back to their face with the
# mass-conservative adjoint R (Vu:2021 convention, as cited by the CRP).
#
# 1D h/p projection building blocks PORTED from specjl SpecCore
# (/data/haiyangw/claude/specjl/specjl/packages/SpecCore/src/spectral/mortars.jl, the
# user's own repo; battery-verified there — audit-speccore: "1-D h+p projection ops
# complete and battery-gated, R∘P = I 1e-12 — never wired into the elliptic operator"),
# adapted to this package's Legendre-GAUSS refelement cache. Key simplification vs the
# Lobatto/Chebyshev reference: on N Gauss nodes the product ℓ_aℓ_b has degree 2N−2 ≤ 2N−1,
# so the N-point Gauss rule integrates it EXACTLY and the diagonal Gauss mass
# M̂ = Diagonal(weights) IS the exact L2 Gram matrix — the reference's `_consistent_ref_mass`
# fine-quadrature Gramian (needed there because the lumped LOBATTO mass under-integrates
# degree-2N products) collapses to the quadrature weights. Hence
#
#   P = Lagrange interpolation face → mortar (lossless: the mortar is the finer grid),
#   R = jac · W_face⁻¹ · Pᵀ · W_mortar,   jac = 1 (full) or 1/2 (half-interval measure),
#
# is the true L2 adjoint w.r.t. the Gauss masses: W_f·R == jac·Pᵀ·W_m, giving R∘P = I on
# the coarse space (full mortars; summed over the two children for a 2:1 split face) and
# exact interface conservation ∫_face R g = jac·∫_mortar g (P preserves constants).
# Restricting instead of interpolating back is load-bearing: interpolation would silently
# drop the high mortar modes' mass (reference docstring, inherited).
#
# Wiring into the wave-4 flux exchange (dg_operator.jl conventions N1-N7 inherited
# VERBATIM, assembly untouched): exchange stays SENDER-convention FaceStates (N3) — the
# sender projects its own face state to the shared mortar(s) (`project_to_mortar`), the
# receiver packs one exterior FaceState per face by concatenating the received mortar
# states in `MortarFlux.mortars` order, and `numerical_flux(::MortarFlux, …)` projects the
# interior trace to each mortar, evaluates the CRP internal penalty flux (eq:numflux,
# main.tex:584-593) THERE, and restricts (Ψ*, (nF)*) back to the receiver's face points.
# The nonconforming penalty σ is the PAIRWISE Vu:2024 law via the exported penalty_sigma
# (the dg-penalty-flux W4-B handoff): σ = C(max(p_int, p_ext)+1)²/min(h_int, h_ext) with p
# the order and h the element size perpendicular to the face, from BOTH sides'
# FaceMatrices. Both sides of a mortar must construct it with identical (σ, Am) — mirrored
# (fm_int, fm_ext) arguments give bitwise-equal σ by the symmetry of the law; the
# receiver-side-σ shortcut breaks SIPG symmetry at O(1) (gated falsifier). `PerFaceFlux`
# routes mixed conforming/nonconforming elements through the single-flux strong_residual!
# seam by dispatching on the face key (fm.dim, fm.side); conforming faces keep the direct
# PenaltyFlux path — the mortar path is BITWISE equal there (gated).
#
# Structural gates (test_mortars.jl, constants calibrated pre-red in scratchpad
# calibrate_dg_mortars.jl): adjoint + independent fine-quadrature integral identity 1e-13;
# R∘P = I 1e-14 (full + half-pair sums); conforming pass-through bitwise; MMS Poisson on a
# 2:1-refined and a p-nonconforming mesh spectral with NO O(1) interface error (interface
# node layer within ×4 of interior, same decay rate); 2:1 interface conservation 1e-12;
# SIPG symmetry 1e-12 + positivity on the affine 2:1 probe; BigFloat(256) exactness.
# Scope: mortar-level coupling only — mesh connectivity containers and the global sweep
# stay test-side until matrix-free-operator; cross-block mortars arrive with the
# domain-blocks-map consumers. DAG node: dg-mortars (elliptic-ssf GLOBAL_DAG row 13,
# TOPOLOGICAL_SCHEDULE W5-A, CRITICAL PATH).

export MortarSize, MORTAR_FULL, MORTAR_LOWER, MORTAR_UPPER,
       projection_prolong, projection_restrict, project_to_mortar,
       Mortar, MortarFlux, PerFaceFlux

"""
    MortarSize

How a mortar covers the face it couples (CRP main.tex:678-686 / spectre convention, ported
from the SpecCore reference): `MORTAR_FULL` — the mortar spans the whole face (conforming
or p-nonconforming interface, and the SMALL side of a 2:1 interface, whose face coincides
with one mortar); `MORTAR_LOWER` / `MORTAR_UPPER` — the mortar spans the lower/upper half
of the face (the BIG side of a 2:1 h-nonconforming interface carries one of each).
"""
@enum MortarSize MORTAR_FULL MORTAR_LOWER MORTAR_UPPER

"""
    projection_prolong(T, Nf, Nm, msize::MortarSize) -> Matrix (Nm × Nf)

Face → mortar prolongation: NODAL Lagrange interpolation of the face polynomial onto the
mortar's Gauss nodes — at the mortar points themselves (`MORTAR_FULL`) or at their images
`(x ± 1)/2` in the face coordinate (`MORTAR_UPPER`/`MORTAR_LOWER` half-interval mortars).
Exact on every polynomial the face represents (degree ≤ Nf−1) because the mortar carries
the finer/equal resolution — `Nm ≥ Nf` is enforced (`ArgumentError`): a coarser mortar
would be lossy, violating the max-resolution rule of CRP main.tex:678-686. Geometry in
`RT = real(float(T))` from the `refelement` cache (BigFloat precision inherited).
"""
function projection_prolong(::Type{T}, Nf::Integer, Nm::Integer,
                            msize::MortarSize) where {T<:Number}
    RT = real(float(T))
    Nf >= 1 || throw(ArgumentError("projection_prolong: face needs Nf ≥ 1 points, got $Nf"))
    Nm >= Nf ||
        throw(ArgumentError("projection_prolong: mortar (N=$Nm) must be at least as fine " *
                            "as the face (N=$Nf) — the mortar carries the LARGER of the " *
                            "two abutting resolutions (CRP main.tex:678-686)"))
    elf = refelement(RT, Int(Nf))
    xm = refelement(RT, Int(Nm)).nodes
    t = msize == MORTAR_FULL ? xm :
        msize == MORTAR_UPPER ? (xm .+ 1) ./ 2 : (xm .- 1) ./ 2
    return interpolation_matrix(elf, t)
end

"""
    projection_restrict(T, Nf, Nm, msize::MortarSize) -> Matrix (Nf × Nm)

Mortar → face restriction: the mass-conservative L2 adjoint of [`projection_prolong`](@ref)
w.r.t. the Gauss mass matrices,

    R = jac · W_face⁻¹ · Pᵀ · W_mortar,    jac = 1 (`MORTAR_FULL`) or 1/2 (half-interval
                                           mortars integrate over half the face measure),

i.e. `W_f·R == jac·Pᵀ·W_m` exactly. On Legendre-Gauss nodes the diagonal quadrature mass
IS the exact L2 Gram matrix (ℓ_aℓ_b has degree 2N−2 ≤ 2N−1), so this is the true L2
projection — algebraically the spectral-truncation/half-interval closed forms of the
SpecCore reference (file header). Satisfies `R∘P = I` on the coarse space (full mortar;
`R_lo·P_lo + R_up·P_up = I` for a 2:1 split face) and `∫_face R g = jac·∫_mortar g`
(interface conservation). PROJECTION, not interpolation: interpolating back would silently
drop the high mortar modes (reference docstring). Gated to 1e-13/1e-14 (1e-70 BigFloat).
"""
function projection_restrict(::Type{T}, Nf::Integer, Nm::Integer,
                             msize::MortarSize) where {T<:Number}
    RT = real(float(T))
    P = projection_prolong(RT, Nf, Nm, msize)
    wf = refelement(RT, Int(Nf)).weights
    wm = refelement(RT, Int(Nm)).weights
    jac = msize == MORTAR_FULL ? one(RT) : one(RT) / 2
    R = Matrix{RT}(undef, Int(Nf), Int(Nm))
    for b in 1:Int(Nm), a in 1:Int(Nf)
        R[a, b] = jac * wm[b] * P[b, a] / wf[a]
    end
    return R
end

"""
    project_to_mortar(P::AbstractMatrix, fs::FaceState) -> FaceState

The SENDER half of the mortar exchange (CRP main.tex:678-686: "both neighboring elements
project their interface data to the mortar"): apply the face → mortar prolongation `P` to
both trace components. The result stays in the sender's own outward-normal convention
(N2/N3 unchanged — projection is linear and normal-agnostic); the receiver packs one
exterior `FaceState` per face by concatenating the received mortar states in its
`MortarFlux.mortars` order.
"""
project_to_mortar(P::AbstractMatrix, fs::FaceState) = FaceState(P * fs.psi, P * fs.nflux)

"""
    Mortar{RT,TC}

One mortar of one element face, as seen from the RECEIVER side: `P`/`R` — that side's
face → mortar prolongation and mass-conservative restriction back; `sigma` — the PAIRWISE
Vu:2024 penalty factor of this interface (identical on both sides by the max/min symmetry
of [`penalty_sigma`](@ref)); `Am` — the flux coefficient `n_i 𝔄^{ij} n_j = 𝔄^{dd}` at the
mortar points (unsigned, like the wave-4 `Aface`; complex for the CRP m-mode operator).
Both sides of an interface must hold identical `(sigma, Am)` for SIPG symmetry — evaluate
the coefficient at the mortar's physical points once per interface (or via a
bitwise-identical composition) and construct with mirrored `(fm_int, fm_ext)`.
Build via the [`Mortar`](@ref Mortar(::FaceMatrices, ::FaceMatrices, ::MortarSize, ::AbstractVector))
constructor; consumed by [`MortarFlux`](@ref).
"""
struct Mortar{RT<:AbstractFloat,TC<:Number}
    P::Matrix{RT}
    R::Matrix{RT}
    sigma::RT
    Am::Vector{TC}
end

"""
    Mortar(fm_int, fm_ext, msize::MortarSize, Am; C = 1.5) -> Mortar

Build the receiver-side mortar coupling face `fm_int` to the neighbor face `fm_ext`
(both wave-3 `FaceMatrices` of the SAME interface): the mortar carries the larger of the
two transverse resolutions, `Nm = max(nf_int, nf_ext)` (CRP main.tex:678-686), spanning
the receiver's face per `msize` (`MORTAR_FULL`, or `MORTAR_LOWER`/`MORTAR_UPPER` on the
big side of a 2:1 interface, where the mortar's span is the small neighbor's face). The
penalty factor is the genuinely two-sided pairing of the exported σ law — the
dg-penalty-flux node's nonconforming handoff:

    σ = penalty_sigma(C, p_int, p_ext, h_int, h_ext),

with `p` = order perpendicular to the face (from each side's face-matrix shape) and
`h` = [`perpendicular_size`](@ref) of each side. On conforming interfaces this reproduces
the receiver-side `PenaltyFlux` σ bitwise. `Am` must hold the `Nm` mortar-point values of
`𝔄^{dd}` (`ArgumentError` otherwise); `C = 1.5` is the CRP default (main.tex:595).
"""
function Mortar(fm_int::FaceMatrices{RT}, fm_ext::FaceMatrices{RT}, msize::MortarSize,
                Am::AbstractVector; C::Real = 1.5) where {RT}
    nf_int = size(fm_int.I, 2)
    nf_ext = size(fm_ext.I, 2)
    Nm = max(nf_int, nf_ext)
    length(Am) == Nm ||
        throw(ArgumentError("Mortar: Am has length $(length(Am)), the mortar has $Nm " *
                            "points (the larger of the two abutting resolutions)"))
    p_int = size(fm_int.I, 1) ÷ nf_int - 1
    p_ext = size(fm_ext.I, 1) ÷ nf_ext - 1
    σ = penalty_sigma(C, p_int, p_ext, perpendicular_size(fm_int), perpendicular_size(fm_ext))
    return Mortar{RT,eltype(Am)}(projection_prolong(RT, nf_int, Nm, msize),
                                 projection_restrict(RT, nf_int, Nm, msize),
                                 RT(σ), collect(Am))
end

"""
    MortarFlux(mortars::Vector{Mortar})
    MortarFlux(m::Mortar)

The mortar instance of the CRP internal penalty flux (eq:numflux, main.tex:584-593) — the
[`NumericalFlux`](@ref) of ONE nonconforming element face, holding that face's mortars in
receiver-face order (`[full]`, or `[lower, upper]` on the big side of a 2:1 interface).
Per mortar `m`, with the interior trace projected in and the exterior trace received AS
mortar data (sender convention N3),

    ψ_i = m.P·ψ_int,  φ_i = m.P·φ_int,
    Ψ*_m    = (ψ_i + ψ_e)/2,
    (nF)*_m = (φ_i − φ_e)/2 − m.sigma · m.Am · (ψ_i − ψ_e)    [receiver normal, N3/N4],

and the face values are the summed restrictions `Ψ* = Σ_m m.R·Ψ*_m` (likewise `(nF)*`) —
"compute the numerical flux there, and restrict the result back" (main.tex:678-686). The
exterior `FaceState` must pack the mortar states concatenated in `mortars` order
(`ArgumentError` on any length mismatch). The receiver's own `fm`/`Aface` arguments of the
flux interface are ignored — σ and 𝔄 live on the mortar. On a conforming face this is
bitwise `PenaltyFlux` (gated); `C = 0` mortars give the central-flux limit.
"""
struct MortarFlux{RT<:AbstractFloat,TC<:Number} <: NumericalFlux
    mortars::Vector{Mortar{RT,TC}}
    function MortarFlux{RT,TC}(ms::Vector{Mortar{RT,TC}}) where {RT<:AbstractFloat,TC<:Number}
        isempty(ms) && throw(ArgumentError("MortarFlux needs at least one mortar"))
        nf = size(ms[1].P, 2)
        all(size(m.P, 2) == nf for m in ms) ||
            throw(ArgumentError("MortarFlux: mortars disagree on the face point count"))
        return new{RT,TC}(ms)
    end
end
MortarFlux(ms::Vector{Mortar{RT,TC}}) where {RT,TC} = MortarFlux{RT,TC}(ms)
MortarFlux(m::Mortar) = MortarFlux([m])

function numerical_flux(mf::MortarFlux{RT,TC}, ::FaceMatrices, ::AbstractVector,
                        int::FaceState, ext::FaceState) where {RT,TC}
    nf = size(mf.mortars[1].P, 2)
    length(int.psi) == nf ||
        throw(ArgumentError("MortarFlux: interior face state has length " *
                            "$(length(int.psi)), the mortars expect $nf"))
    nm_tot = sum(size(m.P, 1) for m in mf.mortars)
    length(ext.psi) == nm_tot ||
        throw(ArgumentError("MortarFlux: exterior state has length $(length(ext.psi)), " *
                            "the mortars carry $nm_tot points (packed in mortars order)"))
    T = promote_type(eltype(int.psi), eltype(int.nflux),
                     eltype(ext.psi), eltype(ext.nflux), TC, RT)
    ψstar = zeros(T, nf)
    φstar = zeros(T, nf)
    row = 0
    for m in mf.mortars
        Nm = size(m.P, 1)
        ψe = view(ext.psi, row+1:row+Nm)
        φe = view(ext.nflux, row+1:row+Nm)
        ψi = m.P * int.psi
        φi = m.P * int.nflux
        ψstar .+= m.R * ((ψi .+ ψe) ./ 2)
        φstar .+= m.R * ((φi .- φe) ./ 2 .- m.sigma .* (m.Am .* (ψi .- ψe)))
        row += Nm
    end
    return ψstar, φstar
end

"""
    PerFaceFlux(faces::NTuple{2,NTuple{2,NumericalFlux}})

Per-face flux routing through the single-flux [`strong_residual!`](@ref) seam:
`faces[d][s]` (s: 1 ↔ side −1, 2 ↔ side +1, the wave-4 face indexing) is the
[`NumericalFlux`](@ref) applied to face (d, s), selected by the face-matrix key
`(fm.dim, fm.side)` the assembly already passes. This is the composition device for mixed
meshes — [`MortarFlux`](@ref) on nonconforming faces, `PenaltyFlux` elsewhere — without
touching the wave-4 assembly (its call site is flux-agnostic by convention N6).
"""
struct PerFaceFlux <: NumericalFlux
    faces::NTuple{2,NTuple{2,NumericalFlux}}
end

numerical_flux(pff::PerFaceFlux, fm::FaceMatrices, Aface::AbstractVector,
               int::FaceState, ext::FaceState) =
    numerical_flux(pff.faces[fm.dim][fm.side == -1 ? 1 : 2], fm, Aface, int, ext)
