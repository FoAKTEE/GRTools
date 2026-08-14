# Generic-spin (s ∈ {0, ±1, ±2}) VACUUM m-mode ℓ-coupled hyperboloidal DG evolution in the
# RT chart (DG-P4; chart decision + projected 1+1 system in rtchart.jl header). Generalizes the
# KerrModes machinery (kerrmodes.jl): same nodal-DG + LLF + RK4 skeleton on NumericsCore.dg,
# but (i) mass block E_ππ = 𝒜_TT is PENTADIAGONAL (ℓ±1 AND ℓ±2) when s·m ≠ 0, (ii) the state
# is COMPLEX (friction/advection/potential carry imaginary parts for a·s ≠ 0 or m ≠ 0),
# (iii) both boundaries are outflow (horizon-penetrating + compactified scri: C_RR(R₊) =
# C_RR(1) = 0, [arXiv:1301.1591 §2.3]) so NO boundary condition is imposed.
#
# FIRST-ORDER REDUCTION (DERIVED-RT-2; repo-pattern reduction variables, cf. [arXiv:2307.01349,
# Eq. (15)] — do not transform them, dossier §3.2): π := −∂_T ψ, χ := ∂_R ψ. Substituting into
# the projected 1+1 system (rtchart.jl header) with ∂_TR ψ = ∂_R(∂_T ψ) = −∂_R π:
#
#   ∂_T ψ⃗ = −π⃗
#   𝒜_TT ∂_T π⃗ + C_TR ∂_R π⃗ − C_RR ∂_R χ⃗ = −𝒜_T π⃗ + C̃_R χ⃗ + diag(C̄₀) ψ⃗
#   ∂_T χ⃗ = −∂_R π⃗
#
# i.e. E U̇ + Â U′ + B̂(U) + V̂(ψ) = 0 with U = [π⃗; χ⃗], E = diag(𝒜_TT, I) (symmetric),
# Â = [[C_TR·I, −C_RR·I],[I, 0]]. The flux matrix A = E⁻¹Â = [[C_TR·W, −C_RR·W],[I, 0]],
# W := 𝒜_TT⁻¹, is REAL (the principal part is s- and m-independent up to the real 𝒜_TT
# coupling; [arXiv:1301.1591, §2.3: "The principal part of the equation is independent of the
# spin weight"]). Characteristic speeds solve λ²𝒜_TT − λC_TR + C_RR = 0; in the decoupled
# (a=0) limit λ± = (C_TR ± √(C_TR² − 4·C_TT·C_RR))/(2C_TT) = the radial coordinate light
# speeds c_R,± of [arXiv:1301.1591, Eq. (19)], gated against their Eqs. (17),(18) values.
# W's regularity on [R₊, 1] (the E-inversion that DIVERGES at scri in the repo layer chart for
# s ≠ 0, dossier §3.4) is asserted at build time, and hyperbolicity (real eigenvalues of A) is
# asserted node-by-node (empirical stand-in for the open s≠0 proof, dossier §8.4).

"""
    SpinModesProblem{T}

Precomputed generic-s RT-chart coupled-mode vacuum problem (DG-P4): per-node real flux
matrices `Aflux` (2n×2n), complex lower-order blocks `Gpi = W𝒜_T`, `Gchi = C̃_R·W`,
`Gpsi = W·diag(C̄₀)`, interface flux matrices + LLF speeds, probe interpolation rows.
Build with [`SpinModesProblem(s, ℓs, m; ...)`](@ref); evolve with [`solve_spin_modes`](@ref).
"""
struct SpinModesProblem{T<:AbstractFloat}
    s::Int
    m::Int
    ℓs::Vector{Int}
    n::Int
    el::RefElement{T}
    mesh::Mesh{T}
    M::T
    a::T
    Rc::Matrix{T}                 # node coordinates R (Np, K)
    Aflux::Array{T,4}             # (2n,2n,Np,K) real flux matrix A = E⁻¹Â
    Gpi::Array{Complex{T},4}      # (n,n,Np,K)  W·𝒜_T
    Gchi::Array{Complex{T},4}     # (n,n,Np,K)  C̃_R·W
    Gpsi::Array{Complex{T},4}     # (n,n,Np,K)  W·diag(C̄₀^{(ℓ)})
    invMc1::Vector{T}
    invMcN::Vector{T}
    twoOverh::Vector{T}
    Aint::Array{T,3}              # (2n,2n,K-1) flux matrix at interior interfaces
    αint::Vector{T}               # LLF dissipation speed per interior interface
    λmax::T                       # global max characteristic speed (CFL)
    Fil::Matrix{T}                # per-element exponential modal filter (identity if disabled)
    use_filter::Bool
    κc::T                         # auxiliary-constraint damping rate (χ − ∂_Rψ sector)
    ψjump_toff::T                 # ψ-jump penalty active for τ < ψjump_toff (transient epoch)
    probe_R::Vector{T}
    probe_k::Vector{Int}
    probe_L::Matrix{T}
end

"""
    SpinModesProblem(s, ℓs, m; N=8, K=48, M=1, a=0, excise=0.05, filter_alpha=36.0,
                     filter_power=16, constraint_kappa=0.3, psi_jump_toff=250.0,
                     extract_R=T[], forced_edges=T[], T=Float64)

Build the generic-s vacuum problem for spin weight `s ∈ {0,±1,±2}`, azimuthal number `m` and
the coupled-mode list `ℓs` (consecutive required when `s·a ≠ 0`: the β± ∝ a·s ℓ±1 couplings
break the parity split, rtchart.jl header). Domain is `R ∈ [R₊ − excise, 1]` — the RT chart
is horizon-penetrating [arXiv:1301.1591, §2.1], and the inner boundary is EXCISED a
coordinate distance `excise` INSIDE the horizon (standard excision):

  - for R < R₊ BOTH radial characteristic speeds are negative (verified at build), so the
    inner boundary is pure outflow exactly like R₊ itself — still no boundary condition;
  - it moves the Δ = 0 point of the radial operator (where the nonphysical branch
    `ψ ~ Δ^{s}e^{−iωu}` at H⁺ [arXiv:1702.04204, Eq. (51); dossier §4] is anchored, and
    where the outgoing characteristic speed vanishes) off the boundary, removing the
    boundary-pinned variant of the noise-seeded growth (MEASURED: s=−2, a=0, boundary AT
    R₊: semi-discrete max Re λ = 0.0099 → 0.034 growing with N,K). The interior
    sonic-point glancing families for s ≠ 0 remain and are handled by the modal filter
    (`Fil` comment below) — together these address the instability class HBB flag as
    "exponentially growing continuum modes triggered by numerical noise" [arXiv:1301.1591,
    §2.3 footnote], which they held below tolerance with quadruple precision.

Numerical-dissipation controls (each measured, see inline comments): `filter_alpha`/
`filter_power` (exponential modal filter; 0 disables), `constraint_kappa` (χ−∂_Rψ damping),
`psi_jump_toff` (transient-epoch ψ-jump penalty cutoff time).

R₊ is a FORCED mesh edge, so horizon extraction stays exact. `extract_R` are extraction
radii appended to the default probes `[R₊, (R₊+1)/2, 1]`. `forced_edges` are additional
forced mesh interfaces (DG-P5 hook: the sourced sector forces the particle radius R₀ onto
the mesh so the jump-conditioned interface flux sits exactly at the orbit; default empty =
bit-identical vacuum meshes).
"""
function SpinModesProblem(s::Integer, ℓs::Vector{Int}, m::Integer; N=8, K=48, M=1.0, a=0.0,
                          excise=0.05, filter_alpha=36.0, filter_power=16, constraint_kappa=0.3,
                          psi_jump_toff=250.0, extract_R=Float64[], forced_edges=Float64[],
                          T::Type{<:AbstractFloat}=Float64)
    abs(s) <= 2 || throw(ArgumentError("spin weight |s| ≤ 2 (got $s)"))
    issorted(ℓs) || throw(ArgumentError("ℓs must be sorted"))
    minimum(ℓs) >= max(abs(s), abs(m)) ||
        throw(ArgumentError("ℓ ≥ max(|s|,|m|) required (got ℓmin=$(minimum(ℓs)))"))
    if s != 0 && a != 0 && length(ℓs) > 1 && any(diff(ℓs) .!= 1)
        throw(ArgumentError("s·a ≠ 0 breaks the parity split (β± ∝ a·s): ℓs must be consecutive"))
    end
    excise >= 0 || throw(ArgumentError("excise ≥ 0 required"))
    n = length(ℓs)
    el = RefElement(N; T=T)
    Rp = T(rt_Rplus(M, a))
    Rin = Rp - T(excise)
    Rin > 0 || throw(ArgumentError("excise too large (inner edge ≤ 0)"))
    mesh = excise == 0 ? build_mesh(Rp, one(T), K; forced=T.(forced_edges), T=T) :
           build_mesh(Rin, one(T), K; forced=[Rp; T.(forced_edges)], T=T)
    Rc = permutedims(node_coords(mesh, el))
    Np, Kk = el.Np, mesh.K
    Aflux = Array{T,4}(undef, 2n, 2n, Np, Kk)
    Gpi = Array{Complex{T},4}(undef, n, n, Np, Kk)
    Gchi = Array{Complex{T},4}(undef, n, n, Np, Kk)
    Gpsi = Array{Complex{T},4}(undef, n, n, Np, Kk)
    λmax = zero(T)
    @inbounds for k in 1:Kk, i in 1:Np
        R = Rc[i, k]
        ATT = rt_ATT(R, s, ℓs, m; M=M, a=a)
        # E-inversion regularity gate (dossier §3.4): 𝒜_TT must be invertible on ALL of
        # [R₊, 1] including scri — this is the property the RT chart buys for s ≠ 0.
        W = inv(ATT)
        all(isfinite, W) || error("singular 𝒜_TT at R=$R (chart regularity violated)")
        ctr = rt_CTR(R; M=M, a=a)
        crr = rt_CRR(R; M=M, a=a)
        A = zeros(T, 2n, 2n)
        A[1:n, 1:n] .= ctr .* W
        A[1:n, n+1:2n] .= -crr .* W
        for q in 1:n
            A[n+q, q] = one(T)
        end
        Aflux[:, :, i, k] = A
        # hyperbolicity (real characteristic speeds) asserted empirically — the s≠0
        # pentadiagonal well-posedness proof is an open item (dossier §8.4).
        ev = eigen(A).values
        maximum(abs.(imag.(ev))) <= 1e-10 * (1 + maximum(abs.(ev))) ||
            error("complex characteristic speeds at R=$R (hyperbolicity lost)")
        # excision outflow check: INSIDE the horizon all characteristic speeds must be ≤ 0
        # (trapped region — the excised inner boundary is then pure outflow, no BC needed)
        R < Rp - 1e-12 && maximum(real.(ev)) > 1e-12 &&
            error("outgoing characteristic inside the horizon at R=$R (excision unsound)")
        λmax = max(λmax, maximum(abs.(real.(ev))))
        AT = rt_AT(R, s, ℓs, m; M=M, a=a)
        Gpi[:, :, i, k] = W * AT
        Gchi[:, :, i, k] = rt_CtildeR(R, m; M=M, a=a, s=s) .* W
        # W·diag(C̄₀): column-c of W scaled by C̄₀^{(ℓ_c)}
        Gpsi[:, :, i, k] = W .* transpose([rt_C0diag(R, s, l, m; M=M, a=a) for l in ℓs])
    end
    Aint = Array{T,3}(undef, 2n, 2n, Kk - 1)
    αint = Vector{T}(undef, Kk - 1)
    for k in 1:(Kk-1)
        Aint[:, :, k] = Aflux[:, :, Np, k]
        αint[k] = maximum(abs.(real.(eigen(Aint[:, :, k]).values)))
    end
    # Exponential modal filter σᵢ = exp(−α((i−1)/N)^power) (Hesthaven–Warburton-standard DG
    # stabilization), applied per element after each RK4 step. NEEDED for s ≠ 0: the RT-chart
    # radial sonic points (λ₊(R₊) = 0, λ₋(1) = 0) carry anti-damped glancing characteristic
    # families for s ≠ 0 — the local counterparts of the nonphysical branches ψ ~ Δˢ at H⁺ and
    # ψ ~ r^{2s} at ℐ⁺ [arXiv:1702.04204, Eqs. (49),(51); dossier §4] — which seed
    # resolution-dependent discrete growth from roundoff (the instability class flagged in
    # [arXiv:1301.1591, §2.3 footnote]; HBB suppressed it below tolerance with quadruple
    # precision, unavailable here). The filter damps only the top modal tail (smooth physical
    # solutions are spectrally unaffected); its effect on the physics is gated by the s=0
    # cross-check (gate A) and the Price-tail gates (gate C). `filter_alpha=0` disables.
    Fil = Matrix{T}(one(T) * Id, Np, Np)
    if filter_alpha > 0
        σv = [exp(-T(filter_alpha) * (T(i - 1) / N)^filter_power) for i in 1:Np]
        Fil = el.V * (σv .* inv(el.V))
    end
    probes = T[Rp, (Rp + 1) / 2, one(T)]
    append!(probes, T.(extract_R))
    pk = [element_of(mesh, x) for x in probes]
    pL = Matrix{T}(undef, length(probes), Np)
    for (ip, x) in enumerate(probes)
        kk = pk[ip]
        ξ = 2 * (x - mesh.xL[kk]) / mesh.h[kk] - 1
        pL[ip, :] = _lagrange_row(el, ξ)
    end
    return SpinModesProblem{T}(Int(s), Int(m), ℓs, n, el, mesh, T(M), T(a), Rc, Aflux,
        Gpi, Gchi, Gpsi, el.invM[:, 1], el.invM[:, el.Np], 2 ./ mesh.h, Aint, αint, λmax,
        Fil, filter_alpha > 0, T(constraint_kappa), T(psi_jump_toff), probes, pk, pL)
end

"apply the per-element modal filter to every component/mode of the state (in place)."
function _apply_filter!(u, p::SpinModesProblem{T}, buf) where {T}
    p.use_filter || return u
    Np, K, n = p.el.Np, p.mesh.K, p.n
    U = reshape(u, Np, K, 3, n)
    F = p.Fil
    @inbounds for q in 1:n, c in 1:3, k in 1:K
        for i in 1:Np
            acc = zero(Complex{T})
            for j in 1:Np
                acc += F[i, j] * U[j, k, c, q]
            end
            buf[i] = acc
        end
        for i in 1:Np
            U[i, k, c, q] = buf[i]
        end
    end
    return u
end

"""
    spin_cfl_dt(p; C=0.6) -> dt

RK4/DG time step `dt = C·Δx_min/λ_max` from the smallest LGL node spacing and the measured
maximum characteristic speed (RT-chart speeds are ≪ 1, [arXiv:1301.1591, Eqs. (17),(18)]).
"""
function spin_cfl_dt(p::SpinModesProblem; C=0.6)
    Δ = Inf
    for k in 1:p.mesh.K
        hk = p.mesh.h[k]
        for i in 1:(p.el.Np-1)
            Δ = min(Δ, hk * (p.el.r[i+1] - p.el.r[i]) / 2)
        end
    end
    return cfl_dt(C, Δ, p.λmax)
end

function _rhs!(du, u, t::Real, p::SpinModesProblem{T}, sc) where {T}
    el = p.el
    Np = el.Np
    K = p.mesh.K
    n = p.n
    U = reshape(u, Np, K, 3, n)
    dU = reshape(du, Np, K, 3, n)
    Dπ, Dχ, Dψ, tmpπ, tmpχ, Fπ, Fχ, FπR, FχR = sc
    D = el.D
    A4 = p.Aflux
    Gpi4 = p.Gpi
    Gch4 = p.Gchi
    Gps4 = p.Gpsi
    κc = p.κc
    @inbounds for k in 1:K
        invJ = p.twoOverh[k]
        for q in 1:n, i in 1:Np           # hand-rolled D·U (avoids per-call mul!/view overhead)
            aπ = zero(Complex{T})
            aχ = zero(Complex{T})
            aψ = zero(Complex{T})
            for j in 1:Np
                aπ += D[i, j] * U[j, k, 2, q]
                aχ += D[i, j] * U[j, k, 3, q]
                aψ += D[i, j] * U[j, k, 1, q]
            end
            Dπ[i, q] = invJ * aπ
            Dχ[i, q] = invJ * aχ
            Dψ[i, q] = invJ * aψ
        end
        for i in 1:Np
            for q in 1:n
                # flux Â U′ row q: Fπ = Σ_c A[q,c]π′ + A[q,n+c]χ′ ; Fχ row = π′ (A[n+q,·]=e_q)
                s1 = zero(Complex{T})
                for c in 1:n
                    s1 += A4[q, c, i, k] * Dπ[i, c] + A4[q, n+c, i, k] * Dχ[i, c]
                end
                # lower order: π̇ ⊃ −(W𝒜_T)π + (C̃_R W)χ + (W diag(C̄₀))ψ  (reduction, header)
                w = zero(Complex{T})
                for c in 1:n
                    w += -Gpi4[q, c, i, k] * U[i, k, 2, c] + Gch4[q, c, i, k] * U[i, k, 3, c] +
                         Gps4[q, c, i, k] * U[i, k, 1, c]
                end
                dU[i, k, 1, q] = -U[i, k, 2, q]
                dU[i, k, 2, q] = -s1 + w
                # χ̇ = −π′ − κc·(χ − ∂_Rψ): the reduction constraint Z := χ − ∂_Rψ satisfies
                # Ż = 0 EXACTLY in the continuum system (frozen sector), so any numerically
                # injected Z (modal filter, discretization noise) would persist as a static
                # residue sourcing π̇ through C̃_R. Standard auxiliary-constraint damping
                # (continuum-consistent: vanishes on Z ≡ 0; strictly dissipative for Z;
                # leaves the principal part/characteristics unchanged — the added ∂_Rψ
                # column has zero characteristic contribution, χ-row cofactor vanishes).
                dU[i, k, 3, q] = -Dπ[i, q] - κc * (U[i, k, 3, q] - Dψ[i, q])
            end
        end
    end
    # interior interface LLF (α = local max characteristic speed, measured at build)
    @inbounds for k in 1:(K-1)
        A = view(p.Aint, :, :, k)
        α = p.αint[k]
        for q in 1:n
            tmpπ[q] = U[Np, k, 2, q]
            tmpχ[q] = U[Np, k, 3, q]
        end
        _AUvec!(Fπ, Fχ, A, tmpπ, tmpχ, n)          # F(U^L)
        for q in 1:n
            tmpπ[q] = U[1, k+1, 2, q]
            tmpχ[q] = U[1, k+1, 3, q]
        end
        _AUvec!(FπR, FχR, A, tmpπ, tmpχ, n)        # F(U^R)
        ck = p.twoOverh[k]
        ck1 = p.twoOverh[k+1]
        for q in 1:n
            uLπ = U[Np, k, 2, q]; uLχ = U[Np, k, 3, q]
            uRπ = U[1, k+1, 2, q]; uRχ = U[1, k+1, 3, q]
            Fsπ = (Fπ[q] + FπR[q]) / 2 - (α / 2) * (uRπ - uLπ)
            Fsχ = (Fχ[q] + FχR[q]) / 2 - (α / 2) * (uRχ - uLχ)
            dRπ = Fπ[q] - Fsπ; dRχ = Fχ[q] - Fsχ
            dLπ = FπR[q] - Fsπ; dLχ = FχR[q] - Fsχ
            # ψ-jump dissipation (LLF penalty of the ZERO flux, same α), TIME-GATED to the
            # transient epoch τ < ψjump_toff: ψ carries no flux of its own (ψ̇ = −π), so
            # interelement ψ-jumps enter NO other term — the scheme otherwise has an EXACT
            # static kernel of elementwise static-ODE solutions stitched by arbitrary ψ-jumps
            # (MEASURED: filter-off step-operator modes with rate ~1e-12, ψ-jump O(1),
            # π ~ 1e-14; excited by the IC transient at ~1e-5 rel, freezing every sub-1e-5
            # tail measurement). The penalty is consistency-preserving ([ψ] = 0 in continuum)
            # and drains that kernel; it is switched OFF (30M linear ramp) after the transient
            # because the LATE-time tail profile is step-like (loc. cit. [arXiv:1301.1591,
            # Fig. 6 discussion]) and a permanent ψ-jump penalty converts its finite
            # resolution into artificial exponential damping (MEASURED: LPI → −20s).
            γψ = t < p.ψjump_toff ? one(T) :
                 t < p.ψjump_toff + 30 ? one(T) - (T(t) - p.ψjump_toff) / 30 : zero(T)
            dψ = γψ * (α / 2) * (U[1, k+1, 1, q] - U[Np, k, 1, q])
            for i in 1:Np
                dU[i, k, 1, q] += ck * p.invMcN[i] * dψ
                dU[i, k+1, 1, q] -= ck1 * p.invMc1[i] * dψ
                dU[i, k, 2, q] += ck * p.invMcN[i] * dRπ
                dU[i, k, 3, q] += ck * p.invMcN[i] * dRχ
                dU[i, k+1, 2, q] -= ck1 * p.invMc1[i] * dLπ
                dU[i, k+1, 3, q] -= ck1 * p.invMc1[i] * dLχ
            end
        end
    end
    # physical boundaries R₊ and 1: pure OUTFLOW (C_RR = 0 there ⇒ no incoming characteristic,
    # [arXiv:1301.1591, §2.3, Eqs. (16)–(19)]): numerical flux = interior flux ⇒ zero penalty,
    # no boundary term (no condition imposed).
    return du
end

"""
    solve_spin_modes(p; tfinal=1500.0, dt=nothing, record_every=100, init_mode=1,
                     id=:ID0, ic=(amp=1.0, R0=0.8, w=3000.0)) -> (τs, Ψ, maxu)

Vacuum evolution with pure-multipole initial data in mode `ℓs[init_mode]`
[arXiv:1301.1591, §4, ID definitions + `G(R) = e^{−(w/2)(R−R₀)²}`]:

  - `id = :ID0` (stationary compact):     ψ = G(R), ∂_Tψ = 0  (π = 0, χ = ∂_R G consistent)
  - `id = :ID1` (non-stationary compact): ∂_Tψ = G(R)         (π = −G, ψ = χ = 0)

`dt = nothing` uses [`spin_cfl_dt`](@ref). Returns recorded times `τs`, complex probe series
`Ψ[probe, mode, sample]` of ψ at `p.probe_R` (defaults: horizon R₊, mid-domain, scri), and
`maxu[sample] = max|U|` over the full state (gate E boundedness monitor).
"""
function solve_spin_modes(p::SpinModesProblem{T}; tfinal=1500.0, dt=nothing,
                          record_every=100, init_mode=1, id::Symbol=:ID0,
                          ic=(amp=1.0, R0=0.8, w=3000.0)) where {T}
    Np, K, n = p.el.Np, p.mesh.K, p.n
    h = dt === nothing ? spin_cfl_dt(p) : dt
    u = zeros(Complex{T}, Np * K * 3 * n)
    U = reshape(u, Np, K, 3, n)
    @inbounds for k in 1:K, i in 1:Np
        R = p.Rc[i, k]
        G = ic.amp * exp(-(ic.w / 2) * (R - ic.R0)^2)
        if id === :ID0
            U[i, k, 1, init_mode] = G
            U[i, k, 3, init_mode] = -ic.w * (R - ic.R0) * G   # χ = ∂_R ψ (consistency)
        elseif id === :ID1
            U[i, k, 2, init_mode] = -G                         # π = −∂_T ψ = −G
        else
            throw(ArgumentError("id must be :ID0 or :ID1"))
        end
    end
    sc = (zeros(Complex{T}, Np, n), zeros(Complex{T}, Np, n), zeros(Complex{T}, Np, n),
          zeros(Complex{T}, n), zeros(Complex{T}, n), zeros(Complex{T}, n),
          zeros(Complex{T}, n), zeros(Complex{T}, n), zeros(Complex{T}, n))
    nprb = length(p.probe_R)
    τs = T[]
    cols = Array{Complex{T},2}[]
    maxu = T[]
    fbuf = zeros(Complex{T}, Np)
    rec(uu, t, step) = begin
        step > 0 && _apply_filter!(uu, p, fbuf)   # modal filter after every accepted step
        if step % record_every == 0
            UU = reshape(uu, Np, K, 3, n)
            push!(τs, t)
            v = Matrix{Complex{T}}(undef, nprb, n)
            for ip in 1:nprb, q in 1:n
                v[ip, q] = dot(view(p.probe_L, ip, :), view(UU, :, p.probe_k[ip], 1, q))
            end
            push!(cols, v)
            push!(maxu, maximum(abs, uu))
        end
    end
    rhs!(du, uu, t) = _rhs!(du, uu, t, p, sc)
    rk4!(rhs!, u, (0.0, tfinal), h; callback=rec)
    Ψ = Array{Complex{T},3}(undef, nprb, n, length(τs))
    for j in eachindex(cols), ip in 1:nprb, q in 1:n
        Ψ[ip, q, j] = cols[j][ip, q]
    end
    return τs, Ψ, maxu
end
