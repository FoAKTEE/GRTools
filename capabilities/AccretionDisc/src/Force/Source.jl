# SOURCE — effective metric force f_μ → 𝒬-projected source field V^μ → full radial forcing 𝔉(r)
# (design §4.10 Force: "fluid-frame effective force → radial source F + grids").
#
# Physics (arXiv:2601.19123): the metric perturbation h_μν exerts an effective force on the disc fluid,
#   f_μ = (ε+p) U^α U^β (½ ∂_μ h_αβ − ∇_β h_μα) + h_μν ∂^ν p ,
# projected with 𝒬 (the disc response operator, Modes.q_operator_z) and the spatial projector
# prj = g^{μν}+U^μU^ν to give the source field V^μ, then the full radial forcing 𝔉(r) of the master
# ODE (eq:source) including the height-averaged θ-divergence and the trace/enthalpy term.
#
# AHA reuse: geometry from Spacetime.{kerr_metric,kerr_metric_inv,christoffel} (the analytic connection —
# the legacy FD kerr_christoffel fork is collapsed, per the Modes.Velocity dedup decision); disc background
# + off-equator structure + 𝒬 + height-average from ..Modes; high-order FD / Chebyshev denoise from
# NumericsCore (the legacy DiscEMRIForce/numerics.jl already lives there — fornberg/highorder_deriv/
# cheb_deriv/cheb_smooth_deriv). TYPE-GENERIC {T}: derives T from the loaded mode; the validated fixed
# FD steps are floored representable at any precision. From DiscEMRIForce/source.jl.

using ..Modes: RelativisticDisc, disc_OmegaK, disc_OmegaK_z, csound, eden0, eden_z, pressure_z,
               Uvec_z, alpha_z, q_operator_z, heightavg, lindblad_locations
import Spacetime: kerr_metric, kerr_metric_inv, christoffel
import NumericsCore: highorder_deriv, cheb_deriv, cheb_smooth_deriv

# ---- finite-difference steps, precision-aware --------------------------------------------------------------
# The validated h1lorenz pipeline uses FIXED FD steps (the metric harmonics carry ~1e-5 FD noise already, so
# a tiny eps(T)-scaled step would AMPLIFY that noise). Validated step sizes kept: 1e-4 metric θ, 1e-6 radial.
_fdr(::Type{T})  where {T} = T(1)/T(10)^6          # radial FD step for background scalars (validated)
_fdθh(::Type{T}) where {T} = T(1)/T(10)^4          # metric θ-FD step in the force (validated)

"Disc background matching the validated a=0 pipeline (a=0, h0=0.1, r_I=1, eden_I=1; working type T)."
disc_a0(::Type{T}=Float64) where {T} = RelativisticDisc(a=T(0), h0=T(1)/T(10), r_I=one(T), eden_I=one(T))

"""
    disc(a; T=Float64, h0=1/10, r_I=1, eden_I=1) -> RelativisticDisc{T}

Spin-parametrized disc background (the spin-general generalization of `disc_a0`): the validated
thermodynamics (constant aspect ratio h0=0.1, ē_I=1 at r_I=1) at arbitrary Kerr spin `a`. The Ω^k,
κ_r, c_s, ISCO and the master coefficients A_m/B_m it feeds are ALL spin-general, so `disc(0)` is
byte-identical to `disc_a0`.
"""
disc(a; T::Type=Float64, h0=T(1)/T(10), r_I=one(T), eden_I=one(T)) =
    RelativisticDisc(a=T(a), h0=T(h0), r_I=T(r_I), eden_I=T(eden_I))

# =============================================================================================================
# EFFECTIVE FORCE f_μ (lower index), parity-aware, off-equator capable.
# =============================================================================================================

"""
    force_mode(d, mode, k, l, m, ω0, z=0; drH_override=nothing) -> Complex{T} vector of length 4

Effective force f_μ from a single (l,m) metric-perturbation mode at radial index `k`, height `z`
(θ=acos z). `ω0 = m Ω(r_p)` is the pattern frequency. `drH_override` supplies a precomputed high-order
∂_r h_μν (the wing-precision array pass); default is 2-point central FD over the grid neighbours k±1.
From DiscEMRIForce/source.jl (h1lorenz force_pa / force_z).
"""
function force_mode(d::RelativisticDisc, mode::H1Mode{T}, k::Integer, l::Integer, m::Integer, ω0, z=zero(T);
                    drH_override=nothing) where {T}
    a = T(d.a); r = mode.r_right[k]; θ = acos(T(z))
    Γ = christoffel(a, r, θ); gi = kerr_metric_inv(a, r, θ)
    U = Uvec_z(d, r, T(z)); cs2 = csound(d, r)^2; εp = eden_z(d, r, T(z))*(1 + cs2)
    # ∂^ν p (raised): ∂_r p, ∂_θ p (via z=cosθ), others 0
    dpr = _fdr(T); hθ = _fdθh(T)
    pof(s, zz) = pressure_z(d, s, zz)
    dpr_ = (pof(r+dpr, T(z)) - pof(r-dpr, T(z)))/(2dpr)
    dpθ  = z == 0 ? zero(T) : (pof(r, cos(θ+hθ)) - pof(r, cos(θ-hθ)))/(2hθ)
    flo_p = Complex{T}[0, dpr_, dpθ, 0]; dpup = gi*flo_p
    # h_μν and its derivatives: ∂_t=−iω₀, ∂_r (2-pt central FD over k±1, edge-clamped, unless a
    # precomputed high-order drH_override is supplied), ∂_θ (central FD), ∂_φ=+im
    H = hmunu(mode, k, l, m, θ; a=a)
    if drH_override === nothing
        kp = min(k+1, length(mode.r_right)); km = max(k-1, 1)   # clamp the FD stencil to the grid (edge-safe)
        rp = mode.r_right[kp]; rm = mode.r_right[km]
        drH = (hmunu(mode, kp, l, m, θ; a=a) .- hmunu(mode, km, l, m, θ; a=a)) ./ (rp - rm)
    else
        drH = drH_override                                     # precomputed high-order ∂_r h_μν (array pass)
    end
    dθH = (hmunu(mode, k, l, m, θ+hθ; a=a) .- hmunu(mode, k, l, m, θ-hθ; a=a)) ./ (2hθ)
    dh = (-im*ω0*H, drH, dθH, im*m*H)
    # covariant deriv pieces C_μ[α,β] = ∂_μ h_αβ − Γ^λ_{μα}h_λβ − Γ^λ_{μβ}h_αλ
    C = [zeros(Complex{T}, 4, 4) for _ in 1:4]
    for μ in 1:4, α in 1:4, β in 1:4
        s = dh[μ][α,β]
        for λ in 1:4; s -= Γ[λ,μ,α]*H[λ,β] + Γ[λ,μ,β]*H[α,λ]; end
        C[μ][α,β] = s
    end
    f = zeros(Complex{T}, 4)
    for μ in 1:4
        t1 = sum(U[α]*U[β]*(C[μ][α,β]/2 - C[β][μ,α]) for α in 1:4, β in 1:4)
        f[μ] = εp*t1 + sum(H[μ,ν]*dpup[ν] for ν in 1:4)
    end
    f
end

# =============================================================================================================
# 𝒬-PROJECTED SOURCE FIELD V^μ (sum over l), off-equator capable.
# =============================================================================================================

"""
    Vvec(d, modes, k, m, r_p, ω0, z=0; lmax, drH_over=nothing) -> Complex{T} vector of length 4

The 𝒬-projected source field V^μ = Σ_l 𝒬[prj·f_l]/(ε+p) at radial index `k`, height `z`.
`modes` is a Dict l=>H1Mode; sums l=m..lmax. The projector prj (=g^{να}+U^νU^α) acts on the LOWER
force f_α (the iter287-fixed index placement, q^ν_α f^α of eq:VelocityReconstruction). `drH_over`
(when supplied) = Dict l => high-order ∂_r h_μν at index k (the wing-precision array pass).
From DiscEMRIForce/source.jl (h1lorenz Vvec_z).
"""
function Vvec(d::RelativisticDisc, modes::Dict{Int,H1Mode{T}}, k::Integer, m::Integer, r_p, ω0, z=zero(T);
             lmax::Integer, drH_over=nothing) where {T}
    a = T(d.a); r = modes[m].r_right[k]; θ = acos(T(z))
    gi = kerr_metric_inv(a, r, θ); U = Uvec_z(d, r, T(z)); prj = gi .+ U*U'
    Q = q_operator_z(d, r, T(z), r_p, m); εp = eden_z(d, r, T(z))*(1 + csound(d, r)^2)
    V = zeros(Complex{T}, 4)
    for l in m:lmax
        haskey(modes, l) || continue
        dro = drH_over === nothing ? nothing : get(drH_over, l, nothing)
        f = force_mode(d, modes[l], k, l, m, ω0, z; drH_override=dro)
        V .+= (Q*(prj*f))/εp
    end
    V
end

# =============================================================================================================
# FULL RADIAL FORCING 𝔉(r) (height-averaged θ-divergence + enthalpy/trace term).
# =============================================================================================================

"""
    source_F(d, modes, k, m, r_p, ω0; lmax, h0=d.h0) -> Complex{T}

Full radial forcing 𝔉(r) of the master ODE at radial index `k` (eq:source): the c_s²-weighted
divergence of the 𝒬-projected source field V^μ (radial + t,φ + height-averaged θ pieces) plus the
R_r V^r drift and the enthalpy trace term i α ω̃ c_s²/2 · g^{μν}h_μν.
From DiscEMRIForce/source.jl (h1lorenz full_F.jl::source_F).
"""
function source_F(d::RelativisticDisc, modes::Dict{Int,H1Mode{T}}, k::Integer, m::Integer, r_p, ω0;
                  lmax::Integer, h0=T(d.h0)) where {T}
    a = T(d.a); r = modes[m].r_right[k]; cs2 = csound(d, r)^2
    sg(z) = (r^2 + a^2*z^2)*sqrt(1 - z^2)                        # √g = Σ sinθ (spin-general; a→0 ⇒ r²)
    Veq = Vvec(d, modes, k, m, r_p, ω0, zero(T); lmax=lmax)
    kp = min(k+1, length(modes[m].r_right)); km = max(k-1, 1)   # clamp the radial-divergence FD stencil
    rp = modes[m].r_right[kp]; rm = modes[m].r_right[km]
    # radial divergence (1/√g) ∂_r(√g V^r), equatorial
    SgVr(kk) = modes[m].r_right[kk]^2 * Vvec(d, modes, kk, m, r_p, ω0, zero(T); lmax=lmax)[2]
    divr = (SgVr(kp) - SgVr(km))/(rp - rm)/r^2
    divtφ = (-im*ω0)*Veq[1] + (im*m)*Veq[4]
    # θ-divergence, height-averaged: ⟨(1/√g)∂_θ(√g V^θ)⟩, ∂_θ = −sinθ ∂_z
    dz = T(2)/T(10)^2
    SgVθ(z) = sg(z)*Vvec(d, modes, k, m, r_p, ω0, z; lmax=lmax)[3]
    gθ(z) = (-sqrt(1 - z^2)*(SgVθ(z+dz) - SgVθ(z-dz))/(2dz))/(r^2*sqrt(1 - z^2))
    divθ_ha = heightavg(z -> gθ(z), h0; hz=dz)
    # trace h_eq + R_r drift
    H = sum(hmunu(modes[l], k, l, m, T(π)/2; a=a) for l in m:lmax if haskey(modes, l))
    gi = kerr_metric_inv(a, r, T(π)/2)
    htr = sum(gi[μ,ν]*H[μ,ν] for μ in 1:4, ν in 1:4)
    α = alpha_z(d, r, zero(T)); ω̃ = ω0 - m*disc_OmegaK(d, r)    # spin-general U^t (a→0 ⇒ 1/√(1−3/r))
    dpr = _fdr(T); pof(s) = csound(d, s)^2*eden0(d, s)
    Rr = (pof(r+dpr) - pof(r-dpr))/(2dpr)/(eden0(d, r)*(1 + cs2))
    cs2*(divr + divtφ + divθ_ha) + Rr*Veq[2] + im*α*ω̃*cs2/2*htr
end

# =============================================================================================================
# WING-PRECISION ARRAY PASS — high-order FD (+ optional Chebyshev spectral) of the source on the WHOLE grid.
# The per-point source_F uses two compounding 2-POINT FD layers over adjacent grid neighbours: (A) drH in
# force_mode, (B) divr. source_F_grid recomputes BOTH with wide high-order Fornberg stencils (one pass per
# grid), optionally Chebyshev-denoising smooth pole-free bands; the ILR/OLR/corotation spike cores stay on
# the original per-point FD (byte-preserved). From DiscEMRIForce/source.jl (see there for the full F5-wing
# noise diagnosis: the wing floor is reconstruction noise in the √g·V^r VALUES, denoised by the :lsq
# least-squares Chebyshev fit; :interp is the genuine spectral derivative for CLEAN BigFloat fields).
# =============================================================================================================

"""
    _highorder_drH_grid(mode, l, m, θ, a; npts=9) -> Vector of 4×4 Complex{T}

Per-grid-index high-order ∂_r h_μν^{lm}(r,θ): assemble `hmunu` on the whole `mode.r_right` grid ONCE,
then apply `NumericsCore.highorder_deriv` (Fornberg, `npts`-point) component-wise. Replaces
force_mode's 2-point k±1 metric FD (layer A).
"""
function _highorder_drH_grid(mode::H1Mode{T}, l::Integer, m::Integer, θ, a; npts::Integer=9) where {T}
    g = mode.r_right; N = length(g)
    Harr = [hmunu(mode, k, l, m, θ; a=a) for k in 1:N]          # 4×4 metric on the whole grid (one pass)
    dH = [zeros(Complex{T}, 4, 4) for _ in 1:N]
    comp = Vector{Complex{T}}(undef, N)
    @inbounds for μ in 1:4, ν in 1:4
        for k in 1:N; comp[k] = Harr[k][μ, ν]; end
        for k in 1:N; dH[k][μ, ν] = highorder_deriv(comp, g, k; npts=npts, order=1); end
    end
    dH
end

"""
    source_F_grid(d, modes, m, r_p, ω0; lmax, kwargs...) -> Vector{Complex{T}} (𝔉 on every grid index)

Array-based, WING-PRECISION evaluation of the full radial forcing 𝔉(r) on the entire
`modes[m].r_right` grid. Identical physics to `source_F`, but the two grid-neighbour radial FD layers
are replaced by high-order Fornberg differentiation (`npts`-point): `drH` via `_highorder_drH_grid`,
`divr` via `highorder_deriv` of the precomputed √g·V^r array — or, when `spectral=true`, a Chebyshev
derivative on each smooth pole-free band (`spectral_method=:lsq` the least-squares DENOISER for noisy
Float64 fields — the original F8-byte-identical behaviour; `:interp` the genuine interpolating
spectral derivative for CLEAN BigFloat fields). The ILR/OLR resonance cores and the corotation r_p
are kept on the original per-point FD via `protect_band`/`protect_band_rp` (spike-faithful,
byte-identical to `source_F` there); `smooth_bands` (a Vector of (lo,hi)) targets SPECIFIC pole-free
valleys for denoising with cosine edge-blends, leaving every point outside the bands byte-identical
to the per-point path. From DiscEMRIForce/source.jl.
"""
function source_F_grid(d::RelativisticDisc, modes::Dict{Int,H1Mode{T}}, m::Integer, r_p, ω0;
                       lmax::Integer, h0=T(d.h0), npts::Integer=9, spectral::Bool=false,
                       spectral_method::Symbol=:lsq,
                       smooth_lo=nothing, smooth_hi=nothing, cheb_deg::Integer=8,
                       protect::Bool=true, protect_band=T(8)/T(10),
                       protect_band_rp=nothing, smooth_bands=nothing, cheb_pad=T(4)/T(10),
                       cheb_blend=T(3)/T(10)) where {T}
    a = T(d.a); g = modes[m].r_right; N = length(g)
    ls = [l for l in m:lmax if haskey(modes, l)]
    # (A) high-order ∂_r h_μν on the whole grid, per l, equatorial (θ=π/2)
    drH_grid = Dict(l => _highorder_drH_grid(modes[l], l, m, T(π)/2, a; npts=npts) for l in ls)
    # equatorial V^μ on the whole grid, using the high-order drH; also √g·V^r = r²V^r
    Veq = Vector{Vector{Complex{T}}}(undef, N)
    SgVr = Vector{Complex{T}}(undef, N)
    @inbounds for k in 1:N
        drk = Dict(l => drH_grid[l][k] for l in ls)
        v = Vvec(d, modes, k, m, r_p, ω0, zero(T); lmax=lmax, drH_over=drk)
        Veq[k] = v
        SgVr[k] = g[k]^2 * v[2]
    end
    # PROTECT the resonance poles (ILR/OLR) + the corotation r_p: there |𝔉| has 1/δr SPIKES whose height
    # is grid-sampling-dependent; a wide stencil rings across them. Within the protect bands keep the
    # ORIGINAL per-point computation (spike structure BYTE-PRESERVED).
    rI, rO = lindblad_locations(d, r_p, m)
    pb_rp = protect_band_rp === nothing ? protect_band : T(protect_band_rp)
    isprotected(r) = protect && ((!isnan(rI) && abs(r - T(rI)) <= protect_band) ||
                                 (!isnan(rO) && abs(r - T(rO)) <= protect_band) ||
                                 abs(r - T(r_p)) <= pb_rp)
    # (B) radial divergence (1/r²)∂_r(√g V^r). When `smooth_bands` is given (targeted mode), points
    # outside both the protect bands and the smooth bands keep the ORIGINAL 2-pt FD (no new edge-rings
    # outside the chosen bands); otherwise the non-protected grid gets the high-order Fornberg FD.
    targeted = smooth_bands !== nothing
    insmooth(r) = !targeted || any(lo <= r <= hi for (lo, hi) in
                  (smooth_bands === nothing ? Tuple{T,T}[] : [(T(b[1]), T(b[2])) for b in smooth_bands]))
    divr = Vector{Complex{T}}(undef, N)
    @inbounds for k in 1:N
        if isprotected(g[k]) || (targeted && !insmooth(g[k]))
            kp = min(k+1, N); km = max(k-1, 1)
            divr[k] = (SgVr[kp] - SgVr[km])/(g[kp] - g[km])/g[k]^2     # original 2-pt central FD
        else
            divr[k] = highorder_deriv(SgVr, g, k; npts=npts, order=1) / g[k]^2
        end
    end
    if spectral
        # smooth pole-free DENOISE bands: default = the single far-wing window [smooth_lo,smooth_hi]
        bands = if smooth_bands !== nothing
            [(T(b[1]), T(b[2])) for b in smooth_bands]
        else
            lo = smooth_lo === nothing ? T(r_p) + 3 : T(smooth_lo)
            hi = smooth_hi === nothing ? g[end]    : T(smooth_hi)
            [(lo, hi)]
        end
        pad = T(cheb_pad)
        for (lo, hi) in bands
            # END-PADDING: FIT over the band PLUS a pad margin each side (still excluding protected spike
            # cores), ASSIGN only inside [lo,hi] — the band edges sit in the fit INTERIOR (no endpoint ring).
            fidx = findall(k -> g[k] >= lo - pad && g[k] <= hi + pad && !isprotected(g[k]), 1:N)
            length(fidx) >= cheb_deg + 2 || continue
            sub_r = g[fidx]; sub_f = SgVr[fidx]
            dsub = spectral_method === :interp ? cheb_deriv(sub_f, sub_r) :
                                                 cheb_smooth_deriv(sub_f, sub_r; deg=cheb_deg)
            # EDGE-BLEND inside the WRITTEN band: cosine-taper back to the pre-spectral divr (C0 join)
            blend = T(cheb_blend)
            for (j, k) in enumerate(fidx)
                r = g[k]
                (lo <= r <= hi) || continue                 # write ONLY the requested band; padding is fit-only
                dnew = dsub[j] / g[k]^2
                wl = clamp((r - lo)/blend, zero(T), one(T))
                wr = clamp((hi - r)/blend, zero(T), one(T))
                w  = min(wl, wr)                            # 0 at the written band ends, 1 in the interior
                w2 = (1 - cos(T(π)*w))/2                    # smooth cosine ramp
                divr[k] = w2*dnew + (1 - w2)*divr[k]
            end
        end
    end
    # assemble 𝔉 per grid index (the remaining pieces are smooth per-point)
    F = Vector{Complex{T}}(undef, N)
    dpr = _fdr(T); pof(s) = csound(d, s)^2*eden0(d, s)
    dz = T(2)/T(10)^2
    @inbounds for k in 1:N
        r = g[k]; cs2 = csound(d, r)^2
        sg(z) = (r^2 + a^2*z^2)*sqrt(1 - z^2)
        v = Veq[k]
        divtφ = (-im*ω0)*v[1] + (im*m)*v[4]
        drk = Dict(l => drH_grid[l][k] for l in ls)
        SgVθ(z) = sg(z)*Vvec(d, modes, k, m, r_p, ω0, z; lmax=lmax, drH_over=drk)[3]
        gθ(z) = (-sqrt(1 - z^2)*(SgVθ(z+dz) - SgVθ(z-dz))/(2dz))/(r^2*sqrt(1 - z^2))
        divθ_ha = heightavg(z -> gθ(z), h0; hz=dz)
        H = sum(hmunu(modes[l], k, l, m, T(π)/2; a=a) for l in ls)
        gi = kerr_metric_inv(a, r, T(π)/2)
        htr = sum(gi[μ,ν]*H[μ,ν] for μ in 1:4, ν in 1:4)
        α = alpha_z(d, r, zero(T)); ω̃ = ω0 - m*disc_OmegaK(d, r)
        Rr = (pof(r+dpr) - pof(r-dpr))/(2dpr)/(eden0(d, r)*(1 + cs2))
        F[k] = cs2*(divr[k] + divtφ + divθ_ha) + Rr*v[2] + im*α*ω̃*cs2/2*htr
    end
    # FIDELITY OVERWRITE with the ORIGINAL per-point source_F: within the protect bands (always), and —
    # in targeted smooth_bands mode — at every point OUTSIDE the denoise bands, so those points are
    # byte-identical to the validated F8-gated per-point path.
    if protect
        @inbounds for k in 3:N-2
            r = g[k]
            if isprotected(r) || (targeted && !insmooth(r))
                F[k] = source_F(d, modes, k, m, r_p, ω0; lmax=lmax, h0=h0)
            end
        end
    end
    F
end
