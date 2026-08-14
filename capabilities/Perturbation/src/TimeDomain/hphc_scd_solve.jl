# hphc_scd_solve.jl — F7-2: the SOURCED second-order Teukolsky solve + TD QQNM-ratio extraction.
#
# PHYSICS / PROVENANCE. This closes the F7 chain: F7-1 (hphc_scd_source.jl) assembles the second-order
# Teukolsky source 𝒮 = d₄⁽⁰⁾𝔰_d + d₃⁽⁰⁾𝔰_t from the F6-A (2,2) QNM tower; here we DRIVE the D3 (f,p)
# reduction of the 2nd-order Teukolsky equation with that source and extract the driven (4,4) QQNM
# amplitude, forming the time-domain quadratic ratio R_TD = A_drive/(A_parent)² for 220⁺×220⁺→(4,4)
# at a=0 — the independent TD-vs-FD cross-validation of the owner's frozen channel-fit anchor.
#
# ── THE SOURCED EVOLUTION (the falloff/rescaling derivation — [PRELIMINARY]→[SOLID] via regularity) ──
# The D3 backend evolves the CODE variable ψ (= scd_f, peeling-stripped Ψ₄; TeukEvolution.jl:208-216
# declares scd_f/scd_p as spin/boost/falloff = −2,−2,1) by the (f,p)=(ψ,∂_Tψ) reduction
#     ∂_T p = pre·(A_pp ∂_R p + A_pq ∂_R²f + B_pp p + B_pq ∂_R f + B_pf f + Δ_S f)          (= kp, hphc.jl)
# which is arXiv:2010.00162 eq:Teukolsky_equation (numerics:1408) solved for ∂_T²ψ, `1/pre` = the
# ∂_T² coefficient bracket [8M(2M−a²R/L²)(1+2MR/L²) − a²sin²ϑ]. That equation is the NP Teukolsky
# operator 𝒯 (spin −2) MULTIPLIED by the peeling factor 2Σ_BL/R (numerics:1399, Σ_BL≡r²+a²cos²ϑ) so
# the leading terms are finite at R=0 (scri). The SECOND-order equation (numerics:1852-1865) is the
# SAME operator with the same 2Σ_BL/R prefactor on BOTH sides:
#     (2Σ_BL/R)·𝒯(R ψ⁽²⁾) = (2Σ_BL/R)·𝒮 ,   coordinate LHS == eq:Teukolsky_equation with ψ⁽¹⁾→ψ⁽²⁾.
# Solving the LHS for ∂_T²ψ⁽²⁾ gives the DRIVEN reduction
#     ∂_T p⁽²⁾ = kp[ψ⁽²⁾] + pre·(2Σ_BL/R)·𝒮_physical .                                       ← the source enters ×pre·(2Σ_BL/R)
# The physical source peels as 𝒮_physical ~ R³ at scri: since (2Σ_BL/R)𝒯(Rψ) = O(1) (eq:Teukolsky_eq
# has an O(1) leading coefficient) and 2Σ_BL/R = 2r²/R ~ 2L⁴/R³ (a=0), we get 𝒯Ψ₄ = 𝒯(Rψ) ~ R³·O(1).
# The F7-1 assembler `scd_source!` returns the PEELING-STRIPPED code source 𝒮̂ = 𝒮_physical/R³ — FINITE
# and nonzero at scri (measured: |𝒮̂(R=0)| ≈ O(1); the recon-form GHP applicators produce code-
# normalized outputs, the F7-1 falloff table's "𝒮 fall 1 == scd_f"). Hence the injected weight is
#     W(R,y) = pre·(2Σ_BL/R)·R³ = pre·2Σ_BL·R² = pre·2(cl⁴ + a²cos²ϑ·R²)      (a=0: W = 2cl⁴·pre)   [scd_source_weight]
# and the driven momentum RHS is  ∂_T p⁽²⁾ = kp[ψ⁽²⁾] + W .* 𝒮̂ .* e^{−iΩt}.
#
# THE POWER q=3 IS DECISIVE (the falloff discriminator, G-F7-2a). W ∝ R^{q−3}: only q=3 makes W finite
# nonzero at scri, so the driven ψ⁽²⁾ inherits the falloff-1 peeling of the free field (Ψ₄⁽²⁾=R·ψ⁽²⁾,
# ψ⁽²⁾ finite at scri). q<3 ⇒ W ∝ R^{<0} DIVERGES at scri (the driven field is not falloff-1 — breaks
# scri regularity/convergence); q>3 ⇒ W vanishes at scri (driven ψ⁽²⁾ falloff>1 — wrong peeling). The
# regularity of the driven scri amplitude under nx-refinement is the pasted evidence.
#
# THE DRIVE IS ANALYTIC [AXIOM-F7-2]. Every source factor ∝ e^{−iωt} (the D4 fact ∂_T X = −iωX), so
# the bilinear 𝒮 ∝ e^{−2iωt} EXACTLY (G-F7-1d, 7.5e-14). We drive with the CLOSED FORM
# 𝒮(t) = 𝒮(0)·e^{−iΩt}, Ω = 2ω₂₂₀ — NO co-evolution of the tower. A cosine RAMP (optional, `tramp`)
# turns the source on smoothly to suppress the sudden-turn-on transient; after `tramp` the drive is the
# pure closed form and the steady driven amplitude is ramp-independent.
#
# EXTRACTION. At a scri-adjacent probe (i=1) the (4,4)=(l=4,m=4) projection of ψ⁽²⁾ is a two-component
# signal: the DRIVEN line at Ω=2ω₂₂₀ (Im −0.1779, decays FAST) + FREE (4,4) QNM transients (ω₄₄₀,
# Im ≈ −0.094, decays SLOWER). We fit A_drive e^{−iΩt} + Σ A_free,k e^{−iω_f,k t} over an honest window
# (Ω FIXED = 2ω₂₂₀; the free ω_f,k are fit — they must land near a plausible (4,4) fundamental). Then
#     R_TD = A_drive / (A_parent)² ,   A_parent = the seeded (2,2) ψ₄ scri amplitude (id.A_pm, same var).
# a = 0 THIS ROUND ([AXIOM-F6A-a0]); {T}-generic (Complex{BigFloat} flows, G-F7-2e).

# ---------------------------------------------------------------------------------------------------
# (1) the injection weight  W = pre·2Σ_BL·R²   (the falloff-derived source rescaling; q=3)
# ---------------------------------------------------------------------------------------------------

"""
    scd_source_weight(prob4::HPHCTeukolskyProblem) -> Matrix{T}

The rescaling `W(R,y)` that injects the F7-1 code source `𝒮̂` (= `scd_source!` output) into the (f,p)
momentum RHS of the second-order solve: `∂_T p⁽²⁾ = kp[ψ⁽²⁾] + W .* 𝒮̂ .* e^{−iΩt}`. Derived from the
RLGP second-order equation `(2Σ_BL/R)𝒯(Rψ⁽²⁾) = (2Σ_BL/R)𝒮` (arXiv:2010.00162 numerics:1852-1865):

    W = pre·(2Σ_BL/R)·R³ = pre·2Σ_BL·R² = 2·pre·(cl⁴ + a²cos²ϑ·R²)      (a=0 ⇒ W = 2·cl⁴·pre),

where `pre` is the D3 ∂_T²-coefficient inverse and `Σ_BL = r²+a²cos²ϑ`, `r = cl²/R` (so `Σ_BL·R² =
cl⁴ + a²cos²ϑ·R²`, finite at scri — no `1/R`). The `R³` power (q=3) is the peeling of `𝒮_physical =
R³·𝒮̂` and is what makes `W` finite nonzero at scri, giving the driven `ψ⁽²⁾` the free field's falloff-1
peeling (G-F7-2a discriminator: q<3 ⇒ `W` diverges at scri; q>3 ⇒ `W` vanishes at scri).
"""
function scd_source_weight(prob4::HPHCTeukolskyProblem{T}) where {T}
    nx, ny = prob4.nx, prob4.ny
    Rv, cl, a = prob4.Rv, prob4.cl, prob4.bhs
    Cv = cos_vals(T, ny)
    W = Matrix{T}(undef, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        SigR2 = cl^4 + (a * Cv[j] * Rv[i])^2          # Σ_BL·R² = cl⁴ + a²cos²ϑ·R²  (finite at R=0)
        W[i, j] = 2 * prob4.pre[i, j] * SigR2
    end
    return W
end

# ---------------------------------------------------------------------------------------------------
# (2) the sourced evolution — the D3 (f,p) RK4 with the analytic drive added to the p-equation
# ---------------------------------------------------------------------------------------------------

"""
    scd_ramp(t, tramp) -> T

Smooth half-cosine turn-on `∈[0,1]`: `0.5(1−cos(π t/tramp))` for `t<tramp`, `1` after. `tramp≤0`
disables the ramp (pure closed-form drive from t=0).
"""
@inline function scd_ramp(t::T, tramp::T) where {T}
    tramp <= zero(T) && return one(T)
    t >= tramp ? one(T) : (one(T) - cos(T(pi) * t / tramp)) / 2
end

"""
    solve_scd_driven(prob4, Ssrc, Ω; nt, ts, tramp=0, weight=scd_source_weight(prob4),
                     l_ext=prob4.mv, m_ext=prob4.mv) -> (; times, signal, f_scri, weight)

Evolve the SOURCED second-order Teukolsky field `ψ⁽²⁾` on the D3 problem `prob4` (spin −2, azimuthal
`prob4.mv`) from REST (`f=p=0`), driven by the analytic closed-form source `W .* Ssrc .* e^{−iΩt}`
added to the momentum RHS ([AXIOM-F7-2]). `Ssrc` is the F7-1 `scd_source!` output in this azimuthal
sector; `Ω = 2ω₂₂₀`. Steps `nt` RK4 steps of `prob4.dt`, recording every `ts` steps the scri (i=1)
field row and its `(l_ext,m_ext)` spin-weighted projection `signal`. Same per-step KO+killtop post-
filter as the homogeneous `solve_hphc`. Returns the snapshot `times`, the complex `signal` time series
(the driven (l,m) QQNM waveform at scri), the recorded scri rows `f_scri`, and the `weight` used.
"""
function solve_scd_driven(prob4::HPHCTeukolskyProblem{T}, Ssrc::AbstractMatrix{Complex{T}}, Ω;
                          nt::Integer, ts::Integer, tramp=zero(T),
                          weight::AbstractMatrix{T}=scd_source_weight(prob4),
                          l_ext::Integer=prob4.mv, m_ext::Integer=prob4.mv) where {T}
    nx, ny = prob4.nx, prob4.ny
    dr, dt = prob4.dr, prob4.dt
    Ωc = Complex{T}(Ω); tr = T(tramp)
    drive0 = weight .* Ssrc                                   # W .* 𝒮̂  (the t=0 drive field)
    A_pp, A_pq = prob4.A_pp, prob4.A_pq
    B_pp, B_pq, B_pf = prob4.B_pp, prob4.B_pq, prob4.B_pf
    pre, laplM, fltrM = prob4.pre, prob4.laplM, prob4.fltrM
    # extraction: the swal(-2,m_ext,l_ext,·) collocation vector + GL weights (Σ w_j f_j swal)
    swl = Complex{T}[swal(prob4.spin, m_ext, l_ext, y) for y in prob4.Yv]

    u = zeros(Complex{T}, 2 * nx * ny)
    f_rd1 = zeros(Complex{T}, nx, ny); f_rd2 = zeros(Complex{T}, nx, ny)
    p_rd1 = zeros(Complex{T}, nx, ny); f_lap = zeros(Complex{T}, nx, ny)
    tmpf = zeros(Complex{T}, nx, ny); tmpp = zeros(Complex{T}, nx, ny)

    function rhs!(du, uu, t::Real)
        UU = reshape(uu, nx, ny, 2); dUU = reshape(du, nx, ny, 2)
        ff = @view UU[:, :, 1]; pp = @view UU[:, :, 2]
        df = @view dUU[:, :, 1]; dp = @view dUU[:, :, 2]
        set_d1!(f_rd1, ff, dr); set_d1!(p_rd1, pp, dr); set_d2!(f_rd2, ff, dr)
        angular_matrix_mult!(f_lap, ff, laplM)
        ph = exp(-im * Ωc * T(t)) * scd_ramp(T(t), tr)
        @inbounds for j in 1:ny, i in 1:nx
            df[i, j] = pp[i, j]
            dp[i, j] = pre[i, j] * (
                A_pp[i, j] * p_rd1[i, j] + A_pq[i, j] * f_rd2[i, j] +
                B_pp[i, j] * pp[i, j] + B_pq[i, j] * f_rd1[i, j] +
                B_pf[i, j] * ff[i, j] + f_lap[i, j]) + drive0[i, j] * ph
        end
        return du
    end
    epsKO = T(1) / 2
    function postfilter!(uu)
        UU = reshape(uu, nx, ny, 2); ff = @view UU[:, :, 1]; pp = @view UU[:, :, 2]
        ko_filter!(ff, tmpf, epsKO); ko_filter!(pp, tmpp, epsKO)
        tmpf .= ff; angular_matrix_mult!(ff, tmpf, fltrM)
        tmpp .= pp; angular_matrix_mult!(pp, tmpp, fltrM)
        return nothing
    end

    times = T[]; signal = Complex{T}[]; f_scri = Vector{Complex{T}}[]
    function cb(uu, t, step)
        step >= 1 && postfilter!(uu)
        if step % ts == 0
            UU = reshape(uu, nx, ny, 2)
            row = copy(@view UU[1, :, 1])
            push!(times, T(t)); push!(f_scri, row)
            # (l_ext,m_ext) projection at scri: Σ_j w_j row_j swal_j (GL orthonormality; swal real).
            push!(signal, inner_product(row, swl))
        end
        return nothing
    end
    rk4!(rhs!, u, (zero(T), nt * dt), dt; callback=cb)
    return (; times, signal, f_scri, weight)
end

"Scri-adjacent `(l,m)` spin-weighted projection of a field row (the swal-orthonormal GL analysis)."
scd_project_lm(row::AbstractVector{Complex{T}}, spin::Integer, m::Integer, l::Integer,
               ny::Integer=length(row)) where {T} =
    inner_product(row, Complex{T}[swal(spin, m, l, y) for y in Y_vals(T, ny)])

# ---------------------------------------------------------------------------------------------------
# (3) the multi-mode ring fit — A_drive at the FIXED drive freq Ω + free QNM modes (freqs fitted)
# ---------------------------------------------------------------------------------------------------

# linear least squares for the mode amplitudes at a FIXED frequency list `ws` (ws[1] = drive Ω),
# restricted to the window [tlo,thi]. Returns (amps, relative-residual, npts).
function _scd_lsq(times::AbstractVector{T}, sig::AbstractVector{Complex{T}},
                  ws::AbstractVector{Complex{T}}; tlo, thi) where {T}
    idx = [k for k in eachindex(times) if tlo <= times[k] <= thi]
    n = length(idx); m = length(ws)
    A = Matrix{Complex{T}}(undef, n, m); b = Vector{Complex{T}}(undef, n)
    @inbounds for (r, k) in enumerate(idx)
        b[r] = sig[k]
        for p in 1:m
            A[r, p] = exp(-im * ws[p] * times[k])
        end
    end
    x = A \ b
    res = norm(A * x - b) / max(norm(b), eps(T))
    return (amps=x, res=res, n=n)
end

"""
    scd_fit_qqnm(times, signal, Ω; window, Kfree=1, wfree0=[Complex(0.809,-0.094)]) -> (; Adrive, wfree, amps, res, n)

Two-(or multi-)mode ring fit of the driven (4,4) signal: `signal ≈ A_drive e^{−iΩt} + Σ_{k=1}^{Kfree}
A_free,k e^{−iω_f,k t}` over the time `window = (tlo,thi)`. The drive frequency `Ω = 2ω₂₂₀` is FIXED
(the source's analytic line, [AXIOM-F7-2]); the `Kfree` free complex frequencies `ω_f,k` are fitted by
Nelder–Mead on the LSQ residual (the amplitudes are the linear least-squares solution at each frequency
guess). Returns the drive amplitude `Adrive`, the fitted free frequencies `wfree` (the first should
land near the (4,4) fundamental ≈ 0.8092−0.0942i), the full amplitude vector, the relative residual,
and the window point count.
"""
function scd_fit_qqnm(times::AbstractVector{T}, signal::AbstractVector{Complex{T}}, Ω;
                      window, Kfree::Integer=1,
                      wfree0::AbstractVector=[Complex{T}(0.80918, -0.09416)]) where {T}
    tlo, thi = window
    Ωc = Complex{T}(Ω)
    length(wfree0) == Kfree || throw(ArgumentError("wfree0 must have length Kfree=$(Kfree)"))
    # objective: residual as a function of the 2·Kfree real free-frequency params
    obj(v) = _scd_lsq(times, signal,
                      vcat(Ωc, [Complex{T}(v[2p-1], v[2p]) for p in 1:Kfree]); tlo=tlo, thi=thi).res
    x0 = T[]; for w in wfree0; push!(x0, T(real(w))); push!(x0, T(imag(w))); end
    xbest = _nelder_mead(obj, x0)
    ws = vcat(Ωc, [Complex{T}(xbest[2p-1], xbest[2p]) for p in 1:Kfree])
    r = _scd_lsq(times, signal, ws; tlo=tlo, thi=thi)
    return (; Adrive=r.amps[1], wfree=ws[2:end], amps=r.amps, res=r.res, n=r.n)
end

# minimal Nelder–Mead (no external deps; the residual surface is smooth near the QNM). {T}-generic.
function _nelder_mead(f, x0::Vector{T}; iters::Integer=500, tol=eps(T)^(3//4)) where {T}
    np = length(x0)
    simplex = [copy(x0)]
    for i in 1:np
        y = copy(x0); y[i] += (y[i] == 0 ? T(3)/100 : abs(y[i]) * T(3)/100 + T(1)/1000); push!(simplex, y)
    end
    fv = [f(s) for s in simplex]
    for _ in 1:iters
        p = sortperm(fv); simplex = simplex[p]; fv = fv[p]
        maximum(abs, fv .- fv[1]) <= tol * (abs(fv[1]) + tol) && break
        c = sum(simplex[1:end-1]) ./ np
        xr = c .+ (c .- simplex[end]); fr = f(xr)
        if fr < fv[1]
            xe = c .+ 2 .* (c .- simplex[end]); fe = f(xe)
            simplex[end], fv[end] = fe < fr ? (xe, fe) : (xr, fr)
        elseif fr < fv[end-1]
            simplex[end], fv[end] = xr, fr
        else
            xk = c .+ (simplex[end] .- c) ./ 2; fk = f(xk)
            if fk < fv[end]
                simplex[end], fv[end] = xk, fk
            else
                for i in 2:np+1
                    simplex[i] = (simplex[1] .+ simplex[i]) ./ 2; fv[i] = f(simplex[i])
                end
            end
        end
    end
    return simplex[argmin(fv)]
end

# ---------------------------------------------------------------------------------------------------
# (4) the full pipeline — assemble source, drive, extract, form R_TD
# ---------------------------------------------------------------------------------------------------

"""
    scd_qqnm_ratio(prob::HPHCReconProblem, qnm; m_pos=2, sector=:pp, l_ext=4,
                   amp=complex(0.1,0.0), tmax=50.0, ts=nothing, tramp=10.0,
                   window=(12.0,48.0), Kfree=2, wfree0=…) -> NamedTuple

The full F7-2 pipeline for the `220⁺×220⁺→(4,4)` (or general `m_pos`, `sector`) time-domain QQNM
ratio at a=0. (1) builds the F6-A constraint-consistent QNM tower `id` from `qnm`; (2) assembles the
second-order source `𝒮̂` in the named `sector` (`:pp` ⇒ the m=2·m_pos sector) via `scd_source_sector`;
(3) drives a fresh spin−2, m=2·m_pos D3 problem from rest with the analytic source (`solve_scd_driven`);
(4) fits the driven (l_ext, 2·m_pos) waveform (`scd_fit_qqnm`); (5) returns
`R_TD = A_drive/(A_parent)²`, `A_parent = id.A_pm` (the seeded (2,2) ψ₄ scri amplitude, SAME variable/
probe convention as A_drive). Returns everything: `R_TD, Adrive, A_parent, wfree, res, omega, times,
signal`. See the file header for the [CONVENTION-DELTA] to the owner's Ma–Yang Γ̄⁴-contour `R_++`.
"""
function scd_qqnm_ratio(prob::HPHCReconProblem{T}, qnm; m_pos::Integer=2, sector::Symbol=:pp,
                        l_ext::Integer=4, amp=complex(0.1, 0.0), tmax=50.0, ts=nothing,
                        tramp=10.0, window=(12.0, 48.0), Kfree::Integer=2,
                        wfree0::AbstractVector=[Complex{T}(0.80918, -0.09416),
                                               Complex{T}(0.40, -0.28)],
                        lam_from::Symbol=:transport) where {T}
    id = hphc_qnm_initialdata(prob, qnm; amp=amp, lam_from=lam_from)
    ops = SCDSourceOps(prob)
    r = scd_source_sector(id, id, m_pos, sector, ops)
    m_ext = r.m_ang                                          # the target azimuthal (= 2·m_pos for :pp)
    prob4 = HPHCTeukolskyProblem(; nx=prob.pm.nx, ny=prob.pm.ny, mv=m_ext, spin=prob.pm.spin,
                                 bhm=prob.pm.bhm, bhs=prob.pm.bhs, cl=prob.pm.cl, T=T)
    dt = prob4.dt
    nt = ceil(Int, T(tmax) / dt)
    tsv = ts === nothing ? max(1, round(Int, T(1)/4 / dt)) : Int(ts)   # ~4 samples per M
    sol = solve_scd_driven(prob4, r.S, r.omega; nt=nt, ts=tsv, tramp=T(tramp),
                           l_ext=l_ext, m_ext=m_ext)
    fit = scd_fit_qqnm(sol.times, sol.signal, r.omega; window=window, Kfree=Kfree,
                       wfree0=Complex{T}.(wfree0))
    A_parent = id.A_pm
    R_TD = fit.Adrive / A_parent^2
    return (; R_TD, Adrive=fit.Adrive, A_parent, wfree=fit.wfree, res=fit.res, n=fit.n,
            omega=r.omega, omega_parent=id.omega, times=sol.times, signal=sol.signal, prob4)
end
