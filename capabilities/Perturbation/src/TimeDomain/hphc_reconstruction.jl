# hphc_reconstruction.jl — DG-D4 (Move D): the TeukEvolution.jl `runtype="reconstruction"` hierarchy
# ported into Perturbation.TimeDomain, DEFECT-FIXED, layered on the landed stack:
#   - Spacetime.NPGHP  (np_spincoeff_hphc NP_0 scalars + hphc_ghp_coeffs pre_* fields + ghp_weights) — Move B
#   - Harmonics.SphereCollocation  (swal raise/lower + laplacian/killtop + angular apply)             — D2
#   - NumericsCore     (set_d1!/set_d2!/ko_filter! + rk4!)                                             — D1
#   - TimeDomain/hphc.jl  (HPHCTeukolskyProblem: grid/coefficients/Gaussian-ID; the psi4 sector)      — D3
#
# PHYSICS (RLGP-N Ripley et al. arXiv:2010.00162 §II.D/§V.C + RLGP-F Loutrel et al. arXiv:2008.11770
# §V; transport chain crosswalk in scratchpad/npghp_moveD/RLGP_equations.md §B). The 2nd-order-in-T
# Teukolsky ψ₄ (reduced to (f,p)=(ψ,∂_Tψ)) drives the ordered vacuum radiation-gauge metric-recon
# hierarchy ψ₄ → {ψ₃, λ} → {ψ₂, h_m̄m̄} → π → h_lm̄ → μh_ll, evolved by ONE method-of-lines RK4
# (arXiv:2010.00162 tex:1703-1705). The μh_ll transport couples the ±m modes (its h_mm/h_lm/π̄
# conjugates map +m↔conj(−m) under e^{imφ}, RLGP_equations.md §B.1(7)), so the state carries BOTH
# ±m and is stepped as one vector.
#
# MANDATORY DEFECT FIXES (cited to the upstream line; grounded in RLGP_equations.md P1/P2 verdicts):
#   P1  — the reconstruction ψ₄ momentum RHS CARRIES the `pre` normalization. Upstream
#         LinearEvolution.jl:50 (`set_psi4_p_k!`) omits it; the paper's `∂_T²ψ` coefficient
#         `1/pre = [8M(2M−a²R/L²)(1+2MR/L²) − a²sin²ϑ]` (RLGP-N eq:Teukolsky_equation, term-by-term
#         verified, RLGP_equations.md §A.2) MUST divide the reduced RHS for the `p ≡ ∂_Tψ` reduction.
#         `recon_set_psi4_p_k!` here = `pre .* bracket` (== Evolution.jl:169, the correct reference).
#   P2  — ONE standard RK4 (1/6,1/3,1/3,1/6) for ALL fields. Upstream Linear_evolution! hand-unrolls
#         a MALFORMED stage-4 (`third` where RK4 needs `sixth` ⇒ Σ=7/6, LinearEvolution.jl:1103-1119)
#         and a `muhll_nm`↔`muhll_pm` index typo (:1007). This port packs the whole 18-field state into
#         ONE `rk4!` (NumericsCore, classical Butcher weights) — the defect is STRUCTURALLY impossible:
#         there is no hand-rolled accumulation to corrupt.
#   P3  — the angular killtop post-step filter IS applied in reconstruction mode (on {ψ₄f,ψ₄p,λ,ψ₃},
#         matching the linear_field path Evolution.jl:299-300). Upstream binds `fltrM_pm/fltrM_nm`
#         (LinearEvolution.jl:692-693) but NEVER applies them — dead. Deliberate defect-fix deviation;
#         it is also what makes the ψ₄ sector match the D3 linear_field backend exactly (G-D4c).
#   quirk-1 — `recon_thorn_prime!` is implemented FRESH from the coefficient fields
#         (pre_thorn_prime_DT/DR of Spacetime.NPGHP.hphc_ghp_coeffs). The upstream body
#         (GHP.jl:281-301, `set_thorn_prime!`) writes an UNDEFINED `thorn_prime_arr` (dead-code typo)
#         and is NEVER transcribed here.
#   deviation (per-m operators) — each mode's ψ₄ sector uses its OWN per-m coefficients (Evo for +m and
#         −m). Upstream passes `Evo_nm` to Linear_evolution! but every set_linear_k! uses `Evo_pm`
#         (LinearEvolution.jl:719,753,…) — so the −m ψ₄ evolves with +m azimuthal coefficients (wrong
#         for s≠0: B_pp's `im·a·m` and `im·a·s·cosϑ` terms do not both flip under conjugation). The
#         port uses the correct per-m operator, required for the −m ψ₄ to match its D3 solve (G-D4c).
#   deviation (fresh cross-mode edth) — the μh_ll ±m coupling consumes ð of the PARTNER mode's h_m̄m̄
#         and h_lm̄, whose ð needs the partner field's ∂_T (its k). Upstream's two set_linear_k! calls
#         make the +m call read a STALE (previous-substep) −m k while the −m call reads a FRESH +m k
#         (an asymmetric ordering artifact). The port computes BOTH chains' k through h_lm̄ first, then
#         the cross-mode ð with FRESH k — symmetric and consistent. This affects only μh_ll (not gated
#         beyond the pointwise G-D4a fixture; ψ₄/ψ₃/λ and res_bianchi3 are independent of it).
#
# {T}-generic throughout (element type flows from the D3 problem's `T`), so Complex{BigFloat} fields
# evolve at inherited precision (G-D4e). NO coefficient closed form is re-hardcoded here — the GHP
# pre_* fields and NP_0 scalars are SOURCED from Spacetime.NPGHP (coefficients-vs-application
# principle; asserted pointwise in G-D4b). Characterized vs the ORIGINAL upstream assemblers in
# test/fixtures/td_hphc_recon_upstream.jl (G-D4a).

# ---------------------------------------------------------------------------------------------------
# (1) GHP grid applicators — thin loops over the NPGHP coefficient fields + D2 raise/lower + NP eps_0
# ---------------------------------------------------------------------------------------------------

"""
    HPHCReconGHP{T}

Per-`m` GHP operator machinery for the HPHC reconstruction backend: the eleven pointwise pre_* GHP
coefficient fields and the six background NP_0 scalars — SOURCED verbatim from
`Spacetime.NPGHP.hphc_ghp_coeffs` / `np_spincoeff_hphc` (no closed form re-hardcoded here) — plus the
`Harmonics.SphereCollocation` spin raise/lower collocation matrices for this `m`. Consumed by the
`recon_{edth,edth_prime,thorn,thorn_prime}!` applicators. Built by [`HPHCReconGHP(::HPHCTeukolskyProblem)`](@ref).
"""
struct HPHCReconGHP{T<:AbstractFloat}
    m::Int
    Rv::Vector{T}
    pre_edth_DT::Matrix{Complex{T}}
    pre_edth_raised::Matrix{Complex{T}}
    pre_edth::Matrix{Complex{T}}
    pre_edth_prime_DT::Matrix{Complex{T}}
    pre_edth_prime_lowered::Matrix{Complex{T}}
    pre_edth_prime::Matrix{Complex{T}}
    pre_thorn_DT::Matrix{Complex{T}}
    pre_thorn_DR::Matrix{Complex{T}}
    pre_thorn::Matrix{Complex{T}}
    pre_thorn_prime_DT::Matrix{Complex{T}}
    pre_thorn_prime_DR::Matrix{Complex{T}}
    mu_0::Matrix{Complex{T}}
    tau_0::Matrix{Complex{T}}
    pi_0::Matrix{Complex{T}}
    rho_0::Matrix{Complex{T}}
    eps_0::Matrix{Complex{T}}
    psi2_0::Matrix{Complex{T}}
    raise::Dict{Int,Matrix{T}}
    lower::Dict{Int,Matrix{T}}
end

"""
    HPHCReconGHP(prob::HPHCTeukolskyProblem) -> HPHCReconGHP

Assemble the per-`m` GHP machinery on the D3 problem's grid: evaluate `Spacetime.NPGHP.hphc_ghp_coeffs`
and `np_spincoeff_hphc` at each `(R, cosϑ, sinϑ)` node (the coefficient-owning kernel — Move B), and
build the `Harmonics` spin raise/lower matrices at `prob.mv` for the spins consumed by the transport
hierarchy (raise: s∈{−2,−1}; lower: s∈{−1}).
"""
function HPHCReconGHP(prob::HPHCTeukolskyProblem{T}) where {T}
    nx, ny, m = prob.nx, prob.ny, prob.mv
    a, cl = prob.bhs, prob.cl
    Rv = prob.Rv
    Cv = cos_vals(T, ny)
    Sv = sin_vals(T, ny)
    C = Matrix{Complex{T}}
    fields = (:pre_edth_DT, :pre_edth_raised, :pre_edth, :pre_edth_prime_DT, :pre_edth_prime_lowered,
              :pre_edth_prime, :pre_thorn_DT, :pre_thorn_DR, :pre_thorn, :pre_thorn_prime_DT,
              :pre_thorn_prime_DR)
    ghp = NamedTuple{fields}(ntuple(_ -> C(undef, nx, ny), length(fields)))
    npf = (:mu_0, :tau_0, :pi_0, :rho_0, :eps_0, :psi2_0)
    np = NamedTuple{npf}(ntuple(_ -> C(undef, nx, ny), length(npf)))
    @inbounds for j in 1:ny, i in 1:nx
        R, cy, sy = Rv[i], Cv[j], Sv[j]
        gc = NPGHP.hphc_ghp_coeffs(R, cy, sy; a=a, cl=cl)
        for f in fields; getfield(ghp, f)[i, j] = getfield(gc, f); end
        nc = NPGHP.np_spincoeff_hphc(R, cy, sy; a=a, cl=cl)
        for f in npf; getfield(np, f)[i, j] = getfield(nc, f); end
    end
    raise = Dict{Int,Matrix{T}}(s => swal_raising_matrix(T, ny, s, m) for s in (-2, -1))
    lower = Dict{Int,Matrix{T}}(s => swal_lowering_matrix(T, ny, s, m) for s in (-1,))
    return HPHCReconGHP{T}(m, Rv, ghp.pre_edth_DT, ghp.pre_edth_raised, ghp.pre_edth,
        ghp.pre_edth_prime_DT, ghp.pre_edth_prime_lowered, ghp.pre_edth_prime, ghp.pre_thorn_DT,
        ghp.pre_thorn_DR, ghp.pre_thorn, ghp.pre_thorn_prime_DT, ghp.pre_thorn_prime_DR,
        np.mu_0, np.tau_0, np.pi_0, np.rho_0, np.eps_0, np.psi2_0, raise, lower)
end

"""
    recon_edth!(edth, raised, f, DT, gh, spin, boost) -> nothing

GHP `ð` on a type-`(spin,boost)` field `f` with time-derivative `DT` (= ∂_T f, needed because the
chart's `m^μ` leg carries a T-component). `edth = pre_edth_DT·DT + pre_edth_raised·(ð_swal f) +
p·pre_edth·f`, `p = ghp_weights(spin,boost)[1]` (Spacetime.NPGHP). Verbatim structure of upstream
`GHP.set_edth!` (GHP.jl:207-217); coefficients from `HPHCReconGHP`, raise matrix `gh.raise[spin]`.
"""
function recon_edth!(edth, raised, f, DT, gh::HPHCReconGHP, spin::Integer, boost::Integer)
    angular_matrix_mult!(raised, f, gh.raise[spin])
    p, _ = NPGHP.ghp_weights(spin, boost)
    nx, ny = size(edth)
    @inbounds for j in 1:ny, i in 1:nx
        edth[i, j] = gh.pre_edth_DT[i, j] * DT[i, j] + gh.pre_edth_raised[i, j] * raised[i, j] +
                     p * gh.pre_edth[i, j] * f[i, j]
    end
    return nothing
end

"""
    recon_edth_prime!(edth_prime, lowered, f, DT, gh, spin, boost) -> nothing

GHP `ð′` on a type-`(spin,boost)` field: `edth_prime = pre_edth_prime_DT·DT +
pre_edth_prime_lowered·(ð′_swal f) + q·pre_edth_prime·f`, `q = ghp_weights(spin,boost)[2]`. Verbatim
structure of upstream `GHP.set_edth_prime!` (GHP.jl:234-244); lower matrix `gh.lower[spin]`.
"""
function recon_edth_prime!(edth_prime, lowered, f, DT, gh::HPHCReconGHP, spin::Integer, boost::Integer)
    angular_matrix_mult!(lowered, f, gh.lower[spin])
    _, q = NPGHP.ghp_weights(spin, boost)
    nx, ny = size(edth_prime)
    @inbounds for j in 1:ny, i in 1:nx
        edth_prime[i, j] = gh.pre_edth_prime_DT[i, j] * DT[i, j] +
                           gh.pre_edth_prime_lowered[i, j] * lowered[i, j] +
                           q * gh.pre_edth_prime[i, j] * f[i, j]
    end
    return nothing
end

"""
    recon_thorn!(thorn, f, DT, DR, gh, spin, boost, falloff) -> nothing

GHP `þ` on a type-`(spin,boost)`, peeling-`falloff` field: `thorn = pre_thorn_DT·DT +
pre_thorn_DR·(R·DR + falloff·f) + m·pre_thorn·f − R·(p·ε₀ + q·conj ε₀)·f`, `{p,q}=ghp_weights(spin,
boost)`, `ε₀` = the NP_0 spin coefficient (encodes the tetrad boost). Verbatim structure of upstream
`GHP.set_thorn!` (GHP.jl:267-272); coefficients + `ε₀` from `HPHCReconGHP` (= Spacetime.NPGHP).
"""
function recon_thorn!(thorn, f, DT, DR, gh::HPHCReconGHP, spin::Integer, boost::Integer, falloff::Integer)
    p, q = NPGHP.ghp_weights(spin, boost)
    m, R = gh.m, gh.Rv
    nx, ny = size(thorn)
    @inbounds for j in 1:ny, i in 1:nx
        thorn[i, j] = gh.pre_thorn_DT[i, j] * DT[i, j] +
                      gh.pre_thorn_DR[i, j] * (R[i] * DR[i, j] + falloff * f[i, j]) +
                      m * gh.pre_thorn[i, j] * f[i, j] -
                      R[i] * (p * gh.eps_0[i, j] + q * conj(gh.eps_0[i, j])) * f[i, j]
    end
    return nothing
end

"""
    recon_thorn_prime!(thorn_prime, f, DT, DR, gh, falloff) -> nothing

GHP `þ′` (the n-leg derivative) FRESH from the coefficient fields (quirk-1 defect fix — the upstream
`GHP.set_thorn_prime!` body writes an undefined `thorn_prime_arr` and is dead): `thorn_prime =
pre_thorn_prime_DT·DT + pre_thorn_prime_DR·(R·DR + falloff·f)`. This is the n^μ∂_μ transport whose
DT-coefficient `pre_thorn_prime_DT = 2+4MR/cl²` (= n^μ T-component) and DR-coefficient
`pre_thorn_prime_DR = R/cl²` are the isolation divisor/numerator used by [`recon_finish!`](@ref).
"""
function recon_thorn_prime!(thorn_prime, f, DT, DR, gh::HPHCReconGHP, falloff::Integer)
    R = gh.Rv
    nx, ny = size(thorn_prime)
    @inbounds for j in 1:ny, i in 1:nx
        thorn_prime[i, j] = gh.pre_thorn_prime_DT[i, j] * DT[i, j] +
                            gh.pre_thorn_prime_DR[i, j] * (R[i] * DR[i, j] + falloff * f[i, j])
    end
    return nothing
end

# ---------------------------------------------------------------------------------------------------
# (2) reconstruction RHS assemblers (verbatim from LinearEvolution.jl:32-238; {T}-generic literals)
# ---------------------------------------------------------------------------------------------------

"""
    recon_set_psi4_p_k!(kp, f_rd1, f_rd2, sph_lap, p_rd1, f, p, A_pp, A_pq, B_pp, B_pq, B_pf, pre)

ψ₄ momentum RHS with the **P1 fix**: `kp = pre .* (A_pp·∂_R p + A_pq·∂_R²f + B_pp·p + B_pq·∂_R f +
B_pf·f + Δ_S f)` (== Evolution.jl:169). Upstream's reconstruction path (LinearEvolution.jl:50) omits
the `pre` factor — the `∂_T²ψ`-coefficient inverse that the `p≡∂_Tψ` reduction requires
(RLGP_equations.md P1). The bracket matches upstream exactly; `pre .*` is the fix.
"""
function recon_set_psi4_p_k!(kp, f_rd1, f_rd2, sph_lap, p_rd1, f, p, A_pp, A_pq, B_pp, B_pq, B_pf, pre)
    nx, ny = size(kp)
    @inbounds for j in 1:ny, i in 1:nx
        kp[i, j] = pre[i, j] * (
            A_pp[i, j] * p_rd1[i, j] + A_pq[i, j] * f_rd2[i, j] +
            B_pp[i, j] * p[i, j] + B_pq[i, j] * f_rd1[i, j] + B_pf[i, j] * f[i, j] + sph_lap[i, j])
    end
    return nothing
end

"ψ₃ transport source (verbatim LinearEvolution.jl:86-88): `k = −4R·μ₀·ψ₃ − R·τ₀·ψ₄ + ðψ₄`."
function recon_set_psi3_k!(k, psi3, psi4, psi4_edth, mu_0, tau_0, R)
    nx, ny = size(k)
    @inbounds for j in 1:ny, i in 1:nx
        k[i, j] = -4 * R[i] * mu_0[i, j] * psi3[i, j] - R[i] * tau_0[i, j] * psi4[i, j] + psi4_edth[i, j]
    end
    return nothing
end

"λ transport source (verbatim LinearEvolution.jl:104): `k = −R·(μ₀+conj μ₀)·λ − ψ₄`."
function recon_set_lam_k!(k, lam, psi4, mu_0, R)
    nx, ny = size(k)
    @inbounds for j in 1:ny, i in 1:nx
        k[i, j] = -R[i] * (mu_0[i, j] + conj(mu_0[i, j])) * lam[i, j] - psi4[i, j]
    end
    return nothing
end

"ψ₂ transport source (verbatim LinearEvolution.jl:122-125): `k = −3R·μ₀·ψ₂ − 2R·τ₀·ψ₃ + ðψ₃`."
function recon_set_psi2_k!(k, psi2, psi3, psi3_edth, mu_0, tau_0, R)
    nx, ny = size(k)
    @inbounds for j in 1:ny, i in 1:nx
        k[i, j] = -3 * R[i] * mu_0[i, j] * psi2[i, j] - 2 * R[i] * tau_0[i, j] * psi3[i, j] +
                  psi3_edth[i, j]
    end
    return nothing
end

"h_m̄m̄ transport source (verbatim LinearEvolution.jl:141-142): `k = R·(μ₀−conj μ₀)·h_m̄m̄ − 2λ`."
function recon_set_hmbmb_k!(k, hmbmb, lam, mu_0, R)
    nx, ny = size(k)
    @inbounds for j in 1:ny, i in 1:nx
        k[i, j] = R[i] * (mu_0[i, j] - conj(mu_0[i, j])) * hmbmb[i, j] - 2 * lam[i, j]
    end
    return nothing
end

"""
π transport source (verbatim LinearEvolution.jl:161-168):
`k = −R·(conj π₀+τ₀)·λ + R²·½·μ₀·(conj π₀+τ₀)·h_m̄m̄ − ψ₃`.
"""
function recon_set_pi_k!(k, lam, hmbmb, psi3, mu_0, tau_0, pi_0, R)
    nx, ny = size(k)
    @inbounds for j in 1:ny, i in 1:nx
        k[i, j] = -R[i] * (conj(pi_0[i, j]) + tau_0[i, j]) * lam[i, j] +
                  (R[i]^2) * (1 // 2) * mu_0[i, j] * (conj(pi_0[i, j]) + tau_0[i, j]) * hmbmb[i, j] -
                  psi3[i, j]
    end
    return nothing
end

"h_lm̄ transport source (verbatim LinearEvolution.jl:186-189): `k = −R·conj μ₀·h_lm̄ − 2π − R·τ₀·h_m̄m̄`."
function recon_set_hlmb_k!(k, hlmb, hmbmb, pi_l, mu_0, tau_0, R)
    nx, ny = size(k)
    @inbounds for j in 1:ny, i in 1:nx
        k[i, j] = -R[i] * conj(mu_0[i, j]) * hlmb[i, j] - 2 * pi_l[i, j] -
                  R[i] * tau_0[i, j] * hmbmb[i, j]
    end
    return nothing
end

"""
μh_ll transport source — the ±m-COUPLED equation (verbatim LinearEvolution.jl:216-234). Consumes the
+m fields (`hlmb_edth, hlmb, pi_edth, pi_l, psi2`) and the −m fields (`pi_nm, hmbmb_edth_nm, hmbmb_nm,
hlmb_edth_nm, hlmb_nm`); the coupling is the `h_mm/h_lm/π̄` conjugates of RLGP-N eq:muhll_recon under
`e^{imφ}` (RLGP_equations.md §B.1(7)).
"""
function recon_set_muhll_k!(k, muhll, hlmb_edth, hlmb, pi_edth, pi_l, psi2, pi_nm, hmbmb_edth_nm,
                            hmbmb_nm, hlmb_edth_nm, hlmb_nm, mu_0, tau_0, pi_0, R)
    nx, ny = size(k)
    @inbounds for j in 1:ny, i in 1:nx
        k[i, j] = -R[i] * conj(mu_0[i, j]) * muhll[i, j] -
                  R[i] * mu_0[i, j] * hlmb_edth[i, j] -
                  (R[i]^2) * mu_0[i, j] * (conj(pi_0[i, j]) + 2 * tau_0[i, j]) * hlmb[i, j] -
                  2 * pi_edth[i, j] -
                  2 * R[i] * conj(pi_0[i, j]) * pi_l[i, j] - 2 * psi2[i, j] -
                  2 * R[i] * pi_0[i, j] * conj(pi_nm[i, j]) -
                  R[i] * pi_0[i, j] * conj(hmbmb_edth_nm[i, j]) +
                  (R[i]^2) * (pi_0[i, j]^2) * conj(hmbmb_nm[i, j]) +
                  R[i] * mu_0[i, j] * conj(hlmb_edth_nm[i, j]) +
                  (R[i]^2) * (
                      -3 * mu_0[i, j] * pi_0[i, j] + 2 * conj(mu_0[i, j]) * pi_0[i, j] -
                      2 * mu_0[i, j] * conj(tau_0[i, j])
                  ) * conj(hlmb_nm[i, j])
    end
    return nothing
end

"""
    recon_finish!(k, f, DR, gh, falloff, dr) -> nothing

The inlined þ′ radial/time isolation (verbatim LinearEvolution.jl:254-256): adds
`−pre_thorn_prime_DR·(R·∂_R f + falloff·f) / pre_thorn_prime_DT` to the already-assembled algebraic
source `k`. The divisor `pre_thorn_prime_DT = 2+4MR/cl²` and numerator coefficient
`pre_thorn_prime_DR = R/cl²` are SOURCED from `HPHCReconGHP` (= Spacetime.NPGHP.hphc_ghp_coeffs);
this is the Δ = n^μ∂_μ transport of RLGP_equations.md §B (the `∂_T(field)` isolation). Computes the
radial derivative `∂_R f` in place into `DR` via `NumericsCore.set_d1!`.
"""
function recon_finish!(k, f, DR, gh::HPHCReconGHP, falloff::Integer, dr)
    set_d1!(DR, f, dr)
    R = gh.Rv
    nx, ny = size(k)
    @inbounds for j in 1:ny, i in 1:nx
        k[i, j] += -gh.pre_thorn_prime_DR[i, j] * (R[i] * DR[i, j] + falloff * f[i, j]) /
                   gh.pre_thorn_prime_DT[i, j]
    end
    return nothing
end

"""
res_bianchi3 diagnostic (verbatim LinearEvolution.jl:1226-1229): the independent Bianchi-3 residual
`= ð′ψ₃ + 4π₀·ψ₃ − þψ₄ + ρ₀·ψ₄ − 3ψ₂₀·λ` (arXiv:2010.00162 §V.H). A quality diagnostic on the
evolved solution — small and convergent under refinement for a consistent scheme (G-D4d).
"""
function recon_set_res_bianchi3!(res, psi4, psi3, lam, psi4_thorn, psi3_edth_prime, pi_0, rho_0, psi2_0)
    nx, ny = size(res)
    @inbounds for j in 1:ny, i in 1:nx
        res[i, j] = psi3_edth_prime[i, j] + 4 * pi_0[i, j] * psi3[i, j] - psi4_thorn[i, j] +
                    rho_0[i, j] * psi4[i, j] - 3 * psi2_0[i, j] * lam[i, j]
    end
    return nothing
end

# ---------------------------------------------------------------------------------------------------
# (3) problem assembly — two D3 psi4 problems (±m) + two GHP-op sets + zero-initialized recon fields
# ---------------------------------------------------------------------------------------------------

# field-slot order in the packed (nx,ny,18) state (9 per mode; +m block then −m block)
const _RECON_FIELDS = (:psi4_f, :psi4_p, :psi3, :lam, :psi2, :hmbmb, :pi, :hlmb, :muhll)
# (spin, boost, falloff) per reconstruction field (R_teuk.md §2; Fields Initialize_Field)
const _RECON_SBF = (psi4_f=(-2, -2, 1), psi4_p=(-2, -2, 1), psi3=(-1, -1, 2), lam=(-2, -1, 1),
                    psi2=(0, 0, 3), hmbmb=(-2, 0, 1), pi=(-1, 0, 2), hlmb=(-1, 1, 2), muhll=(0, 1, 3))

"""
    HPHCReconProblem{T}

Precomputed HPHC metric-reconstruction evolution problem (Move D / D4). Wraps two D3
[`HPHCTeukolskyProblem`](@ref)s (`pm` at +m, `nm` at −m — sharing grid/dt/ID; each with its OWN per-m
ψ₄ coefficients, the deviation), two [`HPHCReconGHP`](@ref) op sets, per-spin killtop matrices for the
P3 post-step filter, and the zero-initialized reconstruction fields (only ψ₄(f,p) get the Gaussian ID;
the transport tower is sourced up the hierarchy). Build with [`HPHCReconProblem(; …)`](@ref); evolve
with [`solve_hphc_recon`](@ref).
"""
struct HPHCReconProblem{T<:AbstractFloat}
    pm::HPHCTeukolskyProblem{T}
    nm::HPHCTeukolskyProblem{T}
    gh_pm::HPHCReconGHP{T}
    gh_nm::HPHCReconGHP{T}
    kt_pm::Dict{Int,Matrix{T}}   # killtop-2 matrices per spin (P3), +m
    kt_nm::Dict{Int,Matrix{T}}
end

"""
    HPHCReconProblem(; nx, ny, m, spin=-2, bhm=1, bhs=0.7, cl=1, cfl=0.5,
                       id_l_ang=2, id_ru=4.0, id_rl=2.5, id_width=1.5,
                       id_amp=complex(0.1,0.0), ingoing=true, T=Float64)

Assemble the ±m reconstruction problem. `m>0` is required (the pair is `±m`). The ψ₄ Gaussian ID,
grid, dt and coefficients come from the D3 [`HPHCTeukolskyProblem`](@ref) constructor (one per sign of
`m`); the GHP machinery from [`HPHCReconGHP`](@ref); killtop-2 matrices for spins {−2,−1} (the
filtered set {ψ₄f,ψ₄p,λ,ψ₃}) from `Harmonics.swal_killtop_matrix`.
"""
function HPHCReconProblem(; nx::Integer, ny::Integer, m::Integer, spin::Integer=-2,
        bhm=1.0, bhs=0.7, cl=1.0, cfl=0.5,
        id_l_ang::Integer=2, id_ru=4.0, id_rl=2.5, id_width=1.5,
        id_amp=complex(0.1, 0.0), ingoing::Bool=true, T::Type{<:AbstractFloat}=Float64)
    m > 0 || throw(ArgumentError("HPHCReconProblem needs m>0 (the ±m pair); got m=$m"))
    kw = (; nx=nx, ny=ny, spin=spin, bhm=bhm, bhs=bhs, cl=cl, cfl=cfl,
          id_l_ang=id_l_ang, id_ru=id_ru, id_rl=id_rl, id_width=id_width,
          id_amp=id_amp, ingoing=ingoing, T=T)
    pm = HPHCTeukolskyProblem(; mv=m, kw...)
    nm = HPHCTeukolskyProblem(; mv=-m, kw...)
    gh_pm = HPHCReconGHP(pm)
    gh_nm = HPHCReconGHP(nm)
    kt_pm = Dict{Int,Matrix{T}}(s => swal_killtop_matrix(T, ny, s, m, 2) for s in (-2, -1))
    kt_nm = Dict{Int,Matrix{T}}(s => swal_killtop_matrix(T, ny, s, -m, 2) for s in (-2, -1))
    return HPHCReconProblem{T}(pm, nm, gh_pm, gh_nm, kt_pm, kt_nm)
end

# ---------------------------------------------------------------------------------------------------
# (4) one-mode RHS chain + packed-state ± m evolution (ONE rk4! — P2 fix) with P3 post-step filter
# ---------------------------------------------------------------------------------------------------

# Per-mode scratch (psi4 radial/angular derivs + the five stored ð outputs + one raised/DR scratch).
struct _ReconScratch{T}
    f_rd1::Matrix{Complex{T}}; f_rd2::Matrix{Complex{T}}; p_rd1::Matrix{Complex{T}}; sph_lap::Matrix{Complex{T}}
    psi4_edth::Matrix{Complex{T}}; psi3_edth::Matrix{Complex{T}}; hlmb_edth::Matrix{Complex{T}}
    pi_edth::Matrix{Complex{T}}; hmbmb_edth::Matrix{Complex{T}}
    raised::Matrix{Complex{T}}; DR::Matrix{Complex{T}}
end
_ReconScratch(::Type{T}, nx, ny) where {T} =
    _ReconScratch{T}(ntuple(_ -> zeros(Complex{T}, nx, ny), 11)...)

# Compute one mode's k for ψ₄(f,p),ψ₃,λ,ψ₂,h_m̄m̄,π,h_lm̄ (with finish) AND the five ð outputs into
# scratch — everything μh_ll needs from THIS mode. `fl` = level values (views), `kk` = k targets
# (views), both NamedTuples over _RECON_FIELDS. Mirrors the upstream set_linear_k! call order
# (LinearEvolution.jl:352-541) MINUS the μh_ll step (done cross-mode afterwards).
function _recon_chain_upto_hlmb!(kk, fl, prob1::HPHCTeukolskyProblem{T}, gh::HPHCReconGHP{T},
                                 sc::_ReconScratch{T}) where {T}
    dr = prob1.dr
    # ψ₄ sector (== D3 rhs!; P1 pre applied)
    set_d1!(sc.f_rd1, fl.psi4_f, dr)
    set_d1!(sc.p_rd1, fl.psi4_p, dr)
    set_d2!(sc.f_rd2, fl.psi4_f, dr)
    angular_matrix_mult!(sc.sph_lap, fl.psi4_f, prob1.laplM)
    recon_set_psi4_p_k!(kk.psi4_p, sc.f_rd1, sc.f_rd2, sc.sph_lap, sc.p_rd1, fl.psi4_f, fl.psi4_p,
                        prob1.A_pp, prob1.A_pq, prob1.B_pp, prob1.B_pq, prob1.B_pf, prob1.pre)
    @inbounds kk.psi4_f .= fl.psi4_p                       # set_psi4_f_k!: ∂_T f = p
    R = gh.Rv
    # ð ψ₄ (DT = ψ₄'s k = p), then ψ₃ and λ
    recon_edth!(sc.psi4_edth, sc.raised, fl.psi4_f, kk.psi4_f, gh, _RECON_SBF.psi4_f[1], _RECON_SBF.psi4_f[2])
    recon_set_psi3_k!(kk.psi3, fl.psi3, fl.psi4_f, sc.psi4_edth, gh.mu_0, gh.tau_0, R)
    recon_finish!(kk.psi3, fl.psi3, sc.DR, gh, _RECON_SBF.psi3[3], dr)
    recon_set_lam_k!(kk.lam, fl.lam, fl.psi4_f, gh.mu_0, R)
    recon_finish!(kk.lam, fl.lam, sc.DR, gh, _RECON_SBF.lam[3], dr)
    # ð ψ₃ (DT = ψ₃'s finished k), then ψ₂
    recon_edth!(sc.psi3_edth, sc.raised, fl.psi3, kk.psi3, gh, _RECON_SBF.psi3[1], _RECON_SBF.psi3[2])
    recon_set_psi2_k!(kk.psi2, fl.psi2, fl.psi3, sc.psi3_edth, gh.mu_0, gh.tau_0, R)
    recon_finish!(kk.psi2, fl.psi2, sc.DR, gh, _RECON_SBF.psi2[3], dr)
    # h_m̄m̄, π, h_lm̄
    recon_set_hmbmb_k!(kk.hmbmb, fl.hmbmb, fl.lam, gh.mu_0, R)
    recon_finish!(kk.hmbmb, fl.hmbmb, sc.DR, gh, _RECON_SBF.hmbmb[3], dr)
    recon_set_pi_k!(kk.pi, fl.lam, fl.hmbmb, fl.psi3, gh.mu_0, gh.tau_0, gh.pi_0, R)
    recon_finish!(kk.pi, fl.pi, sc.DR, gh, _RECON_SBF.pi[3], dr)
    recon_set_hlmb_k!(kk.hlmb, fl.hlmb, fl.hmbmb, fl.pi, gh.mu_0, gh.tau_0, R)
    recon_finish!(kk.hlmb, fl.hlmb, sc.DR, gh, _RECON_SBF.hlmb[3], dr)
    # the ð outputs μh_ll consumes (own + cross): ð h_lm̄, ð π, ð h_m̄m̄  (DT = finished k's)
    recon_edth!(sc.hlmb_edth, sc.raised, fl.hlmb, kk.hlmb, gh, _RECON_SBF.hlmb[1], _RECON_SBF.hlmb[2])
    recon_edth!(sc.pi_edth, sc.raised, fl.pi, kk.pi, gh, _RECON_SBF.pi[1], _RECON_SBF.pi[2])
    recon_edth!(sc.hmbmb_edth, sc.raised, fl.hmbmb, kk.hmbmb, gh, _RECON_SBF.hmbmb[1], _RECON_SBF.hmbmb[2])
    return nothing
end

# μh_ll k for one mode, using own scratch (this mode) + partner scratch/level (the ±m coupling), then finish.
function _recon_muhll!(kk, fl, sc, fl_p, sc_p, prob1::HPHCTeukolskyProblem{T}, gh::HPHCReconGHP{T}) where {T}
    recon_set_muhll_k!(kk.muhll, fl.muhll, sc.hlmb_edth, fl.hlmb, sc.pi_edth, fl.pi, fl.psi2,
                       fl_p.pi, sc_p.hmbmb_edth, fl_p.hmbmb, sc_p.hlmb_edth, fl_p.hlmb,
                       gh.mu_0, gh.tau_0, gh.pi_0, gh.Rv)
    recon_finish!(kk.muhll, fl.muhll, sc.DR, gh, _RECON_SBF.muhll[3], prob1.dr)
    return nothing
end

# NamedTuple of @views onto the 9 field slots for one mode (offset 0 for +m, 9 for −m) of the (nx,ny,18) state.
@inline function _mode_views(A3, off)
    NamedTuple{_RECON_FIELDS}(ntuple(k -> (@view A3[:, :, off+k]), 9))
end

"""
    solve_hphc_recon(prob; nt, ts) -> (; times, pm, nm)

Evolve the HPHC reconstruction problem for `nt` RK4 steps (step `prob.pm.dt`), recording all nine
fields for both modes every `ts` steps. The 18-field state `(ψ₄f,ψ₄p,ψ₃,λ,ψ₂,h_m̄m̄,π,h_lm̄,μh_ll)×(±m)`
is packed into one `Vector{Complex{T}}` and driven by ONE `NumericsCore.rk4!` (the P2 fix — no
hand-rolled per-field weights). The RHS computes both mode chains through h_lm̄ (fresh cross-mode ð),
then both μh_ll (the ±m coupling). The P3 post-step filter (KO ε=0.5 + killtop-2 angular) is applied
to {ψ₄f,ψ₄p,λ,ψ₃} after each step. `pm`/`nm` are NamedTuples of snapshot vectors per field.
"""
function solve_hphc_recon(prob::HPHCReconProblem{T}; nt::Integer, ts::Integer) where {T}
    nx, ny = prob.pm.nx, prob.pm.ny
    dt = prob.pm.dt
    u = zeros(Complex{T}, nx * ny * 18)
    U = reshape(u, nx, ny, 18)
    # ID: ψ₄(f,p) from the D3 Gaussian; recon tower zero
    @views U[:, :, 1] .= prob.pm.f0
    @views U[:, :, 2] .= prob.pm.p0
    @views U[:, :, 10] .= prob.nm.f0
    @views U[:, :, 11] .= prob.nm.p0
    sc_pm = _ReconScratch(T, nx, ny)
    sc_nm = _ReconScratch(T, nx, ny)
    epsKO = T(1) / 2
    tmpf = zeros(Complex{T}, nx, ny)

    function rhs!(du, uu, ::Real)
        UU = reshape(uu, nx, ny, 18); dUU = reshape(du, nx, ny, 18)
        fl_pm = _mode_views(UU, 0); kk_pm = _mode_views(dUU, 0)
        fl_nm = _mode_views(UU, 9); kk_nm = _mode_views(dUU, 9)
        _recon_chain_upto_hlmb!(kk_pm, fl_pm, prob.pm, prob.gh_pm, sc_pm)
        _recon_chain_upto_hlmb!(kk_nm, fl_nm, prob.nm, prob.gh_nm, sc_nm)
        _recon_muhll!(kk_pm, fl_pm, sc_pm, fl_nm, sc_nm, prob.pm, prob.gh_pm)
        _recon_muhll!(kk_nm, fl_nm, sc_nm, fl_pm, sc_pm, prob.nm, prob.gh_nm)
        return du
    end

    # P3 post-step: KO(ε=0.5) then killtop-2 angular, on {ψ₄f,ψ₄p,λ,ψ₃} of both modes.
    function postfilter!(uu)
        UU = reshape(uu, nx, ny, 18)
        for (off, kt) in ((0, prob.kt_pm), (9, prob.kt_nm))
            for (slot, s) in ((1, -2), (2, -2), (4, -2), (3, -1))   # ψ₄f,ψ₄p,λ : s=−2 ; ψ₃ : s=−1
                v = @view UU[:, :, off + slot]
                ko_filter!(v, tmpf, epsKO)
                tmpf .= v; angular_matrix_mult!(v, tmpf, kt[s])
            end
        end
        return nothing
    end

    snaps = Dict(f => Matrix{Complex{T}}[] for f in _RECON_FIELDS)
    snaps_nm = Dict(f => Matrix{Complex{T}}[] for f in _RECON_FIELDS)
    times = T[]
    function cb(uu, t, step)
        step >= 1 && postfilter!(uu)
        if step % ts == 0
            UU = reshape(uu, nx, ny, 18)
            for (k, f) in enumerate(_RECON_FIELDS)
                push!(snaps[f], copy(@view UU[:, :, k]))
                push!(snaps_nm[f], copy(@view UU[:, :, 9 + k]))
            end
            push!(times, T(t))
        end
        return nothing
    end

    rk4!(rhs!, u, (zero(T), nt * dt), dt; callback=cb)
    pm = NamedTuple{_RECON_FIELDS}(Tuple(snaps[f] for f in _RECON_FIELDS))
    nm = NamedTuple{_RECON_FIELDS}(Tuple(snaps_nm[f] for f in _RECON_FIELDS))
    return (; times, pm, nm)
end

"""
    hphc_recon_residual3(prob, snap_pm) -> Matrix

Independent Bianchi-3 residual (arXiv:2010.00162 §V.H) evaluated on a +m snapshot `snap_pm`
(a NamedTuple of `nx×ny` field matrices, e.g. `(; ψ₄f, ψ₄p, ψ₃, λ, …)` — one recorded step of
`solve_hphc_recon(...).pm`). Recomputes the ψ₃ RHS at this state (for `ð′ψ₃`'s DT) and `þψ₄`
(DT = ψ₄'s k = ψ₄p), then `recon_set_res_bianchi3!`. Returns the residual field; its `maximum(abs, …)`
is the G-D4d quality diagnostic.
"""
function hphc_recon_residual3(prob::HPHCReconProblem{T}, snap_pm) where {T}
    nx, ny = prob.pm.nx, prob.pm.ny
    gh = prob.gh_pm; dr = prob.pm.dr
    psi4 = snap_pm.psi4_f; psi4_p = snap_pm.psi4_p
    psi3 = snap_pm.psi3; lam = snap_pm.lam
    # ψ₃'s k at this state (algebraic ðψ₄ source + finish) — needed as DT for ð′ψ₃
    f_rd1 = zeros(Complex{T}, nx, ny); f_rd2 = zeros(Complex{T}, nx, ny)
    p_rd1 = zeros(Complex{T}, nx, ny); sph_lap = zeros(Complex{T}, nx, ny)
    raised = zeros(Complex{T}, nx, ny); lowered = zeros(Complex{T}, nx, ny)
    DR = zeros(Complex{T}, nx, ny); psi4_edth = zeros(Complex{T}, nx, ny)
    psi4_k = copy(psi4_p)                                  # ∂_T ψ₄ = ψ₄p
    recon_edth!(psi4_edth, raised, psi4, psi4_k, gh, -2, -2)
    psi3_k = zeros(Complex{T}, nx, ny)
    recon_set_psi3_k!(psi3_k, psi3, psi4, psi4_edth, gh.mu_0, gh.tau_0, gh.Rv)
    recon_finish!(psi3_k, psi3, DR, gh, _RECON_SBF.psi3[3], dr)
    # ð′ψ₃ (DT = ψ₃'s k) and þψ₄ (DT = ψ₄p, DR = ∂_R ψ₄)
    psi3_edth_prime = zeros(Complex{T}, nx, ny)
    recon_edth_prime!(psi3_edth_prime, lowered, psi3, psi3_k, gh, _RECON_SBF.psi3[1], _RECON_SBF.psi3[2])
    psi4_DR = zeros(Complex{T}, nx, ny); set_d1!(psi4_DR, psi4, dr)
    psi4_thorn = zeros(Complex{T}, nx, ny)
    recon_thorn!(psi4_thorn, psi4, psi4_k, psi4_DR, gh, _RECON_SBF.psi4_f[1], _RECON_SBF.psi4_f[2],
                 _RECON_SBF.psi4_f[3])
    res = zeros(Complex{T}, nx, ny)
    recon_set_res_bianchi3!(res, psi4, psi3, lam, psi4_thorn, psi3_edth_prime, gh.pi_0, gh.rho_0, gh.psi2_0)
    return res
end
