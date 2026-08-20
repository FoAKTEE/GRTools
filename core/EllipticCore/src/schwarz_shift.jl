# schwarz_shift.jl — the complex-shift POLICY of the Schwarz smoother and the
# Schwarz-preconditioned-GMRES glue (CRP main.tex:680, advancement (2)): Helmholtz-type
# m-mode problems in t-slicing are indefinite at the scale (mΩ)², so noise excites
# spurious oscillations that render the subdomain solves ineffective; a small complex
# shift on the order of m²Ω² added to the subdomain matrices dampens them
# (Erlangga2006, vanGijzen2007 — complex-shifted-Laplacian preconditioning). The
# literal-diagonal shift MECHANISM (−i·s·I on the subdomain matrices, cache rekeying on
# s — the spectre ExplicitInverse ComplexShift convention) is owned by schwarz.jl
# (W6-B, gated bitwise there); this file owns the shift VALUE from problem data, its
# slicing policy, the mass-weighted APPLICATION route, and the right-preconditioner
# adapter the row-18 preconditioned-GMRES oracles run through.
#
# MASS-WEIGHTED APPLICATION (BINDING for this operator family; established by the
# pre-red calibration falsifiers, evidence/schwarz-complex-shift/):
# EllipticCore's operator is the MASSIVE Gauss formulation (dg_operator.jl, matching
# the CRP production config Massive: True + Quadrature: Gauss, anc yaml:167-171), so
# the OPERATOR-level damping perturbation β → β − i·s discretizes to −i·s·M (M = the
# diagonal Gauss mass), NOT −i·s·I. The policy value s = (mΩ)² is an operator-level
# (β-units) quantity: apply it by building the SMOOTHER's operator with the shifted
# β − i·s coefficient — its subdomain matrices are then A_ss − i·s·M_ss (gated) — while
# GMRES runs on the ORIGINAL operator. Round-1 calibration FALSIFIED the literal
# route at this value: −i·s·I with s = (mΩ)² ~ 1e2 dwarfs the massive matrices (mass
# diagonal ~1e-3 on the unit-square models), degrades every patch solve to ≈ (i/s)·r,
# and Schwarz-pGMRES fails on ALL shifted runs while s = 0 converges — the literal
# schwarz.jl seam remains correct spectre-convention MECHANISM, but its useful values
# live on the matrix scale, not the operator scale.
#
# WHAT THE SHIFT BUYS (round-3 calibration, the honest t-sliced pathology): away from
# subdomain resonances the unshifted smoother is fine; when k² = (mΩ)²·profile hits a
# patch eigenvalue (inevitable SOMEWHERE as m sweeps 0..20), the unshifted patch
# matrix is near-singular (σ_min ~ 1e-13 measured ON resonance) and its solves amplify
# the resonant mode — "noise ... excites oscillations, rendering preconditioners
# ineffective" (main.tex:680): pGMRES stagnates ~8 orders above tolerance. The
# mass-weighted shift bounds every patch inverse by ~1/(s·mass) and restores UNIFORM
# convergence, insensitive to the distance from resonance (its ≈ const across the
# scanned δ ∈ [−0.05, 0.05] window) at ~2x the generic-k cost.
#
# Sign convention (documented against the operator's time convention): the locked
# m-mode ansatz phase is e^{+i m (φ̄ − Ω t)} (ScalarSelfForce mode-ansatz docstring;
# CRP eq:mode_sum), i.e. time dependence e^{−iωt} with ω = mΩ. A friction term +2ε∂_t
# of the time-domain wave operator lands in the frequency-domain operator as −2iεω:
# for ω > 0 physical absorption enters on the −i side — the same side as the spectre
# shift A_ss − i·s·I (ExplicitInverse.hpp:249-254), the W6-B damped-Helmholtz test
# convention β = −k² − iε, and this policy's β − i·s. For the REAL indefinite
# principal part the two sides are spectrally equivalent (|λ ∓ i·s| ≥ s·mass either
# way), so the policy returns s = (mΩ)² ≥ 0, even in mΩ (retrograde orbits included).
#
# Slicing policy (CRP main.tex:680): the shift is ESSENTIAL for oscillating problems
# in t-slicing, and "can be reduced or disabled altogether when hyperboloidal or null
# slicing is used, as the solution loses its global wave-like behavior". The shipped
# production configuration — a vtu-slicing run — sets ComplexShift: 0.0 at both stack
# levels (anc/ScalarSelfForce.yaml:194,203), so the policy returns 0 for
# :hyperboloidal and :vtu (the conservative end of "reduced or disabled"; callers may
# still shift explicitly).
#
# DAG node: schwarz-complex-shift (elliptic-ssf GLOBAL_DAG row 18, TOPOLOGICAL_SCHEDULE
# W7-A). Predecessor scope cited: schwarz-smoother owns the −i·s·I seam + cache key;
# complex-gmres owns the right-preconditioning hook. The FULL preconditioner stack
# (multigrid V-cycle + bottom solve + Schwarz) is the solver-stack-core node's scope —
# the adapter here composes gmres with the Schwarz smoother ALONE.

export helmholtz_shift, schwarz_preconditioner

"""
    helmholtz_shift(m, Ω; slicing = :t) -> s

The CRP complex-shift policy (main.tex:680): the operator-level damping shift `s`
(β-units) for the Schwarz subdomain solves, from the problem data — azimuthal mode
number `m` and orbital frequency `Ω`:

  - `slicing = :t` — `s = (mΩ)²` (the CRP "small complex shift on the order of m²Ω²";
    the complex-shifted-Laplacian scale of Erlangga2006/vanGijzen2007, shift ~ k² for
    the oscillation wavenumber k = mΩ). Even in `mΩ`: retrograde orbits (Ω < 0) get
    the same nonnegative `s`.
  - `slicing = :hyperboloidal` or `:vtu` — `s = 0`: the solution loses its global
    wave-like behavior and the shift "can be reduced or disabled altogether"
    (main.tex:680); the CRP's shipped vtu production config sets `ComplexShift: 0.0`
    (anc/ScalarSelfForce.yaml:194,203).

APPLICATION on this package's MASSIVE Gauss operator (file header, calibration-gated):
build the smoother's operator with the damped coefficient `β − i·s` — its subdomain
matrices are the mass-weighted `A_ss − i·s·M_ss` — and precondition GMRES on the
ORIGINAL operator with it:

    s   = helmholtz_shift(m, Ω)
    ops = MeshOperator(T, mesh, A, (x, y) -> β(x, y) - im * s, γ; C, hooks, bcs)
    sm  = SchwarzSmoother(ops; overlap = 2)
    gmres(op, b; Minv! = schwarz_preconditioner(sm))

Do NOT pass this operator-level `s` to the literal `SchwarzSmoother(op; shift = s)`
seam: that spectre-convention mechanism adds the UNWEIGHTED `−i·s·I`, which at
operator scale dwarfs the massive matrices and destroys the preconditioner (the
round-1 calibration falsifier, evidence/schwarz-complex-shift/) — its useful values
live on the matrix scale.

Sign: with the locked ansatz phase `e^{+im(φ̄ − Ωt)}` (time dependence `e^{−iωt}`,
ω = mΩ), physical absorption enters the frequency-domain operator on the −i side
(`+2ε∂_t → −2iεω` for ω > 0) — `β − i·s` damps on that same side (the `β = −k² − iε`
Helmholtz convention of test_schwarz.jl). The shift enters ONLY the preconditioner,
never the solved operator, so the converged GMRES answer is independent of `s`
(gated in test_schwarz_shift.jl).

Returns `abs2(m*Ω)` in the type of `(m*Ω)²` (BigFloat `Ω` stays BigFloat).
"""
function helmholtz_shift(m::Integer, Ω::Real; slicing::Symbol = :t)
    slicing === :t && return abs2(m * Ω)
    (slicing === :hyperboloidal || slicing === :vtu) && return zero(abs2(m * Ω))
    throw(ArgumentError("helmholtz_shift: unknown slicing $slicing " *
                        "(expected :t, :hyperboloidal, or :vtu)"))
end

"""
    schwarz_preconditioner(sm::SchwarzSmoother; passes = 2) -> Minv!

Right-preconditioner adapter over `sm`: a callable `Minv!(z, v)` running `passes`
Schwarz smoothing passes on `A_sm·z = v` from `z = 0` (the pure LINEARIZED iteration —
[`smooth!`](@ref) from zero needs no operator-on-zero correction on `v`), for the
`Minv!` hook of [`gmres`](@ref). Defaults to 2 passes (CRP main.tex:676). Build `sm`
on the DAMPED operator (see [`helmholtz_shift`](@ref)) to shift-precondition: the
smoothed operator `A_sm` then differs from the solved one, but the GMRES residual,
convergence test, and converged answer remain those of the ORIGINAL system (right
preconditioning).

This composes GMRES with the Schwarz smoother ALONE (the row-18 oracle
configuration); the production multigrid-V-cycle stack is the solver-stack-core node.
"""
function schwarz_preconditioner(sm::SchwarzSmoother; passes::Integer = 2)
    passes ≥ 1 ||
        throw(ArgumentError("schwarz_preconditioner: passes must be ≥ 1, got $passes"))
    return function (z::AbstractVector, v::AbstractVector)
        fill!(z, zero(eltype(z)))
        return smooth!(z, sm, v; passes = passes)
    end
end
