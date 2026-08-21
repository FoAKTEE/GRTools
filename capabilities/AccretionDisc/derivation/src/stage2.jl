# stage2.jl — the admission gate of the route-R3 master system (the CRIT node EGD::jd-perturbation-derive).
#
# The gate is INTERNAL to route R3: nothing here reads the notebook (R1) or the Mathematica (R2)
# expressions.  It answers "do the derived δû^t, δû^φ, δε̂ and dy/dr = C·y + S actually solve the
# perturbed system of the paper?", by putting them back into the linearised equations of
# perturbation.jl and demanding an EXACT zero at exact-rational points.
#
# The t-component of the Euler equation is the non-tautological one: it is never used in the
# reduction (u_μℰ^μ ≡ 0 makes it dependent ON the background Euler equation), so its residual
# vanishing on-shell certifies the whole elimination chain.

"""
    backsub(e)

Put an O(λ) equation on the derived solution: eliminate δû^t, δû^φ, δε̂, change to the master
variables, then replace y0′, y1′ by the derived right-hand sides C·y + S.
"""
backsub(e) = lin_subs(lin_subs(eliminate(e), MASTER_SUB),
                      Dict(y0_r => YPRIME[1], y1_r => YPRIME[2]))

"Names of the six master entries, in the order of `MASTER_RAW` / `MASTER_ONSHELL`."
const MASTER_NAMES = ("C00", "C01", "C10", "C11", "S0", "S1")

"The metric-mode symbols ĥ_μν and ĥ′_μν."
const HSYMS = [htt, htr, htphi, hrr, hrphi, hphiphi, hzz,
               htt_r, htr_r, htphi_r, hrr_r, hrphi_r, hphiphi_r, hzz_r]

_names(e) = Set(string.(Symbolics.get_variables(e isa Complex ? real(e) + imag(e) : e)))

"Symbols of `list` that occur in `e`."
occurring(e, list) = [string(s) for s in list if string(Symbolics.value(s)) in _names(e)]

"""
    stage2_certificates(; npts=50, seed=2608, verbose=true) -> Vector{Certificate}

The stage-2 certificate suite.
"""
function stage2_certificates(; npts::Int = 50, seed::Int = 2608, verbose::Bool = true)
    C = Certificate[]

    # --- the displacement/velocity relation (16) that defines ξ^r -------------------------------
    push!(C, certify_symbolic("(16) inversion  𝓛_uξ^r − δû^r", "(16)",
                              lie_u_xi_r(XI_R) - dur; verbose))
    push!(C, certify_symbolic("z-Euler  δℰ^z ≡ 0", "(15c)", DEULER_UP[IZ]; verbose))

    # --- structural facts that force (and validate) the reduction order ------------------------
    st = reduction_structure()
    push!(C, certify_bool("φ-Euler algebraic in δû^φ", "(M:604)",
                          !("dp_r" in st.euler_phi) && !("dur_r" in st.euler_phi) &&
                          !("dut" in st.euler_phi) && ("duphi" in st.euler_phi),
                          "δℰ^φ carries $(st.euler_phi)"; verbose))
    push!(C, certify_bool("r-Euler carries δp̂′ but not δû^r′", "(M:606)",
                          ("dp_r" in st.euler_r) && !("dur_r" in st.euler_r),
                          "δℰ^r carries $(st.euler_r)"; verbose))
    push!(C, certify_bool("energy eq. carries δû^r′ but not δp̂′", "(M:606)",
                          ("dur_r" in st.energy) && !("dp_r" in st.energy),
                          "δ(energy) carries $(st.energy)"; verbose))

    # --- back-substitution: the derived solution solves the perturbed system --------------------
    push!(C, certify_exact("back-substitution  δℰ^r", "(15c)", backsub(DEULER_UP[IR]);
                           npts, seed, verbose))
    push!(C, certify_exact("back-substitution  energy equation", "(15b)", backsub(DENERGY);
                           npts, seed, verbose))
    push!(C, certify_exact("back-substitution  δℰ^φ", "(15c)", backsub(DEULER_UP[IPHI]);
                           npts, seed, verbose))
    push!(C, certify_exact("back-substitution  normalisation", "(M:604)", backsub(DNORM);
                           npts, seed, verbose))
    push!(C, certify_exact("back-substitution  adiabatic closure", "(19)-(21)",
                           backsub(CLOSURE_RESIDUAL); npts, seed, verbose))
    push!(C, certify_exact("back-substitution  δℰ^t (on-shell, NOT used in the reduction)",
                           "(15c),(11)", onshell(backsub(DEULER_UP[IT]));
                           npts, seed, onshell = true, verbose))
    push!(C, certify_exact("back-substitution  δℰ^r, energy — Landau-shifted ω_d ∈ ℂ", "(62)",
                           [backsub(DEULER_UP[IR]), backsub(DENERGY)];
                           npts = max(10, npts ÷ 5), seed = seed + 1, complex_wd = true, verbose))

    # --- shape of the answer --------------------------------------------------------------------
    hfree = all(isempty(occurring(e, HSYMS)) for e in MASTER_ONSHELL[1:4])
    push!(C, certify_bool("C is ĥ-free", "(47)", hfree,
                          hfree ? "C00…C11 contain no ĥ" :
                          "ĥ leaked into C: $(occurring.(MASTER_ONSHELL[1:4], Ref(HSYMS)))"; verbose))
    zeroh = Dict(s => Num(0) for s in HSYMS)
    slin = all(_is_zero_expr(subs(e, zeroh)) for e in MASTER_ONSHELL[5:6])
    push!(C, certify_bool("S vanishes at ĥ = 0 (purely sourced)", "(47)", slin,
                          "S0, S1 → 0 when every ĥ_μν = 0"; verbose))
    nosec = all(no_second_derivatives(e) for e in MASTER_ONSHELL) &&
            all(isempty(occurring(e, ELIMINATED_DERIVS)) for e in MASTER_ONSHELL)
    push!(C, certify_bool("no second derivatives, no eliminated-velocity derivatives", "(47)", nosec,
                          "entries depend on first background derivatives only"; verbose))
    push!(C, certify_bool("S source content", "(47)",
                          true, "S0 ← $(occurring(MASTER_ONSHELL[5], HSYMS)); " *
                                "S1 ← $(occurring(MASTER_ONSHELL[6], HSYMS))"; verbose))

    # --- corotation singular structure -----------------------------------------------------------
    push!(C, certify_pole_order("corotation pole ≤ 1 in C, S", "(M:648)", collect(MASTER_ONSHELL);
                                maxorder = 1, npts = 8, seed, verbose))
    push!(C, certify_pole_order("corotation pole ≤ 1 in δû^t, δû^φ, δε̂, δû_φ", "(M:648)",
                                collect(ELIMINATED_ONSHELL); maxorder = 1, npts = 8, seed, verbose))
    C
end
