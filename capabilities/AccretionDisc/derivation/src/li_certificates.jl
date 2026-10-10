# Certificates of the pressure law and reduction only. No LI Newtonian limit,
# stress/flux diagnostics, boundary condition, or runtime target is certified here.

"Θ′ as a derived expression in the original R3 alphabet."
const THETA_R = R3.ddr(p / eps)

"Constant Θ, with Γ₁ matched to the G26 energy-density convention."
const LI_OVERLAP = Dict(Gamma1 => 1 + p / eps, p_r => p * eps_r / eps)

"C/S and reconstruction expressions, without circular metric zeroing."
function closure_exports(reduction)
    vcat([name => ex for (name, ex) in zip(MASTER_NAMES, reduction.MASTER_ONSHELL)],
         [string(name) => ex for (name, ex) in pairs(reduction.ELIMINATED_ONSHELL)])
end

"""
    closure_certificates(outdir; npts=50, seed=191101428, verbose=true)

Run symbolic/structural gates and the existing exact Rational{BigInt} modality.
Sampling uses rational_point's nonuniform measure, including the rational conic
v=2t/(1+t²), W=(1+t²)/(1-t²); no Schwartz–Zippel probability is claimed. Poles
are excluded: p, ε, ε+p, e^ν, e^ψ, e^μ, W and ω_d must be nonzero, as must every
algebraic pivot. The 2×2 determinant is checked explicitly. The certified degree
bound is total degree ≤1 in response/metric amplitudes (background held fixed).
No global numerator-degree bound in the background jets is asserted.

The emitted expressions used for C-LI-0 retain h_tr and h_rφ. Coefficientwise
symbolic identities are holomorphic in the amplitude symbols, so they also hold
for general complex modes. These small text exports are certificate artifacts,
not the separate LI runtime-emission target.
"""
function closure_certificates(outdir; npts::Int = 50, seed::Int = 191101428,
                              verbose::Bool = true)
    npts > 0 || throw(ArgumentError("npts must be positive"))
    C = Certificate[]
    push!(C, certify_symbolic("LI Lagrangian temperature displacement", "Option A",
        ((dp + XI_R * p_r) - (p / eps) * (deps + XI_R * eps_r) -
         eps * THETA_R * XI_R) / p - CLOSURE_RESIDUAL; verbose))
    push!(C, certify_symbolic("LI density reconstruction", "Option A",
        DEPS_SOL - eps * dp / p; verbose))
    push!(C, certify_symbolic("Gamma-only substitution misses Θ′ξʳ", "Option A",
        subs(R3.CLOSURE_RESIDUAL, Dict(Gamma1 => 1 + p / eps)) -
        CLOSURE_RESIDUAL - eps * THETA_R * XI_R / p; verbose))
    push!(C, certify_symbolic("C-LI-0 pressure residual", "Theta'=0, Gamma1=1+Theta",
        subs(R3.CLOSURE_RESIDUAL - CLOSURE_RESIDUAL, LI_OVERLAP); verbose))

    A = [dsym(eq, x) for eq in (EQ_R, EQ_E), x in (y0_r, y1_r)]
    detA = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
    push!(C, certify_symbolic("LI first-order 2x2 derivative matrix", "(47)-(49)",
        [A[1, 1] - (eps + p) / R3.emu^2, A[1, 2], A[2, 1],
         A[2, 2] - im * (eps + p), detA - im * (eps + p)^2 / R3.emu^2]; verbose))
    push!(C, certify_exact("LI backsub Euler/energy/norm/pressure", "(15), Option A",
        [backsub(DEULER_UP[IR]), backsub(DENERGY), backsub(DEULER_UP[IPHI]),
         backsub(DNORM), backsub(CLOSURE_RESIDUAL)]; npts, seed, complex_wd = true, verbose))
    push!(C, certify_exact("LI unused t-Euler on shell", "(15c),(11)",
        backsub(DEULER_UP[IT]); npts, seed = seed + 1, onshell = true,
        complex_wd = true, verbose))

    amplitudes = vcat([y0, y1], HSYMS)
    ex = closure_exports(@__MODULE__)
    all_exprs = Any[MASTER_RAW...; last.(ex)...]
    affine = all(_is_zero_expr(dsym(dsym(e, x), z))
                 for e in Any[YPRIME...; all_exprs...] for x in amplitudes for z in amplitudes)
    hfree = all(isempty(occurring(e, HSYMS)) for e in MASTER_ONSHELL[1:4])
    source = all(_is_zero_expr(subs(e, Dict(h => Num(0) for h in HSYMS)))
                 for e in MASTER_ONSHELL[5:6])
    noextra = all(no_second_derivatives(e) && isempty(occurring(e, ELIMINATED_DERIVS)) &&
                  isempty(occurring(e, [y0_r, y1_r, dut, duphi, deps, dp, dp_r, dur, dur_r]))
                  for e in all_exprs)
    push!(C, certify_bool("LI y'=C*y+S, first order and affine", "(47)-(49)",
        affine && hfree && source && noextra,
        "amplitude degree <=1; C metric-free; S vanishes at h=0; no remaining derivatives"; verbose))
    known = Set(string.(R3.ALL_SYMBOLS))
    used = union((_names(e) for e in Any[all_exprs...; THETA_R; CLOSURE_RESIDUAL])...)
    push!(C, certify_bool("LI no new primitive symbols", "Theta=p/eps",
        issubset(used, known) && !("Gamma1" in used) && !("Gamma1_r" in used),
        "Theta and Theta' are expressions; LI coefficients do not contain Gamma1"; verbose))

    # Emit each closure BEFORE imposing the overlap. Compare what was actually
    # read back, including all four reconstructed quantities and complex parts.
    mkpath(outdir)
    paths = [joinpath(outdir, "overlap_$(label).jls.txt") for label in ("adiabatic", "li")]
    for (path, reduction) in zip(paths, (R3, @__MODULE__))
        R3.write_exports(path, closure_exports(reduction);
            header = "# C-LI-0 certificate coefficients; general metric modes, no circular zeroing")
    end
    ad, li = R3.read_artifact.(paths)
    @assert first.(ad) == first.(li)
    residuals = [subs(a - b, LI_OVERLAP)
                 for ((_, av), (_, bv)) in zip(ad, li) for (a, b) in zip(av, bv)]
    # Exact-rational crosscheck FIRST: it rejects a non-identity in seconds, whereas
    # simplify on a nonzero emitted coefficient ran >10 min without terminating
    # (review mutants M2/M3). The symbolic gate runs only once this one passes.
    cross = certify_exact("C-LI-0 emitted exact-rational crosscheck", "Theta'=0, Gamma1=1+Theta",
        residuals; npts, seed = seed + 2, complex_wd = true, verbose)
    push!(C, cross)
    cross.pass || return C
    # Both reductions have the SAME geometry. This identity holds even with W
    # free, so avoid the conic pullback in certify_symbolic: it needlessly swells
    # the rational polynomials. No numeric sampling enters this symbolic gate.
    for (j, (name, _)) in enumerate(ad)
        verbose && (println("C-LI-0 symbolic emitted coefficient: ", name); flush(stdout))
        exact_zero = all(R3._zero_num(residuals[k]) for k in (2j - 1, 2j))
        push!(C, certify_bool("C-LI-0 emitted $name", "Theta'=0, Gamma1=1+Theta",
            exact_zero, "both complex parts simplify to zero with W free";
            verbose, modality = :symbolic))
    end
    C
end
