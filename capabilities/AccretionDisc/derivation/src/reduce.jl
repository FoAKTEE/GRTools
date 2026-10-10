# reduce.jl — route R3 stage 2, part 2: the reduction scheme of arXiv:2608.17003 M:604–608 executed
# in Julia, from the linearised equations of perturbation.jl.
#
# The four steps, in the order forced by the variable content of the equations (printed by
# `reduction_structure()`; it matches the route-R2 report exactly):
#
#   (ii)  pressure closure: adiabatic (19)–(21) (default) or Option A locally isothermal,
#         ξ^r = iδû^r/(u^tω_d)                                              →  δε̂   [needs δp̂, δû^r]
#   (iii) φ-Euler δℰ^φ = 0 — algebraic in δû^φ, no δp̂′/δû^r′/δû^t          →  δû^φ  [needs δp̂, δû^r, ĥ]
#   (i)   perturbed normalisation 2u_αδû^α + ĥ_{αβ}u^αu^β = 0              →  δû^t  [needs δû^φ, ĥ]
#   (iv)  r-Euler (carries δp̂′, not δû^r′) + energy equation (carries δû^r′,
#         not δp̂′) → a 2×2 linear solve for (y0′, y1′)                     →  C, S
#
# with the master variables (48),(49)  y0 = δÛ = δp̂/(ε+p),  y1 = −iδû^r.
#
# Nothing is "collected by hand": every elimination is `solve_linear`, i.e. x = −E|_{x=0}/(∂E/∂x)
# together with an assertion that E really is affine in x.  Blanket `simplify` is deliberately NOT
# applied to the entries — the admission gate is the exact back-substitution certificate below.

# ---------------------------------------------------------------------------------------------
# linear-algebra helpers over Num / Complex{Num}
# ---------------------------------------------------------------------------------------------

"∂e/∂x, componentwise on a complex expression."
dsym(e::Num, x) = Symbolics.derivative(e, x)
dsym(e::Complex{Num}, x) = Symbolics.derivative(real(e), x) + im * Symbolics.derivative(imag(e), x)

"Structural zero test after `simplify`, componentwise on a complex expression."
_is_zero_expr(e::Num) = isequal(Symbolics.simplify(e), Num(0))
_is_zero_expr(e::Complex{Num}) = _is_zero_expr(real(e)) && _is_zero_expr(imag(e))

"Substitute a dictionary into a real or complex expression."
subs(e::Num, d) = Symbolics.substitute(e, d)
subs(e::Complex{Num}, d) = subs(real(e), d) + im * subs(imag(e), d)
subs(e::AbstractArray, d) = map(x -> subs(x, d), e)

"""
    lin_subs(E, rules)

Replace symbols by possibly COMPLEX expressions inside an expression that is AFFINE in them, via
coefficient extraction: `Symbolics.substitute` refuses to bind a complex value to a real-typed symbol
(`promote_symtype(complex, …)` asserts `S <: Real`), so every perturbation-amplitude replacement of
this reduction goes through here.  Affineness in each replaced symbol is asserted, never assumed —
which is also a live check that the linearisation really is linear.
"""
function lin_subs(E, rules::AbstractDict)
    acc = E
    for (x, val) in rules
        a = dsym(acc, x)
        iszero_entry(a) && continue
        _is_zero_expr(dsym(a, x)) || error("lin_subs: expression is not affine in $x")
        acc = subs(acc, Dict(x => Num(0))) + a * val
    end
    acc
end

"""
    solve_linear(E, x) -> x_sol

Solve `E = 0` for `x`, given that `E` is affine in `x`: x = −E|_{x=0} / (∂E/∂x).  Asserts affineness
(∂²E/∂x² = 0) so that a hidden nonlinearity can never pass silently.
"""
function solve_linear(E, x)
    a = dsym(E, x)
    iszero_entry(a) && error("solve_linear: expression does not contain $x")
    _is_zero_expr(dsym(a, x)) || error("solve_linear: expression is not affine in $x")
    b = subs(E, Dict(x => Num(0)))
    -b / a
end

"""
    solve_linear2(E1, E2, x1, x2) -> (x1_sol, x2_sol)

Solve the 2×2 affine system `E1 = E2 = 0` for `(x1, x2)` by Cramer's rule.
"""
function solve_linear2(E1, E2, x1, x2)
    a11, a12 = dsym(E1, x1), dsym(E1, x2)
    a21, a22 = dsym(E2, x1), dsym(E2, x2)
    z = Dict(x1 => Num(0), x2 => Num(0))
    b1, b2 = subs(E1, z), subs(E2, z)
    det = a11 * a22 - a12 * a21
    ((-b1 * a22 + b2 * a12) / det, (-a11 * b2 + a21 * b1) / det)
end

# ---------------------------------------------------------------------------------------------
# step (ii) — pressure closure  ⇒  δε̂   (adiabatic (19)-(21) by default; Option A LI selectable)
# ---------------------------------------------------------------------------------------------

"""
    pressure_closure([Val(:adiabatic)])
    pressure_closure(Val(:locally_isothermal))

Pressure residual for the reduction. The default is G26 (19)-(21). Option A
`Val(:locally_isothermal)` imposes p = Θ(r)ε, δΘ = 0 in the G26 singular-field
gauge in Boyer–Lindquist coordinates, with Θ = p/ε (no extra primitive).
Energy/Euler equations and conserved baryons are unchanged; there is no bath.

Since Δ = δ + ξʳ∂r, the LI residual is equivalently
    [Δp - ΘΔε - εΘ′ξʳ]/p,
where εΘ′ = p′ - (p/ε)ε′. Thus the displacement term is retained even though
it cancels from the Eulerian residual. Changing Γ₁ alone would omit it.
Together with retained (19), Δp = -Θ(ε+p)q^{αβ}Δg_{αβ}/2 + εΘ′ξʳ.
This pressure law is gauge-fixed when Θ′ ≠ 0; it is not a cooling model.
"""
pressure_closure(::Val{:adiabatic} = Val(:adiabatic)) =
    (dp + XI_R * p_r) / p - Gamma1 * (deps + XI_R * eps_r) / (eps + p)
pressure_closure(::Val{:locally_isothermal}) = dp / p - deps / eps

"""
    CLOSURE_RESIDUAL

`pressure_closure(PRESSURE_CLOSURE)` when the including module defines `PRESSURE_CLOSURE`
(the `LocallyIsothermal` submodule), otherwise the adiabatic default Δp/p − Γ₁Δε/(ε+p) with
Δ = δ + 𝓛_ξ acting on the scalars p, ε: Δp̂ = δp̂ + ξ^r p′, Δε̂ = δε̂ + ξ^r ε′ (p and ε depend on
r only) — Eqs. (19)+(20) of the paper.
"""
const CLOSURE_RESIDUAL = pressure_closure(
    isdefined(@__MODULE__, :PRESSURE_CLOSURE) ? PRESSURE_CLOSURE : Val(:adiabatic))

"δε̂ from the closure — step (ii)."
const DEPS_SOL = solve_linear(CLOSURE_RESIDUAL, deps)

# ---------------------------------------------------------------------------------------------
# step (iii) — φ-Euler  ⇒  δû^φ ;  step (i) — normalisation  ⇒  δû^t
# ---------------------------------------------------------------------------------------------

"δû^φ from the φ-component of the linearised Euler equation — step (iii)."
const DUPHI_SOL = solve_linear(lin_subs(DEULER_UP[IPHI], Dict(deps => DEPS_SOL)), duphi)

"δû^t from the perturbed normalisation — step (i)."
const DUT_SOL = solve_linear(lin_subs(DNORM, Dict(duphi => DUPHI_SOL)), dut)

"""
    ALG_SUB

The three algebraic eliminations as one substitution, applied in dependency order
(δε̂ → δû^φ → δû^t).  After it, an expression contains only δp̂, δp̂′, δû^r, δû^r′ and ĥ.
"""
const ALG_SUB = [dut => DUT_SOL, duphi => DUPHI_SOL, deps => DEPS_SOL]

"Apply the algebraic eliminations of the reduction, in dependency order (δû^t ← δû^φ ← δε̂ ← δp̂,δû^r)."
eliminate(e) = foldl((acc, r) -> lin_subs(acc, Dict(r)), ALG_SUB; init = e)

# ---------------------------------------------------------------------------------------------
# step (iv) — the two ODEs  ⇒  C, S
# ---------------------------------------------------------------------------------------------

"""
    MASTER_SUB

δp̂ = (ε+p)y0, δp̂′ = (ε′+p′)y0 + (ε+p)y0′, δû^r = i y1, δû^r′ = i y1′ — the change of variables to
the master state vector (49).
"""
const MASTER_SUB = Dict(dp => (eps + p) * y0,
                        dp_r => (eps_r + p_r) * y0 + (eps + p) * y0_r,
                        dur => im * y1,
                        dur_r => im * y1_r)

"r-Euler equation in master variables (the row carrying y0′)."
const EQ_R = lin_subs(eliminate(DEULER_UP[IR]), MASTER_SUB)
"Energy equation in master variables (the row carrying y1′)."
const EQ_E = lin_subs(eliminate(DENERGY), MASTER_SUB)

"(y0′, y1′) as explicit expressions — the solution of the 2×2 system of step (iv)."
const YPRIME = solve_linear2(EQ_R, EQ_E, y0_r, y1_r)

"""
    master_entries(x0p, x1p) -> (C00, C01, C10, C11, S0, S1)

Split y′ = C·y + S: C_{ij} = ∂y_i′/∂y_j and S_i = y_i′|_{y=0}.  Linearity in (y0, y1) is the content
of the perturbative expansion, and is asserted by `stage2_certificates`.
"""
function master_entries(x0p, x1p)
    z = Dict(y0 => Num(0), y1 => Num(0))
    (dsym(x0p, y0), dsym(x0p, y1), dsym(x1p, y0), dsym(x1p, y1), subs(x0p, z), subs(x1p, z))
end

"The RAW master entries (ν′ and p′ independent, i.e. the background NOT put on (11)), in ω_d."
const MASTER_RAW = to_omega_d.(master_entries(YPRIME...))

"""
    ELIM_NU_R

The ON-SHELL substitution that eliminates ν′ through the hydrostatic balance (11):
ν′ = v²ψ′ − e^{ψ−ν}vω̃′ − p′/[(ε+p)W²].  Route R2 used exactly this direction ("on-shell form
eliminates ν′ by (11)"), which is why the emitted entries carry p′ but not ν′.
"""
const ELIM_NU_R = Dict(nu_r => v^2 * psi_r - epsi * v * omt_r / enu - p_r / ((eps + p) * W^2))

"Put an expression on the background Euler equation (11) by eliminating ν′."
onshell_nu(e) = subs(e, ELIM_NU_R)

"The ON-SHELL master entries (C00, C01, C10, C11, S0, S1) — the route-R2 convention."
const MASTER_ONSHELL = onshell_nu.(MASTER_RAW)

"δû^t, δû^φ, δε̂ and the COVARIANT δû_φ in master variables (y0, y1, ĥ), on-shell, in ω_d."
const ELIMINATED_ONSHELL = let
    quantities = (dut = DUT_SOL, duphi_up = DUPHI_SOL, deps = DEPS_SOL,
                  duphi_low = eliminate(DUPHI_LOW))
    map(q -> onshell_nu(to_omega_d(lin_subs(q, MASTER_SUB))), quantities)
end

"""
    reduction_structure()

The variable content of each linearised equation — the structural facts that FORCE the reduction
order (and that route R2 certified independently).  Returned as a `NamedTuple` of sorted name lists.
"""
function reduction_structure()
    names(e) = sort(string.(collect(Symbolics.get_variables(e isa Complex ? real(e) + imag(e) : e))))
    keep = Set(["dp", "dp_r", "dur", "dur_r", "dut", "duphi", "deps",
                "htt", "htr", "htphi", "hrr", "hrphi", "hphiphi", "hzz",
                "htt_r", "htr_r", "htphi_r", "hrr_r", "hrphi_r", "hphiphi_r", "hzz_r"])
    f(e) = sort(collect(intersect(Set(names(e)), keep)))
    (euler_t = f(DEULER_UP[IT]), euler_r = f(DEULER_UP[IR]), euler_phi = f(DEULER_UP[IPHI]),
     euler_z = f(DEULER_UP[IZ]), energy = f(DENERGY), norm = f(DNORM))
end
