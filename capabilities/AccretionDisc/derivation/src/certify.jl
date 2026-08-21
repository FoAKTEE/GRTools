# certify.jl — the two admission modalities of this project.
#
#   :symbolic        `Symbolics.simplify` of the residual returns an exact structural 0.
#   :exact_rational  the residual is evaluated at ≥N random points whose coordinates are EXACT
#                    `Rational{BigInt}` (ω_d optionally `Complex{Rational{BigInt}}`), and every
#                    residual is EXACTLY `0//1`.  This is an EXACT modality — rational arithmetic has
#                    no rounding — and it is the admission evidence wherever `simplify` stalls on
#                    expression swell (which is the documented risk of Symbolics for these entries).
#
# Points respect the two algebraic constraints of the alphabet exactly:
#   • W²(1−v²) = 1 via the rational parametrisation v = 2t/(1+t²), W = (1+t²)/(1−t²), t ∈ ℚ, |t| < 1;
#   • optionally the hydrostatic balance (11), by COMPUTING p′ from the other coordinates
#     (`onshell = true`) rather than substituting into the expression.

using Random
using Printf
using SymbolicUtils

"Every symbol of the alphabet that a certificate point must bind."
const ALL_SYMBOLS = [enu, epsi, emu, omt, nu_r, psi_r, mu_r, omt_r,
                     nu_rr, psi_rr, mu_rr, omt_rr,
                     v, v_r, v_rr, W,
                     eps, eps_r, eps_rr, p, p_r, p_rr, Gamma1, Gamma1_r,
                     m, omega, omega_d,
                     dp, dp_r, dp_rr, deps, deps_r, dut, dut_r, dur, dur_r, dur_rr, duphi, duphi_r,
                     y0, y0_r, y1, y1_r,
                     htt, htr, htphi, hrr, hrphi, hphiphi, hzz,
                     htt_r, htr_r, htphi_r, hrr_r, hrphi_r, hphiphi_r, hzz_r,
                     htt_rr, htr_rr, htphi_rr, hrr_rr, hrphi_rr, hphiphi_rr, hzz_rr]

const RatB = Rational{BigInt}

_rand_rat(rng; lo = 1, hi = 37) = RatB(rand(rng, lo:hi), rand(rng, lo:hi))
_rand_signed(rng; kw...) = (rand(rng, Bool) ? 1 : -1) * _rand_rat(rng; kw...)

"""
    rational_point(rng; onshell=false, complex_wd=false) -> Dict

One exact-rational certificate point.  `v`, `W` come from the rational parametrisation so that
W²(1−v²) = 1 holds EXACTLY; with `onshell=true`, `p′` is computed from (11) so the point sits on the
background Euler equation; with `complex_wd=true`, ω_d (and hence ω) carries a rational imaginary
part — the Landau shift (62) — and every residual is a `Complex{Rational{BigInt}}`.
"""
function rational_point(rng::AbstractRNG; onshell::Bool = false, complex_wd::Bool = false)
    a = rand(rng, 1:29)
    b = rand(rng, (a + 1):30)
    t = RatB(a, b) * (rand(rng, Bool) ? 1 : -1)      # |t| < 1  ⇒  |v| < 1, W > 1
    vval = 2t / (1 + t^2)
    Wval = (1 + t^2) / (1 - t^2)

    pt = Dict{Num,Any}()
    for s in ALL_SYMBOLS
        pt[s] = _rand_signed(rng)
    end
    pt[enu] = _rand_rat(rng)                          # e^ν, e^ψ, e^μ > 0
    pt[epsi] = _rand_rat(rng)
    pt[emu] = _rand_rat(rng)
    pt[v] = vval
    pt[W] = Wval
    pt[eps] = _rand_rat(rng)
    pt[p] = _rand_rat(rng)
    pt[Gamma1] = _rand_rat(rng)
    pt[m] = RatB(rand(rng, 1:9))                      # azimuthal number: a positive integer
    pt[omega_d] = complex_wd ? Complex(_rand_signed(rng), _rand_rat(rng)) : _rand_signed(rng)
    # ω = ω_d + mΩ (51) keeps ω and ω_d mutually consistent for expressions written in either.
    Ωval = exact_eval(Ω_expr, pt)
    pt[omega] = pt[omega_d] + pt[m] * Ωval
    if onshell
        pt[p_r] = (pt[eps] + pt[p]) * exact_eval(HYDRO_RHS, pt)
    end
    pt
end

"""
    exact_eval(x, pt)

Substitute the point `pt` into `x` and return the EXACT number (`Rational{BigInt}` or
`Complex{Rational{BigInt}}`).  Errors if the result is not fully numeric or is a float.
"""
function exact_eval(x::Num, pt)
    val = Symbolics.value(Symbolics.substitute(x, pt))
    _assert_exact(val, x)
    val
end
function exact_eval(x::Complex{Num}, pt)
    r = Symbolics.substitute(x, pt)
    re = Symbolics.value(real(r))
    imv = Symbolics.value(imag(r))
    _assert_exact(re, x)
    _assert_exact(imv, x)
    re + im * imv
end
exact_eval(x::Number, pt) = x

function _assert_exact(val, x)
    val isa Number ||
        error("exact_eval: expression not fully substituted (residual symbols): $(val)")
    val isa AbstractFloat && error("exact_eval: point produced a FLOAT — modality would not be exact")
    val isa Complex && real(val) isa AbstractFloat &&
        error("exact_eval: point produced a FLOAT — modality would not be exact")
    nothing
end

"""
    Certificate

One admitted (or refused) identity: `name`, `eqref` (the paper equation), `modality`
(`:symbolic` | `:exact_rational`), number of points, verdict, and a human note.
"""
struct Certificate
    name::String
    eqref::String
    modality::Symbol
    npoints::Int
    pass::Bool
    note::String
end

"""
    certify_exact(name, eqref, expr; npts=50, seed=2608, onshell=false, complex_wd=false)

Exact-rational certificate: `expr` (a scalar, or any iterable of scalars) must evaluate to exactly
zero at `npts` independent random rational points.
"""
function certify_exact(name, eqref, expr; npts::Int = 50, seed::Int = 2608,
                       onshell::Bool = false, complex_wd::Bool = false, verbose::Bool = true)
    rng = MersenneTwister(seed)
    exprs = expr isa AbstractArray ? collect(expr) : [expr]
    bad = 0
    firstbad = ""
    for k in 1:npts
        pt = rational_point(rng; onshell, complex_wd)
        for (j, e) in enumerate(exprs)
            val = exact_eval(e, pt)
            if !iszero(val)
                bad += 1
                isempty(firstbad) && (firstbad = "point $k, component $j: residual = $(val)")
            end
        end
    end
    pass = bad == 0
    note = pass ? "$(length(exprs)) component(s) × $npts points, all residuals exactly 0" :
           "$bad nonzero residual(s); first: $firstbad"
    verbose && @printf("[%s] %-46s %-16s %s\n", pass ? "PASS" : "FAIL", name,
                       "exact-rational", note)
    Certificate(name, eqref, :exact_rational, npts, pass, note)
end

"""
    certify_symbolic(name, eqref, expr; kwargs...)

Symbolic certificate: `Symbolics.simplify` of the residual (after eliminating W by
W → 1/√(1−v²), the only algebraic constraint of the alphabet) must be structurally zero.
"""
function certify_symbolic(name, eqref, expr; verbose::Bool = true)
    exprs = expr isa AbstractArray ? collect(expr) : [expr]
    bad = 0
    firstbad = ""
    for (j, e) in enumerate(exprs)
        s = _simplify_zero(e)
        if !s
            bad += 1
            isempty(firstbad) && (firstbad = "component $j did not simplify to 0")
        end
    end
    pass = bad == 0
    note = pass ? "$(length(exprs)) component(s) simplified to exact 0" :
           "$bad component(s) not zero; first: $firstbad"
    verbose && @printf("[%s] %-46s %-16s %s\n", pass ? "PASS" : "FAIL", name, "symbolic", note)
    Certificate(name, eqref, :symbolic, 0, pass, note)
end

# The alphabet carries exactly ONE algebraic constraint, W²(1−v²) = 1 (8).  Substituting the standard
# birational parametrisation of that conic,  v = 2t/(1+t²),  W = (1+t²)/(1−t²),  removes it and turns
# every expression into an element of the RATIONAL FUNCTION FIELD ℚ(t, e^ν, e^ψ, e^μ, ω̃, ν′, …), where
# `simplify_fractions` is a decision procedure for zero.  The map is dominant onto the W > 0 branch of
# the conic, so a pullback that vanishes identically in t proves the original identity modulo the
# constraint.  (Substituting W → 1/√(1−v²) instead leaves uncollapsed (1−v²)^{n/2} powers and STALLS —
# recorded as a finding in the node report.)
@variables tpar
const _CONIC_PARAM = Dict(v => 2tpar / (1 + tpar^2), W => (1 + tpar^2) / (1 - tpar^2))

_dropW(x::Num) = Symbolics.substitute(x, _CONIC_PARAM)
_dropW(x::Complex{Num}) = _dropW(real(x)) + im * _dropW(imag(x))

_zero_num(s::Num) = isequal(Symbolics.simplify_fractions(s), Num(0)) ||
                    isequal(Symbolics.simplify(s; expand = true), Num(0))

function _simplify_zero(e)
    s = _dropW(e)
    s isa Complex{Num} && return _zero_num(real(s)) && _zero_num(imag(s))
    _zero_num(s)
end

"""
    certify_bool(name, eqref, ok, note)

Structural certificate (a variable-content / linearity fact established by inspection of the
expression tree, not by evaluation).
"""
function certify_bool(name, eqref, ok::Bool, note::AbstractString; verbose::Bool = true,
                      modality::Symbol = :structural, npoints::Int = 0)
    verbose && @printf("[%s] %-46s %-16s %s\n", ok ? "PASS" : "FAIL", name, String(modality), note)
    Certificate(name, eqref, modality, npoints, ok, note)
end

"""
    certify_pole_order(name, eqref, exprs; maxorder=1, npts=8, seed=2608)

Corotation pole structure, EXACTLY.  All symbols except ω_d are bound to exact rationals, leaving a
univariate rational function of ω_d; `expand` + `simplify_fractions` puts it in single-fraction
normal form N(ω_d)/D(ω_d) (plain `simplify_fractions` alone leaves a SUM of fractions and is not a
normal form — recorded as a finding); the pole order at corotation is then
ord₀(D) − ord₀(N), each multiplicity obtained by exact differentiation at ω_d = 0//1.  The paper's
singular structure at M:648 (and the Newtonian (50)) demands order ≤ `maxorder` = 1.
"""
function certify_pole_order(name, eqref, exprs; maxorder::Int = 1, npts::Int = 8,
                            seed::Int = 2608, verbose::Bool = true)
    rng = MersenneTwister(seed)
    exprs = exprs isa AbstractArray ? collect(exprs) : [exprs]
    bad = 0
    firstbad = ""
    orders = fill(-99, length(exprs))
    for k in 1:npts
        full = rational_point(rng)
        partial = Dict(kv for kv in full if !isequal(kv.first, omega_d))
        for (j, e) in enumerate(exprs)
            o = pole_order(subs_partial(e, partial), omega_d)
            orders[j] = max(orders[j], o)
            if o > maxorder
                bad += 1
                isempty(firstbad) && (firstbad = "point $k, component $j: pole order $o > $maxorder")
            end
        end
    end
    pass = bad == 0
    note = pass ? "$(length(exprs)) entrie(s) × $npts points; pole order at ω_d = 0: $(orders)" :
           "$bad failure(s); first: $firstbad"
    verbose && @printf("[%s] %-46s %-16s %s\n", pass ? "PASS" : "FAIL", name, "pole-order", note)
    Certificate(name, eqref, :exact_rational, npts, pass, note)
end

"""
    pole_order(e, x) -> Int

Order of the pole of the univariate rational expression `e` at `x = 0` (negative ⇒ a zero).
"""
const NO_POLE = -99   # returned for an identically vanishing component

function pole_order(e::Num, x)
    f = Symbolics.simplify_fractions(Symbolics.expand(e))
    isequal(f, Num(0)) && return NO_POLE
    n, d = _num_den(f)
    isequal(Symbolics.expand(n), Num(0)) && return NO_POLE
    _zero_order(d, x) - _zero_order(n, x)
end
pole_order(e::Complex{Num}, x) = max(pole_order(real(e), x), pole_order(imag(e), x))

function _num_den(f::Num)
    v = Symbolics.value(f)
    if SymbolicUtils.iscall(v) && SymbolicUtils.operation(v) === (/)
        a = SymbolicUtils.arguments(v)
        return (Num(a[1]), Num(a[2]))
    end
    (f, Num(1))
end

function _zero_order(e::Num, x; maxord::Int = 12)
    q = e
    for j in 0:maxord
        val = Symbolics.value(Symbolics.substitute(q, Dict(x => RatB(0))))
        val isa Number || error("_zero_order: not univariate in $x")
        iszero(val) || return j
        q = Symbolics.derivative(q, x)
    end
    error("_zero_order: multiplicity exceeds $maxord")
end

"Partial substitution (leaves some symbols free)."
subs_partial(e::Num, d) = Symbolics.substitute(e, d)
subs_partial(e::Complex{Num}, d) = subs_partial(real(e), d) + im * subs_partial(imag(e), d)

"Pretty markdown table of a certificate list (the log body of a node's evidence)."
function certificate_table(certs::AbstractVector{Certificate})
    io = IOBuffer()
    println(io, "| identity | eq. | modality | points | result |")
    println(io, "|---|---|---|---|---|")
    for c in certs
        println(io, "| `", c.name, "` | ", c.eqref, " | ", c.modality, " | ",
                c.npoints == 0 ? "—" : string(c.npoints), " | ", c.pass ? "PASS" : "**FAIL**", " |")
    end
    String(take!(io))
end
