# codegen.jl — route R3 stage 4 (EGD::jd-codegen): the Julia printer that turns the serialized R3
# expressions into the GENERATED sources of `AccretionDisc.SpiralWaves`.
#
# INPUT  : `out/R3_master.jls.txt`, `out/R3_diagnostics.jls.txt` — the hand-off artefacts of
#          `serialize.jl` (one `name = <expr>` line per exported quantity, in Symbolics' printed
#          form).  The emitter reads the ARTEFACT, never the in-memory expressions, so the emitted
#          file is traceable to a committed, content-hashed input.
# OUTPUT : `src/SpiralWaves/master_derived.jl`  — `_master_entries(inp::MasterInputs)`
#          `src/SpiralWaves/r2_derived.jl`      — `_r2_delta_u_phi` / `_r2_amf` / `_r2_torque_density`
#          with the SAME names and signatures as the route-R1/R2 generated files they replace, so the
#          swap (EGD::jd-swap) is a one-line include change at every call site.
#
# WHY NOT `Symbolics.build_function` (REUSE.md item 4 — evaluated first, recorded here):
#   `build_function(exprs, args...; expression = Val{true})` does emit a correct Julia function, but
#   its output is a `RuntimeGeneratedFunction`-shaped `:(function (ˍ₋arg1, …) …)` with generated
#   argument names, no docstring, no `@inline`, no `MasterInputs` destructuring and no control over
#   the emission conventions (it prints the same expression tree, so the n//d rationals and the
#   `enu^a` hoisting come for free either way — those are properties of the ARTEFACT alphabet, not of
#   the printer).  What is actually needed here is the WRAPPER (signature, preamble, docstring,
#   header provenance) plus the two alphabet conversions below; the expression body itself is
#   `string(::Num)`, i.e. exactly the printer that produced the artefact.  So: `build_function` is
#   used for nothing, and `string` + a 40-line wrapper emitter does the whole job.  Measured: the
#   emitted bodies parse back with `Meta.parse` and evaluate EXACTLY equal to the Symbolics
#   expressions at random `Rational{BigInt}` points (`codegen_selfcheck`).
#
# TWO ALPHABET CONVERSIONS (both are consumer conventions, not physics):
#   1. `circmetriczero` — ĥ_tr = ĥ_rφ = 0 for the paper's singular field (52)-(54) (see
#      `CIRC_METRIC_ZERO`: a property of that SOURCE, not of circular orbits in general — the reason
#      `MasterInputs` has no fields for them).  Route R1 applies it in the notebook; here it is
#      applied at emission time and is what makes every imaginary part of the master system vanish
#      identically ("Complex only through omega_d").
#   2. `MasterInputs` carries Ω and Ω′, route R3 carries ω̃ and v′.  The two are related by (10) and
#      its r-derivative — Ω = ω̃ + v e^{ν−ψ}, Ω′ = ω̃′ + e^{ν−ψ}[v′ + v(ν′−ψ′)] — so
#      ω̃ = Ω − v e^{ν−ψ} and v′ = (Ω′−ω̃′)e^{ψ−ν} − v(ν′−ψ′).  Substituting these into the master
#      entries makes the derived `_master_entries` a DROP-IN replacement for the route-R1 one (same
#      26 inputs, same signature).  The diagnostics keep the route-R2 signature
#      `(inp, omt, vp, y0, y1)`, which passes ω̃ and v′ explicitly, so they are NOT converted.

using SHA

# ---------------------------------------------------------------------------------------------
# the consumer alphabet
# ---------------------------------------------------------------------------------------------

# Ω and Ω′ are NOT symbols of the derivation (Ω is the expression (10)); they are the two fields
# `MasterInputs` actually carries, so the emitter needs them as symbols to print into.
@variables Omega Omega_r

"Every symbol the artefact may mention, by printed name."
const ARTIFACT_ALPHABET = let d = Dict{Symbol,Num}()
    for s in vcat(ALL_SYMBOLS, [Omega, Omega_r])
        d[Symbol(string(Symbolics.value(s)))] = s
    end
    d
end

"Operators the artefact's printed form may use (anything else is a symbol that must be known)."
const ARTIFACT_OPS = Set([:+, :-, :*, :/, :^, ://])

"""
    CIRC_METRIC_ZERO

ĥ_tr = ĥ_rφ = 0 — `circmetriczero` of the notebook.  This is a property of THE PAPER'S SOURCE, not a
theorem about circular orbits (adversarial review, EGD::jd-adversarial-review, finding 3): the
singular field (52)-(53) is ĥ_μν ∝ g_μν(w) + 2v_μv_ν, the circular worldline (54) has v^r = 0, and the
line element (4) has g_tr = g_rφ = 0 — hence those two components of THIS field vanish.  A different
source or gauge need not satisfy it, which is why the derivation itself keeps `htr`, `hrphi`
SYMBOLIC and the zeroing happens only here, at the boundary with `MasterInputs` (which has no fields
for them).  Consuming code must supply modes in the same gauge/source class; `SpiralWaves.source.jl`
does (`singular_field_mode` returns exactly ĥ_tt, ĥ_tφ, ĥ_φφ, ĥ_rr, ĥ_zz)."""
const CIRC_METRIC_ZERO = Dict{Num,Num}(htr => Num(0), hrphi => Num(0),
                                       htr_r => Num(0), hrphi_r => Num(0),
                                       htrc => Num(0), hrphic => Num(0),
                                       htrc_r => Num(0), hrphic_r => Num(0))

"ω̃ and v′ in terms of the Ω, Ω′ that `MasterInputs` carries — Eq. (10) and its r-derivative."
const TO_MASTERINPUTS = Dict{Num,Num}(omt => Omega - v * enu / epsi,
                                      v_r => (Omega_r - omt_r) * epsi / enu - v * (nu_r - psi_r))

"`MasterInputs` fields, in the order of the struct (= the preamble order of the generated files)."
const MI_FIELDS = [:m, :nu, :nu_r, :psi, :psi_r, :mu, :mu_r, :omt_r, :W, :v, :Omega, :Omega_r,
                   :omega_d, :p, :p_r, :eps, :eps_r, :Gamma1,
                   :htt, :htt_r, :htphi, :htphi_r, :hphiphi, :hphiphi_r, :hrr, :hzz]

"Barred symbol → the expression that binds it in a generated preamble."
const CONJ_BINDINGS = [:y0c => "conj(y0)", :y1c => "conj(y1)", :httc => "conj(htt)",
                       :htphic => "conj(htphi)", :hphiphic => "conj(hphiphi)",
                       :hrrc => "conj(hrr)", :hzzc => "conj(hzz)"]

# ---------------------------------------------------------------------------------------------
# reading the artefact
# ---------------------------------------------------------------------------------------------

_artifact_node(x, tbl) = x
_artifact_node(s::Symbol, tbl) =
    s in ARTIFACT_OPS ? s :
    (haskey(tbl, s) ? tbl[s] : error("codegen: artefact mentions the unknown symbol `$s`"))
_artifact_node(e::Expr, tbl) = Expr(e.head, Any[_artifact_node(a, tbl) for a in e.args]...)

"""
    parse_artifact_expr(rhs) -> Num

Read one printed expression of the artefact back into a `Num`: `Meta.parse` it, replace the symbol
leaves by the symbols of [`ARTIFACT_ALPHABET`](@ref) and evaluate the resulting arithmetic tree.
An unknown identifier is an ERROR, never a fresh symbol.
"""
parse_artifact_expr(rhs::AbstractString) = eval(_artifact_node(Meta.parse(rhs), ARTIFACT_ALPHABET))

"""
    read_artifact(path) -> Vector{Pair{String,Tuple{Num,Num}}}

The exported quantities of an artefact file, IN FILE ORDER, as `name => (re, im)` (a real export has
`im = 0`).  Comment and blank lines are skipped; every other line must be `name[.re|.im] = <expr>`.
"""
function read_artifact(path::AbstractString)
    out = Pair{String,Tuple{Num,Num}}[]
    pos = Dict{String,Int}()
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        occursin(" = ", s) || error("codegen: malformed artefact line: $s")
        lhs, rhs = split(s, " = ", limit = 2)
        lhs = strip(lhs)
        base, slot = endswith(lhs, ".re") ? (lhs[1:end-3], 1) :
                     endswith(lhs, ".im") ? (lhs[1:end-3], 2) : (lhs, 1)
        e = parse_artifact_expr(rhs)
        if !haskey(pos, base)
            push!(out, String(base) => (Num(0), Num(0)))
            pos[base] = length(out)
        end
        cur = out[pos[base]].second
        out[pos[base]] = out[pos[base]].first => (slot == 1 ? (e, cur[2]) : (cur[1], e))
    end
    out
end

# ---------------------------------------------------------------------------------------------
# conversions + printing
# ---------------------------------------------------------------------------------------------

_sub(e::Num, d) = Symbolics.substitute(e, d)

"Apply `circmetriczero` (both parts of a complex export)."
circmetriczero(p::Tuple{Num,Num}) = (_sub(p[1], CIRC_METRIC_ZERO), _sub(p[2], CIRC_METRIC_ZERO))

"Apply `circmetriczero` AND the ω̃, v′ → Ω, Ω′ change of variables of `MasterInputs`."
to_masterinputs(p::Tuple{Num,Num}) =
    (_sub(_sub(p[1], CIRC_METRIC_ZERO), TO_MASTERINPUTS),
     _sub(_sub(p[2], CIRC_METRIC_ZERO), TO_MASTERINPUTS))

"Printed names of the symbols occurring in `e`."
symbol_names(e::Num) = Set(Symbol(string(s)) for s in Symbolics.get_variables(e))
symbol_names(p::Tuple{Num,Num}) = union(symbol_names(p[1]), symbol_names(p[2]))

"""
    print_scalar(e) -> String

The emitted form of one expression: `string(::Num)`, i.e. the same printer that produced the
artefact, with the emission conventions ASSERTED rather than assumed — no float literal, no
transcendental call (`enu`, `epsi`, `emu`, `W` are symbols), and no identifier outside the alphabet.
"""
function print_scalar(e::Num)
    s = string(e)
    occursin(r"[0-9]\.[0-9]|[0-9]e[-+][0-9]|Inf|NaN", s) &&
        error("codegen: emitted body contains a float literal — precision would not be inherited")
    occursin(r"\b(exp|sqrt|log|sin|cos|tan|abs)\(", s) &&
        error("codegen: emitted body contains a transcendental call — the alphabet must stay rational")
    for id in Set(Symbol(m.match) for m in eachmatch(r"[A-Za-z_][A-Za-z_0-9]*", s))
        haskey(ARTIFACT_ALPHABET, id) ||
            error("codegen: emitted body mentions the unknown identifier `$id`")
    end
    s
end

"""
    print_value(p) -> String

`re` for a real export, `(re) + im*(imag)` for a complex one (the route-R2 convention: the only
`Complex` inputs of the generated bodies are `omega_d` and the master amplitudes `y0, y1`).
"""
function print_value(p::Tuple{Num,Num})
    is_structural_zero(p[2]) && return print_scalar(p[1])
    is_structural_zero(p[1]) && return string("im * (", print_scalar(p[2]), ")")
    string("(", print_scalar(p[1]), ") + im * (", print_scalar(p[2]), ")")
end

"Structural zero test that never runs `simplify` on a swollen tree (which is the documented stall)."
is_structural_zero(e::Num) = isequal(Symbolics.value(e), 0) ||
    (length(string(e)) < 200 && isequal(Symbolics.simplify(e), Num(0)))

# ---------------------------------------------------------------------------------------------
# the emitted files
# ---------------------------------------------------------------------------------------------

_sha(path) = bytes2hex(SHA.sha256(read(path)))

"sha256 over the concatenated route-R3 derivation sources, in sorted file order."
function derivation_sha(dir::AbstractString = joinpath(@__DIR__))
    ctx = SHA.SHA256_CTX()
    for f in sort(filter(endswith(".jl"), readdir(dir)))
        SHA.update!(ctx, read(joinpath(dir, f)))
    end
    bytes2hex(SHA.digest!(ctx))
end

function _header(io, file, artifact, artifact_path, content)
    println(io, "# $file -- GENERATED FILE, DO NOT EDIT BY HAND (AccretionDisc.SpiralWaves, arXiv:2608.17003).")
    println(io, "# route     : R3 -- the JULIA-NATIVE symbolic derivation of capabilities/AccretionDisc/derivation/")
    println(io, "#             (Symbolics.jl, from the paper's equations only). NO external CAS is in this build path:")
    println(io, "#             the supplementary notebook (route R1) and the Mathematica derivation (route R2) are")
    println(io, "#             ORACLES that gate this file in test/SpiralWaves/test_derived_equivalence.jl.")
    println(io, "# generator : capabilities/AccretionDisc/derivation/src/codegen.jl")
    println(io, "#             (rerun: julia --project=capabilities/AccretionDisc/derivation \\")
    println(io, "#                          capabilities/AccretionDisc/derivation/scripts/emit_sources.jl)")
    println(io, "# input     : capabilities/AccretionDisc/derivation/out/$artifact  sha256 $(_sha(artifact_path))")
    println(io, "# derivation: sha256 $(derivation_sha()) over derivation/src/*.jl (sorted, concatenated)")
    println(io, "# julia     : $(VERSION)   Symbolics 7.36.0")
    for line in content
        println(io, "# $line")
    end
end

const _CONVENTIONS = [
    "conventions: Rationals are emitted as n//d and e^{a nu + b psi + c mu} as enu^a*epsi^b*emu^c (hoisted",
    "            into the preamble), so every float type T (Float64, BigFloat, ...) inherits its own",
    "            precision; there is no transcendental in a body (W = 1/sqrt(1-v^2) is an INPUT).",
    "            htr = hrphi = 0 (`circmetriczero`) is applied at emission: it is a property of the",
    "            paper's singular field (52)-(54) (h ∝ g(w) + 2 v v, v^r = 0, g_tr = g_rphi = 0), NOT a",
    "            theorem about circular orbits; the derivation keeps both components symbolic.",
]

"""
    emit_master(path, artifact_path) -> path

Write `master_derived.jl`: `_master_entries(inp::MasterInputs)` returning `(C00, C01, C10, C11, S0, S1)`
of the master system (47), with the SAME name, signature and 26 inputs as the route-R1
`master_generated.jl` it replaces.
"""
function emit_master(path::AbstractString, artifact_path::AbstractString)
    ex = Dict(read_artifact(artifact_path))
    entries = [n => to_masterinputs(ex[n]) for n in ("C00", "C01", "C10", "C11", "S0", "S1")]
    for (n, p) in entries
        for s in symbol_names(p)
            (s in MI_FIELDS || s in (:enu, :epsi, :emu)) ||
                error("codegen: $n needs `$s`, which MasterInputs does not carry")
        end
    end
    open(path, "w") do io
        _header(io, "master_derived.jl", "R3_master.jls.txt", artifact_path,
                vcat(["content   : dy/dr = C.y + S with y = (deltaU, urbyi) = (delta p-hat/(eps+p), -i delta u-hat^r);",
                      "            C = [C00 C01; C10 C11], S = [S0, S1]  (arXiv:2608.17003 eq 47).",
                      "            The R3 alphabet carries omega-tilde and v'; MasterInputs carries Omega and Omega'",
                      "            -- Omega = omt + v e^{nu-psi} (10) -- so omt and v_r are substituted out at emission.",
                      "            Entries are the ON-SHELL form (nu' eliminated by the hydrostatic balance (11)), which",
                      "            is what the R1/R2 routes also emit; off-shell the three routes differ by multiples of",
                      "            the (11) residual and must NOT be compared (see test_derived_equivalence.jl)."],
                     _CONVENTIONS))
        println(io)
        println(io, """
\"\"\"
    _master_entries(inp::MasterInputs) -> (C00, C01, C10, C11, S0, S1)

Entries of the relativistic master system dy/dr = C y + S (arXiv:2608.17003 eq 47), DERIVED in Julia
by route R3 (see the file header). Pure, allocation-free, type-generic: returns `Complex` entries whose
real type is the promotion of the `MasterInputs` field types (the Landau-shifted `omega_d` is Complex).
\"\"\"""")
        println(io, "@inline function _master_entries(inp::MasterInputs)")
        println(io, "    ", join(("$f = inp.$f" for f in MI_FIELDS), "; "))
        println(io, "    enu = exp(nu); epsi = exp(psi); emu = exp(mu)")
        for (n, p) in entries
            println(io, "    $n = ", print_value(p))
        end
        println(io, "    return (C00, C01, C10, C11, S0, S1)")
        println(io, "end")
    end
    path
end

"One emitted route-R2-signature diagnostics kernel."
struct DiagKernel
    fname::String
    export_name::String
    artifact::Symbol          # :master or :diagnostics
    doc::String
end

const DIAG_KERNELS = [
    DiagKernel("_r2_delta_u_phi", "duphi_low", :master, """
covariant azimuthal velocity perturbation delta-u-hat_phi = g_{phi mu} delta-u-hat^mu + h_{phi mu} u^mu
(the eliminated phi-velocity of the (47) reduction), in the master variables y = (y0, y1). Newtonian
limit: r * eq (50a) (ERRATUM E-M1: the printed (50a) is the ORTHONORMAL component).

DERIVED by route R3 (`duphi_low` of R3_master.jls.txt; see the file header). Pure, allocation-free,
type-generic."""),
    DiagKernel("_r2_amf", "amf_mode", :diagnostics, """
per-mode angular-momentum-flux integrand of (72a) via the on-shell identity (69):
AMF_m = u_phi delta-p-hat delta-u-hat^{r*} + (eps+p) delta-u-hat^r delta-u-hat_phi^*, so that
F_J(r) = sqrt(-g) * 4pi * sum_{m>=1} Re[AMF_m] with sqrt(-g) = e^{nu+psi+mu} (73).

DERIVED by route R3 (`amf_mode` of R3_diagnostics.jls.txt; see the file header). Pure,
allocation-free, type-generic."""),
    DiagKernel("_r2_torque_density", "torque_mode", :diagnostics, """
per-mode torque-density integrand of (43b)/(72b): 1/2 sum_{ab} delta_xi-T-hat^{ab} (i m h-hat_{ab})^*,
so that <T>(r) = 4pi * sum_{m>=1} Re[torque_m]; the integrand of the (72b) torque 𝔗(r).

DERIVED by route R3 (`torque_mode` of R3_diagnostics.jls.txt; see the file header). Pure,
allocation-free, type-generic."""),
]

"""
    emit_diagnostics(path, master_artifact, diag_artifact) -> path

Write `r2_derived.jl`: `_r2_delta_u_phi`, `_r2_amf`, `_r2_torque_density` with the SAME names and the
SAME `(inp::MasterInputs, omt, vp, y0, y1)` signature as the route-R2 `r2_generated.jl` they replace
(ω̃ and v′ are passed explicitly there, so no change of variables is needed).
"""
function emit_diagnostics(path::AbstractString, master_artifact::AbstractString,
                          diag_artifact::AbstractString)
    exm = Dict(read_artifact(master_artifact))
    exd = Dict(read_artifact(diag_artifact))
    bodies = [(k, circmetriczero(k.artifact === :master ? exm[k.export_name] : exd[k.export_name]))
              for k in DIAG_KERNELS]
    open(path, "w") do io
        _header(io, "r2_derived.jl", "R3_diagnostics.jls.txt", diag_artifact,
                vcat(["input2    : capabilities/AccretionDisc/derivation/out/R3_master.jls.txt  sha256 $(_sha(master_artifact))",
                      "content   : the quantities of the (47) reduction and of the diagnostics (69)-(73) in the master",
                      "            variables y = (y0, y1) = (delta-U-hat, -i delta-u-hat^r): the COVARIANT delta-u-hat_phi,",
                      "            the per-mode AMF integrand (72a) and the per-mode torque density (72b).",
                      "signature : (inp::MasterInputs, omt, vp, y0, y1) -- omt = omega-tilde(r) and vp = dv/dr are the two",
                      "            background scalars route R1 does not use; `r2_context` of master.jl assembles them.",
                      "            Conjugate amplitudes are bound in each preamble (the h-modes are real, conj is identity)."],
                     _CONVENTIONS))
        for (k, p) in bodies
            names = symbol_names(p)
            conjs = [(a, b) for (a, b) in CONJ_BINDINGS if a in names]
            # a conjugate binding `xc = conj(x)` needs its base field bound even when the expression
            # mentions only the barred symbol (the AMF bilinear does exactly that)
            bases = Set(Symbol(match(r"conj\((\w+)\)", b)[1]) for (a, b) in conjs)
            fields = [f for f in MI_FIELDS if f in names || f in bases]
            println(io)
            println(io, "\"\"\"")
            println(io, "    $(k.fname)(inp::MasterInputs, omt, vp, y0, y1)")
            println(io)
            println(io, k.doc)
            println(io, "\"\"\"")
            println(io, "@inline function $(k.fname)(inp::MasterInputs, omt, vp, y0, y1)")
            for f in fields
                println(io, "    $f = inp.$f")
            end
            :v_r in names && println(io, "    v_r = vp")
            for (a, b) in conjs
                println(io, "    $a = $b")
            end
            for (sym, fld) in ((:enu, :nu), (:epsi, :psi), (:emu, :mu))
                sym in names && println(io, "    $sym = exp(inp.$fld)")
            end
            println(io, "    return ", print_value(p))
            println(io, "end")
        end
    end
    path
end

# ---------------------------------------------------------------------------------------------
# round-trip self-check: the EMITTED TEXT, parsed back, evaluates exactly like the R3 expression
# ---------------------------------------------------------------------------------------------
#
# This is the admission gate of the printer itself and it is EXACT: the emitted Julia is re-read with
# `Meta.parse`, interpreted over `Rational{BigInt}` (no compilation, no floats), and compared with
# `exact_eval` of the Symbolics expression at the same random rational point.  Residual must be
# EXACTLY zero.  It certifies, in one shot: the artefact round-trip, `circmetriczero`, the
# ω̃, v′ → Ω, Ω′ change of variables and the printer.

"Every `Expr(:function, …)` node of a parsed file, innermost-first."
function _find_functions(ex, acc = Any[])
    if ex isa Expr
        ex.head === :function && push!(acc, ex)
        for a in ex.args
            _find_functions(a, acc)
        end
    end
    acc
end

_fname(sig) = sig isa Expr && sig.head === :call ? String(sig.args[1]) : error("codegen: odd signature $sig")

"Collect `name = <expr>` assignments (LHS in `keep`) and the `return` expression of a function body."
function _body_exprs(body, keep::Set{Symbol}, acc)
    body isa Expr || return acc
    if body.head === :(=) && body.args[1] isa Symbol && body.args[1] in keep
        push!(acc, String(body.args[1]) => body.args[2])
    elseif body.head === :return
        push!(acc, "return" => body.args[1])
    else
        for a in body.args
            _body_exprs(a, keep, acc)
        end
    end
    acc
end

"""
    emitted_expressions(path) -> Vector{Pair{String,Any}}

`"<function>/<entry>" => <parsed Julia expression>` for every emitted body of a generated file
(`entry` = the assigned name, or `return`).  Read back from the FILE, not from memory.
"""
function emitted_expressions(path::AbstractString)
    keep = Set([:C00, :C01, :C10, :C11, :S0, :S1])
    out = Pair{String,Any}[]
    for f in _find_functions(Meta.parseall(read(path, String)))
        n = _fname(f.args[1])
        got = _body_exprs(f.args[2], keep, Pair{String,Any}[])
        # a function that assigns named entries returns the TUPLE of them — not a body of its own
        any(p -> p.first != "return", got) && filter!(p -> p.first != "return", got)
        for (k, e) in got
            push!(out, string(n, "/", k) => e)
        end
    end
    out
end

"Interpret a parsed arithmetic tree over exact numbers (no `eval`, no compilation, no floats)."
function eval_tree(ex, tbl::Dict{Symbol,Any})
    ex isa Number && return ex
    ex isa Symbol && return (haskey(tbl, ex) ? tbl[ex] : error("eval_tree: unbound `$ex`"))
    ex isa Expr || error("eval_tree: cannot interpret $(typeof(ex))")
    ex.head === :call || error("eval_tree: unexpected head $(ex.head)")
    op = ex.args[1]
    as = [eval_tree(a, tbl) for a in ex.args[2:end]]
    op === :+ && return sum(as)
    op === :* && return prod(as)
    op === :- && return length(as) == 1 ? -as[1] : as[1] - as[2]
    op === :/ && return as[1] / as[2]
    op === :^ && return as[1]^as[2]
    op === :// && return as[1] // as[2]
    error("eval_tree: unsupported operator $op")
end

"""
    codegen_selfcheck(master_path, r2_path, master_artifact, diag_artifact; npts=20, seed=2608)

Exact round-trip certificate of the emitted sources.  For `npts` random exact-rational points:
re-read the emitted Julia, interpret it over `Rational{BigInt}`, and demand EXACT equality with
`exact_eval` of the corresponding route-R3 expression (`circmetriczero` applied, Ω and Ω′ computed
from ω̃ and v′ by (10)).
"""
function codegen_selfcheck(master_path, r2_path, master_artifact, diag_artifact;
                           npts::Int = 20, seed::Int = 2608, verbose::Bool = true)
    exm = Dict(read_artifact(master_artifact))
    exd = Dict(read_artifact(diag_artifact))
    refs = Dict{String,Tuple{Num,Num}}()
    for n in ("C00", "C01", "C10", "C11", "S0", "S1")
        refs["_master_entries/$n"] = to_masterinputs(exm[n])
    end
    for k in DIAG_KERNELS
        refs["$(k.fname)/return"] =
            circmetriczero(k.artifact === :master ? exm[k.export_name] : exd[k.export_name])
    end
    emitted = vcat(emitted_expressions(master_path), emitted_expressions(r2_path))
    length(emitted) == length(refs) ||
        error("codegen: emitted $(length(emitted)) bodies, expected $(length(refs))")

    rng = MersenneTwister(seed)
    bad = 0
    firstbad = ""
    for k in 1:npts
        pt = rational_point(rng)
        for s in keys(CIRC_METRIC_ZERO)
            pt[s] = RatB(0)
        end
        pt[Omega] = pt[omt] + pt[v] * pt[enu] / pt[epsi]
        pt[Omega_r] = pt[omt_r] + pt[enu] * (pt[v_r] + pt[v] * (pt[nu_r] - pt[psi_r])) / pt[epsi]
        tbl = Dict{Symbol,Any}(:im => im)
        for (s, val) in pt
            tbl[Symbol(string(Symbolics.value(s)))] = val
        end
        for (name, ex) in emitted
            haskey(refs, name) || error("codegen: emitted body `$name` has no reference expression")
            got = eval_tree(ex, tbl)
            r = refs[name]
            want = exact_eval(r[1], pt) + im * exact_eval(r[2], pt)
            got isa AbstractFloat && error("codegen: round-trip produced a FLOAT")
            if !iszero(got - want)
                bad += 1
                isempty(firstbad) && (firstbad = "point $k, body $name: residual = $(got - want)")
            end
        end
    end
    c = Certificate("emitted sources round-trip to the R3 expressions", "(47),(72)",
                    :exact_rational, npts, bad == 0,
                    bad == 0 ? "$(length(emitted)) emitted bodies, exact over Rational{BigInt}" :
                               "$bad exact-zero failures; first: $firstbad")
    verbose && @printf("[%s] %-46s %-16s %s\n", c.pass ? "PASS" : "FAIL", c.name,
                       "exact-rational", c.note)
    c
end
