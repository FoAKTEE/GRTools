"""
    PostNewtonian.MMAInputForm

The unified MMA-InputForm→`Symbolics.Num` symbolic walker — the symbolic side of the PN catalog
(DISTINCT from `Data.read_mma_inputform`, which reads the NUMERIC `.wl` dialect; design §10 format
separation). Rewrites a stored `.m` `Association`/`SeriesData`/`Module` text into a Julia `Expr` with
`Meta.parse`, then walks that tree (NO `eval`, NO world-age) into a `Dict`/`Vector`/[`PNSeries`](@ref)
value. Fuses the two legacy walker forks (PostNewtonianSelfForce `_eval_pn_ast` + TeukolskyPN
`_eval_mma_ast`) behind one entry.

Alphabet covers the circular/spherical/generic store: series variables `p`/`y`/`v`, spin `a`/`q`,
inclination `x`/`Y`, eccentricity `e`, `Infinity` expansion point, `Pi`/`EulerGamma`/`Logy`/`Logp`/
`Logv` transcendentals, opaque `Log[k]`/`Zeta[k]`/`Psi[..]`, the imaginary unit (`I` and lowercase
`i`), the complex digamma `ExtraFunctions` (`PolyGamma`/`RePsi`/`ImPsi`/`PsiA`/`PsiB`), the `Psibar`
formal quantity, and scalar `HoldForm` unwrap. Giant (~35k-term) coefficient sums materialize via
balanced reduction (`_accumulate_balanced`) over an iteratively-flattened `±` chain, so the
Spherical flux giants no longer overflow. The full ECCENTRIC `e`-graded store deserializes with ZERO
unparsed heads: InputForm implicit multiplication (`(341/16) e^2`, `240 (1-e^2)^(3/2)(1+…)`), the
`Rational[p,q]` head, nested `SeriesData[e,…]` coefficients, `ResummedSeriesData[sf,SeriesData[e,…]]`,
the held generating-function heads (`HoldForm[oneP5PN[e,…]]`/`chi3PN`/`chi4PNConv`/… → a FORMAL
Symbolics atom, since their `:=` programs live in the brace-balanced-dropped `"ExtraFunctions"`
block), the Kerr `TagBox[RowBox[…],HoldForm]` box-language, MMA numeric forms (`mantissa`prec`, `*^`
scientific), the `O[e]^n` truncation term, and the `Dot` (`.`) scalar product.
"""
module MMAInputForm

using ..Series: PNSeries, Normal
import Symbolics
import SpecialFunctions

export read_mma_inputform, parse_pn_data_file

# The imaginary unit symbol (shared convention with MST.imagunit: `𝕚`). MMA InputForm writes `I`;
# kept as a symbol so a coefficient stays a single `Num` (substitute `=> im` to evaluate).
const _IMAG = Symbolics.variable(:𝕚)

# Symbol alphabet for the circular/spherical store. Built once so repeated parses reuse the same
# Symbolics variables (substitution-compatible). `Log[k]`/`Zeta[k]` are minted on demand as opaque
# atoms.
const _PN_SYMS = Dict{Symbol,Symbolics.Num}(
    :a => Symbolics.variable(:a),
    :q => Symbolics.variable(:q),   # alternate Kerr-spin symbol (Kerr/Circular/Local store)
    :x => Symbolics.variable(:x),
    :p => Symbolics.variable(:p),
    :Infinity => Symbolics.variable(:Infinity),
    :Pi => Symbolics.variable(:Pi),
    :EulerGamma => Symbolics.variable(:EulerGamma),
    :Logy => Symbolics.variable(:Logy),
    :Logp => Symbolics.variable(:Logp),   # Log[p] expansion-parameter atom (Spherical-flux store)
    :Logv => Symbolics.variable(:Logv),   # Log[v] expansion-parameter atom (Generic-flux store)
    :I => _IMAG,                          # MMA imaginary unit
    :i => _IMAG,                          # lowercase imaginary unit (Kerr/Circular/Local/Spin store)
    :e => Symbolics.variable(:e),         # eccentricity formal PN symbol (Generic/Eccentric store)
    :Y => Symbolics.variable(:Y),         # inclination cos-θ formal PN symbol (Generic-flux store)
    :v => Symbolics.variable(:v),         # PN velocity formal symbol (Generic/Kerr-Eccentric store)
    :vv => Symbolics.variable(:vv),       # Kerr Eccentric Redshift log-argument formal symbol (Log[vv])
    :ka => Symbolics.variable(:ka),       # Kerr Eccentric Redshift log-argument formal symbol (Log[ka])
    # MMA named-character free parameters (`\[Sigma]` etc., rewritten to bare identifiers).
    :Sigma => Symbolics.variable(:Sigma),
    :CapitalDelta => Symbolics.variable(:CapitalDelta),
    :Kappa => Symbolics.variable(:Kappa),
    :Lambda => Symbolics.variable(:Lambda),
    :Psi => Symbolics.variable(:Psi),
)

"""
    read_mma_inputform(src::AbstractString) -> Any

Deserialize one MMA InputForm text fragment (`Association`/`SeriesData`/`Module`/list/number) into a
Julia value: numbers/rationals → `Num`, strings → `String`, `{...}` → `Vector`, `<|...|>` → `Dict`,
`SeriesData[...]` → [`PNSeries`](@ref). No `eval`/world-age. Unknown heads or out-of-alphabet atoms
throw (an unexpected token surfaces loudly, never silently mis-parsed).
"""
function read_mma_inputform(src::AbstractString)
    jl = _mma_data_to_julia(String(src))
    ast = Meta.parse(jl)
    return _eval_pn_ast(ast, _PN_SYMS, Dict{Symbol,Any}())
end

"""
    parse_pn_data_file(path) -> Dict{String,Any}

Read and deserialize a stored PN `.m` data file (a literal MMA Association) into a
`Dict{String,Any}` with `"Series"::PNSeries`. The `SelfForce` layer wraps the result in
`PostNewtonianData`. Upstream analogue: `Get[file]`.
"""
function parse_pn_data_file(path::AbstractString)
    val = read_mma_inputform(read(path, String))
    val isa Dict{String,Any} || throw(ArgumentError(
        "PN data file did not deserialize to an Association; got $(typeof(val))"))
    return val
end

# Rewrite MMA InputForm tokens into Julia-parseable syntax:
#   (* ... *) -> removed ; <| |> -> Dict( ) ; -> => ; { } -> [ ]
# `Module[{var}, var = SeriesData[...]; <|...|>]` wrappers (86 of 91 files) become a Julia
# `begin ... end` block (the local-var declaration dropped; the body's `;`-statements become block
# statements, the final Association the block's value).
function _mma_data_to_julia(s::String)
    s = replace(s, r"\(\*.*?\*\)"s => "")              # strip MMA comments
    s = replace(s, r"\\\[([A-Za-z]+)\]" => s"\1")      # \[Name] -> Name
    # Drop the `"ExtraFunctions" -> { ... }` DEFINITION entry: a held pattern-rule/program block
    # (`_` patterns, `:>`/`:=`, `D[..]`, `Sum[..,{n,1,..}]`, `Series[..,{e,0,..}]`) that is not data
    # and not Julia-parseable. Its RePsi/ImPsi digamma semantics are reproduced by the walker atoms;
    # its eccentric gen-func programs are kept FORMAL (the held head is minted as a Symbolics atom).
    # The eccentric blocks carry NESTED `{n,1,..}`/`{e,0,..}` braces, so a brace-BALANCED scan is
    # required (the old `[^{}]*` regex stopped at the first inner `}` and left a dangling tail).
    s = _drop_extrafunctions(s)
    s = _unwrap_tagbox(s)
    # Protect Mathematica numeric tokens that Julia cannot represent directly.  Explicit
    # precision/accuracy marks and `*^` scientific exponents become a private data head whose
    # evaluator preserves BigFloat metadata (or machine Float64 for an unmarked decimal).
    s = _rewrite_mma_numbers(s)
    # Insert MMA InputForm implicit multiplication (juxtaposition = Times) as explicit `*`, e.g.
    # `(341/16) e^2`, `5 a (1-e^2)^2`, `240 (1-e^2)^(3/2)(1+..)`, `Pi {list}` — Meta.parse rejects the
    # bare juxtaposition. Runs while `{}`=list / `[]`=call / `()`=group are still MMA-shaped so a `{`
    # is a value-start (a list factor) but a call bracket `Head[` is not. Strings are opaque.
    s = _mma_implicit_mult(s)
    m = match(r"^\s*Module\[\{[^}]*\},(.*)\]\s*$"s, s)
    if m !== nothing
        s = "begin " * m.captures[1] * " end"
    end
    s = replace(s, "<|" => "Dict(", "|>" => ")")
    s = replace(s, "->" => "=>")
    s = replace(s, "{" => "[", "}" => "]")
    return s
end

# Rewrite Mathematica numeric tokens carrying a backtick precision/accuracy mark or a `*^`
# scientific exponent to `MMAReal["verbatim-token"]`.  The scan is string-aware so number-looking
# text in Association values is untouched.  Ordinary decimal tokens are copied byte-for-byte and
# retain the historical Meta.parse -> Float64 behavior.
function _rewrite_mma_numbers(s::String)
    out = IOBuffer()
    n = lastindex(s)
    i = firstindex(s)
    while i <= n
        c = s[i]
        if c == '"'
            j = nextind(s, i)
            write(out, c)
            while j <= n
                cj = s[j]
                write(out, cj)
                if cj == '\\' && j < n
                    j = nextind(s, j)
                    write(out, s[j])
                elseif cj == '"'
                    break
                end
                j = nextind(s, j)
            end
            i = j <= n ? nextind(s, j) : nextind(s, n)
            continue
        end

        nexti = nextind(s, i)
        number_start = isdigit(c) || (c == '.' && nexti <= n && isdigit(s[nexti]))
        if !number_start
            write(out, c)
            i = nexti
            continue
        end

        start = i
        while i <= n && isdigit(s[i])
            i = nextind(s, i)
        end
        if i <= n && s[i] == '.'
            i = nextind(s, i)
            while i <= n && isdigit(s[i])
                i = nextind(s, i)
            end
        end

        marked = false
        double_tick = false
        if i <= n && s[i] == '`'
            marked = true
            i = nextind(s, i)
            if i <= n && s[i] == '`'
                double_tick = true
                i = nextind(s, i)
            end
            if i <= n && (s[i] == '+' || s[i] == '-')
                i = nextind(s, i)
            end
            while i <= n && isdigit(s[i])
                i = nextind(s, i)
            end
            if i <= n && s[i] == '.'
                i = nextind(s, i)
                while i <= n && isdigit(s[i])
                    i = nextind(s, i)
                end
            end
        end

        scientific = false
        if i <= n && s[i] == '*' && nextind(s, i) <= n && s[nextind(s, i)] == '^'
            scientific = true
            i = nextind(s, nextind(s, i))
            if i <= n && (s[i] == '+' || s[i] == '-')
                i = nextind(s, i)
            end
            expstart = i
            while i <= n && isdigit(s[i])
                i = nextind(s, i)
            end
            expstart == i && throw(ArgumentError(
                "Mathematica scientific literal needs an integer exponent near index $start"))
        end

        stop = prevind(s, i)
        token = s[start:stop]
        if marked || scientific
            double_tick && !occursin(r"``[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)", token) &&
                throw(ArgumentError("Mathematica accuracy literal is missing its accuracy: $token"))
            write(out, "MMAReal[")
            print(out, repr(token))
            write(out, ']')
        else
            write(out, token)
        end
    end
    return String(take!(out))
end

# Brace-balanced removal of the `"ExtraFunctions" -> { ... }` entry (and its trailing comma). Scans
# from the first `{` after the key to its MATCHING `}` (depth counter), so nested `{n,1,..}`/`{e,0,..}`
# in the eccentric held programs are spanned. No-op when the key is absent.
function _drop_extrafunctions(s::String)
    i = findfirst("\"ExtraFunctions\"", s)
    i === nothing && return s
    istart = first(i)
    j = findnext('{', s, nextind(s, last(i)))
    j === nothing && return s
    depth = 0
    stop = 0
    k = j
    while k <= lastindex(s)
        c = s[k]
        if c == '{'
            depth += 1
        elseif c == '}'
            depth -= 1
            if depth == 0
                stop = k
                break
            end
        end
        k = nextind(s, k)
    end
    stop == 0 && return s
    e = nextind(s, stop)
    while e <= lastindex(s) && (s[e] == ',' || isspace(s[e]))
        e = nextind(s, e)
    end
    return s[1:prevind(s, istart)] * (e <= lastindex(s) ? s[e:end] : "")
end

# Rewrite the MMA StandardForm box-language `\!\(\* TagBox[RowBox[{"head","[",RowBox[{a1,",",a2}],
# "]"}], HoldForm] \)` (how the Kerr eccentric Redshift store typesets a held gen-func call) into the
# plain InputForm `HoldForm[head[a1,a2]]`, so the held head flows through the same formal-atom path as
# the InputForm `HoldForm[oneP5PN[…]]` occurrences. The RowBox string tokens carry the head and its
# comma-separated arguments; the interleaved `","` separator tokens are dropped.
function _unwrap_tagbox(s::String)
    box = r"\\!\\\(\\\*\s*TagBox\[\s*RowBox\[\{\s*\"(\w+)\"\s*,\s*\"\[\"\s*,\s*RowBox\[\{([^}]*)\}\]\s*,\s*\"\]\"\s*\}\]\s*,\s*HoldForm\]\\\)"s
    return replace(s, box => function (m)
        cap = match(box, m)
        head = cap.captures[1]
        argtoks = [x.captures[1] for x in eachmatch(r"\"(\w+)\"", cap.captures[2])]
        return string("HoldForm[", head, "[", join(argtoks, ","), "]]")
    end)
end

# Insert MMA InputForm implicit multiplication (juxtaposition = Times) as explicit `*`. A single
# character scan classifies each token as a product LEFT factor (`prev_end`: a number, identifier,
# `)`, `]`, `}`) or a RIGHT factor start (a number, identifier, `(`, `{`); `*` is emitted between an
# end and a start. Whitespace is dropped (Julia self-delimits); string literals are copied verbatim
# (escapes honoured) and never juxtapose. Every operator/bracket char is passed through unchanged, so
# `Head[` stays a call and `,`/`+`/`-`/`^`/`/` suppress the insertion.
function _mma_implicit_mult(s::String)
    out = IOBuffer()
    n = lastindex(s)
    i = firstindex(s)
    prev_end = false
    while i <= n
        c = s[i]
        if c == '"'
            j = nextind(s, i)
            write(out, '"')
            while j <= n
                cj = s[j]
                write(out, cj)
                if cj == '\\' && j < n
                    j = nextind(s, j)
                    write(out, s[j])
                elseif cj == '"'
                    break
                end
                j = nextind(s, j)
            end
            i = j <= n ? nextind(s, j) : nextind(s, n)
            prev_end = false
        elseif isdigit(c)
            prev_end && write(out, '*')
            while i <= n && isdigit(s[i])
                write(out, s[i])
                i = nextind(s, i)
            end
            prev_end = true
        elseif c == '.'
            # A `.` between two DIGITS is a decimal point (kept in the number); otherwise it is MMA
            # `Dot`, which on the scalar resummed series here is Times → emit `*`.
            prevc = i > firstindex(s) ? s[prevind(s, i)] : ' '
            nextc = nextind(s, i) <= n ? s[nextind(s, i)] : ' '
            if isdigit(prevc) && isdigit(nextc)
                write(out, '.')
            else
                write(out, '*')
            end
            i = nextind(s, i)
            prev_end = false
        elseif c == '*' && nextind(s, i) <= n && s[nextind(s, i)] == '^'
            # MMA scientific notation `mantissa*^exp` = mantissa·10^exp → `mantissa*10^exp`.
            write(out, "*10^")
            i = nextind(s, nextind(s, i))
            prev_end = false
        elseif c == '`'
            # `_rewrite_mma_numbers` consumes every valid number mark before this pass.
            throw(ArgumentError("unparsed Mathematica number mark near index $i"))
        elseif isletter(c)
            prev_end && write(out, '*')
            while i <= n && (isletter(s[i]) || isdigit(s[i]))
                write(out, s[i])
                i = nextind(s, i)
            end
            prev_end = true
        elseif c == '(' || c == '{'
            prev_end && write(out, '*')
            write(out, c)
            i = nextind(s, i)
            prev_end = false
        elseif c == ')' || c == ']' || c == '}'
            write(out, c)
            i = nextind(s, i)
            prev_end = true
        elseif isspace(c)
            i = nextind(s, i)
        else
            write(out, c)
            i = nextind(s, i)
            prev_end = false
        end
    end
    return String(take!(out))
end

# Walk a parsed Julia Expr (from the rewritten MMA data) into a Julia value.
function _eval_pn_ast(node, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    if node isa Bool
        throw(ArgumentError("_eval_pn_ast: unexpected boolean $node"))
    elseif node isa Integer
        return Symbolics.Num(node)
    elseif node isa AbstractFloat
        return Symbolics.Num(node)
    elseif node isa Rational
        return Symbolics.Num(node)
    elseif node isa AbstractString
        return String(node)
    elseif node isa LineNumberNode
        return nothing
    elseif node isa Symbol
        haskey(locals, node) && return locals[node]
        node === :Infinity && return env[:Infinity]
        haskey(env, node) || throw(ArgumentError("_eval_pn_ast: unbound symbol $node"))
        return env[node]
    elseif node isa Expr
        if node.head === :block
            val = nothing
            for st in node.args
                v = _eval_pn_ast(st, env, locals)
                v === nothing || (val = v)
            end
            return val
        elseif node.head === :(=)
            sym = node.args[1]
            sym isa Symbol || throw(ArgumentError("_eval_pn_ast: assignment LHS not a symbol: $(repr(sym))"))
            locals[sym] = _eval_pn_ast(node.args[2], env, locals)
            return nothing
        elseif node.head === :call
            op = node.args[1]
            if op === :Dict
                d = Dict{String,Any}()
                for pr in node.args[2:end]
                    k, v = _eval_pn_pair(pr, env, locals)
                    d[k] = v
                end
                return d
            elseif op === :(=>)
                return Pair(_eval_pn_ast(node.args[2], env, locals), _eval_pn_ast(node.args[3], env, locals))
            elseif op === :+ || op === :-
                # `+`/`-` chains: the Spherical-giant coefficients are ~35k-term signed sums that
                # Meta.parse renders as a LEFT-NESTED binary `±` tree (mixed `+`/`-` breaks Julia's
                # flat-`+` folding). Walking that tree with plain recursion overflows the stack, so
                # flatten it ITERATIVELY into signed leaves and rebuild by BALANCED reduction — see
                # `_eval_addsub_chain`. Small/medium chains take the same path (identical value).
                return _eval_addsub_chain(node, env, locals)
            else
                as = [_eval_pn_ast(a, env, locals) for a in node.args[2:end]]
                op === :* && return prod(as)
                op === :/ && return _pn_div(as[1], as[2])
                op === :^ && return _pn_pow(as[1], as[2])
                throw(ArgumentError("_eval_pn_ast: unsupported call op $op"))
            end
        elseif node.head === :macrocall
            # Oversized integer literals: Meta.parse emits `@big_str`/`@int128_str` rather than a
            # folded Integer. Recover the exact value from the string arg.
            mac = node.args[1]
            macname = mac isa GlobalRef ? mac.name : mac
            endswith(string(macname), "_str") ||
                throw(ArgumentError("_eval_pn_ast: unsupported macrocall $mac"))
            str = node.args[end]
            str isa AbstractString || throw(ArgumentError("_eval_pn_ast: macrocall arg not a string: $(repr(str))"))
            return Symbolics.Num(parse(BigInt, str))
        elseif node.head === :vect
            return Any[_eval_pn_ast(a, env, locals) for a in node.args]
        elseif node.head === :ref
            head = node.args[1]
            head === :SeriesData && return _build_pnseries(node.args[2:end], env, locals)
            head === :MMAReal && return _atom_mma_real(node.args[2:end])
            head === :Rational && return _atom_rational(node.args[2:end], env, locals)
            head === :ResummedSeriesData && return _atom_resummed(node.args[2:end], env, locals)
            head === :O && return Symbolics.Num(0)   # MMA `O[e]^n` truncation term → 0 on materialize
            head === :Sqrt && return sqrt(_eval_pn_ast(node.args[2], env, locals))
            head === :Log && return _parse_log(node.args[2:end], env, locals)
            head === :PolyGamma && return _atom_polygamma(node.args[2:end], env, locals)
            head === :RePsi && return _atom_repsi(node.args[2:end], env, locals; imag = false)
            head === :ImPsi && return _atom_repsi(node.args[2:end], env, locals; imag = true)
            head === :PsiA && return _atom_psi_combo(node.args[2:end], env, locals; imag = false)
            head === :PsiB && return _atom_psi_combo(node.args[2:end], env, locals; imag = true)
            head === :Psibar && return _atom_psibar(node.args[2:end], env, locals)
            head === :HoldForm && return _unwrap_holdform(node.args[2:end], env, locals)
            head === :Re && return _atom_reim(real, node.args[2:end], env, locals)
            head === :Im && return _atom_reim(imag, node.args[2:end], env, locals)
            (head === :Zeta || head === :Psi) &&
                return _opaque_atom(head, node.args[2:end], env, locals)
            throw(ArgumentError("_eval_pn_ast: unsupported bracketed head $head"))
        else
            throw(ArgumentError("_eval_pn_ast: unexpected expr head $(node.head)"))
        end
    else
        throw(ArgumentError("_eval_pn_ast: unexpected node $(repr(node))"))
    end
end

# Parse the private `MMAReal["..."]` data head introduced by `_rewrite_mma_numbers`.
# Mathematica `number`p` specifies p decimal digits of precision, while `number``a` specifies a
# decimal accuracy a, so the corresponding precision is a + log10(abs(number)).  MPFR/BigFloat
# stores binary precision but not a separate accuracy tag; choosing at least the converted bit count
# preserves the stated uncertainty and the numeric value.  A bare single tick and an unmarked
# decimal `*^` token remain machine Float64, matching Mathematica and the parser's prior behavior.
const _MMA_REAL_TOKEN =
    r"^([0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:(`{1,2})([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+))?)?(?:\*\^([+-]?[0-9]+))?$"

function _atom_mma_real(args)
    length(args) == 1 || throw(ArgumentError("MMAReal expects one literal, got $(length(args))"))
    token = args[1]
    token isa AbstractString || throw(ArgumentError("MMAReal literal must be a string"))
    m = match(_MMA_REAL_TOKEN, String(token))
    m === nothing && throw(ArgumentError("unsupported Mathematica real literal: $token"))
    mantissa, ticks, digits, exponent = m.captures
    normalized = exponent === nothing ? mantissa : string(mantissa, 'e', exponent)

    if ticks === nothing
        if occursin('.', mantissa)
            return Symbolics.Num(parse(Float64, normalized))
        end
        e = try
            parse(Int, exponent)
        catch
            throw(ArgumentError("Mathematica scientific exponent is outside the supported Int range: $token"))
        end
        base = parse(BigInt, mantissa)
        return Symbolics.Num(e >= 0 ? base * big(10)^e : base // big(10)^(-e))
    elseif ticks == "`" && digits === nothing
        return Symbolics.Num(parse(Float64, normalized))
    elseif digits === nothing
        throw(ArgumentError("Mathematica precision/accuracy literal is missing its digit count: $token"))
    end

    stated = try
        parse(Float64, digits)
    catch
        throw(ArgumentError("invalid Mathematica precision/accuracy in literal: $token"))
    end
    isfinite(stated) || throw(ArgumentError("non-finite Mathematica precision/accuracy in literal: $token"))

    decimal_precision = if ticks == "`"
        stated
    else
        rough = setprecision(BigFloat, max(precision(BigFloat), 256)) do
            parse(BigFloat, normalized)
        end
        iszero(rough) ? stated : stated + Float64(log10(abs(rough)))
    end
    rawbits = (max(decimal_precision, 0.0) * log2(10.0)) + 2
    isfinite(rawbits) && rawbits <= typemax(Int) ||
        throw(ArgumentError("Mathematica precision is outside the supported BigFloat range: $token"))
    bits = max(2, ceil(Int, rawbits))
    value = setprecision(BigFloat, bits) do
        parse(BigFloat, normalized)
    end
    return Symbolics.Num(value)
end

# Flatten a nested `+`/`-` `Expr` chain into `(sign, leaf-expr)` pairs WITHOUT recursion (an explicit
# worklist), so a 35k-deep left-nested `±` tree — the shape Meta.parse gives a giant mixed-sign
# coefficient — never blows the call stack. A leaf is any node that is not itself a `+`/`-` call
# (numbers, atoms, products, powers, nested parenthesised sub-sums handled later by `_eval_pn_ast`).
function _flatten_addsub(root)
    terms = Tuple{Int,Any}[]
    work = Tuple{Int,Any}[(1, root)]
    while !isempty(work)
        sgn, node = pop!(work)
        if node isa Expr && node.head === :call && (node.args[1] === :+ || node.args[1] === :-)
            a = node.args
            if node.args[1] === :+
                for i in 2:length(a)
                    push!(work, (sgn, a[i]))
                end
            elseif length(a) == 2                     # unary minus
                push!(work, (-sgn, a[2]))
            else                                       # binary/n-ary minus: first `+`, rest `-`
                push!(work, (sgn, a[2]))
                for i in 3:length(a)
                    push!(work, (-sgn, a[i]))
                end
            end
        else
            push!(terms, (sgn, node))
        end
    end
    return terms
end

# Fold a length-`N` `Num` vector into one `Num` by pairwise divide-and-conquer instead of a
# left-associated sum, bounding the built expression-tree depth to `log2(N)`. Sum is associative and
# commutative, so the value is invariant to this reordering (exact on the rational/symbolic leaves);
# the reordering only avoids the deep nesting that overflows Symbolics' recursive `+` builder.
function _accumulate_balanced(v::AbstractVector{Symbolics.Num})
    n = length(v)
    n == 0 && return Symbolics.Num(0)
    n == 1 && return @inbounds v[1]
    mid = (n + 1) >> 1
    return _accumulate_balanced(@view v[1:mid]) + _accumulate_balanced(@view v[mid+1:n])
end

# Evaluate a `+`/`-` chain: flatten to signed leaves (iterative), walk each leaf (shallow), then
# balanced-reduce. Each leaf walk is bounded-depth, so the only unbounded dimension (chain length)
# is handled without recursion.
function _eval_addsub_chain(root, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    signed = _flatten_addsub(root)
    terms = Vector{Symbolics.Num}(undef, length(signed))
    for i in eachindex(signed)
        sgn, leaf = signed[i]
        # A leaf may be a nested `SeriesData[e,…]` (an eccentricity sub-series added to a poly-in-`e`
        # coefficient) — materialize it to its `Num` before combining. `_coeff_to_num` is the identity
        # on ordinary `Num` leaves.
        v = _coeff_to_num(_eval_pn_ast(leaf, env, locals))
        terms[i] = sgn == 1 ? v : -v
    end
    return _accumulate_balanced(terms)
end

# `Log[arg]`: `Log[k]` for a concrete integer `k` is an opaque constant atom (`Log2`); `Log[expr]` of
# a symbolic argument is a genuine logarithm.
function _parse_log(args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    length(args) == 1 || throw(ArgumentError("_parse_log: Log expects 1 arg, got $(length(args))"))
    v = _eval_pn_ast(args[1], env, locals)
    vv = Symbolics.value(v)
    (vv isa Integer) && return Symbolics.variable(Symbol("Log", Int(vv)))
    return log(v)
end

# `PolyGamma[m, z]` -> SpecialFunctions.polygamma(m, z) as a Symbolics term. The order `m` and the
# (possibly complex, symbolic-in-`a`) argument `z` are walked recursively; the term folds to a number
# once `a` and `𝕚` are substituted to concrete values (`polygamma(0,1)=-γ`, `polygamma(1,1)=π²/6`).
function _atom_polygamma(args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    length(args) == 2 || throw(ArgumentError("_atom_polygamma: PolyGamma expects 2 args, got $(length(args))"))
    m = _eval_pn_ast(args[1], env, locals)
    z = _eval_pn_ast(args[2], env, locals)
    return Symbolics.Num(Symbolics.term(SpecialFunctions.polygamma, m, z))
end

# The digamma-family real/imaginary combination at the complex Kerr argument
# `base ± 𝕚 a·mult/√(1-a²)`, the shape shared by every stored `"ExtraFunctions"` rule:
#   sum-branch  = ψ^(order)(base + arg) + ψ^(order)(base - arg)   (real part, 𝕚→im)
#   diff-branch = (ψ^(order)(base + arg) - ψ^(order)(base - arg)) / 𝕚   (imaginary part)
# `half` divides the branch by 2 (the `PsiA`/`PsiB` normalisation); `RePsi`/`ImPsi` use the raw
# (un-halved) branch. `a` is the Kerr-spin symbol from the alphabet.
function _psi_digamma_combo(order, mult, env::Dict{Symbol,Symbolics.Num}; base, imag::Bool, half::Bool)
    a = env[:a]
    arg = _IMAG * a * mult / sqrt(1 - a^2)
    pg(sgn) = Symbolics.Num(Symbolics.term(SpecialFunctions.polygamma, order, base + sgn * arg))
    combo = imag ? (pg(+1) - pg(-1)) / _IMAG : pg(+1) + pg(-1)
    return half ? combo / 2 : combo
end

# `RePsi[m,n]` / `ImPsi[m,n]` — digamma parts at base-point 1 (Kavanagh–Ottewill–Wardell
# convention), the un-halved sum/difference of `ψ^(m)(1 ± 𝕚 a n/√(1-a²))`.
function _atom_repsi(args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any}; imag::Bool)
    length(args) == 2 || throw(ArgumentError("_atom_repsi: RePsi/ImPsi expect 2 args, got $(length(args))"))
    m = _eval_pn_ast(args[1], env, locals)
    n = _eval_pn_ast(args[2], env, locals)
    return _psi_digamma_combo(m, n, env; base = 1, imag = imag, half = false)
end

# `PsiA[n,m]` / `PsiB[n,m]` — the `"ExtraFunctions"` digamma combinations at base-point 3, exactly
# as the stored def gives them (see Kerr/Circular/Flux/Energy/Infinity.m):
#   PsiA[n,m] = (ψ^(n)(3 + m a 𝕚/√(1-a²)) + ψ^(n)(3 - m a 𝕚/√(1-a²))) / 2
#   PsiB[n,m] = (ψ^(n)(3 + m a 𝕚/√(1-a²)) - ψ^(n)(3 - m a 𝕚/√(1-a²))) / 2 / 𝕚
# The `.m` def entry itself is dropped by the tokenizer; the semantics are reproduced here.
function _atom_psi_combo(args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any}; imag::Bool)
    length(args) == 2 || throw(ArgumentError("_atom_psi_combo: PsiA/PsiB expect 2 args, got $(length(args))"))
    n = _eval_pn_ast(args[1], env, locals)
    m = _eval_pn_ast(args[2], env, locals)
    return _psi_digamma_combo(n, m, env; base = 3, imag = imag, half = true)
end

# `Psibar[n,m]` (`\[Psi]bar[n,m]`) — a named formal PN quantity (no in-file definition; its symbol
# identity is what the RedShift series needs, per the round-note). Minted as `Psibar_n_m`.
function _atom_psibar(args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    ks = Int[_to_int(_eval_pn_ast(a, env, locals)) for a in args]
    return Symbolics.variable(Symbol("Psibar_" * join(ks, "_")))
end

# Bracketed heads the walker deserializes to a value; anything else wrapped in `HoldForm[·[…]]` is a
# held eccentric GENERATING-FUNCTION head (`oneP5PN`, `chi3PN`, `chi4PNConv`, `fiveP5PN`, …) whose
# `:=` program was stripped with the `"ExtraFunctions"` block — it stays a FORMAL atom.
const _KNOWN_REF_HEADS = Set{Symbol}([
    :SeriesData, :Rational, :ResummedSeriesData, :Sqrt, :Log, :PolyGamma, :RePsi, :ImPsi,
    :PsiA, :PsiB, :Psibar, :HoldForm, :Re, :Im, :Zeta, :Psi,
])

# `HoldForm[atom]` — the MMA "keep unevaluated" wrapper. A held SCALAR atom (`PsiA`/`PsiB`/`PolyGamma`
# …) unwraps to its held value. A held eccentric generating-function CALL (`oneP5PN[e,highestPower]`,
# `chi3PN[…]`, …) is NOT a walker-known head; its `:=` definition was stripped with `"ExtraFunctions"`,
# so it is minted as a FORMAL Symbolics atom (`oneP5PN_e_highestPower`) rather than raising.
function _unwrap_holdform(args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    length(args) == 1 || throw(ArgumentError("_unwrap_holdform: HoldForm expects 1 arg, got $(length(args))"))
    node = args[1]
    if node isa Expr && node.head === :ref && node.args[1] isa Symbol && !(node.args[1] in _KNOWN_REF_HEADS)
        return _atom_genfunc(node)
    end
    return _eval_pn_ast(node, env, locals)
end

# A held eccentric generating-function call `head[a1,a2,…]` -> the formal Symbolics atom
# `head_a1_a2_…`. The arguments (`e`, `highestPower`, …) are taken as their raw symbol names — NOT
# evaluated (`highestPower` is a formal, unbound in the alphabet) — so the atom's identity is stable
# and substitution-compatible. The `:=` generating program lives only in the stripped
# `"ExtraFunctions"` block, so the held head is the strongest thing the stored data determines.
function _atom_genfunc(node::Expr)
    head = node.args[1]::Symbol
    parts = String[string(head)]
    for a in node.args[2:end]
        push!(parts, string(a))
    end
    return Symbolics.variable(Symbol(join(parts, "_")))
end

# `Re[expr]` / `Im[expr]` -> Symbolics real/imag term (folds once `𝕚 => im` and `a` are numeric).
function _atom_reim(f, args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    length(args) == 1 || throw(ArgumentError("_atom_reim: Re/Im expect 1 arg, got $(length(args))"))
    z = _eval_pn_ast(args[1], env, locals)
    return Symbolics.Num(Symbolics.term(f, z))
end

# An opaque transcendental atom `Head[k...]` -> a single Symbolics variable `Head` joined with its
# integer args (`Zeta3`, `Log2`). A symbolic arg surfaces loudly for a later module to handle.
function _opaque_atom(head::Symbol, args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    ks = Int[]
    for a in args
        v = _eval_pn_ast(a, env, locals)
        push!(ks, _to_int(v))
    end
    name = Symbol(string(head, join(ks, "_")))
    return Symbolics.variable(name)
end

# A `key => value` rule inside an Association.
function _eval_pn_pair(pr, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    (pr isa Expr && pr.head === :call && pr.args[1] === :(=>)) ||
        throw(ArgumentError("_eval_pn_ast: Association entry is not a rule: $(repr(pr))"))
    k = _eval_pn_ast(pr.args[2], env, locals)
    k isa AbstractString || throw(ArgumentError("_eval_pn_ast: non-string Association key $(repr(k))"))
    return (String(k), _eval_pn_ast(pr.args[3], env, locals))
end

# Build a PNSeries from the six SeriesData arguments [var, x0, {coeffs}, nmin, nmax, den]. `var` is
# the raw symbol (NOT resolved through the coeff alphabet).
function _build_pnseries(args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    length(args) == 6 || throw(ArgumentError("SeriesData expects 6 args, got $(length(args))"))
    var = args[1]
    var isa Symbol || throw(ArgumentError("SeriesData variable must be a symbol, got $(repr(var))"))
    x0 = _eval_pn_ast(args[2], env, locals)
    coeffs = _eval_pn_ast(args[3], env, locals)
    coeffs isa AbstractVector || throw(ArgumentError("SeriesData coefficients must be a list"))
    nmin = _to_int(_eval_pn_ast(args[4], env, locals))
    nmax = _to_int(_eval_pn_ast(args[5], env, locals))
    den = _to_int(_eval_pn_ast(args[6], env, locals))
    return PNSeries(var, x0, Symbolics.Num[_coeff_to_num(c) for c in coeffs], nmin, nmax, den)
end

# A stored `SeriesData` coefficient may itself be a `SeriesData[e,0,{..}]` (an inner eccentricity
# series carried as the coeff of an outer y/p/v slot) — materialize such a nested `PNSeries` to its
# `Num` polynomial in `e` via `Normal`. Any other coefficient is already a `Num`.
_coeff_to_num(c) = c isa PNSeries ? Normal(c) : Symbolics.Num(c)

# `Rational[p, q]` (MMA `Rational` head) -> the EXACT rational `p//q` as a `Num` (never a Float64).
function _atom_rational(args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    length(args) == 2 || throw(ArgumentError("_atom_rational: Rational expects 2 args, got $(length(args))"))
    p = Symbolics.value(_eval_pn_ast(args[1], env, locals))
    q = Symbolics.value(_eval_pn_ast(args[2], env, locals))
    (p isa Integer && q isa Integer) ||
        throw(ArgumentError("_atom_rational: non-integer Rational args ($(repr(p)), $(repr(q)))"))
    return Symbolics.Num(p // q)
end

# `ResummedSeriesData[pref, series]` -> the materialized `pref * Normal(series)` (`Num`). Mirrors the
# upstream collapse: when `series` is a `SeriesData` (a `PNSeries` here) it is `Normal`-expanded; the
# `!MatchQ[_,_SeriesData]` fallback keeps `pref*series` for a scalar second argument.
function _atom_resummed(args, env::Dict{Symbol,Symbolics.Num}, locals::Dict{Symbol,Any})
    length(args) == 2 || throw(ArgumentError("_atom_resummed: ResummedSeriesData expects 2 args, got $(length(args))"))
    pref = Symbolics.Num(_eval_pn_ast(args[1], env, locals))
    series = _eval_pn_ast(args[2], env, locals)
    return series isa PNSeries ? pref * Normal(series) : pref * Symbolics.Num(series)
end

_to_int(num) = Int(Symbolics.value(num))

# Power that keeps an integer base to a NEGATIVE integer exponent exact (and crash-free): plain
# `Int^negInt` throws `DomainError`, and the `mantissa*10^exp` scientific rewrite produces `10^-6`
# for `*^-6`. Promote such a base to `Rational{BigInt}` so `10^-6 == 1//1000000`. All other cases
# (symbolic base, non-negative or symbolic exponent) fall through to the generic `^`.
function _pn_pow(x::Symbolics.Num, y::Symbolics.Num)
    vx = Symbolics.value(x)
    vy = Symbolics.value(y)
    if vx isa Integer && vy isa Integer && vy < 0
        return Symbolics.Num(Rational{BigInt}(vx)^vy)
    end
    return x^y
end

# Division that keeps integer/rational literals EXACT. A bare `Num(p)/Num(q)` folds to Float64, which
# then contaminates any symbolic sum it joins; when both operands are concrete integers/rationals,
# divide with `//` to stay exact.
function _pn_div(x::Symbolics.Num, y::Symbolics.Num)
    vx = Symbolics.value(x)
    vy = Symbolics.value(y)
    if (vx isa Integer || vx isa Rational) && (vy isa Integer || vy isa Rational)
        return Symbolics.Num(vx // vy)
    elseif vy isa Integer || vy isa Rational
        return x * Symbolics.Num(1 // vy)
    end
    return x / y
end

end # module MMAInputForm
