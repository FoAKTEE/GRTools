# Mathematica `.wl`/`.m` InputForm NUMERIC-literal reader — the unified `read_mma_inputform` target.
# Self-contained recursive descent (no `eval`, no world-age) over the WaSABI asset subset: `(* *)`
# comments (nesting), `{}` lists, `<||>` associations, `"..."` strings, `Around[c,e]` → `c`,
# `` `precision `` markers and `*^` exponents. De-dups WaSABI/wl_parser.jl.
#
# FORMAT SEPARATION (design §10 — Data owns the NUMERIC dialect ONLY). The SYMBOLIC MMA dialect — the
# PN `.m` `SeriesData`/`Association`/`Module` store walked into `Symbolics.Num` (`_eval_pn_ast` /
# `parse_pn_data_file`) — is DELEGATED to `PostNewtonian.MMAInputForm` (the Catalog's symbolic side,
# behind the Symbolics quarantine). That is a sibling reader, NOT a Data target: Data must stay
# Symbolics-free on the numeric path. So there is deliberately NO `series_data` / symbolic-`Num` walker
# in Data — such a walker here would be a mis-filed responsibility, not a gap.

# --- comment stripping ----------------------------------------------------

# Remove InputForm line-continuations (`\` before a newline — they split long high-precision numbers),
# then strip `(* ... *)` comments (which may nest).
function _strip_wl_comments(s::AbstractString)
    cs = collect(replace(s, r"\\\r?\n" => ""))
    out = Char[]
    depth = 0
    i = 1
    n = length(cs)
    while i <= n
        if i < n && cs[i] == '(' && cs[i+1] == '*'
            depth += 1
            i += 2
        elseif depth > 0 && i < n && cs[i] == '*' && cs[i+1] == ')'
            depth -= 1
            i += 2
        elseif depth > 0
            i += 1
        else
            push!(out, cs[i])
            i += 1
        end
    end
    return out
end

# --- numeric token --------------------------------------------------------

# Reduce a Mathematica numeric token to a plain decimal string: drop the `` `precision `` marker and
# rewrite `*^exp` as `eexp`. E.g. "2.486*^-7" -> "2.486e-7", "5.469`96.37*^6" -> "5.469e6".
function _wl_number_base(tok::AbstractString)
    base = tok
    if occursin('`', tok)
        parts = split(tok, '`')
        base = parts[1]                # single backtick -> precision; `` -> accuracy
        rest = parts[end]              # may still carry a trailing *^exp
        if occursin("*^", rest)
            _, e = split(rest, "*^")
            base = base * "e" * e
        end
    elseif occursin("*^", tok)
        m, e = split(tok, "*^")
        base = m * "e" * e
    end
    return base
end

# Parse one numeric token. `want_big` keeps full file precision as `BigFloat` (default 256-bit — above
# any `` `precision `` marker in these assets), so all file digits are preserved.
function _parse_number_token(tok::AbstractString; want_big::Bool)
    base = _wl_number_base(tok)
    T = want_big ? BigFloat : Float64
    if occursin('/', base)                     # rational literal `a/b`
        num, den = split(base, '/'; limit = 2)
        return parse(T, num) / parse(T, den)
    end
    return parse(T, base)
end

# Parse a real WL token to `Float64`, honouring rational literals `a/b`.
function _parse_real_token(tok::AbstractString)
    base = _wl_number_base(strip(tok))
    if occursin('/', base)
        num, den = split(base, '/'; limit = 2)
        return parse(Float64, num) / parse(Float64, den)
    end
    return parse(Float64, base)
end

# --- recursive descent over Vector{Char} ----------------------------------

_isws(c::Char) = c == ' ' || c == '\t' || c == '\n' || c == '\r'

function _skipws(cs, i)
    n = length(cs)
    while i <= n && _isws(cs[i])
        i += 1
    end
    return i
end

_matches(cs, i, lit) = i + length(lit) - 1 <= length(cs) &&
    all(cs[i+k-1] == lit[k] for k in 1:length(lit))

function _parse_value(cs, i; want_big::Bool)
    i = _skipws(cs, i)
    c = cs[i]
    if c == '{'
        return _parse_list(cs, i; want_big)
    elseif c == '<' && _matches(cs, i, "<|")
        return _parse_assoc(cs, i; want_big)
    elseif c == '"'
        return _parse_string(cs, i)
    elseif _matches(cs, i, "Around[")
        return _parse_around(cs, i; want_big)
    else
        return _parse_number(cs, i; want_big)
    end
end

function _parse_list(cs, i; want_big::Bool)
    @assert cs[i] == '{'
    i += 1
    items = Any[]
    while true
        i = _skipws(cs, i)
        if cs[i] == '}'
            return items, i + 1
        end
        v, i = _parse_value(cs, i; want_big)
        push!(items, v)
        i = _skipws(cs, i)
        if cs[i] == ','
            i += 1
        end
    end
end

function _parse_assoc(cs, i; want_big::Bool)
    @assert _matches(cs, i, "<|")
    i += 2
    d = Dict{String,Any}()
    while true
        i = _skipws(cs, i)
        if _matches(cs, i, "|>")
            return d, i + 2
        end
        key, i = _parse_string(cs, i)
        i = _skipws(cs, i)
        @assert _matches(cs, i, "->") "expected -> after association key $(repr(key))"
        i += 2
        v, i = _parse_value(cs, i; want_big)
        d[key] = v
        i = _skipws(cs, i)
        if cs[i] == ','
            i += 1
        end
    end
end

function _parse_string(cs, i)
    @assert cs[i] == '"'
    i += 1
    buf = Char[]
    n = length(cs)
    while i <= n
        c = cs[i]
        if c == '\\' && i < n && (cs[i+1] == '"' || cs[i+1] == '\\')
            push!(buf, cs[i+1])     # true escape: \" -> " and \\ -> \
            i += 2
        elseif c == '\\'
            push!(buf, c)           # keep named-character escapes verbatim, e.g. \[CapitalOmega]
            i += 1
        elseif c == '"'
            return String(buf), i + 1
        else
            push!(buf, c)
            i += 1
        end
    end
    error("Unterminated string literal.")
end

function _parse_around(cs, i; want_big::Bool)
    i += length("Around[")
    central, i = _parse_value(cs, i; want_big)   # keep central value only ("Around-stripping")
    i = _skipws(cs, i)
    @assert cs[i] == ',' "Around[ expects a comma"
    i += 1
    _err, i = _parse_value(cs, i; want_big)       # error term parsed then discarded
    i = _skipws(cs, i)
    @assert cs[i] == ']' "Around[ expects a closing bracket"
    return central, i + 1
end

function _parse_number(cs, i; want_big::Bool)
    n = length(cs)
    j = i
    while j <= n && !(cs[j] in (',', '}', ']', '|', '<', '"')) && !_isws(cs[j])
        j += 1
    end
    tok = String(cs[i:j-1])
    return _parse_number_token(tok; want_big), j
end

# --- public API -----------------------------------------------------------

"""
    read_mma_inputform(path; want_big=false) -> parsed

The canonical umbrella reader: parse a Mathematica InputForm `.wl`/`.m` file (numeric dialect) into
nested `Vector{Any}` lists / `Dict{String,Any}` associations / strings / numbers, with `Around[c,e]`
reduced to `c`, `` `precision `` and `*^` honoured. `want_big=true` keeps full precision as `BigFloat`.
No `eval`, no world-age — pure string parsing.
"""
function read_mma_inputform(path::AbstractString; want_big::Bool = false)
    cs = _strip_wl_comments(read(path, String))
    val, _ = _parse_value(cs, 1; want_big)
    return val
end

"""
    read_wl_table(path) -> Vector{Tuple{Float64,Float64}}

Parse a `.wl` flat list of numeric `{x, y}` pairs into `Float64` tuples (e.g. the 1SF flux tables).
"""
function read_wl_table(path::AbstractString)
    val = read_mma_inputform(path; want_big = false)
    return [(Float64(p[1]), Float64(p[2])) for p in val]
end

"""
    read_wl_data(path) -> Vector{Dict{String,Any}}

Parse a `.data` list of associations into one `Dict` per `<|...|>` record, preserving full precision
as `BigFloat` (so a later JLD2 round-trip is bit-for-bit). `Around[c,e]` reduced to `c`.
"""
function read_wl_data(path::AbstractString)
    val = read_mma_inputform(path; want_big = true)
    return Vector{Dict{String,Any}}(val)
end

"""
    read_wl_grid2d(path) -> Vector{Tuple{Tuple{Float64,Float64},Float64}}

Parse a `.wl` 2D grid table `{{x, y}, z}` into `((x, y), z)` `Float64` tuples.
"""
function read_wl_grid2d(path::AbstractString)
    val = read_mma_inputform(path; want_big = false)
    return [((Float64(p[1][1]), Float64(p[1][2])), Float64(p[2])) for p in val]
end

"""
    read_wl_grid3d(path) -> Vector{Tuple{NTuple{3,Float64},Float64}}

Parse a `.wl` 3D grid table `{{x, y, z}, f}` into `((x, y, z), f)` `Float64` tuples.
"""
function read_wl_grid3d(path::AbstractString)
    val = read_mma_inputform(path; want_big = false)
    return [((Float64(p[1][1]), Float64(p[1][2]), Float64(p[1][3])), Float64(p[2])) for p in val]
end

# Parse one Wolfram complex literal (`re`, `im*I`, `re ± im*I`, with `*^` exponents) into ComplexF64.
function _wl_complex(s::AbstractString)
    e = replace(replace(strip(s), "*^" => "e"), " " => "")
    im = 0.0
    re_str = e
    mim = match(r"([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\*I", e)
    if mim !== nothing
        im = parse(Float64, mim.captures[1])
        re_str = replace(e, mim.match => "")
    elseif occursin("I", e)
        mb = match(r"([+-]?)I", e)
        im = (mb !== nothing && mb.captures[1] == "-") ? -1.0 : 1.0
        re_str = mb === nothing ? e : replace(e, mb.match => "")
    end
    re = (isempty(re_str) || re_str == "+" || re_str == "-") ? 0.0 : parse(Float64, re_str)
    return ComplexF64(re, im)
end

"""
    read_wl_table_complex(path) -> Vector{Tuple{Float64,ComplexF64}}

Parse a `.wl` flat list of `{x, y}` pairs whose `y` is a complex literal `re ± im*I` (e.g. the
Teukolsky amplitude `Z` tables, the 2SF `{Log[r0], Log[flux]}` table). `*^` exponents honoured.
"""
function read_wl_table_complex(path::AbstractString)
    s = String(_strip_wl_comments(read(path, String)))
    out = Tuple{Float64,ComplexF64}[]
    for m in eachmatch(r"\{([^{}]+)\}", s)
        inner = m.captures[1]
        ci = findfirst(',', inner)
        ci === nothing && continue
        x = _parse_real_token(inner[1:prevind(inner, ci)])
        y = _wl_complex(inner[nextind(inner, ci):end])
        push!(out, (x, y))
    end
    return out
end

export read_mma_inputform, read_wl_table, read_wl_data, read_wl_grid2d, read_wl_grid3d,
    read_wl_table_complex
