"""
    ReproHarness

Shared harness for the thin reproduction drivers under `reproduce/`.

A driver `include`s this file, registers its numbers with `@repro_result`, and ends with
`exit(ReproHarness.conclude!(opts; paper, id, script=@__FILE__, rtol))`:

    include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl")); using .ReproHarness
    opts = ReproHarness.parse_args()          # --check (default) | --regen | --tier | --paper | --only
    ReproHarness.setup!()                      # BLAS single-thread + run-log header
    using GrTools.Spacetime
    @repro_result "isco_a0.9" kerr_isco(0.9)
    exit(ReproHarness.conclude!(opts; paper = "_selftest", id = "harness-selftest",
                                script = @__FILE__, rtol = 1e-12))

Goldens are JSON files `reproduce/<paper-dir>/expected/<id>.json` with a header
`{paper, id, julia, generated_at_commit, script_sha256}` and a flat `values` object
(scalars and small arrays; the file must stay ≤ 200 KB). JSON is read and written
by the minimal serializer in this file (stdlib only), so `--check` needs no JSON
package. Encoded Julia types: Bool, Int, Float64 (non-finite values and BigFloat,
BigInt, Complex, N-d arrays as tagged `{"__type__": ...}` objects), String, Vector.

Exit codes (used by `run_all.jl`): 0 pass · 1 mismatch · 2 missing golden or harness
error · 3 HOLE/skip (the driver called `hole(reason)`; reason printed).
"""
module ReproHarness

using LinearAlgebra: BLAS, norm
using Printf
using Dates
using SHA

export @repro_result, register!, results, reset_results!,
       parse_args, setup!, expected_path, paper_dir, write_expected, read_expected,
       compare_expected, conclude!, hole, require_tier, Verdict, RowVerdict, Options,
       json_write, json_read, encode_value, decode_value

const REPRO_ROOT = normpath(joinpath(@__DIR__, ".."))
const GOLDEN_MAX_BYTES = 200 * 1024
const MAX_THREADS = 8

# ───────────────────────────── CLI options ──────────────────────────────────────────────

struct Options
    mode::Symbol                     # :check | :regen
    tier::Symbol                     # :fast | :heavy | :all  (requested tier ceiling)
    paper::Union{Nothing,String}
    only::Union{Nothing,String}
    threads::Int
    extra::Vector{String}            # unrecognised arguments, passed through to the driver
end

function _tier(s::AbstractString)
    t = Symbol(lowercase(s))
    t in (:fast, :heavy, :all) || error("--tier must be fast|heavy|all, got $(s)")
    return t
end

"""
    parse_args(args = ARGS) -> Options

Accepts `--check` (default), `--regen`, `--tier fast|heavy|all` (also `--tier=...`),
`--paper <id>`, `--only <id>`, `--threads <n>` (capped at $(MAX_THREADS)); anything else is
kept in `extra`.
"""
function parse_args(args::AbstractVector{<:AbstractString} = ARGS)
    mode = :check; tier = :all; paper = nothing; only = nothing; threads = 4
    extra = String[]
    i = 1
    takeval(flag) = begin
        i + 1 <= length(args) || error("$(flag) needs a value")
        i += 1
        String(args[i])
    end
    while i <= length(args)
        a = String(args[i])
        if a == "--check"
            mode = :check
        elseif a == "--regen"
            mode = :regen
        elseif a == "--tier"
            tier = _tier(takeval(a))
        elseif startswith(a, "--tier=")
            tier = _tier(a[8:end])
        elseif a == "--paper"
            paper = takeval(a)
        elseif startswith(a, "--paper=")
            paper = a[9:end]
        elseif a == "--only"
            only = takeval(a)
        elseif startswith(a, "--only=")
            only = a[8:end]
        elseif a == "--threads"
            threads = parse(Int, takeval(a))
        elseif startswith(a, "--threads=")
            threads = parse(Int, a[11:end])
        else
            push!(extra, a)
        end
        i += 1
    end
    threads = clamp(threads, 1, MAX_THREADS)
    return Options(mode, tier, paper, only, threads, extra)
end

"Rank of a tier for ceiling comparisons (fast < heavy; `all` admits everything)."
tier_rank(t::Symbol) = t === :fast ? 1 : t === :heavy ? 2 : 3

"""
    require_tier(opts, tier)

Exit with code 3 (HOLE/skip) when the driver's declared `tier` (`:fast`/`:heavy`) exceeds
the requested ceiling `opts.tier`. Call before any heavy computation.
"""
function require_tier(opts::Options, tier::Symbol)
    if tier_rank(tier) > tier_rank(opts.tier)
        println("SKIP: result tier '$(tier)' exceeds requested ceiling '--tier $(opts.tier)'")
        exit(3)
    end
    return nothing
end

"Print `HOLE: reason` and exit with code 3 (documented skip — never fabricate a result)."
function hole(reason::AbstractString)
    println("HOLE: ", reason)
    exit(3)
end

# ───────────────────────────── run-log header ───────────────────────────────────────────

function git_head()
    try
        return readchomp(`git -C $(REPRO_ROOT) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

"""
    setup!(; blas_threads = 1, io = stdout) -> NamedTuple

Pin BLAS to `blas_threads` (default 1: results must not depend on the thread count) and
print the run-log header (Julia version, repository commit, Julia/BLAS threads, CPU count,
load average, UTC time). Returns the recorded values.
"""
function setup!(; blas_threads::Int = 1, io::IO = stdout)
    BLAS.set_num_threads(blas_threads)
    info = (julia = string(VERSION), commit = git_head(), threads = Threads.nthreads(),
            blas_threads = BLAS.get_num_threads(), blas = string(BLAS.get_config()),
            cpus = Sys.CPU_THREADS, loadavg = round.(Sys.loadavg(); digits = 2),
            time_utc = string(Dates.now(Dates.UTC)))
    println(io, "[repro] julia $(info.julia)  commit $(info.commit)")
    println(io, "[repro] threads=$(info.threads)  blas_threads=$(info.blas_threads)  cpus=$(info.cpus)  loadavg=$(info.loadavg)")
    println(io, "[repro] blas: $(info.blas)")
    println(io, "[repro] time: $(info.time_utc) UTC")
    return info
end

# ───────────────────────────── result registry ──────────────────────────────────────────

const RESULTS = Vector{Pair{String,Any}}()

"Register a named result value (names must be unique within a driver)."
function register!(name::AbstractString, value)
    any(p -> first(p) == name, RESULTS) && error("repro result '$(name)' registered twice")
    push!(RESULTS, String(name) => value)
    return value
end

"Registered results in registration order, as an ordered `Dict`-like vector of pairs."
results() = copy(RESULTS)
reset_results!() = (empty!(RESULTS); nothing)

"""
    @repro_result "name" expr

Evaluate `expr`, register its value under `name`, and return the value.
"""
macro repro_result(name, expr)
    quote
        local v = $(esc(expr))
        register!($(esc(name)), v)
        v
    end
end

# ───────────────────────────── paths ────────────────────────────────────────────────────

"Directory under `reproduce/` for a paper id: `_x` → `_x`; `arXiv-…` → as is; else `arXiv-<id>`."
function paper_dir(paper::AbstractString)
    (startswith(paper, "_") || startswith(paper, "arXiv-")) && return String(paper)
    return "arXiv-" * paper
end

"`reproduce/<paper-dir>/expected/<id>.json` (absolute path rooted at this checkout)."
expected_path(paper::AbstractString, id::AbstractString) =
    joinpath(REPRO_ROOT, paper_dir(paper), "expected", String(id) * ".json")

# ───────────────────────────── minimal JSON ─────────────────────────────────────────────
# Subset sufficient for goldens: objects (String keys), arrays, numbers, strings, true/false/null.

function _json_escape(io::IO, s::AbstractString)
    print(io, '"')
    for c in s
        if c == '"'; print(io, "\\\"")
        elseif c == '\\'; print(io, "\\\\")
        elseif c == '\n'; print(io, "\\n")
        elseif c == '\r'; print(io, "\\r")
        elseif c == '\t'; print(io, "\\t")
        elseif c == '\b'; print(io, "\\b")
        elseif c == '\f'; print(io, "\\f")
        elseif UInt32(c) < 0x20; @printf(io, "\\u%04x", UInt32(c))
        else; print(io, c)
        end
    end
    print(io, '"')
end

_json_scalar(x::Bool) = x ? "true" : "false"
_json_scalar(::Nothing) = "null"
_json_scalar(x::Integer) = string(x)
function _json_scalar(x::AbstractFloat)
    isfinite(x) || error("non-finite float reached the JSON writer; encode_value tags them")
    return repr(Float64(x))
end

_is_flat(v::AbstractVector) = all(x -> !(x isa AbstractVector || x isa AbstractDict), v)

function _json_write(io::IO, x, indent::Int)
    pad = repeat(' ', indent)
    if x isa AbstractDict
        ks = sort!(collect(String.(keys(x))))
        if isempty(ks); print(io, "{}"); return; end
        print(io, "{\n")
        for (k, key) in enumerate(ks)
            print(io, pad, "  "); _json_escape(io, key); print(io, ": ")
            _json_write(io, x[key], indent + 2)
            k < length(ks) && print(io, ",")
            print(io, "\n")
        end
        print(io, pad, "}")
    elseif x isa AbstractVector
        if isempty(x); print(io, "[]"); return; end
        if _is_flat(x)
            print(io, "[")
            for (k, v) in enumerate(x)
                k > 1 && print(io, ", ")
                _json_write(io, v, indent + 2)
            end
            print(io, "]")
        else
            print(io, "[\n")
            for (k, v) in enumerate(x)
                print(io, pad, "  "); _json_write(io, v, indent + 2)
                k < length(x) && print(io, ",")
                print(io, "\n")
            end
            print(io, pad, "]")
        end
    elseif x isa AbstractString
        _json_escape(io, x)
    elseif x isa Bool || x isa Integer || x isa AbstractFloat || x === nothing
        print(io, _json_scalar(x))
    else
        error("json_write: unsupported type $(typeof(x)) — pass values through encode_value first")
    end
end

"Serialise `x` (Dict/Vector/String/Real/Bool/nothing tree) as indented JSON text."
function json_write(x)
    io = IOBuffer(); _json_write(io, x, 0); print(io, "\n")
    return String(take!(io))
end

struct _P; s::String; n::Int; end   # parser state (immutable; positions passed explicitly)

_skipws(p::_P, i::Int) = (while i <= p.n && (p.s[i] in (' ', '\t', '\n', '\r')); i = nextind(p.s, i); end; i)

function _perror(p::_P, i::Int, msg)
    error("json_read: $(msg) at byte $(i)/$(p.n)")
end

function _parse_value(p::_P, i::Int)
    i = _skipws(p, i)
    i <= p.n || _perror(p, i, "unexpected end of input")
    c = p.s[i]
    if c == '{'; return _parse_object(p, i)
    elseif c == '['; return _parse_array(p, i)
    elseif c == '"'; return _parse_string(p, i)
    elseif startswith(SubString(p.s, i), "true"); return (true, i + 4)
    elseif startswith(SubString(p.s, i), "false"); return (false, i + 5)
    elseif startswith(SubString(p.s, i), "null"); return (nothing, i + 4)
    else; return _parse_number(p, i)
    end
end

function _parse_object(p::_P, i::Int)
    d = Dict{String,Any}()
    i = _skipws(p, i + 1)
    i <= p.n || _perror(p, i, "unterminated object")
    if p.s[i] == '}'; return (d, i + 1); end
    while true
        i = _skipws(p, i)
        (i <= p.n && p.s[i] == '"') || _perror(p, i, "expected string key")
        key, i = _parse_string(p, i)
        i = _skipws(p, i)
        (i <= p.n && p.s[i] == ':') || _perror(p, i, "expected ':'")
        val, i = _parse_value(p, i + 1)
        d[key] = val
        i = _skipws(p, i)
        i <= p.n || _perror(p, i, "unterminated object")
        if p.s[i] == ','; i += 1; continue
        elseif p.s[i] == '}'; return (d, i + 1)
        else; _perror(p, i, "expected ',' or '}'")
        end
    end
end

function _parse_array(p::_P, i::Int)
    a = Any[]
    i = _skipws(p, i + 1)
    i <= p.n || _perror(p, i, "unterminated array")
    if p.s[i] == ']'; return (a, i + 1); end
    while true
        val, i = _parse_value(p, i)
        push!(a, val)
        i = _skipws(p, i)
        i <= p.n || _perror(p, i, "unterminated array")
        if p.s[i] == ','; i += 1; continue
        elseif p.s[i] == ']'; return (a, i + 1)
        else; _perror(p, i, "expected ',' or ']'")
        end
    end
end

function _parse_string(p::_P, i::Int)
    io = IOBuffer()
    i += 1   # opening quote
    while i <= p.n
        c = p.s[i]
        if c == '"'
            return (String(take!(io)), i + 1)
        elseif c == '\\'
            i += 1
            i <= p.n || _perror(p, i, "bad escape")
            e = p.s[i]
            if e == 'n'; print(io, '\n')
            elseif e == 't'; print(io, '\t')
            elseif e == 'r'; print(io, '\r')
            elseif e == 'b'; print(io, '\b')
            elseif e == 'f'; print(io, '\f')
            elseif e == '"'; print(io, '"')
            elseif e == '\\'; print(io, '\\')
            elseif e == '/'; print(io, '/')
            elseif e == 'u'
                i + 4 <= p.n || _perror(p, i, "bad \\u escape")
                cp = parse(UInt32, p.s[i+1:i+4]; base = 16)
                i += 4
                if 0xD800 <= cp <= 0xDBFF && i + 6 <= p.n && p.s[i+1] == '\\' && p.s[i+2] == 'u'
                    lo = parse(UInt32, p.s[i+3:i+6]; base = 16)
                    cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00)
                    i += 6
                end
                print(io, Char(cp))
            else
                _perror(p, i, "unknown escape \\$(e)")
            end
            i += 1
        else
            print(io, c)
            i = nextind(p.s, i)
        end
    end
    _perror(p, i, "unterminated string")
end

function _parse_number(p::_P, i::Int)
    j = i
    while j <= p.n && (p.s[j] in "+-0123456789.eE"); j += 1; end
    tok = p.s[i:j-1]
    isempty(tok) && _perror(p, i, "unexpected character '$(p.s[i])'")
    if occursin(r"[.eE]", tok)
        v = tryparse(Float64, tok)
        v === nothing && _perror(p, i, "bad number '$(tok)'")
        return (v, j)
    else
        v = tryparse(Int, tok)
        if v === nothing
            b = tryparse(BigInt, tok)
            b === nothing && _perror(p, i, "bad number '$(tok)'")
            return (b, j)
        end
        return (v, j)
    end
end

"Parse JSON text into Dict{String,Any} / Vector{Any} / String / Int / Float64 / Bool / nothing."
function json_read(s::AbstractString)
    p = _P(String(s), ncodeunits(s))
    v, i = _parse_value(p, 1)
    i = _skipws(p, i)
    i > p.n || _perror(p, i, "trailing characters")
    return v
end

# ───────────────────────────── Julia ⇄ JSON value encoding ─────────────────────────────

_tagged(t, kv...) = Dict{String,Any}("__type__" => t, kv...)

"Encode a Julia value into the JSON-serialisable tree used inside `values`."
encode_value(x::Bool) = x
encode_value(x::Integer) = (typemin(Int) <= x <= typemax(Int)) ? Int(x) : _tagged("BigInt", "value" => string(x))
function encode_value(x::AbstractFloat)
    if x isa BigFloat
        return _tagged("BigFloat", "value" => string(x), "precision" => precision(x))
    end
    isfinite(x) && return Float64(x)
    return _tagged("Float64", "value" => string(Float64(x)))
end
encode_value(x::Complex) = _tagged("Complex", "re" => encode_value(real(x)), "im" => encode_value(imag(x)))
encode_value(x::AbstractString) = String(x)
encode_value(x::Symbol) = String(x)
encode_value(::Nothing) = nothing
encode_value(::Missing) = nothing
encode_value(x::AbstractVector) = Any[encode_value(v) for v in x]
encode_value(x::AbstractArray) = _tagged("Array", "size" => collect(size(x)),
                                         "data" => Any[encode_value(v) for v in vec(x)])
encode_value(x::Tuple) = Any[encode_value(v) for v in x]
encode_value(x::NamedTuple) = Dict{String,Any}(String(k) => encode_value(v) for (k, v) in pairs(x))
encode_value(x::AbstractDict) = Dict{String,Any}(String(k) => encode_value(v) for (k, v) in x)
encode_value(x) = error("encode_value: unsupported result type $(typeof(x)); reduce it to numbers/arrays first")

"Decode the JSON tree back into Julia values (inverse of `encode_value`)."
function decode_value(x)
    if x isa AbstractDict
        t = get(x, "__type__", nothing)
        if t === nothing
            return Dict{String,Any}(String(k) => decode_value(v) for (k, v) in x)
        elseif t == "BigFloat"
            prec = Int(get(x, "precision", 256))
            return setprecision(BigFloat, prec) do
                parse(BigFloat, x["value"])
            end
        elseif t == "BigInt"
            return parse(BigInt, x["value"])
        elseif t == "Float64"
            return parse(Float64, x["value"])
        elseif t == "Complex"
            return complex(decode_value(x["re"]), decode_value(x["im"]))
        elseif t == "Array"
            sz = Tuple(Int.(x["size"]))
            data = Any[decode_value(v) for v in x["data"]]
            return reshape(_narrow(data), sz)
        else
            error("decode_value: unknown __type__ '$(t)'")
        end
    elseif x isa AbstractVector
        return _narrow(Any[decode_value(v) for v in x])
    else
        return x
    end
end

"Narrow a Vector{Any} of homogeneous scalars to a concrete element type where possible."
function _narrow(v::Vector{Any})
    isempty(v) && return v
    T = mapreduce(typeof, promote_type, v)
    (T === Any || T isa Union) && return v
    T <: Number || T === String || T === Bool || return v
    return T[T(x) for x in v]
end

# ───────────────────────────── golden I/O ───────────────────────────────────────────────

"sha256 hex digest of a file's bytes (used for `script_sha256`)."
file_sha256(path::AbstractString) = bytes2hex(open(sha256, path))

"""
    write_expected(path, values::AbstractDict; paper, id, script) -> nbytes

Write the golden JSON `{paper, id, julia, generated_at_commit, script_sha256, values}`.
Values are encoded with `encode_value`; the file must be ≤ $(GOLDEN_MAX_BYTES ÷ 1024) KB
(reduce big arrays to the plotted quantities before registering them).
"""
function write_expected(path::AbstractString, values::AbstractDict; paper::AbstractString,
                        id::AbstractString, script::AbstractString = "")
    doc = Dict{String,Any}(
        "paper" => String(paper), "id" => String(id), "julia" => string(VERSION),
        "generated_at_commit" => git_head(),
        "script_sha256" => (isempty(script) || !isfile(script)) ? "" : file_sha256(script),
        "values" => Dict{String,Any}(String(k) => encode_value(v) for (k, v) in values))
    txt = json_write(doc)
    n = ncodeunits(txt)
    n <= GOLDEN_MAX_BYTES || error("golden $(path) would be $(n) bytes > $(GOLDEN_MAX_BYTES); reduce the registered arrays")
    mkpath(dirname(path))
    write(path, txt)
    return n
end

write_expected(path::AbstractString, values::AbstractVector{<:Pair}; kw...) =
    write_expected(path, Dict{String,Any}(values); kw...)

"Read a golden; returns `(header::Dict, values::Dict{String,Any})` with decoded values."
function read_expected(path::AbstractString)
    doc = json_read(read(path, String))
    doc isa AbstractDict || error("golden $(path) is not a JSON object")
    haskey(doc, "values") || error("golden $(path) has no 'values' entry")
    header = Dict{String,Any}(k => v for (k, v) in doc if k != "values")
    vals = Dict{String,Any}(String(k) => decode_value(v) for (k, v) in doc["values"])
    return header, vals
end

# ───────────────────────────── comparison ───────────────────────────────────────────────

struct RowVerdict
    name::String
    expected::String
    computed::String
    relerr::Float64
    pass::Bool
    note::String
end

struct Verdict
    path::String
    rows::Vector{RowVerdict}
    pass::Bool
    n_pass::Int
    n_fail::Int
    notes::Vector{String}
end

_fmt(x::Bool) = string(x)
_fmt(x::Integer) = string(x)
_fmt(x::AbstractFloat) = isfinite(x) ? (@sprintf("%.12g", x)) : string(x)
_fmt(x::BigFloat) = begin s = @sprintf("%.20g", x); s end
_fmt(x::Complex) = _fmt(real(x)) * (imag(x) < 0 ? "-" : "+") * _fmt(abs(imag(x))) * "im"
_fmt(x::AbstractString) = repr(x)
_fmt(::Nothing) = "null"
function _fmt(x::AbstractArray)
    n = length(x)
    if n <= 3
        return "[" * join(_fmt.(x), ", ") * "]"
    end
    shape = ndims(x) == 1 ? "Vector($(n))" : "Array$(size(x))"
    return shape * " [" * _fmt(first(x)) * ", …, " * _fmt(last(x)) * "]"
end
_fmt(x) = string(x)

_relerr(d, e) = abs(e) > 0 ? d / abs(e) : (d == 0 ? 0.0 : Inf)

"Compare one expected/computed pair → (pass, relerr, note)."
function compare_values(e, c; rtol::Real, atol::Real)
    if e isa Bool || c isa Bool || e isa AbstractString || c isa AbstractString || e === nothing || c === nothing
        ok = isequal(e, c)
        return (ok, ok ? 0.0 : NaN, ok ? "" : "value mismatch (exact compare)")
    elseif e isa Number && c isa Number
        if !(isfinite(e) && isfinite(c))
            ok = isequal(e, c) || (isnan(e) && isnan(c))
            return (ok, ok ? 0.0 : NaN, ok ? "non-finite, exact" : "non-finite mismatch")
        end
        d = Float64(abs(c - e))
        ok = isapprox(c, e; rtol = rtol, atol = atol)
        return (ok, _relerr(d, Float64(abs(e))), ok ? "" : "|Δ|=" * _fmt(d))
    elseif e isa AbstractArray && c isa AbstractArray
        size(e) == size(c) || return (false, NaN, "size $(size(c)) ≠ expected $(size(e))")
        isempty(e) && return (true, 0.0, "")
        if !(eltype(e) <: Number && eltype(c) <: Number)
            ok = isequal(e, c)
            return (ok, ok ? 0.0 : NaN, ok ? "" : "array of non-numbers mismatch")
        end
        fe = isfinite.(e); fc = isfinite.(c)
        fe == fc || return (false, NaN, "non-finite pattern mismatch")
        all(isequal.(e[.!fe], c[.!fc])) || return (false, NaN, "non-finite entries mismatch")
        ef = e[fe]; cf = c[fc]
        isempty(ef) && return (true, 0.0, "all non-finite, exact")
        d = norm(cf .- ef); ne = norm(ef)
        ok = isapprox(cf, ef; rtol = rtol, atol = atol)        # norm-based, Julia semantics
        maxabs = maximum(abs.(cf .- ef))
        return (ok, _relerr(d, ne), (ok ? "" : "FAIL ") * "norm-rel; max|Δ|=" * _fmt(maxabs))
    else
        return (false, NaN, "type mismatch: expected $(typeof(e)), computed $(typeof(c))")
    end
end

"""
    compare_expected(path, computed; rtol = 1e-10, atol = 0.0, tolerances = Dict(), io = stdout) -> Verdict

Compare a golden file against the computed values (`Dict` or vector of pairs), print the
table `name | expected | computed | rel err | pass`, and return the `Verdict`.
`tolerances` may override `(rtol, atol)` per name (`name => (rtol = …, atol = …)`). Names
present on one side only are failures (run `--regen` after adding a result).
Arrays are compared norm-wise (`isapprox` semantics); the reported error is
`‖c−e‖/‖e‖`.
"""
function compare_expected(path::AbstractString, computed; rtol::Real = 1e-10, atol::Real = 0.0,
                          tolerances::AbstractDict = Dict{String,Any}(), io::IO = stdout)
    comp = collect(computed)
    header, expected = read_expected(path)
    rows = RowVerdict[]
    notes = String[]
    seen = Set{String}()
    for (name, cval) in comp
        name = String(name)
        push!(seen, name)
        tol = get(tolerances, name, nothing)
        rt = tol === nothing ? rtol : get(tol, :rtol, rtol)
        at = tol === nothing ? atol : get(tol, :atol, atol)
        if !haskey(expected, name)
            push!(rows, RowVerdict(name, "—", _fmt(cval), NaN, false, "not in golden (run --regen)"))
            continue
        end
        eval_ = expected[name]
        ok, rel, note = compare_values(eval_, cval; rtol = rt, atol = at)
        push!(rows, RowVerdict(name, _fmt(eval_), _fmt(cval), rel, ok, note))
    end
    for name in sort!(collect(setdiff(keys(expected), seen)))
        push!(rows, RowVerdict(name, _fmt(expected[name]), "—", NaN, false, "in golden but not computed"))
    end
    npass = count(r -> r.pass, rows); nfail = length(rows) - npass
    # header provenance notes (informational, never a failure)
    gj = get(header, "julia", ""); gj == string(VERSION) || push!(notes, "golden generated with julia $(gj) (running $(VERSION))")
    gc = get(header, "generated_at_commit", ""); isempty(gc) || push!(notes, "golden generated at commit $(gc)")
    _print_table(io, path, rows, notes)
    return Verdict(String(path), rows, nfail == 0 && !isempty(rows), npass, nfail, notes)
end

function _print_table(io::IO, path, rows, notes)
    w1 = max(4, maximum(length(r.name) for r in rows; init = 4))
    w2 = max(8, min(34, maximum(length(r.expected) for r in rows; init = 8)))
    w3 = max(8, min(34, maximum(length(r.computed) for r in rows; init = 8)))
    clip(s, w) = length(s) > w ? first(s, w - 1) * "…" : s
    println(io, "[repro] golden: ", relpath(path, REPRO_ROOT))
    println(io, rpad("name", w1), " | ", rpad("expected", w2), " | ", rpad("computed", w3), " | ", rpad("rel err", 10), " | pass")
    println(io, repeat("-", w1), "-+-", repeat("-", w2), "-+-", repeat("-", w3), "-+-", repeat("-", 10), "-+-----")
    for r in rows
        rel = isnan(r.relerr) ? "n/a" : @sprintf("%.2e", r.relerr)
        tail = isempty(r.note) ? "" : "   # " * r.note
        println(io, rpad(r.name, w1), " | ", rpad(clip(r.expected, w2), w2), " | ", rpad(clip(r.computed, w3), w3),
                " | ", rpad(rel, 10), " | ", r.pass ? "PASS" : "FAIL", tail)
    end
    for n in notes
        println(io, "[repro] note: ", n)
    end
end

# ───────────────────────────── driver epilogue ──────────────────────────────────────────

"""
    conclude!(opts; paper, id, script = "", rtol = 1e-10, atol = 0.0, tolerances = Dict()) -> Int

End of a driver: on `--regen` write the golden from the registered results; on `--check`
compare and print the table. Returns the process exit code (0 pass, 1 mismatch,
2 missing golden). Drivers end with `exit(conclude!(...))`.
"""
function conclude!(opts::Options; paper::AbstractString, id::AbstractString, script::AbstractString = "",
                   rtol::Real = 1e-10, atol::Real = 0.0, tolerances::AbstractDict = Dict{String,Any}(), io::IO = stdout)
    path = expected_path(paper, id)
    vals = results()
    isempty(vals) && (println(io, "[repro] ERROR: no results registered (use @repro_result)"); return 2)
    if opts.mode === :regen
        n = write_expected(path, vals; paper = paper, id = id, script = script)
        println(io, "[repro] wrote ", relpath(path, REPRO_ROOT), " ($(length(vals)) values, $(n) bytes)")
        for (k, v) in vals
            println(io, "  ", rpad(k, 28), " = ", _fmt(v))
        end
        return 0
    end
    if !isfile(path)
        println(io, "[repro] ERROR: no golden at ", relpath(path, REPRO_ROOT), " — run with --regen first")
        return 2
    end
    v = compare_expected(path, vals; rtol = rtol, atol = atol, tolerances = tolerances, io = io)
    if !isempty(script) && isfile(script)
        h, _ = read_expected(path)
        sh = get(h, "script_sha256", "")
        if !isempty(sh) && sh != file_sha256(script)
            println(io, "[repro] note: driver script changed since the golden was generated (sha256 differs)")
        end
    end
    println(io, "[repro] ", v.pass ? "PASS" : "FAIL", " $(v.n_pass)/$(length(v.rows)) values within tolerance (rtol=$(rtol), atol=$(atol))")
    return v.pass ? 0 : 1
end

end # module
