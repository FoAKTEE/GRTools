# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/s4c_capture_fixtures.jl (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = b89a396763cd3db1de0818f1eef69cbc3084ce4daf9919f483508329ba76f730
# fixture(s)  : r2_coupling_ref.jl, aak_bessel_ref.jl
# oracle      : PRE-HOIST core/Kernels @ 51efa9f (Wigner + SpecialFunctions still Kernels deps); Wigner.coupling
# run         : julia --project=core/Kernels test/fixtures/capture/s4c_capture_fixtures.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const GRTOOLS_ROOT         = get(ENV, "GRTOOLS_ROOT",                 normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..")))
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

# S4-c fixture capture (ONE-SHOT, run against the PRE-CHANGE committed Kernels @ 51efa9f state).
# Run: julia --project=grtools/core/Kernels scratchpad/s4c_capture_fixtures.jl
#
# Bakes, from the CURRENT committed code (Wigner + SpecialFunctions still Kernels deps):
#   (1) test/fixtures/r2_coupling_ref.jl  — the host-side coupling arrays c the committed r2_N!
#       kernels computed internally (Wigner.coupling) for every (n, branch) of test_r2.jl, plus the
#       (s3,s1,s2) spec PARSED from the committed r2.jl source (no hand transcription).
#   (2) test/fixtures/aak_bessel_ref.jl   — ak_bessel_coeffs(n,e) 5-tuples at every (n,e) the aak
#       tests drive (keys bit-exact to the kernel's spline-resampled requests).
#   (3) stdout: pre-hoist Kernels.aak_summation output literals on the EMRIWaveforms test_aak inputs
#       (to be baked as the EMRIWaveforms behavior-identity gate).

using Kernels
using Wigner: coupling
using NumericsCore: CubicSplineInterpolant

const KROOT = joinpath(GRTOOLS_ROOT, "core", "Kernels")

# ---------- (1) parse the committed r2.jl coupling spec ----------
spec = Dict{Int,Vector{NTuple{3,Int}}}()
let cur = 0
    for line in eachline(joinpath(KROOT, "src", "r2.jl"))
        m = match(r"^function r2_(\d+)!", line)
        m !== nothing && (cur = parse(Int, m.captures[1]); spec[cur] = NTuple{3,Int}[])
        c = match(r"^\s*c(\d+) = coupling\(l3,m3,(-?\d+),l1,m1,(-?\d+),l2,m2,(-?\d+)\);", line)
        if c !== nothing
            @assert cur >= 1
            k = parse(Int, c.captures[1])
            @assert k == length(spec[cur]) + (cur in (4, 5, 8, 9) ? 8 : cur in (7, 10) ? 18 : 1) "c-index order drift in r2_$cur"
            push!(spec[cur], (parse(Int, c.captures[2]), parse(Int, c.captures[3]), parse(Int, c.captures[4])))
        end
    end
end
@assert sort(collect(keys(spec))) == 1:10
@assert map(n -> length(spec[n]), 1:10) == [5, 5, 5, 10, 10, 7, 10, 10, 10, 10]

# branches: verbatim from test/test_r2.jl R2_BRANCHES
branches = [
    (l1=2, l2=2, l3=2, m3=2, m1=1, m2=1),
    (l1=0, l2=0, l3=0, m3=0, m1=0, m2=0),
    (l1=1, l2=0, l3=1, m3=1, m1=1, m2=0),
    (l1=0, l2=1, l3=1, m3=1, m1=0, m2=1),
    (l1=1, l2=1, l3=2, m3=2, m1=1, m2=1),
    (l1=0, l2=2, l3=2, m3=1, m1=0, m2=1),
    (l1=2, l2=0, l3=2, m3=1, m1=1, m2=0),
    (l1=1, l2=2, l3=3, m3=2, m1=1, m2=1),
    (l1=2, l2=1, l3=3, m3=2, m1=1, m2=1),
]

open(joinpath(FIXTURE_DIR, "r2_coupling_ref.jl"), "w") do io
    println(io, "# Baked host-side coupling arrays for the 10 r2 kernels x 9 (l1,l2) oracle branches.")
    println(io, "# PROVENANCE: generated ONE-SHOT (2026-07-02, scratchpad/s4c_capture_fixtures.jl) from the LAST")
    println(io, "# Wigner-coupled committed Kernels (51efa9f: r2_N! evaluated Wigner.coupling host-side): each entry")
    println(io, "# is exactly the c-vector the committed kernel computed internally for (n, branch), i.e.")
    println(io, "# coupling(l3,m3,s3, l1,m1,s1, l2,m2,s2) over the (s3,s1,s2) spec PARSED from the committed r2.jl")
    println(io, "# (grounded on PRE-change behavior; the spec is also baked so Kernels.R2_COUPLING_SPEC is gated")
    println(io, "# against transcription drift). Kernels' r2 tests must stay Wigner-free (design §5.5 purity:")
    println(io, "# deps == {NumericsCore, KernelAbstractions}), so the couplings are baked, not recomputed.")
    println(io, "# Regenerate against grtools/core/Wigner if R2_BRANCHES changes.")
    println(io)
    println(io, "const R2_COUPLING_SPEC_REF = (")
    for n in 1:10
        println(io, "    (", join(["(" * join(t, ",") * ")" for t in spec[n]], ","), "),")
    end
    println(io, ")")
    println(io)
    println(io, "const R2_COUPLING_REF = [")
    for n in 1:10
        println(io, "    [ # n = $n")
        for b in branches
            c = Float64[coupling(b.l3, b.m3, s3, b.l1, b.m1, s1, b.l2, b.m2, s2) for (s3, s1, s2) in spec[n]]
            println(io, "        Float64[", join(repr.(c), ", "), "],")
        end
        println(io, "    ],")
    end
    println(io, "]")
end
println("wrote r2_coupling_ref.jl")

# ---------- (2) aak bessel fixture ----------
# (a) _aak_point test: e=0.3, nmodes=3
pairs = Tuple{Int,Float64}[(n, 0.3) for n in 1:3]
# (b) aak_summation test: replicate the committed aak_summation stacking + spline exactly.
N = 5
t = collect(0.0:10.0:10.0 * (N - 1))
e = [0.3 + 0.001 * i for i in 1:N]
Phivec = [0.1 * i for i in 1:N]
gimvec = [0.2 + 0.01 * i for i in 1:N]
alpvec = [0.4 + 0.02 * i for i in 1:N]
nuvec = [1e-3 + 1e-5 * i for i in 1:N]
gimdotvec = [1e-4 * i for i in 1:N]
OmegaPhi = [2e-3 + 1e-5 * i for i in 1:N]
lam = [0.6 + 0.01 * i for i in 1:N]
Y = Matrix{Float64}(undef, 8, N)
Y[1, :] .= e;        Y[2, :] .= Phivec
Y[3, :] .= gimvec;   Y[4, :] .= alpvec
Y[5, :] .= nuvec;    Y[6, :] .= gimdotvec
Y[7, :] .= OmegaPhi; Y[8, :] .= lam
spline = CubicSplineInterpolant(collect(Float64, t), Y)
for i in 1:4                       # out_len = 4, dt = 10.0
    td = 10.0 * (i - 1)
    ed = spline(td, 1)
    for n in 1:3                   # nmodes = 3
        push!(pairs, (n, ed))
    end
end
open(joinpath(FIXTURE_DIR, "aak_bessel_ref.jl"), "w") do io
    println(io, "# Baked Peters-Mathews Bessel 5-tuples (J_{n-2..n+2}(n*e), J_{-1}=-J_1 at n=1) for the aak tests.")
    println(io, "# PROVENANCE: generated ONE-SHOT (2026-07-02, scratchpad/s4c_capture_fixtures.jl) from the LAST")
    println(io, "# SpecialFunctions-coupled committed Kernels (51efa9f state): values are the committed")
    println(io, "# Kernels.ak_bessel_coeffs(n,e) (SpecialFunctions.besselj) at every (n,e) the tests drive:")
    println(io, "#   * _aak_point fixture: e=0.3, n=1:3;")
    println(io, "#   * aak_summation fixture: n=1:3 at the dense-node eccentricities spline(td,1), td=0,10,20,30,")
    println(io, "#     reproduced by the SAME CubicSplineInterpolant stacking aak_summation performs, so the Float64")
    println(io, "#     keys match the kernel's provider requests BIT-EXACTLY (a provider KeyError = input drift).")
    println(io, "# Kernels' aak tests must stay SpecialFunctions-free (design §5.5 purity). ak_bessel_coeffs now")
    println(io, "# lives in EMRIWaveforms.AAK (its besselj grounding moved there too).")
    println(io)
    println(io, "const AAK_BESSEL_REF = Dict{Tuple{Int,Float64},NTuple{5,Float64}}(")
    for (n, ev) in sort(unique(pairs))
        J = ak_bessel_coeffs(n, ev)
        println(io, "    (", n, ", ", repr(ev), ") => (", join(repr.(J), ", "), "),")
    end
    println(io, ")")
end
println("wrote aak_bessel_ref.jl")

# ---------- (3) pre-hoist EMRIWaveforms aak_summation behavior literals ----------
te = [0.0, 100.0, 200.0, 300.0, 400.0]
ee = [0.3, 0.29, 0.28, 0.27, 0.26]
Ph = [0.0, 1.0, 2.0, 3.0, 4.0]
gi = [0.0, 0.5, 1.0, 1.5, 2.0]
al = [0.0, 0.2, 0.4, 0.6, 0.8]
nu = [2.0e-3, 2.01e-3, 2.02e-3, 2.03e-3, 2.04e-3]
gd = fill(1.0e-5, 5)
Om = [1.3e-2, 1.31e-2, 1.32e-2, 1.33e-2, 1.34e-2]
la = fill(0.6, 5)
h = aak_summation(te, ee, Ph, gi, al, nu, gd, Om, la,
                  1.0e6, 0.5, 10.0, 0.7, 0.8, 0.9, 1.0, 1.5, 4; dt = 100.0, out_len = 5)
println("PRE-HOIST EMRIWaveforms aak_summation literals (bake into EMRIWaveforms test_aak.jl):")
for x in h
    println("            ", repr(x), ",")
end
