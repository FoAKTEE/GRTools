# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/capture_d2hdr2_oracle.jl (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = 46f31b113a62d4fc9123dc078ef4a8411708a5a909b08cfd8e748cec29bbab14
# fixture(s)  : d2hdr2_oracle.jl
# oracle      : ref-code/SecondOrderRicci/src/h1.cc L59-122 @ abeb8c0 compiled with g++ -O2 (Boost multi_array)
# run         : julia --project=capabilities/SelfForce test/fixtures/capture/capture_d2hdr2_oracle.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const REFCODE_ROOT         = get(ENV, "GRTOOLS_REFCODE_ROOT",         "/data/haiyangw/claude/emri-disk/ref-code")
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

# ONE-SHOT capture script (S5-3): bake the live g++ d2hdr2 oracle of
# grtools/capabilities/SelfForce/test/second_order_analytic.jl into a fixture.
#
# Reproduces EXACTLY what the live test did:
#   * slices d2hdr2 verbatim from ref-code/SecondOrderRicci/src/h1.cc L59-122 into a
#     Boost-multi_array g++ TU (same wrapper source, same `g++ -O2 -I/usr/include`),
#   * generates the SAME 40 input cases via MersenneTwister(20240617) with the identical
#     call sequence (i=1..10 × (l,m) ∈ {(2,0),(3,1),(4,2),(5,-2)}, r0=8.0),
#   * runs the oracle binary once, parses the (re,im) outputs,
# then writes test/fixtures/d2hdr2_oracle.jl with a full provenance header and the
# inputs+outputs baked as round-trip-exact Float64 literals.
#
# This script is the ONLY thing that touches ref-code/; the committed test includes the
# fixture only. Run from anywhere: julia --project=<SelfForce> capture_d2hdr2_oracle.jl
# (needs only Random + SHA stdlib and a g++/boost toolchain, TODAY, once).

using Random, SHA, Dates

const H1CC = "$(REFCODE_ROOT)/SecondOrderRicci/src/h1.cc"
const FIXTURE = "$(FIXTURE_DIR)/d2hdr2_oracle.jl"

lines = readlines(H1CC)
# The live test's "source slice anchor" assertions — enforced HERE at capture time:
@assert occursin("d2hdr2", lines[59]) "h1.cc line-drift: L59 lost d2hdr2 signature"
@assert occursin("return ddhbar", lines[121]) "h1.cc line-drift: L121 lost return ddhbar"
body = join(lines[59:122], "\n")

# ── build the oracle exactly as the live test did ────────────────────────────────────
src = """
#include <boost/multi_array.hpp>
#include <complex>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;
using boost::multi_array;
typedef boost::multi_array_types::extent_range range;
static const double M = 1.0;

$body

int main() {
  cout << setprecision(17);
  int i, l, m; double r, f, fp, r0;
  while (cin >> i >> l >> m >> r >> f >> fp >> r0) {
    multi_array<complex<double>,1> hbar(boost::extents[range(1,11)]);
    multi_array<complex<double>,1> dhbar(boost::extents[range(1,11)]);
    for (int k=1;k<=10;++k){double re,im; cin>>re>>im; hbar[k]=complex<double>(re,im);}
    for (int k=1;k<=10;++k){double re,im; cin>>re>>im; dhbar[k]=complex<double>(re,im);}
    complex<double> out = d2hdr2(i,l,m,r,f,fp,r0,hbar,dhbar);
    cout << out.real() << " " << out.imag() << "\\n";
  }
  return 0;
}
"""
cpp = tempname() * ".cpp"; write(cpp, src)
bin = tempname()
run(`g++ -O2 -I/usr/include -o $bin $cpp`)

# ── generate the SAME inputs as the live test (identical RNG call sequence) ──────────
rng = MersenneTwister(20240617)
seeds = Vector{Tuple{Int,Int,Int,Float64,Float64,Float64,Vector{ComplexF64},Vector{ComplexF64}}}()
io = IOBuffer()
r0 = 8.0
for i in 1:10, (l, m) in ((2, 0), (3, 1), (4, 2), (5, -2))
    r = 4.0 + 3.0 * rand(rng)
    f = 1 - 2 / r; fp = 2 / r^2
    hb = [complex(2rand(rng) - 1, 2rand(rng) - 1) for _ in 1:10]
    db = [complex(2rand(rng) - 1, 2rand(rng) - 1) for _ in 1:10]
    push!(seeds, (i, l, m, r, f, fp, hb, db))
    print(io, "$i $l $m $r $f $fp $r0")
    for c in hb; print(io, " $(real(c)) $(imag(c))"); end
    for c in db; print(io, " $(real(c)) $(imag(c))"); end
    print(io, "\n")
end
out = read(pipeline(IOBuffer(String(take!(io))), `$bin`), String)
rows = [parse.(Float64, split(line)) for line in split(strip(out), "\n")]
@assert length(rows) == length(seeds) == 40

# ── provenance ───────────────────────────────────────────────────────────────────────
gxx = first(split(read(`g++ --version`, String), "\n"))
boost_ver = let m = match(r"BOOST_LIB_VERSION \"([^\"]+)\"", read("/usr/include/boost/version.hpp", String))
    m === nothing ? "unknown" : m.captures[1]
end
h1cc_sha = bytes2hex(sha256(read(H1CC)))
slice_sha = bytes2hex(sha256(body))
git_state = strip(read(setenv(`git log -1 "--format=%H %ci" -- src/h1.cc`;
                              dir="$(REFCODE_ROOT)/SecondOrderRicci"), String))

# ── write fixture ────────────────────────────────────────────────────────────────────
cplx(c) = "complex($(repr(real(c))), $(repr(imag(c))))"
open(FIXTURE, "w") do fio
    println(fio, """
# BAKED ORACLE FIXTURE — d2hdr2 (second-order analytic first-order fields).
# Generated ONE-SHOT by scratchpad/capture_d2hdr2_oracle.jl on $(Dates.today()).
#
# Provenance:
#   legacy fn : d2hdr2(i,l,m,r,f,fp,r0,hbar,dhbar), sliced VERBATIM from
#               ref-code/SecondOrderRicci/src/h1.cc lines 59-122
#   h1.cc git : $git_state
#   h1.cc sha256      : $h1cc_sha
#   slice(59-122) sha256: $slice_sha
#   compiler  : $gxx
#   boost     : $boost_ver (BOOST_LIB_VERSION), header-only multi_array; wrapped in a tiny
#               TU (no HDF5/GSL) with M=1.0, cout<<setprecision(17), invoked as
#               `g++ -O2 -I/usr/include -o bin slice.cpp`
#   inputs    : the EXACT 40 cases of the pre-fixture live test — MersenneTwister(20240617),
#               i=1..10 × (l,m) ∈ {(2,0),(3,1),(4,2),(5,-2)}, r ∈ 4+3*rand, f=1-2/r,
#               fp=2/r², r0=8.0, hbar/dhbar ∈ ([-1,1]+[-1,1]im)^10 (ALL i hit; m≠0 hits
#               dt=-imΩ and Ω²m²)
#   tolerance : rtol = 1e-12 (UNCHANGED from the live test; setprecision(17) round-trips
#               Float64 exactly, so the only slack needed is port-vs-g++ FP ordering)
#
# Slice anchors verified at capture time (h1.cc line-drift guard of the old live test):""")
    println(fio, "const D2HDR2_SLICE_FIRST_LINE = ", repr(lines[59]))
    println(fio, "const D2HDR2_SLICE_LAST_LINE = ", repr(lines[121]))
    println(fio, "\nconst D2HDR2_ORACLE_R0 = 8.0")
    println(fio, "\n# (i, l, m, r, f, fp, hbar[1:10], dhbar[1:10], oracle d2hdr2)")
    println(fio, "const D2HDR2_ORACLE_CASES = [")
    for (k, (i, l, m, r, f, fp, hb, db)) in enumerate(seeds)
        println(fio, "    ($i, $l, $m, $(repr(r)), $(repr(f)), $(repr(fp)),")
        println(fio, "     ComplexF64[", join(cplx.(hb), ", "), "],")
        println(fio, "     ComplexF64[", join(cplx.(db), ", "), "],")
        println(fio, "     ", cplx(complex(rows[k][1], rows[k][2])), "),")
    end
    println(fio, "]")
end
println("wrote $FIXTURE with $(length(seeds)) cases")
