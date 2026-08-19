# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/s4b_capture_hdf5_fixture.jl (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = bbb58e0d7c95bce483d63d17069357632268cd1df1db76c4dedbaaab4771c7ba
# fixture(s)  : hdf5_readers_ref.jl
# oracle      : PRE-MOVE EMRIWaveforms.AmplitudeModels._load_roman_nn / KerrAmpSpline.read_kerrampgrid + FEW corpus ref-code/FastEMRIWaveforms/.fewdata @ ce6e56ae
# run         : julia --project=capabilities/EMRIWaveforms -t4 test/fixtures/capture/s4b_capture_hdf5_fixture.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const EMRIDISK_ROOT        = get(ENV, "GRTOOLS_EMRIDISK_ROOT",        "/data/haiyangw/claude/emri-disk")
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

# S4-b ONE-SHOT PRE-CHANGE fixture capture (tests-first grounding).
# Runs against the COMMITTED EMRIWaveforms HDF5 loader bodies BEFORE any S4-b code move:
#   (1) EMRIWaveforms.AmplitudeModels._load_roman_nn          (src/AmplitudeModels.jl:55, ROMAN NN weights)
#   (2) EMRIWaveforms.AmplitudeModels.KerrAmpSpline.read_kerrampgrid (src/KerrAmpSpline.jl:178, ZNAmps grid)
#   (3) KerrAmpSpline.kerr_mode_table                          (moves to Data alongside the loader)
# Output: grtools/core/Data/test/fixtures/hdf5_readers_ref.jl — repr-round-trip literals that the NEW
# Data.read_roman_nn / Data.read_kerrampgrid must reproduce BIT-EXACTLY (pure file reads, Float64).
#
# Run: julia --project=grtools/capabilities/EMRIWaveforms -t4 scratchpad/s4b_capture_hdf5_fixture.jl

using EMRIWaveforms
using EMRIWaveforms.AmplitudeModels: _load_roman_nn
using EMRIWaveforms.AmplitudeModels.KerrAmpSpline: read_kerrampgrid, kerr_mode_table

const ROOT = EMRIDISK_ROOT
const FEWDATA = joinpath(ROOT, "ref-code", "FastEMRIWaveforms", ".fewdata")
const NN_FILE = joinpath(FEWDATA, "SchwarzschildEccentricInput.hdf5")
const ZNAMPS = joinpath(FEWDATA, "ZNAmps_l10_m10_n55_DS2Outer.h5")
const OUT = joinpath(FIXTURE_DIR, "hdf5_readers_ref.jl")

isfile(NN_FILE) || error("missing $NN_FILE")
isfile(ZNAMPS) || error("missing $ZNAMPS")
mkpath(dirname(OUT))

commit = strip(read(`git -C $ROOT rev-parse --short HEAD`, String))

io = IOBuffer()
println(io, "# PRE-CHANGE characterization fixture for the S4-b HDF5-readers-to-Data-ext move.")
println(io, "# Captured by scratchpad/s4b_capture_hdf5_fixture.jl from the COMMITTED EMRIWaveforms loaders")
println(io, "# (_load_roman_nn @ src/AmplitudeModels.jl:55; KerrAmpSpline.read_kerrampgrid @ src/KerrAmpSpline.jl:178;")
println(io, "# kerr_mode_table @ src/KerrAmpSpline.jl:137) at repo commit $commit, BEFORE any code move.")
println(io, "# repr-round-trip literals: the Data-ext readers must reproduce these BIT-EXACTLY (pure Float64 file reads).")
println(io)

# ---------------- (1) ROMAN NN loader ----------------
weights, biases, rb, tfi, bi, ntm = _load_roman_nn(NN_FILE)
nl = length(weights)
println(io, "const ROMAN_REF = (")
println(io, "    num_layers = $nl,")
println(io, "    weight_sizes = $(repr([size(w) for w in weights])),")
println(io, "    bias_lengths = $(repr([length(b) for b in biases])),")
println(io, "    rb_size = $(repr(size(rb))),")
println(io, "    transform_factor_inv = $(repr(tfi)),")
println(io, "    break_index = $bi,")
println(io, "    num_teuk_modes = $ntm,")
# weight probes: (layer, i, j, value) — corners + interior of every layer
wp = Tuple{Int,Int,Int,Float64}[]
for (L, w) in enumerate(weights)
    nr, nc = size(w)
    for (i, j) in ((1, 1), (nr, nc), (max(1, nr ÷ 2), max(1, nc ÷ 2)))
        push!(wp, (L, i, j, w[i, j]))
    end
end
println(io, "    weight_probes = $(repr(wp)),")
bp = Tuple{Int,Int,Float64}[]
for (L, b) in enumerate(biases)
    n = length(b)
    for i in (1, n, max(1, n ÷ 2))
        push!(bp, (L, i, b[i]))
    end
end
println(io, "    bias_probes = $(repr(bp)),")
rp = Tuple{Int,Int,ComplexF64}[]
for (i, j) in ((1, 1), (size(rb, 1), size(rb, 2)), (size(rb, 1) ÷ 2, size(rb, 2) ÷ 2), (10, 3))
    push!(rp, (i, j, rb[i, j]))
end
println(io, "    rb_probes = $(repr(rp)),")
println(io, ")")
println(io)

# ---------------- (2) kerr_mode_table ----------------
modes_all, index = kerr_mode_table(; lmax = 10, nmax = 55)
println(io, "const KERR_MODE_TABLE_REF = (")
println(io, "    n_modes = $(length(modes_all)),")
println(io, "    first = $(repr(modes_all[1])),")
println(io, "    m550 = $(repr(modes_all[550])),")
println(io, "    last = $(repr(modes_all[end])),")
println(io, "    idx_220 = $(index[(2, 2, 0)]),")
println(io, "    idx_321 = $(index[(3, 2, 1)]),")
println(io, "    idx_40m3 = $(index[(4, 0, -3)]),")
println(io, "    idx_101055 = $(index[(10, 10, 55)]),")
println(io, ")")
println(io)

# ---------------- (3) ZNAmps read_kerrampgrid (same mode subset the EMRIW test drives) ----------------
req = [(2, 2, 0), (3, 2, 1), (4, 0, -3)]
g = read_kerrampgrid(ZNAMPS; modes = req)
println(io, "const KERRAMP_REF = (")
println(io, "    modes = $(repr(g.modes)),")
println(io, "    lmax = $(g.lmax), mmax = $(g.mmax), nmax = $(g.nmax),")
println(io, "    num_total_modes = $(g.num_total_modes),")
println(io, "    z_knots = $(repr(g.z_knots)),")
println(io, "    twA = $(repr(g.twA)),")
println(io, "    tuA = $(repr(g.tuA)),")
println(io, "    twB = $(repr(g.twB)),")
println(io, "    tuB = $(repr(g.tuB)),")
println(io, "    sizeA = $(repr(size(g.coeffsA))),")
println(io, "    sizeB = $(repr(size(g.coeffsB))),")
cpA = Tuple{Int,Int,Int,Int,Float64}[]
for idx in ((1, 1, 1, 1), (545, 2, 2, 17), (1089, 1, 3, 33), (700, 2, 1, 5))
    push!(cpA, (idx..., g.coeffsA[idx...]))
end
println(io, "    coeffsA_probes = $(repr(cpA)),")
cpB = Tuple{Int,Int,Int,Int,Float64}[]
for idx in ((1, 1, 1, 1), (144, 2, 2, 17), (289, 1, 3, 33), (200, 2, 1, 5))
    push!(cpB, (idx..., g.coeffsB[idx...]))
end
println(io, "    coeffsB_probes = $(repr(cpB)),")
println(io, ")")

write(OUT, String(take!(io)))
println("WROTE $OUT")
