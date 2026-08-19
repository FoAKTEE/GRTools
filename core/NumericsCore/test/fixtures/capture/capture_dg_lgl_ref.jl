# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/capture_dg_lgl_ref.jl (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = de7990201b2167662412ff55a401d705c82a95dc79ca646f4cdfb89113f00d1d
# fixture(s)  : dg_lgl_ref.jl
# oracle      : /data/haiyangw/claude/DG-Teulolsky/src/DGCore.jl @ 0a88967 (read-only source repo)
# run         : julia --project=$DGTEUK_ROOT test/fixtures/capture/capture_dg_lgl_ref.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const DGTEUK_ROOT          = get(ENV, "GRTOOLS_DGTEUK_ROOT",          "/data/haiyangw/claude/DG-Teulolsky")
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

# One-shot characterization capture for DG-P1 (run ONCE with the ORIGINAL repo's env):
#
#   cd /data/haiyangw/claude/DG-Teulolsky && julia --project=. \
#       /data/haiyangw/claude/emri-disk/scratchpad/capture_dg_lgl_ref.jl
#
# Includes the ORIGINAL /data/haiyangw/claude/DG-Teulolsky/src/DGCore.jl (which uses
# FastGaussQuadrature.gausslobatto) and bakes RefElement(N) reference values for N=2..16 into
# grtools/core/NumericsCore/test/fixtures/dg_lgl_ref.jl:
#   - nodes r and weights w for ALL N=2..16 (full-precision Float64 round-trip literals),
#   - full D and M matrices for the representative subset N in (2, 5, 9, 16),
#   - Frobenius norms of D and M for ALL N=2..16 (cheap whole-matrix hashes).
# The DG-Teulolsky repo is READ-ONLY source; this script only reads it and writes the fixture.

include("$(DGTEUK_ROOT)/src/DGCore.jl")
using .DGCore
using LinearAlgebra: norm

const NS = 2:16
const NS_FULL = (2, 5, 9, 16)
const OUT = "$(FIXTURE_DIR)/dg_lgl_ref.jl"

veclit(v) = "[" * join((repr(x) for x in v), ", ") * "]"
matlit(A) = "[" * join((join((repr(A[i, j]) for j in 1:size(A, 2)), " ") for i in 1:size(A, 1)), "; ") * "]"

open(OUT, "w") do io
    println(io, "# AUTO-GENERATED one-shot fixture — DO NOT EDIT.")
    println(io, "# Baked ORIGINAL DG-Teulolsky DGCore.RefElement(N) reference values (N = 2..16), captured once")
    println(io, "# by scratchpad/capture_dg_lgl_ref.jl running $(DGTEUK_ROOT)/src/DGCore.jl")
    println(io, "# in that repo's own environment (LGL nodes/weights from FastGaussQuadrature.gausslobatto;")
    println(io, "# D = Vr/V, M = inv(V*V'), orthonormal-Legendre Vandermonde). Full-precision Float64 repr")
    println(io, "# literals (exact round-trip). Characterizes NumericsCore's self-contained Newton-on-P'_N LGL")
    println(io, "# rule + RefElement against the FastGaussQuadrature-based original at <=1e-14 rel.")
    println(io, "#   - DG_LGL_REF[N]: nodes r, weights w for every N in 2..16")
    println(io, "#   - DG_DM_REF[N]:  full D and M matrices for N in $(NS_FULL)")
    println(io, "#   - DG_DM_FRO[N]:  (Dfro=norm(D), Mfro=norm(M)) Frobenius hashes for every N in 2..16")
    println(io)
    println(io, "const DG_LGL_REF = Dict{Int,NamedTuple{(:r, :w),Tuple{Vector{Float64},Vector{Float64}}}}(")
    for N in NS
        el = RefElement(N)
        println(io, "    $N => (r = $(veclit(el.r)),")
        println(io, "          w = $(veclit(el.w))),")
    end
    println(io, ")")
    println(io)
    println(io, "const DG_DM_REF = Dict{Int,NamedTuple{(:D, :M),Tuple{Matrix{Float64},Matrix{Float64}}}}(")
    for N in NS_FULL
        el = RefElement(N)
        println(io, "    $N => (D = $(matlit(el.D)),")
        println(io, "          M = $(matlit(el.M))),")
    end
    println(io, ")")
    println(io)
    println(io, "const DG_DM_FRO = Dict{Int,NTuple{2,Float64}}(")
    for N in NS
        el = RefElement(N)
        println(io, "    $N => ($(repr(norm(el.D))), $(repr(norm(el.M)))),")
    end
    println(io, ")")
end
println("wrote $OUT")
