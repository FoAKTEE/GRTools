# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/npghp_moveD/TeukEvolution_run/minrun/bake_numericscore_fd.jl (untracked one-shot, mtime 2026-07-04)
#               sha256(original) = e503ec0cb85d9f31a1478794883c7c739c3ebb910085375b1503b38edb71feaa
# fixture(s)  : numericscore_fd_uniform_ref.jl
# oracle      : ref-code/TeukEvolution.jl/src/Radial.jl @ 3c8e47f + FastGaussQuadrature.gausslegendre(24)
# run         : julia --project=<env with FastGaussQuadrature> test/fixtures/capture/bake_numericscore_fd.jl
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
# upstream oracle tree (TeukEvolution.jl); the capture ran from a working copy of it
const TEUK_ROOT            = get(ENV, "GRTOOLS_TEUKEVOLUTION_ROOT",   joinpath(REFCODE_ROOT, "TeukEvolution.jl"))
# ----------------------------------------------------------------------------------------------

# Move-D D1 bake: characterization fixtures for grtools NumericsCore fd_uniform.jl + gauss_legendre.jl.
# Runs the ORIGINAL TeukEvolution.jl/src/Radial.jl (set_d1!/set_d2!/filter!) and
# FastGaussQuadrature.gausslegendre(24) in the fixture-faithful upstream environment and writes
# full-precision Float64 literals into NumericsCore's test fixture. Never modifies src/.
#
# Run:  julia --project=. minrun/bake_numericscore_fd.jl   (from TeukEvolution_run/)

cd(TEUK_ROOT)
include(joinpath(TEUK_ROOT, "src", "Radial.jl"))
using .Radial
import FastGaussQuadrature

const OUT = "$(FIXTURE_DIR)/numericscore_fd_uniform_ref.jl"

# --- grid spacing dr from the minrun config (TeukEvolution.jl:63-65; a=0.7, M=1, nx=48) ---
const NX = 48
const NY = 4
bhs = 0.7; bhm = 1.0
minr = bhm * (1 + sqrt(1 - (bhs / bhm)^2))   # uncompactified horizon
maxR = 1 / minr
dr = maxR / (NX - 1)
const EPS_KO = 0.5

# --- seeded pseudo-random field (nx=48, ny=4), literal-int seed, no RNG-from-time ---
# self-contained 64-bit LCG (Knuth MMIX constants); high 53 bits -> [-1,1).
function lcg_field(nx::Int, ny::Int, seed::UInt64)
    a = 0x5851f42d4c957f2d
    c = 0x14057b7ef767814f
    state = seed
    f = Array{Float64}(undef, nx, ny)
    @inbounds for idx in eachindex(f)
        state = a * state + c              # UInt64 wraps mod 2^64
        u = Float64(state >> 11) / Float64(UInt64(1) << 53)   # [0,1)
        f[idx] = 2u - 1
    end
    return f
end

v = lcg_field(NX, NY, UInt64(20240704))

# --- upstream oracle outputs ---
d1 = similar(v); Radial.set_d1!(d1, v, dr)
d2 = similar(v); Radial.set_d2!(d2, v, dr)
filt = copy(v); tmp = similar(v); Radial.filter!(filt, tmp, EPS_KO)

# --- Gauss-Legendre(24) from FastGaussQuadrature ---
gl_nodes, gl_weights = FastGaussQuadrature.gausslegendre(24)

# --- writer: full-precision Float64 literals (exact round-trip via Julia show) ---
mat_lit(M) = "[" * join([join([repr(M[i, j]) for j in 1:size(M, 2)], " ") for i in 1:size(M, 1)], ";\n     ") * "]"
vec_lit(x) = "[" * join([repr(v) for v in x], ", ") * "]"

open(OUT, "w") do io
    println(io, "# AUTO-GENERATED one-shot fixture — DO NOT EDIT.")
    println(io, "# Baked from the ORIGINAL TeukEvolution.jl (ref-code) Radial.jl set_d1!/set_d2!/filter!")
    println(io, "# and FastGaussQuadrature.gausslegendre(24), captured once by")
    println(io, "# scratchpad/npghp_moveD/TeukEvolution_run/minrun/bake_numericscore_fd.jl in the")
    println(io, "# upstream env (julia $(VERSION), FastGaussQuadrature). Full-precision Float64 repr")
    println(io, "# literals (exact round-trip). Characterizes NumericsCore fd_uniform.jl + gauss_legendre.jl.")
    println(io, "#   FD_UNIFORM_REF: seeded (nx=$NX,ny=$NY) field v + Radial oracle outputs d1,d2,filt; dr,eps_KO.")
    println(io, "#   GL24_REF:       FastGaussQuadrature.gausslegendre(24) nodes+weights.")
    println(io)
    println(io, "const FD_UNIFORM_REF = (")
    println(io, "  nx = $NX, ny = $NY,")
    println(io, "  dr = $(repr(dr)),")
    println(io, "  eps_KO = $(repr(EPS_KO)),")
    println(io, "  v = ", mat_lit(v), ",")
    println(io, "  d1 = ", mat_lit(d1), ",")
    println(io, "  d2 = ", mat_lit(d2), ",")
    println(io, "  filt = ", mat_lit(filt), ",")
    println(io, ")")
    println(io)
    println(io, "const GL24_REF = (")
    println(io, "  n = 24,")
    println(io, "  nodes = ", vec_lit(gl_nodes), ",")
    println(io, "  weights = ", vec_lit(gl_weights), ",")
    println(io, ")")
end

println("dr = ", repr(dr))
println("v[1,1]=", v[1,1], " v[48,4]=", v[NX,NY])
println("d1[24,2]=", d1[24,2], " d2[24,2]=", d2[24,2], " filt[24,2]=", filt[24,2])
println("gl_nodes[1]=", gl_nodes[1], " gl_nodes[24]=", gl_nodes[24], " gl_weights[1]=", gl_weights[1])
println("WROTE ", OUT)
