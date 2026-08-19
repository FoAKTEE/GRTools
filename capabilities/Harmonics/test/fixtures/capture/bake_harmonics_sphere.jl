# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/npghp_moveD/TeukEvolution_run/minrun/bake_harmonics_sphere.jl (untracked one-shot, mtime 2026-07-04)
#               sha256(original) = 4357d0372e2552e0e9865a448ca5ed1b4010d57a6c19cd6db8959ff425910ad2
# fixture(s)  : sphere_collocation_ref.jl
# oracle      : ref-code/TeukEvolution.jl/src/Sphere.jl @ 3c8e47f + FastGaussQuadrature.gausslegendre(24)
# run         : julia --project=<env with FastGaussQuadrature> test/fixtures/capture/bake_harmonics_sphere.jl
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

# Move-D D2 bake: characterization fixtures for grtools Harmonics SphereCollocation (swal + GL
# collocation matrices). Runs the ORIGINAL TeukEvolution.jl/src/Sphere.jl (swal_vals,
# swal_laplacian_matrix, swal_raising_matrix, swal_lowering_matrix, swal_killtop_matrix(...,2),
# swal_filter_matrix, Y_vals/cos_vals/sin_vals + GL weights) in the fixture-faithful upstream env
# and writes full-precision Float64 literals into the Harmonics test fixture. Never modifies src/.
#
# Run:  julia --project=. minrun/bake_harmonics_sphere.jl   (from TeukEvolution_run/)

cd(TEUK_ROOT)
include(joinpath(TEUK_ROOT, "src", "Sphere.jl"))
using .Sphere
import FastGaussQuadrature

const OUT = "$(FIXTURE_DIR)/sphere_collocation_ref.jl"

const NY = 24
const SM = [(-2, 2), (-2, -2), (0, 2), (2, 2)]   # production spins (s,m)
const KILLTOP_N = 2

# GL grid + weights (upstream FastGaussQuadrature path)
gl_nodes, gl_weights = FastGaussQuadrature.gausslegendre(NY)
yv   = Sphere.Y_vals(NY)      # = roots
cosv = Sphere.cos_vals(NY)    # = -roots = cosθ
sinv = Sphere.sin_vals(NY)    # = sinθ

# --- literal writers: full-precision Float64 repr (exact round-trip) ---
mat_lit(M) = "[" * join([join([repr(M[i, j]) for j in 1:size(M, 2)], " ") for i in 1:size(M, 1)], ";\n     ") * "]"
vec_lit(x) = "[" * join([repr(v) for v in x], ", ") * "]"

open(OUT, "w") do io
    println(io, "# AUTO-GENERATED one-shot fixture — DO NOT EDIT.")
    println(io, "# Baked from the ORIGINAL TeukEvolution.jl (ref-code) Sphere.jl at ny=$NY:")
    println(io, "#   Y_vals/cos_vals/sin_vals + FastGaussQuadrature.gausslegendre($NY) weights, and per")
    println(io, "#   (s,m) the swal_vals / laplacian / raising / lowering / killtop(N=$KILLTOP_N) / filter matrices.")
    println(io, "# Captured once by scratchpad/npghp_moveD/TeukEvolution_run/minrun/bake_harmonics_sphere.jl")
    println(io, "# in the upstream env (julia $(VERSION), Jacobi + FastGaussQuadrature). Full-precision")
    println(io, "# Float64 repr literals. Characterizes Harmonics SphereCollocation (swal adapter + matrices).")
    println(io)
    println(io, "const SPHERE_NY = $NY")
    println(io, "const SPHERE_KILLTOP_N = $KILLTOP_N")
    println(io, "const SPHERE_YV  = ", vec_lit(yv))
    println(io, "const SPHERE_COS = ", vec_lit(cosv))
    println(io, "const SPHERE_SIN = ", vec_lit(sinv))
    println(io, "const SPHERE_W   = ", vec_lit(gl_weights))
    println(io)
    println(io, "# SPHERE_MATS[(s,m)] = (vals, lap, raise, lower, killtop2, filt)")
    println(io, "const SPHERE_MATS = Dict{Tuple{Int,Int},NamedTuple}(")
    for (s, m) in SM
        vals    = Sphere.swal_vals(NY, s, m)
        lap     = Sphere.swal_laplacian_matrix(NY, s, m)
        raise   = Sphere.swal_raising_matrix(NY, s, m)
        lower   = Sphere.swal_lowering_matrix(NY, s, m)
        killtop = Sphere.swal_killtop_matrix(NY, s, m, KILLTOP_N)
        filt    = Sphere.swal_filter_matrix(NY, s, m)
        println(io, "  ($s, $m) => (")
        println(io, "    vals = ", mat_lit(vals), ",")
        println(io, "    lap = ", mat_lit(lap), ",")
        println(io, "    raise = ", mat_lit(raise), ",")
        println(io, "    lower = ", mat_lit(lower), ",")
        println(io, "    killtop2 = ", mat_lit(killtop), ",")
        println(io, "    filt = ", mat_lit(filt), ",")
        println(io, "  ),")
    end
    println(io, ")")
end

# console sanity
println("ny=$NY  num_l=", Sphere.num_l(NY))
println("yv[1]=", yv[1], "  cos[1]=", cosv[1], "  sin[1]=", sinv[1], "  w[1]=", gl_weights[1])
for (s, m) in SM
    v = Sphere.swal_vals(NY, s, m)
    println("(s=$s,m=$m) swal_vals size=", size(v), "  vals[1,1]=", v[1, 1])
end
println("WROTE ", OUT)
