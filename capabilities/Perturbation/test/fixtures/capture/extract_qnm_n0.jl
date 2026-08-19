# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/npghp_phase2/extract_qnm_n0.jl (untracked one-shot, mtime 2026-07-05)
#               sha256(original) = 11f815e4a3b1e3cef002e1cf2654f929b632b46bb5c51a897ce771b9703a4e3a
# fixture(s)  : td_qnm_a0_l2_m2_n0.jl
# oracle      : ref-code/TeukEvolution.jl/qnm/a0.0_l2_m2.h5 @ 3c8e47f (group [n=0]); the final .jl bake step (f6_gen_fixture.jl) is LOST
# run         : julia --project=<env with HDF5> test/fixtures/capture/extract_qnm_n0.jl
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

import HDF5
import SHA
import Printf: @printf

const H5 = "$(REFCODE_ROOT)/TeukEvolution.jl/qnm/a0.0_l2_m2.h5"
const KEY = "[n=0]"

# sha256 of the source file
sha = bytes2hex(SHA.sha256(read(H5)))

# read the n=0 group via HDF5.jl (same call Id.jl:137 uses: h5read(file, "[n=0]"))
g = HDF5.h5read(H5, KEY)

println("=== keys in group ", KEY, " ===")
for k in sort(collect(keys(g)))
    v = g[k]
    println("  ", rpad(k, 14), " :: ", typeof(v), "  size=", (isa(v, AbstractArray) ? size(v) : "scalar"))
end

om   = g["omega"]
lam  = g["lambda"]
nr   = g["nr"]
rc   = g["radial_coef"]
ac   = g["angular_coef"]
rf   = g["radial_func"]
rv   = g["rvals"]
amp  = g["amp"]
damp = g["damp"]

println("\n=== scalars ===")
println("omega         = ", repr(om))
println("lambda        = ", repr(lam))
println("nr            = ", repr(nr))

println("\n=== array lengths ===")
println("radial_coef   len=", length(rc), " eltype=", eltype(rc))
println("angular_coef  len=", length(ac), " eltype=", eltype(ac))
println("radial_func   len=", length(rf), " eltype=", eltype(rf))
println("rvals         len=", length(rv), " eltype=", eltype(rv))
println("amp           len=", length(amp), " eltype=", eltype(amp))
println("damp          len=", length(damp), " eltype=", eltype(damp))

println("\n=== radial_coef first 5 ===")
for i in 1:min(5,length(rc)); println("  [", i, "] ", repr(rc[i])); end
println("=== radial_coef last 3 ===")
for i in max(1,length(rc)-2):length(rc); println("  [", i, "] ", repr(rc[i])); end

println("\n=== angular_coef (all ", length(ac), ") ===")
for i in 1:length(ac); println("  [", i, "] ", repr(ac[i])); end

println("\n=== rvals first 3 / last 3 ===")
println("  first: ", rv[1], ", ", rv[2], ", ", rv[3])
println("  last:  ", rv[end-2], ", ", rv[end-1], ", ", rv[end])

# helper: full-precision complex literal
fpc(z) = string("(", repr(real(z)), ", ", repr(imag(z)), ")")

# ----- write full-precision plain-text extraction -----
OUT = get(ENV, "GRTOOLS_QNM_TEXT_OUT", joinpath(FIXTURE_DIR, "capture", "qnm_a0_l2_m2_n0.txt"))
open(OUT, "w") do io
    println(io, "# QNM initial data extraction — TeukEvolution.jl qnm/a0.0_l2_m2.h5")
    println(io, "# PROVENANCE")
    println(io, "#   source h5    : ref-code/TeukEvolution.jl/qnm/a0.0_l2_m2.h5")
    println(io, "#   sha256       : ", sha)
    println(io, "#   group key    : ", KEY)
    println(io, "#   generator    : TeukolskyQNMFunctions.jl (JLRipley314), per qnm/README.md")
    println(io, "#   extractor    : scratchpad/npghp_phase2/extract_qnm_n0.jl (HDF5.jl h5read)")
    println(io, "#   julia        : ", VERSION)
    println(io, "#   config       : a=0.0 (Schwarzschild), s=-2, l=2, m=2, overtone n=0")
    println(io, "#   shapes       : radial_coef[", length(rc), "] angular_coef[", length(ac),
                    "] radial_func[", length(rf), "] rvals[", length(rv), "] amp[", length(amp),
                    "] damp[", length(damp), "]")
    println(io, "#   Cheby domain : radial_coef are ChebyshevT coeffs (Polynomials.jl); evaluated as")
    println(io, "#                  P(x), x = 2*R/maximum(Rv) - 1 in [-1,1]  (Id.jl:138,148); R is the")
    println(io, "#                  compactified evolution radial coord (r = cl^2/R, Id.jl:61).")
    println(io, "#   angular use  : swal reprojection sum_{l} (-1)^l angular_coef[l+1] * swal(s,m,l+lmin,y),")
    println(io, "#                  lmin = max(|s|,|m|) = 2  (Id.jl:140,149-151).")
    println(io, "#   p (time der) : p = f,t = -i*omega*f  (Id.jl:163-168).")
    println(io)
    println(io, "omega_real = ", repr(real(om)))
    println(io, "omega_imag = ", repr(imag(om)))
    println(io, "lambda_real = ", repr(real(lam)))
    println(io, "lambda_imag = ", repr(imag(lam)))
    println(io, "nr = ", repr(nr))
    println(io, "rmax_grid_convention = 2*R/maximum(Rv) - 1   # x in [-1,1]")
    println(io)
    println(io, "# radial_coef : ChebyshevT coefficients (index 0-based degree), full ComplexF64 precision")
    println(io, "# format: index  real  imag")
    println(io, "radial_coef_len = ", length(rc))
    for i in 1:length(rc)
        @printf(io, "radial_coef[%d] = %s %s\n", i-1, repr(real(rc[i])), repr(imag(rc[i])))
    end
    println(io)
    println(io, "# angular_coef : spherical-harmonic reprojection coefficients (l index 0-based, l+lmin physical)")
    println(io, "angular_coef_len = ", length(ac))
    for i in 1:length(ac)
        @printf(io, "angular_coef[%d] = %s %s   # swal l = %d\n", i-1, repr(real(ac[i])), repr(imag(ac[i])), (i-1)+2)
    end
end
println("\nWROTE ", OUT)
println("SHA256(src) = ", sha)
