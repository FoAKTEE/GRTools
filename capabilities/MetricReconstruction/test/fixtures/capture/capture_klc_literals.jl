# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/cap12.jl (untracked one-shot, mtime 2026-06-30)
#               original scratchpad name: cap12.jl
#               sha256(original) = ab8fd23809f6e9ca44b191ae3daa1011d94ae16b541fc732b995fb6cac9b75a8
# fixture(s)  : (stdout literals baked in test/lorenz_background.jl, lorenz_operators.jl)
# oracle      : legacy packages/KerrLorenzCirc @ b24f8356
# run         : julia --project=$GRTOOLS_LEGACY_PACKAGES_ROOT/KerrLorenzCirc test/fixtures/capture/capture_klc_literals.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
# ----------------------------------------------------------------------------------------------

using KerrLorenzCirc
const K = KerrLorenzCirc
# --- Module 1: angular CG mixing + jet algebra fixtures at sample (a,r0,ω,m) ---
m = 2
println("# frak fixtures (s1,s2,l,lp,m) and (s,l,lp,m)")
for (s1,s2,l,lp) in ((1,0,2,3),(-1,0,3,2),(2,1,2,2),(-2,-1,3,3),(-1,1,2,3))
    println("frak_s($s1,$s2,$l,$lp,$m) = ", repr(K.frak_s(s1,s2,l,lp,m)))
end
for (s1,s2,l,lp) in ((2,0,2,4),(-2,0,4,2))
    println("frak_s2($s1,$s2,$l,$lp,$m) = ", repr(K.frak_s2(s1,s2,l,lp,m)))
end
for (s,l,lp) in ((2,2,3),(0,3,4),(-1,4,5))
    println("frak_c($s,$l,$lp,$m) = ", repr(K.frak_c(s,l,lp,m)))
end
for (s,l,lp) in ((0,2,2),(0,2,4),(2,3,3))
    println("frak_c2($s,$l,$lp,$m) = ", repr(K.frak_c2(s,l,lp,m)))
end
# jet_ops sample at (r0,a,ω,m)
r0=6.0; a=0.5; ω=0.4
ops = K.jet_ops(r0,a,ω,m; N=6)
println("# jet_ops(6.0,0.5,0.4,2;N=6)")
println("Delta = ", repr(ops.Δjet))
println("Kjet = ", repr(ops.Kjet))
println("invDelta = ", repr(ops.invΔjet))
println("cjet = ", repr(ops.cjet))
# jet mul/inv sample
f = ComplexF64[1.0+0im, 2.0, -1.0, 0.5, 0,0,0]
g = ComplexF64[2.0+1im, -1.0, 0.3, 0,0,0,0]
println("jmul_f_g = ", repr(K.jmul(f,g)))
println("jinv_f = ", repr(K.jinv(f)))
# --- Module 2: orbit_background calA/calB at a=0.99 r0=6 m=2 ---
bg = K.orbit_background(0.99, 6.0, 2)
println("# orbit_background(0.99,6.0,2)")
for fld in (:E,:L,:Ω,:ut,:Δ0,:calA,:calB,:ω,:aω)
    println("bg.$fld = ", repr(getfield(bg, fld)))
end
# BigFloat check
setprecision(BigFloat,256) do
    bgb = K.orbit_background(big"0.99", big"6.0", 2)
    println("bg_big.calA typeof = ", typeof(bgb.calA))
end
