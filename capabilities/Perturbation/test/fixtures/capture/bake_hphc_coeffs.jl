# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/npghp_moveD/TeukEvolution_run/minrun/bake_hphc_coeffs.jl (untracked one-shot, mtime 2026-07-04)
#               sha256(original) = e4643a560b0464f6bbea0a02b160fdb03503cf9ba2f655e9292a244d57ba7b15
# fixture(s)  : (literals in test/test_td_hphc.jl — Evo_lin_f coefficient arrays A_pp/A_pq/B_pp/B_pq/B_pf/pre)
# oracle      : ref-code/TeukEvolution.jl/src/Evolution.jl @ 3c8e47f (a=0.7, s=-2, m=2, nx=48, ny=24, cl=1)
# run         : julia --project=<TeukEvolution env> test/fixtures/capture/bake_hphc_coeffs.jl
#
# Paths rewritten to the ENV-overridable GRTOOLS_TEUKEVOLUTION_ROOT / GRTOOLS_REFCODE_ROOT
# (defaults reproduce the capture host).  Body otherwise VERBATIM from the original one-shot.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const TEUK_ROOT = get(ENV, "GRTOOLS_TEUKEVOLUTION_ROOT",
                      joinpath(get(ENV, "GRTOOLS_REFCODE_ROOT", "/data/haiyangw/claude/emri-disk/ref-code"),
                               "TeukEvolution.jl"))
# ----------------------------------------------------------------------------------------------

# Bake the ORIGINAL Evo_lin_f coefficient arrays (A_pp, A_pq, B_pp, B_pq, B_pf, pre) at subsample
# (i,j) indices, for the Move-D fixture config (a=0.7, s=-2, m=2, nx=48, ny=24, cl=1). Also bake the
# exact (R, cy, sy) inputs so the port's {T}-generic closed forms can be unit-tested pointwise
# WITHOUT depending on node reproduction. Read-only src (../src). Output -> stdout as Julia literals.
cd(TEUK_ROOT)
using Printf
include(joinpath(TEUK_ROOT, "src") * "/Evolution.jl")
import .Evolution
include(joinpath(TEUK_ROOT, "src") * "/Radial.jl");  import .Radial
include(joinpath(TEUK_ROOT, "src") * "/Sphere.jl");  import .Sphere

nx=48; ny=24; cl=1.0; cfl=0.5; bhs=0.7; bhm=1.0; spin=-2; mv=2
minr = bhm*(1 + sqrt(1 - (bhs/bhm)^2)); maxR = 1/minr; dr = maxR/(nx-1)
Rv = Radial.R_vals(nx, dr)
Cv = Sphere.cos_vals(ny)
Sv = Sphere.sin_vals(ny)

evo = Evolution.Evo_lin_f(Rvals=Rv, Cvals=Cv, Svals=Sv, bhm=bhm, bhs=bhs, mv=mv, cl=cl, spin=spin)

isub = [1,12,24,36,48]
jsub = [6,12,18,24]

println("# AUTO-BAKED from ORIGINAL Evolution.Evo_lin_f — config a=0.7,s=-2,m=2,nx=48,ny=24,cl=1,bhm=1")
println("const HPHC_COEFF = (")
println("    bhm=$bhm, bhs=$bhs, cl=$cl, spin=$spin, mv=$mv,")
print("    rows=[")
first=true
for i in isub, j in jsub
    global first
    R=Rv[i]; cy=Cv[j]; sy=Sv[j]
    App=evo.A_pp[i,j]; Apq=evo.A_pq[i,j]
    Bpp=evo.B_pp[i,j]; Bpq=evo.B_pq[i,j]; Bpf=evo.B_pf[i,j]; pre=evo.pre[i,j]
    sep = first ? "" : ",\n         "
    first=false
    print(sep, "(i=$i, j=$j, R=", @sprintf("%.17e",R), ", cy=", @sprintf("% .17e",cy),
          ", sy=", @sprintf("%.17e",sy),
          ", App=", @sprintf("% .17e",App), ", Apq=", @sprintf("% .17e",Apq),
          ", Bpp=(", @sprintf("% .17e",real(Bpp)), ",", @sprintf("% .17e",imag(Bpp)), ")",
          ", Bpq=(", @sprintf("% .17e",real(Bpq)), ",", @sprintf("% .17e",imag(Bpq)), ")",
          ", Bpf=(", @sprintf("% .17e",real(Bpf)), ",", @sprintf("% .17e",imag(Bpf)), ")",
          ", pre=", @sprintf("% .17e",pre), ")")
end
println("],")
println(")")
