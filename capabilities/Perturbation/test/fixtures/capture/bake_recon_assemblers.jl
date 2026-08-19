# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/npghp_moveD/TeukEvolution_run/minrun/bake_recon_assemblers.jl (untracked one-shot, mtime 2026-07-04)
#               sha256(original) = 672a935b16665c210347c62895d905a1c6aa8825a6d952d9e004270683001132
# fixture(s)  : td_hphc_recon_upstream.jl
# oracle      : ref-code/TeukEvolution.jl/src/LinearEvolution.jl @ 3c8e47f (set_psi3_k!/set_lam_k!/set_psi2_k!/set_hmbmb_k!/set_pi_k!/set_hlmb_k!/set_muhll_k!/set_psi4_p_k!)
# run         : julia --project=<TeukEvolution env> test/fixtures/capture/bake_recon_assemblers.jl > <fixture>
#
# Paths rewritten to the ENV-overridable GRTOOLS_TEUKEVOLUTION_ROOT / GRTOOLS_REFCODE_ROOT
# (defaults reproduce the capture host).  Body otherwise VERBATIM from the original one-shot.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const TEUK_ROOT = get(ENV, "GRTOOLS_TEUKEVOLUTION_ROOT",
                      joinpath(get(ENV, "GRTOOLS_REFCODE_ROOT", "/data/haiyangw/claude/emri-disk/ref-code"),
                               "TeukEvolution.jl"))
# ----------------------------------------------------------------------------------------------

# Bake per-stage k outputs of the ORIGINAL TeukEvolution reconstruction RHS assemblers
# (LinearEvolution.jl set_psi3_k! / set_lam_k! / set_psi2_k! / set_hmbmb_k! / set_pi_k! /
#  set_hlmb_k! / set_muhll_k! + the P1-relevant raw set_psi4_p_k!) on SEEDED pseudo-random field
# data, for the Move-D (D4) reconstruction port fixture. Read-only src (../src). Output -> stdout as
# Julia literals (redirect into grtools .../test/fixtures/td_hphc_recon_upstream.jl).
#
# Design: every assembler is POINTWISE in (i,j) (verified by reading LinearEvolution.jl:32-238), so we
# bake, per subsampled (i,j): the grid (R,cy,sy), the (upstream) background NP scalars mu_0/tau_0/pi_0,
# ALL seeded complex input values, and each assembler's output k. The D4 test rebuilds column vectors
# from these rows, computes NP via Spacetime.NPGHP (bit-identical to upstream BackgroundNP at Float64 —
# the √2 quirk-4 fix only diverges at BigFloat), and calls the PORTED assemblers, comparing k pointwise.
# For psi4 we bake the RAW upstream set_psi4_p_k! (the P1-buggy no-`pre` bracket) at m=+2 and m=-2 so
# the port can gate `kp_port == pre .* raw` (the P1 fix — pre~0.03-0.06, so raw is ~16-30x too large).
cd(TEUK_ROOT)
using Printf
include(joinpath(TEUK_ROOT, "src") * "/LinearEvolution.jl"); import .LinearEvolution
include(joinpath(TEUK_ROOT, "src") * "/Radial.jl");          import .Radial
include(joinpath(TEUK_ROOT, "src") * "/Sphere.jl");          import .Sphere
include(joinpath(TEUK_ROOT, "src") * "/BackgroundNP.jl");    import .BackgroundNP
include(joinpath(TEUK_ROOT, "src") * "/Evolution.jl");       import .Evolution

const LE = LinearEvolution

# ---- config + grid (a=0.7, nx=32, ny=24, cl=1) ---------------------------------------------------
nx = 32; ny = 24; cl = 1.0; bhs = 0.7; bhm = 1.0
minr = bhm * (1 + sqrt(1 - (bhs / bhm)^2)); maxR = 1 / minr; dr = maxR / (nx - 1)
Rv = Radial.R_vals(nx, dr)
Yv = Sphere.Y_vals(ny)
Cv = Sphere.cos_vals(ny)
Sv = Sphere.sin_vals(ny)
NP = BackgroundNP.NP_0{Float64}(Rvals = Rv, Yvals = Yv, Cvals = Cv, Svals = Sv,
                                bhm = bhm, bhs = bhs, cl = cl)

# ---- deterministic LCG fill (values baked verbatim; test never reproduces the LCG) ---------------
mutable struct LCG; s::UInt64; end
function nf!(g::LCG)
    g.s = g.s * 0x5851f42d4c957f2d + 0x14057b7ef767814f
    return Float64(g.s >> 11) / Float64(UInt64(1) << 53) * 2 - 1   # ∈ (-1,1)
end
fc!(A, g) = (for idx in eachindex(A); A[idx] = complex(nf!(g), nf!(g)); end; A)
fr!(A, g) = (for idx in eachindex(A); A[idx] = nf!(g); end; A)
newc() = Array{ComplexF64}(undef, nx, ny)
newr() = Array{Float64}(undef, nx, ny)

g = LCG(0x00000000d4c95709)

# distinct seeded complex input arrays (recon assemblers)
psi4  = fc!(newc(), g); psi3  = fc!(newc(), g); psi2  = fc!(newc(), g); lam = fc!(newc(), g)
pi_l  = fc!(newc(), g); hmbmb = fc!(newc(), g); hlmb  = fc!(newc(), g); muhll = fc!(newc(), g)
psi4_edth = fc!(newc(), g); psi3_edth = fc!(newc(), g); hlmb_edth = fc!(newc(), g); pi_edth = fc!(newc(), g)
pi_nm = fc!(newc(), g); hmbmb_nm = fc!(newc(), g); hlmb_nm = fc!(newc(), g)
hmbmb_edth_nm = fc!(newc(), g); hlmb_edth_nm = fc!(newc(), g)

# ---- run the ORIGINAL recon assemblers -----------------------------------------------------------
k_psi3 = newc(); k_lam = newc(); k_psi2 = newc(); k_hmbmb = newc(); k_pi = newc(); k_hlmb = newc(); k_muhll = newc()

LE.set_psi3_k!(k = k_psi3, psi3 = psi3, psi4 = psi4, psi4_edth = psi4_edth,
               mu_0 = NP.mu_0, tau_0 = NP.tau_0, R = Rv)
LE.set_lam_k!(k = k_lam, lam = lam, psi4 = psi4, mu_0 = NP.mu_0, R = Rv)
LE.set_psi2_k!(k = k_psi2, psi2 = psi2, psi3 = psi3, psi3_edth = psi3_edth,
               mu_0 = NP.mu_0, tau_0 = NP.tau_0, R = Rv)
LE.set_hmbmb_k!(k = k_hmbmb, hmbmb = hmbmb, lam = lam, mu_0 = NP.mu_0, R = Rv)
LE.set_pi_k!(k = k_pi, lam = lam, hmbmb = hmbmb, psi3 = psi3,
             mu_0 = NP.mu_0, tau_0 = NP.tau_0, pi_0 = NP.pi_0, R = Rv)
LE.set_hlmb_k!(k = k_hlmb, hlmb = hlmb, hmbmb = hmbmb, pi_l = pi_l,
               mu_0 = NP.mu_0, tau_0 = NP.tau_0, R = Rv)
LE.set_muhll_k!(k = k_muhll, muhll = muhll, hlmb_edth = hlmb_edth, hlmb = hlmb,
                pi_edth = pi_edth, pi_l = pi_l, psi2 = psi2, pi_nm = pi_nm,
                hmbmb_edth_nm = hmbmb_edth_nm, hmbmb_nm = hmbmb_nm,
                hlmb_edth_nm = hlmb_edth_nm, hlmb_nm = hlmb_nm,
                mu_0 = NP.mu_0, tau_0 = NP.tau_0, pi_0 = NP.pi_0, R = Rv)

# ---- psi4 P1: seeded bracket + coeffs, raw upstream set_psi4_p_k! (no `pre`) at m=+2 and m=-2 ------
f4  = fc!(newc(), g); p4  = fc!(newc(), g)
f_rd1 = fc!(newc(), g); f_rd2 = fc!(newc(), g); p_rd1 = fc!(newc(), g); sph_lap = fc!(newc(), g)
A_pp = fr!(newr(), g); A_pq = fr!(newr(), g); pre = fr!(newr(), g)
B_pp = fc!(newc(), g); B_pq = fc!(newc(), g); B_pf = fc!(newc(), g)
raw_kp = newc()
LE.set_psi4_p_k!(kp = raw_kp, f_rd1 = f_rd1, f_rd2 = f_rd2, sph_lap = sph_lap, p_rd1 = p_rd1,
                 f = f4, p = p4, A_pp = A_pp, A_pq = A_pq, B_pp = B_pp, B_pq = B_pq, B_pf = B_pf)

# ---- subsample + emit ----------------------------------------------------------------------------
isub = [1, 8, 16, 24, 32]; jsub = [6, 12, 18, 24]
c(z) = @sprintf("(% .17e,% .17e)", real(z), imag(z))
r(x) = @sprintf("% .17e", x)

println("# AUTO-BAKED — DO NOT EDIT. ORIGINAL TeukEvolution reconstruction RHS assemblers (Move D / D4).")
println("#   ref-code : ref-code/TeukEvolution.jl  git 3c8e47f4  (LinearEvolution.jl:32-238)")
println("#   config   : a(bhs)=$bhs bhm=$bhm cl=$cl nx=$nx ny=$ny ; seeded pseudo-random field data (LCG).")
println("#   NP (mu_0,tau_0,pi_0) = upstream BackgroundNP; bit-identical to Spacetime.NPGHP at Float64.")
println("#   Each recon assembler is pointwise in (i,j); psi4 raw = set_psi4_p_k! WITHOUT `pre` (P1 gate).")
println("const RECON_ASM = (")
println("  bhs=$bhs, bhm=$bhm, cl=$cl, nx=$nx, ny=$ny,")
print("  pts=[")
first = true
for i in isub, j in jsub
    global first
    sep = first ? "" : ",\n       "; first = false
    print(sep, "(i=$i, j=$j, R=", r(Rv[i]), ", cy=", r(Cv[j]), ", sy=", r(Sv[j]),
          ", mu0=", c(NP.mu_0[i,j]), ", tau0=", c(NP.tau_0[i,j]), ", pi0=", c(NP.pi_0[i,j]),
          # recon inputs
          ", psi4=", c(psi4[i,j]), ", psi3=", c(psi3[i,j]), ", psi2=", c(psi2[i,j]),
          ", lam=", c(lam[i,j]), ", pi_l=", c(pi_l[i,j]), ", hmbmb=", c(hmbmb[i,j]),
          ", hlmb=", c(hlmb[i,j]), ", muhll=", c(muhll[i,j]),
          ", psi4_edth=", c(psi4_edth[i,j]), ", psi3_edth=", c(psi3_edth[i,j]),
          ", hlmb_edth=", c(hlmb_edth[i,j]), ", pi_edth=", c(pi_edth[i,j]),
          ", pi_nm=", c(pi_nm[i,j]), ", hmbmb_nm=", c(hmbmb_nm[i,j]), ", hlmb_nm=", c(hlmb_nm[i,j]),
          ", hmbmb_edth_nm=", c(hmbmb_edth_nm[i,j]), ", hlmb_edth_nm=", c(hlmb_edth_nm[i,j]),
          # recon outputs
          ", k_psi3=", c(k_psi3[i,j]), ", k_lam=", c(k_lam[i,j]), ", k_psi2=", c(k_psi2[i,j]),
          ", k_hmbmb=", c(k_hmbmb[i,j]), ", k_pi=", c(k_pi[i,j]), ", k_hlmb=", c(k_hlmb[i,j]),
          ", k_muhll=", c(k_muhll[i,j]),
          # psi4 P1 inputs + raw output
          ", f4=", c(f4[i,j]), ", p4=", c(p4[i,j]), ", f_rd1=", c(f_rd1[i,j]), ", f_rd2=", c(f_rd2[i,j]),
          ", p_rd1=", c(p_rd1[i,j]), ", sph_lap=", c(sph_lap[i,j]),
          ", A_pp=", r(A_pp[i,j]), ", A_pq=", r(A_pq[i,j]), ", pre=", r(pre[i,j]),
          ", B_pp=", c(B_pp[i,j]), ", B_pq=", c(B_pq[i,j]), ", B_pf=", c(B_pf[i,j]),
          ", raw_kp=", c(raw_kp[i,j]), ")")
end
println("],")
println(")")
