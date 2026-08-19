# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/capture_spinforce_s51.jl (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = d559b1c2f63c2f0e32cf8bad7aec4bb0a4c8f91617fed776e5936214a77f2102
# fixture(s)  : spacetime_spinforce_ref.jl
# oracle      : legacy packages/SpinningOrbit @ 394b7844 (its own Mathematica-verified green suite inputs)
# run         : julia --project=$GRTOOLS_LEGACY_PACKAGES_ROOT/SpinningOrbit test/fixtures/capture/capture_spinforce_s51.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

# ============================================================================================
# ONE-SHOT capture (S5-1): legacy oracle packages/SpinningOrbit spin-curvature force terms +
# numeric-correction engine + orbit-correction assembly → baked fixture literals for
# grtools/capabilities/Spacetime test/fixtures/spacetime_spinforce_ref.jl.
#
# Evaluation points are the EXACT input tuples of the oracle's own green test suite
# (packages/SpinningOrbit/test/ref_{spin_dev3,numcorr_force*,numcorr_chi,numcorr_dft,
# numcorr_delu,numcorr_solve,orbcorr_helpers,orbcorr}.jl), so the baked values inherit the
# oracle's Mathematica-verified grounding at the oracle's own tolerances.
#
# Run ONCE:  julia --project=packages/SpinningOrbit scratchpad/capture_spinforce_s51.jl
# The committed grtools test depends ONLY on the baked fixture file, never on packages/*.
# ============================================================================================
using SpinningOrbit, KerrGeodesics
const SO = SpinningOrbit

r6(x) = repr(x)   # full round-trip Float64 literal

io = IOBuffer()
println(io, "# CITATION: baked ONCE from the legacy oracle packages/SpinningOrbit (green suite) by")
println(io, "# scratchpad/capture_spinforce_s51.jl on $(read(`date -u +%Y-%m-%dT%H:%MZ`, String) |> strip).")
println(io, "# Evaluation points = the EXACT input tuples of the oracle's own Mathematica-verified test")
println(io, "# refs (packages/SpinningOrbit/test/ref_*.jl); functions: SpinningOrbit._δt3/_δr3/_δz3/_δφ3,")
println(io, "# the 47 numeric-correction force/coefficient terms, _chir_hat/_chith_hat, _numeric_force_dft,")
println(io, "# _delta_u_st, _numeric_correction, the 11 _oc_* orbit-correction helpers, and")
println(io, "# KerrSpinOrbitCorrection/KerrSpinOrbitNum. Tolerances mirror the oracle suite: 1e-12 on the")
println(io, "# closed-form force terms (pure algebra, machine-precision reproducible), 1e-8/1e-9 on the")
println(io, "# χ-maps (trajectory-mediated), 1e-7/1e-6 on DFT/δu/solve/assembly (DFT-truncation-limited,")
println(io, "# matches the oracle's own gates).")
println(io)

# ---------------- 1. spin-curvature worldline deviations (REF_SPIN_DEV3 inputs) ----------------
dev_in = [
    (r=10.2, z=0.3,   Ur=0.5,  Uz=-0.4, En=0.95, Jz=2.8, K=13.0, a=0.9),
    (r=12.5, z=-0.25, Ur=-0.3, Uz=0.6,  En=0.96, Jz=3.2, K=12.0, a=0.5),
    (r=8.1,  z=0.45,  Ur=0.2,  Uz=-0.5, En=0.94, Jz=2.0, K=11.5, a=0.9),
    (r=14.0, z=0.1,   Ur=-0.6, Uz=0.3,  En=0.97, Jz=3.5, K=15.0, a=0.7),
]
println(io, "const SPINFORCE_DEV3 = [")
for t in dev_in
    dt3 = SO._δt3(t.r, t.z, t.Ur, t.Uz, t.K, t.a)
    dr3 = SO._δr3(t.r, t.z, t.En, t.Jz, t.K, t.a)
    dz3 = SO._δz3(t.r, t.z, t.En, t.Jz, t.K, t.a)
    dphi3 = SO._δφ3(t.r, t.z, t.Ur, t.Uz, t.K, t.a)
    println(io, "    (r=$(r6(t.r)), z=$(r6(t.z)), Ur=$(r6(t.Ur)), Uz=$(r6(t.Uz)), En=$(r6(t.En)), Jz=$(r6(t.Jz)), K=$(r6(t.K)), a=$(r6(t.a)), dt3=$(r6(dt3)), dr3=$(r6(dr3)), dz3=$(r6(dz3)), dphi3=$(r6(dphi3))),")
end
println(io, "]")

# ------------- 2. main force layer (REF_NUMCORR_FORCE/2/3 inputs; + _MR_force capture) -------------
main_in = [
    (r=10.2, z=0.3,   dr=0.5,  dz=-0.4, En=0.95, Lz=2.8, K=13.0, a=0.9, dut=0.07,  dphi=-0.03, d2r=-0.012, d2z=-0.08, Q=9.216975000000001,  Ur=3.1, Uth=3.9),
    (r=12.5, z=-0.25, dr=-0.3, dz=0.6,  En=0.96, Lz=3.2, K=12.0, a=0.5, dut=-0.05, dphi=0.04,  d2r=-0.009, d2z=-0.06, Q=4.601599999999999,  Ur=2.7, Uth=3.4),
    (r=8.1,  z=0.45,  dr=0.2,  dz=-0.5, En=0.94, Lz=2.0, K=11.5, a=0.9, dut=0.02,  dphi=0.06,  d2r=-0.02,  d2z=-0.1,  Q=10.168284,          Ur=3.5, Uth=4.2),
]
println(io, "const SPINFORCE_MAIN = [")
for t in main_in
    vals = (
        Rt   = SO._Rt_force(t.r, t.z, t.dr, t.dz, t.K, t.a),
        Rphi = SO._Rphi_force(t.r, t.z, t.dr, t.dz, t.K, t.a),
        I2   = SO._I2_coef(t.r, t.z, t.En, t.Lz, t.K, t.a),
        I3   = SO._I3_coef(t.r, t.z, t.En, t.Lz, t.K, t.a),
        J    = SO._J_source(t.r, t.z, t.dut, t.dphi, t.En, t.Lz, t.K, t.a),
        d2rf = SO._d2r_dλ2(t.r, t.En, t.Lz, t.K, t.a),
        d2zf = SO._d2z_dλ2(t.z, t.En, t.Lz, t.Q, t.a),
        GR   = SO._GR_force(t.r, t.z, t.dr, t.a),
        GZ   = SO._GZ_force(t.r, t.z, t.dz, t.a),
        I1r  = SO._I1r_coef(t.r, t.z, t.dr, t.d2r, t.Ur, t.a),
        I1th = SO._I1θ_coef(t.r, t.z, t.dz, t.Uth, t.a),
        KR   = SO._KR_force(t.r, t.z, t.dr, t.a),
        KZ   = SO._KZ_force(t.r, t.z, t.dz, t.a),
        HR   = SO._HR_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a),
        HZ   = SO._HZ_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a),
        N1r  = SO._N1r_coef(t.r, t.z, t.dr, t.Ur, t.a),
        N1th = SO._N1th_coef(t.r, t.z, t.dz, t.Uth, t.a),
        N2   = SO._N2_coef(t.r, t.z, t.En, t.Lz, t.a),
        N3   = SO._N3_coef(t.r, t.z, t.En, t.Lz, t.a),
        U2   = SO._U2_coef(t.r, t.z, t.En, t.Lz, t.a),
        U3   = SO._U3_coef(t.r, t.z, t.En, t.Lz, t.a),
        QZ   = SO._QZ_force(t.z),
        MZ   = SO._MZ_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a),
        MR   = SO._MR_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a),
        Vgeo = SO._Vgeo(t.r, t.z, t.En, t.Lz, t.K, t.a),
        V    = SO._V_source(t.r, t.z, t.dut, t.dphi, t.En, t.Lz, t.K, t.a),
        P    = SO._P_source(t.r, t.z, t.dut, t.dphi, t.En, t.Lz, t.a),
    )
    ins = "r=$(r6(t.r)), z=$(r6(t.z)), dr=$(r6(t.dr)), dz=$(r6(t.dz)), En=$(r6(t.En)), Lz=$(r6(t.Lz)), K=$(r6(t.K)), a=$(r6(t.a)), dut=$(r6(t.dut)), dphi=$(r6(t.dphi)), d2r=$(r6(t.d2r)), d2z=$(r6(t.d2z)), Q=$(r6(t.Q)), Ur=$(r6(t.Ur)), Uth=$(r6(t.Uth))"
    outs = join(["$k=$(r6(v))" for (k, v) in pairs(vals)], ", ")
    println(io, "    ($ins, $outs),")
end
println(io, "]")

# ---------- 3. χ/sin(I)-dependent forces, first block (REF_NUMCORR_FORCE4 inputs) ----------
chi4_in = [
    (r=10.2, z=0.3,   dr=0.5,  dz=-0.4, a=0.9, e=0.2,  p=10.0, chir=0.8, chith=2.1, incl=1.0472, Ur=3.1),
    (r=12.5, z=-0.25, dr=-0.3, dz=0.6,  a=0.5, e=0.3,  p=12.0, chir=1.7, chith=0.9, incl=0.7754, Ur=2.7),
    (r=8.1,  z=0.45,  dr=0.2,  dz=-0.5, a=0.9, e=0.25, p=8.5,  chir=2.5, chith=3.4, incl=1.231,  Ur=3.5),
]
println(io, "const SPINFORCE_CHI4 = [")
for t in chi4_in
    vals = (
        Fr  = SO._Fr_force(t.chir, t.e, t.p),
        Gr  = SO._Gr_force(t.r, t.z, t.dr, t.chir, t.a, t.e, t.p),
        Gth = SO._Gth_force(t.r, t.z, t.dz, t.chith, t.incl, t.a),
        Kr  = SO._Kr_force(t.r, t.z, t.dr, t.chir, t.a, t.e, t.p),
        Kth = SO._Kth_force(t.r, t.z, t.dz, t.chith, t.incl, t.a),
        Qth = SO._Qth_force(t.z, t.chith, t.incl),
        Sr  = SO._Sr_force(t.r, t.z, t.dr, t.chir, t.a, t.e, t.p),
        SR  = SO._SR_force(t.r, t.z, t.dr, t.a),
        Sth = SO._Sth_force(t.r, t.z, t.dz, t.chith, t.incl, t.a),
        SZ  = SO._SZ_force(t.r, t.z, t.dz, t.a),
        U1r = SO._U1r_coef(t.r, t.z, t.dr, t.Ur, t.a),
    )
    ins = "r=$(r6(t.r)), z=$(r6(t.z)), dr=$(r6(t.dr)), dz=$(r6(t.dz)), a=$(r6(t.a)), e=$(r6(t.e)), p=$(r6(t.p)), chir=$(r6(t.chir)), chith=$(r6(t.chith)), incl=$(r6(t.incl)), Ur=$(r6(t.Ur))"
    outs = join(["$k=$(r6(v))" for (k, v) in pairs(vals)], ", ")
    println(io, "    ($ins, $outs),")
end
println(io, "]")

# ---------- 4. heaviest 2nd-order χ forces (REF_NUMCORR_FORCE5 inputs) ----------
chi5_in = [
    (r=10.2, z=0.3,  dr=0.5, dz=-0.4, d2r=-0.012, d2z=-0.08, En=0.95, Lz=2.8, a=0.9, e=0.2,  p=10.0, chir=0.8, chith=2.1, incl=1.0472, Uth=3.9),
    (r=8.1,  z=0.45, dr=0.2, dz=-0.5, d2r=-0.02,  d2z=-0.1,  En=0.94, Lz=2.0, a=0.9, e=0.25, p=8.5,  chir=2.5, chith=3.4, incl=1.231,  Uth=4.2),
]
println(io, "const SPINFORCE_CHI5 = [")
for t in chi5_in
    vals = (
        Hr   = SO._Hr_force(t.r, t.z, t.dr, t.dz, t.d2r, t.chir, t.En, t.Lz, t.a, t.e, t.p),
        Hth  = SO._Hth_force(t.r, t.z, t.dr, t.dz, t.chith, t.incl, t.En, t.Lz, t.a),
        Mr   = SO._Mr_force(t.r, t.z, t.dr, t.dz, t.chir, t.En, t.Lz, t.a, t.e, t.p),
        Mth  = SO._Mth_force(t.r, t.z, t.dr, t.dz, t.chith, t.incl, t.En, t.Lz, t.a),
        Tr   = SO._Tr_force(t.r, t.z, t.dr, t.dz, t.chir, t.En, t.Lz, t.a, t.e, t.p),
        TR   = SO._TR_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a),
        Tth  = SO._Tth_force(t.r, t.z, t.dr, t.dz, t.d2z, t.chith, t.incl, t.En, t.Lz, t.a),
        TZ   = SO._TZ_force(t.r, t.z, t.dr, t.dz, t.d2z, t.incl, t.En, t.Lz, t.a),
        U1th = SO._U1th_coef(t.r, t.z, t.dz, t.d2z, t.chith, t.En, t.Lz, t.a, t.Uth),
    )
    ins = "r=$(r6(t.r)), z=$(r6(t.z)), dr=$(r6(t.dr)), dz=$(r6(t.dz)), d2r=$(r6(t.d2r)), d2z=$(r6(t.d2z)), En=$(r6(t.En)), Lz=$(r6(t.Lz)), a=$(r6(t.a)), e=$(r6(t.e)), p=$(r6(t.p)), chir=$(r6(t.chir)), chith=$(r6(t.chith)), incl=$(r6(t.incl)), Uth=$(r6(t.Uth))"
    outs = join(["$k=$(r6(v))" for (k, v) in pairs(vals)], ", ")
    println(io, "    ($ins, $outs),")
end
println(io, "]")

# ---------- 5. Darwin-angle maps χrhat/χθhat (REF_NUMCORR_CHI inputs) ----------
chih_in = [(a=0.9, p=12.0, e=0.2, x=0.5, lam=l) for l in (0.05, 0.18, 0.4, 0.65)]
append!(chih_in, [(a=0.5, p=10.0, e=0.2, x=0.7, lam=l) for l in (0.05, 0.18, 0.4, 0.65)])
println(io, "const SPINFORCE_CHIHAT = [")
for t in chih_in
    chir = SO._chir_hat(t.a, t.p, t.e, t.x, t.lam)
    chith = SO._chith_hat(t.a, t.p, t.e, t.x, t.lam)
    println(io, "    (a=$(r6(t.a)), p=$(r6(t.p)), e=$(r6(t.e)), x=$(r6(t.x)), lam=$(r6(t.lam)), chir=$(r6(chir)), chith=$(r6(chith))),")
end
println(io, "]")

# ---------- 6. 2-torus force DFT (Im parts, nmax=kmax=2, REF_NUMCORR_DFT orbit) ----------
let a=0.9, p=12.0, e=0.2, x=0.5, nmax=2, kmax=2
    Rtf, Rpf = SO._numeric_force_dft(a, p, e, x, nmax, kmax)
    println(io, "const SPINFORCE_DFT = (a=$(r6(a)), p=$(r6(p)), e=$(r6(e)), x=$(r6(x)), nmax=$nmax, kmax=$kmax,")
    println(io, "    RtIm = $(repr(imag.(Rtf))),")
    println(io, "    RphiIm = $(repr(imag.(Rpf))))")
end

# ---------- 7. spin-shift potentials δu_{ds,t/φ} (REF_NUMCORR_DELU phase points) ----------
let a=0.9, p=12.0, e=0.2, x=0.5, nmax=2, kmax=2
    Rtf, Rpf = SO._numeric_force_dft(a, p, e, x, nmax, kmax)
    fr = KerrGeoFrequencies(a, p, e, x; Time="Mino")
    println(io, "const SPINFORCE_DELU = [")
    for (lr, lt) in ((0.1, 0.05), (0.3, 0.18), (0.55, 0.4), (0.2, 0.6))
        dut = SO._delta_u_st(Rtf, fr.Ur, fr.Uth, nmax, kmax, lr, lt)
        dup = SO._delta_u_st(Rpf, fr.Ur, fr.Uth, nmax, kmax, lr, lt)
        println(io, "    (lr=$(r6(lr)), lt=$(r6(lt)), delut=$(r6(dut)), deluph=$(r6(dup))),")
    end
    println(io, "]")
end

# ---------- 8. FINAL numeric-correction solve (REF_NUMCORR_SOLVE orbit) ----------
let a=0.9, p=12.0, e=0.2, x=0.5, nmax=2, kmax=2
    c = SO._numeric_correction(a, p, e, x, nmax, kmax)
    println(io, "const SPINFORCE_SOLVE = (a=$(r6(a)), p=$(r6(p)), e=$(r6(e)), x=$(r6(x)), nmax=$nmax, kmax=$kmax, YrS=$(r6(c.ΥrS)), YthS=$(r6(c.ΥθS)), udt0S=$(r6(c.udt0S)), udphi0S=$(r6(c.udφ0S)), err=$(r6(c.error)))")
end

# ---------- 9. orbit-correction four-velocity/metric helpers (REF_ORBCORR_HELPERS inputs) ----------
oc_in = [
    (r=10.2, z=0.3,   K=12.3, En=0.96, Lz=1.9, a=0.9),
    (r=12.5, z=-0.25, K=11.8, En=0.95, Lz=2.1, a=0.5),
    (r=8.1,  z=0.45,  K=13.0, En=0.94, Lz=1.7, a=0.9),
]
println(io, "const SPINFORCE_OC_HELPERS = [")
for t in oc_in
    vals = (
        h1 = SO._oc_h1(t.r, t.z, t.K, t.a),
        h2 = SO._oc_h2(t.r, t.z, t.K, t.En, t.Lz, t.a),
        h3 = SO._oc_h3(t.r, t.z, t.K, t.a),
        Sgtt = SO._oc_Sgtt(t.r, t.z, t.a),
        Sgphph = SO._oc_Sgphph(t.r, t.z, t.a),
        Sgtph = SO._oc_Sgtph(t.r, t.a),
        dSgttdr = SO._oc_dSgttdr(t.r, t.a),
        dSgttdz = SO._oc_dSgttdz(t.z, t.a),
        dSgtphdr = SO._oc_dSgtphdr(t.r, t.a),
        dSgphphdr = SO._oc_dSgphphdr(t.r, t.a),
        dSgphphdz = SO._oc_dSgphphdz(t.z, t.a),
    )
    ins = "r=$(r6(t.r)), z=$(r6(t.z)), K=$(r6(t.K)), En=$(r6(t.En)), Lz=$(r6(t.Lz)), a=$(r6(t.a))"
    outs = join(["$k=$(r6(v))" for (k, v) in pairs(vals)], ", ")
    println(io, "    ($ins, $outs),")
end
println(io, "]")

# ---------- 10. END-TO-END KerrSpinOrbitCorrection + KerrSpinOrbitNum (REF_ORBCORR orbit) ----------
let a=0.9, p=12.0, e=0.2, x=0.5, nmax=2, kmax=2, spar=0.1
    corr = KerrSpinOrbitCorrection(a, p, e, x, nmax, kmax)
    num = KerrSpinOrbitNum(corr, spar)
    sampWr = [0.7, 2.1, 3.5, 1.0]; sampWz = [1.3, 0.4, 5.0, 2.5]
    DtS   = [corr.TrajectoryCorrection[1](wr, wz) for (wr, wz) in zip(sampWr, sampWz)]
    rS    = [corr.TrajectoryCorrection[2](wr, wz) for (wr, wz) in zip(sampWr, sampWz)]
    zS    = [corr.TrajectoryCorrection[3](wr, wz) for (wr, wz) in zip(sampWr, sampWz)]
    DphiS = [corr.TrajectoryCorrection[4](wr, wz) for (wr, wz) in zip(sampWr, sampWz)]
    oc    = [corr.OrbitCorrection(wr, wz) for (wr, wz) in zip(sampWr, sampWz)]
    println(io, "const SPINFORCE_ORBCORR = (a=$(r6(a)), p=$(r6(p)), e=$(r6(e)), x=$(r6(x)), nmax=$nmax, kmax=$kmax, spar=$(r6(spar)),")
    println(io, "    ScriptI=$(r6(corr.ScriptI)), Ehat=$(r6(corr.Ehat)), Lzhat=$(r6(corr.Lzhat)), Khat=$(r6(corr.Khat)),")
    println(io, "    GammaHat=$(r6(corr.MinoFrequenciesGeo[1])), UrHat=$(r6(corr.MinoFrequenciesGeo[2])), UthHat=$(r6(corr.MinoFrequenciesGeo[3])), UphHat=$(r6(corr.MinoFrequenciesGeo[4])),")
    println(io, "    GammaS=$(r6(corr.MinoFrequenciesCorrection[1])), UrS=$(r6(corr.MinoFrequenciesCorrection[2])), UthS=$(r6(corr.MinoFrequenciesCorrection[3])), UphS=$(r6(corr.MinoFrequenciesCorrection[4])),")
    println(io, "    ES=$(r6(corr.ES)), LS=$(r6(corr.LS)), err=$(r6(corr.error)),")
    println(io, "    sampWr=$(repr(sampWr)), sampWz=$(repr(sampWz)),")
    println(io, "    DtS=$(repr(DtS)), rS=$(repr(rS)), zS=$(repr(zS)), DphiS=$(repr(DphiS)),")
    println(io, "    ocDtS=$(repr([o.ΔtS for o in oc])), ocrS=$(repr([o.rS for o in oc])), oczS=$(repr([o.zS for o in oc])), ocDphiS=$(repr([o.ΔφS for o in oc])),")
    println(io, "    ocUrS=$(repr([o.UrS for o in oc])), ocUzS=$(repr([o.UzS for o in oc])),")
    println(io, "    numGamma=$(r6(num.MinoFrequencies[1])), numUr=$(r6(num.MinoFrequencies[2])), numUth=$(r6(num.MinoFrequencies[3])), numUph=$(r6(num.MinoFrequencies[4])),")
    println(io, "    numBLr=$(r6(num.BLFrequencies[1])), numBLth=$(r6(num.BLFrequencies[2])), numBLph=$(r6(num.BLFrequencies[3])),")
    println(io, "    numES=$(r6(num.ES)), numLS=$(r6(num.LS)), numEhat=$(r6(num.Ehat)), numLzhat=$(r6(num.Lzhat)), numKhat=$(r6(num.Khat)))")
end

# ---------- 11. tests.nb PUBLISHED golden (KerrSpinOrbitCorrection[0.95,10,0.2,0.8,8,8]) ----------
# The oracle suite reproduces the AUTHORS' own tests.nb Out[] literals; re-capture + time it here
# to decide whether the 8,8 case is affordable in the grtools suite.
let t0 = time()
    ocg = KerrSpinOrbitCorrection(0.95, 10.0, 0.2, 0.8, 8, 8)
    dt = time() - t0
    println(io, "# 8,8 golden captured in $(round(dt; digits=1)) s (oracle-side wall time)")
    println(io, "const SPINFORCE_GOLDEN88 = (a=0.95, p=10.0, e=0.2, x=0.8, nmax=8, kmax=8,")
    println(io, "    Ehat=$(r6(ocg.Ehat)), GammaHat=$(r6(ocg.MinoFrequenciesGeo[1])),")
    println(io, "    GammaS=$(r6(ocg.MinoFrequenciesCorrection[1])), BLc1=$(r6(ocg.BLFrequenciesCorrection[1])))")
    @assert isapprox(ocg.Ehat, 0.9544746770348815028245737919294451657838972804860443023072; rtol=1e-12)
    @assert isapprox(ocg.MinoFrequenciesGeo[1], 124.7629736005981287023759422670878097995901511945922024627606; rtol=1e-12)
    @assert isapprox(ocg.MinoFrequenciesCorrection[1], -0.30438197634945600312974191602274806466; rtol=1e-10)
    @assert isapprox(ocg.BLFrequenciesCorrection[1], 0.00128402591293930362965445382249630736; rtol=1e-10)
    println(stderr, "GOLDEN88 reproduced the published tests.nb literals; wall $(round(dt; digits=1)) s")
end

out = "$(FIXTURE_DIR)/spacetime_spinforce_ref.jl"
write(out, String(take!(io)))
println("WROTE $out")
