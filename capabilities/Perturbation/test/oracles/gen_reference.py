# TRACKED ORACLE SCRIPT  (GRR::g14-oracle-scripts, repro-audit gap G14)
# 
# origin : emri-disk/scratchpad/qqnm_phase3/gen_reference.py (untracked one-shot)
# sha256 : c369ac7828f097aacbd07d93ccbd07ad95d785f482f47aea2712e188bab97bfb
# backs  : test/fixtures/qnm_fast_reference.jl (via scratchpad/qqnm_phase3/qnm_package_reference.txt, sha256 157b38c6... recorded in the fixture header)
# engine : python3 + the `qnm` package (ref-code/qnm @ f3abd18e)
# 
# Tracked for provenance/regeneration only: nothing in the Julia test suite executes it
# (no Wolfram/py dependency is added).  Re-running it requires the engine above.

import os, sys, time, datetime
WORK="/data/haiyangw/claude/emri-disk/scratchpad/qqnm_phase3/qnm_run"
sys.path.insert(0, WORK)
os.environ["QNMCACHEDIR"]="/data/haiyangw/claude/emri-disk/scratchpad/qqnm_phase3/qnm_cache"
os.environ.setdefault("NUMBA_CACHE_DIR", "/data/haiyangw/claude/emri-disk/scratchpad/qqnm_phase3/numba_cache")
import numpy as np
import qnm
from qnm.spinsequence import KerrSpinSeq
from qnm.schwarzschild.overtonesequence import SchwOvertoneSeq
from qnm.schwarzschild.tabulated import QNMDict

def cfmt(z):
    z=complex(z)
    return "%+.17e %+.17ej" % (z.real, z.imag)

OUT=[]
def emit(line=""):
    OUT.append(line)

emit("="*100)
emit("QNM PACKAGE REFERENCE TABLES  (Q1c-PREP-A cross-fixtures)")
emit("="*100)
emit("Provenance:")
emit("  package : qnm  (Leo C. Stein), method: Cook-Zalutskiy 1410.7698 spectral-angular + Leaver-CF radial")
emit("  source  : /data/haiyangw/claude/emri-disk/ref-code/qnm  @ git f3abd18e (2026-06-22)")
emit("            run from READ-ONLY copy at %s" % WORK)
emit("  python  : %s" % sys.version.split()[0])
emit("  numpy %s | scipy %s | numba %s" % (np.__version__, __import__('scipy').__version__, __import__('numba').__version__))
emit("  generated: %s UTC" % datetime.datetime.utcnow().isoformat())
emit("  convention: s=-2 gravitational; omega,A in units M=1; a = chi (dimensionless spin), prograde m>0.")
emit("  A = qnm angular separation constant = eigenvalue of Cook-Zalutskiy spherical-spheroidal M-matrix")
emit("      (at a=0: A = l(l+1)-s(s+1); l_max=20 truncation, cf_tol=1e-10, tol=sqrt(eps)).")
emit("  Each omega(a) is a GENUINE Leaver root-polish at that exact a (not bare interpolation).")
emit("="*100)
emit()

# ---------- ANCHOR verification ----------
seq220 = qnm.modes_cache(s=-2,l=2,m=2,n=0)
om0,A0,_ = seq220(a=0.0)
ref = 0.37367168441804177-0.08896231568893546j
emit("ANCHOR CHECK  (2,2,0) a=0  vs in-repo TeukEvolution.jl/qnm/a0.0_l2_m2.h5")
emit("  qnm omega = %s" % cfmt(om0))
emit("  h5  omega = %s" % cfmt(ref))
emit("  diff      = %s   |diff| = %.3e" % (cfmt(om0-ref), abs(om0-ref)))
emit()

# ---------- 5 campaign multiplets over chi grid ----------
chis=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]
multiplets=[("(2,2,0)",2,2,0),("(3,3,0)",3,3,0),("(4,4,0)",4,4,0),("(2,0,0)",2,0,0),("(2,-2,0)",2,-2,0)]
emit("#"*100)
emit("# TABLE 1: campaign multiplets, s=-2, over chi grid {0,0.1,...,0.9,0.95}")
emit("#"*100)
seqs={}
for lbl,l,m,n in multiplets:
    t0=time.time()
    seq = qnm.modes_cache(s=-2,l=l,m=m,n=n)  # builds full a=0..0.99 sequence from scratch (offline)
    dt=time.time()-t0
    seqs[lbl]=seq
    emit()
    emit("--- multiplet %s  (l=%d, m=%d, n=%d)   [sequence build: %.2f s, %d nodes] ---"%(lbl,l,m,n,dt,len(seq.a)))
    emit("  %-5s  %-42s  %-42s  %-10s"%("chi","omega (Re +Im j)","A_sep (Re +Im j)","cf_err"))
    for chi in chis:
        om,A,C = seq(a=chi)
        cf_err,n_frac = seq.solver.get_cf_err()
        emit("  %-5.2f  %s   %s   %.2e"%(chi, cfmt(om), cfmt(A), cf_err))
emit()

# ---------- a=0 (2,2,n) n=0..10 overtone ladder ----------
emit("#"*100)
emit("# TABLE 2: a=0 Schwarzschild (2,2,n) overtone ladder, n=0..10, s=-2")
emit("#"*100)
t0=time.time()
so = SchwOvertoneSeq(s=-2,l=2,n_max=10)
so.find_sequence()
dt=time.time()-t0
emit("  [SchwOvertoneSeq build n=0..10: %.2f s]   A_sep = %s (real, = l(l+1)-s(s+1)=4)"%(dt, cfmt(so.A)))
qd=QNMDict(init=True)
emit("  %-3s  %-42s  %-12s  %-6s  %-42s"%("n","omega (Re +Im j) [SchwOvertoneSeq]","cf_err","n_frac","omega [pickle QNMDict, cross-check]"))
for n in range(11):
    om=so.omega[n]; cf=so.cf_err[n]; nf=int(so.n_frac[n])
    omp,cfp,nfp = qd(-2,2,n)
    emit("  %-3d  %s   %.2e   %-6d  %s"%(n, cfmt(om), cf, nf, cfmt(omp)))
emit()

txt="\n".join(OUT)
with open("/data/haiyangw/claude/emri-disk/scratchpad/qqnm_phase3/qnm_package_reference.txt","w") as f:
    f.write(txt+"\n")
print(txt)
