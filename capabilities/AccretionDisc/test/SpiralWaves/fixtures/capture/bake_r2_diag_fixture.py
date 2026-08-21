#!/usr/bin/env python3
"""bake_r2_diag_fixture.py -- bake the route-R2 DIAGNOSTICS fixture for test_r2_generated.jl.

Source: the R1-vs-R2 gate fixture of the emri-disk mission repository (EGD::R1-vs-R2-gate),
    results/emri-gasdisk/paper_EMRI-GASDISK/codes/cas/deriv-master-system-R2/master_fixture.json
    (generator master_fixture.wls, seed 20260819, 30 significant digits from 40-digit arithmetic).
Where bake_master_r2_fixture.py keeps the master ENTRIES (C, S) of eq (47), this script keeps the
columns the generated route-R2 diagnostics need (test/SpiralWaves/fixtures/oracle/r2_generated.jl,
emitted by emri-disk codes/cas/R2_to_julia.wls from R2_master_system.m):

  * the full R2 input set at the field point -- nu, psi, mu, psip, v, vp, omt, eps, epsp, p, pp, Gam1,
    omega_d and the COVARIANT test-profile modes htt, htp, hpp, hrr, hzz (+ the derivatives http, htpp,
    hppp, hrrp, hzzp and nup, mup, omtp, Omega, Omp for the end-to-end background gate);
  * sqrt(-g) = e^{nu+psi+mu} and the background u_phi (72a normalisation);
  * the sample master vector (sample_y0, sample_y1) and Mathematica's values of
      deltaUphi  (as the three coefficients [c_y0, c_y1, c_h]: deltaUphi = c_y0 y0 + c_y1 y1 + c_h),
      AMF               (per-mode AMF integrand of (72a) at the sample y),
      TorqueDensity     (per-mode torque density of (72b) at the sample y).

All numbers stay decimal STRINGS (parsed at Float64 or BigFloat precision in the test); complex values
are (re, im) string pairs.  Usage: python3 bake_r2_diag_fixture.py [master_fixture.json] [out .jl]
"""
import hashlib, json, os, sys

here = os.path.dirname(os.path.abspath(__file__))
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, *([".."] * 7), "emri-disk", "results", "emri-gasdisk",
        "paper_EMRI-GASDISK", "codes", "cas", "deriv-master-system-R2", "master_fixture.json")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "r2_diag_fixture.jl")
raw = open(src, "rb").read()
sha = hashlib.sha256(raw).hexdigest()
d = json.loads(raw)
pts = d["points"]
assert len(pts) == 50, len(pts)

BG = ["nu", "nup", "psi", "psip", "mu", "mup", "omt", "omtp", "v", "vp", "eps", "epsp", "p", "pp",
      "Gam1", "omd", "htt", "http", "htp", "htpp", "hpp", "hppp", "hrr", "hrrp", "hzz", "hzzp", "Omp"]


def om_key(bg):
    for k in bg:
        if "omt + E" in k.replace("\n", " "):
            return k
    raise KeyError("Omega key")


rows = []
for p in pts:
    bg = p["background"]
    fields = ["index = %d" % p["index"]]
    for k in ("chi", "aEMRI", "h0", "pdisk", "qdisk"):
        fields.append('%s = "%s"' % (k, p[k]))
    fields.append("m = %d" % p["m"])
    for k in ("r", "OmegaEMRI", "sqrtmg", "u_phi"):
        fields.append('%s = "%s"' % (k, p[k]))
    fields.append('Omega = "%s"' % bg[om_key(bg)])
    for k in BG:
        fields.append('%s = "%s"' % (k, bg[k]))
    for k in ("sample_y0", "sample_y1", "R2_AMF_at_sample", "R2_TorqueDensity_at_sample"):
        fields.append('%s = ("%s", "%s")' % ((k,) + tuple(p[k])))
    c = p["R2_deltaUphi_coeffs_y0_y1_h"]
    fields.append('deltaUphi_coeffs = (("%s", "%s"), ("%s", "%s"), ("%s", "%s"))'
                  % tuple(x for pair in c for x in pair))
    # 3 fields per line keeps the generated file diff-friendly
    body = ",\n     ".join(", ".join(fields[i:i + 4]) for i in range(0, len(fields), 4))
    rows.append("    (" + body + ")")

hdr = ("# r2_diag_fixture.jl -- GENERATED FILE, DO NOT EDIT BY HAND (AccretionDisc.SpiralWaves test fixture).\n"
       "# baked by : test/SpiralWaves/fixtures/capture/bake_r2_diag_fixture.py\n"
       "# source   : emri-disk results/emri-gasdisk/paper_EMRI-GASDISK/codes/cas/deriv-master-system-R2/master_fixture.json\n"
       "#            sha256 %s  (generator master_fixture.wls, seed %s, %s digits; EGD::R1-vs-R2-gate PASS)\n"
       "# content  : 50 on-shell points -- the full route-R2 input set at the field point (nu, psi, mu, psip, v, vp, omt,\n"
       "#            eps, epsp, p, pp, Gam1, omega_d + the COVARIANT test-profile modes htt, htp, hpp, hrr, hzz and their\n"
       "#            r-derivatives), sqrt(-g) and u_phi, the sample master vector y = (sample_y0, sample_y1) and\n"
       "#            Mathematica's deltaUphi (as coefficients of y0, y1 and the h-term), AMF (72a) and TorqueDensity (72b).\n"
       "#            All numbers are decimal STRINGS; complex values are (re, im) string pairs.\n"
       "const R2_DIAG_FIXTURE = [\n") % (sha, d["meta"]["seed"], d["meta"]["digits"])
text = hdr + ",\n".join(rows) + "\n]\n"
open(out, "w").write(text)
print("source sha256", sha)
print("wrote", os.path.abspath(out), "bytes", len(text.encode()))
