#!/usr/bin/env python3
"""bake_master_onshell_fixture.py -- bake the ON-SHELL master-system fixture for test_master.jl.

Source: the R1-vs-R2 gate fixture of the emri-disk mission repository (EGD::R1-vs-R2-gate),
    results/emri-gasdisk/paper_EMRI-GASDISK/codes/cas/deriv-master-system-R2/master_fixture.json
    (generator master_fixture.wls, seed 20260819, 30 significant digits from 40-digit arithmetic;
    sha256 8db949d8a48ec4fdb880decd2fb0383ee202f7135a002006c762c6bcf3ff2c56, 207 KB).
50 physical points (chi, a_EMRI, h0, p, q, m, r): Kerr (A1) background, power-law disc (12)-(13), (22),
v(r) = prograde root of the hydrostatic balance (11), covariant test profile
    h_tt = exp(-((r - a)/(a/2))^2), h_tphi = 0.3 r^2 Omega_EMRI h_tt, h_phiphi = 0.7 r^2 h_tt,
    h_rr = 0.5 h_tt, h_zz = 0.2 h_tt,
and the master entries C00, C01, C10, C11, S0, S1 of arXiv:2608.17003 eq (47) by BOTH routes (R1 notebook
transcription == R2 independent derivation, exact zero difference on-shell). This script keeps only what
test_master.jl needs (parameters, Omega, Omega', v, omega_d, and the six R1 entries) as decimal STRINGS
(parse at Float64 or BigFloat precision) -> master_onshell_fixture.jl (< 50 KB, tracked beside this script).
Usage: python3 bake_master_onshell_fixture.py [path/to/master_fixture.json] [output .jl]
"""
import hashlib, json, os, sys

here = os.path.dirname(os.path.abspath(__file__))
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, *([".."] * 7), "emri-disk", "results", "emri-gasdisk",
        "paper_EMRI-GASDISK", "codes", "cas", "deriv-master-system-R2", "master_fixture.json")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "master_onshell_fixture.jl")
raw = open(src, "rb").read()
sha = hashlib.sha256(raw).hexdigest()
d = json.loads(raw)
pts = d["points"]
assert len(pts) == 50, len(pts)

def om_key(bg):
    # the fixture stores Omega under a typeset key ("omt + E^(nu-psi) v" with a line break); find it robustly
    for k in bg:
        if "omt + E" in k.replace("\n", " "):
            return k
    raise KeyError("Omega key")

rows = []
for p in pts:
    bg = p["background"]
    okey = om_key(bg)
    r1 = p["R1"]
    rows.append(
        "    (index = %d, chi = \"%s\", aEMRI = \"%s\", h0 = \"%s\", pdisk = \"%s\", qdisk = \"%s\", m = %d, r = \"%s\",\n"
        "     Omega = \"%s\", Omega_r = \"%s\", v = \"%s\", omega_d = \"%s\", OmegaEMRI = \"%s\",\n"
        "     C00 = (\"%s\", \"%s\"), C01 = (\"%s\", \"%s\"), C10 = (\"%s\", \"%s\"), C11 = (\"%s\", \"%s\"), S0 = (\"%s\", \"%s\"), S1 = (\"%s\", \"%s\"))"
        % (p["index"], p["chi"], p["aEMRI"], p["h0"], p["pdisk"], p["qdisk"], p["m"], p["r"],
           bg[okey], bg["Omp"], bg["v"], bg["omd"], p["OmegaEMRI"],
           *r1["C00"], *r1["C01"], *r1["C10"], *r1["C11"], *r1["S0"], *r1["S1"]))

hdr = ("# master_onshell_fixture.jl -- GENERATED FILE, DO NOT EDIT BY HAND (AccretionDisc.SpiralWaves test fixture).\n"
       "# baked by : test/SpiralWaves/fixtures/capture/bake_master_onshell_fixture.py\n"
       "# source   : emri-disk results/emri-gasdisk/paper_EMRI-GASDISK/codes/cas/deriv-master-system-R2/master_fixture.json\n"
       "#            sha256 %s  (generator master_fixture.wls, seed %s, %s digits; EGD::R1-vs-R2-gate PASS: R1 == R2 on-shell)\n"
       "# content  : 50 on-shell points (chi, aEMRI, h0, pdisk, qdisk, m, r) of the arXiv:2608.17003 disc: Omega(r), dOmega/dr, v from the\n"
       "#            exact hydrostatic balance (11), omega_d = m(Omega_EMRI(aEMRI) - Omega(r)), and the master entries C00..S1 of eq (47)\n"
       "#            for the COVARIANT test profile h_tt = exp(-((r-a)/(a/2))^2), h_tphi = 0.3 r^2 Omega_EMRI h_tt, h_phiphi = 0.7 r^2 h_tt,\n"
       "#            h_rr = 0.5 h_tt, h_zz = 0.2 h_tt (derivatives analytic). All numbers are decimal STRINGS; complex as (re, im).\n"
       "const MASTER_ONSHELL_FIXTURE = [\n") % (sha, d["meta"]["seed"], d["meta"]["digits"])
text = hdr + ",\n".join(rows) + "\n]\n"
open(out, "w").write(text)
print("source sha256", sha)
print("wrote", os.path.abspath(out), "bytes", len(text.encode()))
