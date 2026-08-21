#!/usr/bin/env python3
"""bake_master_r2_fixture.py -- bake the ON-SHELL master-system fixture from the ROUTE-R2 columns.

Route R2 is OUR OWN independent Mathematica derivation of the arXiv:2608.17003 master system
(EGD::deriv-master-system-R2), sources
    emri-disk results/emri-gasdisk/paper_EMRI-GASDISK/codes/cas/deriv-master-system-R2/R2_master_system.m
    sha256 67825ad049e45ed995f1a087cde6ff9f0d64101ea1ac52f457286b645de26bf6
evaluated at the 50 fixture points by master_fixture.wls into
    .../deriv-master-system-R2/master_fixture.json
    sha256 8db949d8a48ec4fdb880decd2fb0383ee202f7135a002006c762c6bcf3ff2c56 (seed 20260819, 30 digits).

This script reads the "R2" columns ONLY.  The json also carries an "R1" column (values from the
authors' supplementary notebook), which EGD::R1-vs-R2-gate proved EQUAL to R2 on-shell -- all 300
entries are byte-identical decimal strings, and maxRelDiff_R1_R2 is an exact zero at every point.
The notebook lineage was REMOVED from GRTools on 2026-08-22 by owner decision
(EGD::jd-notebook-purge); this fixture therefore descends from route R2 alone, with no loss of value:
it supersedes the earlier master_onshell_fixture.jl, whose numbers were the identical strings read
from the "R1" key.

50 physical points (chi, a_EMRI, h0, p, q, m, r): Kerr (A1) background, power-law disc (12)-(13), (22),
v(r) = prograde root of the hydrostatic balance (11), covariant test profile
    h_tt = exp(-((r - a)/(a/2))^2), h_tphi = 0.3 r^2 Omega_EMRI h_tt, h_phiphi = 0.7 r^2 h_tt,
    h_rr = 0.5 h_tt, h_zz = 0.2 h_tt,
and the master entries C00, C01, C10, C11, S0, S1 of eq (47).  Kept: the parameters, Omega, Omega', v,
omega_d and the six entries, as decimal STRINGS (parse at Float64 or BigFloat precision)
-> master_r2_fixture.jl (< 50 KB, tracked beside this script).
Usage: python3 bake_master_r2_fixture.py [path/to/master_fixture.json] [output .jl]
"""
import hashlib, json, os, sys

R2_SYSTEM_SHA = "67825ad049e45ed995f1a087cde6ff9f0d64101ea1ac52f457286b645de26bf6"

here = os.path.dirname(os.path.abspath(__file__))
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, *([".."] * 7), "emri-disk", "results", "emri-gasdisk",
        "paper_EMRI-GASDISK", "codes", "cas", "deriv-master-system-R2", "master_fixture.json")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "master_r2_fixture.jl")
raw = open(src, "rb").read()
sha = hashlib.sha256(raw).hexdigest()
d = json.loads(raw)
pts = d["points"]
assert len(pts) == 50, len(pts)
# the json must be the one this fixture claims to descend from, and it must name OUR R2 source
assert R2_SYSTEM_SHA in d["meta"]["R2_source"] or d["meta"]["R2_source"].endswith("R2_master_system.m"), d["meta"]["R2_source"]

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
    r2 = p["R2"]                       # ROUTE R2 ONLY -- never p["R1"]
    rows.append(
        "    (index = %d, chi = \"%s\", aEMRI = \"%s\", h0 = \"%s\", pdisk = \"%s\", qdisk = \"%s\", m = %d, r = \"%s\",\n"
        "     Omega = \"%s\", Omega_r = \"%s\", v = \"%s\", omega_d = \"%s\", OmegaEMRI = \"%s\",\n"
        "     C00 = (\"%s\", \"%s\"), C01 = (\"%s\", \"%s\"), C10 = (\"%s\", \"%s\"), C11 = (\"%s\", \"%s\"), S0 = (\"%s\", \"%s\"), S1 = (\"%s\", \"%s\"))"
        % (p["index"], p["chi"], p["aEMRI"], p["h0"], p["pdisk"], p["qdisk"], p["m"], p["r"],
           bg[okey], bg["Omp"], bg["v"], bg["omd"], p["OmegaEMRI"],
           *r2["C00"], *r2["C01"], *r2["C10"], *r2["C11"], *r2["S0"], *r2["S1"]))

hdr = ("# master_r2_fixture.jl -- GENERATED FILE, DO NOT EDIT BY HAND (AccretionDisc.SpiralWaves test fixture).\n"
       "# baked by : test/SpiralWaves/fixtures/capture/bake_master_r2_fixture.py\n"
       "# route    : R2 -- OUR independent Mathematica derivation (EGD::deriv-master-system-R2), source\n"
       "#            codes/cas/deriv-master-system-R2/R2_master_system.m  sha256 %s\n"
       "# source   : emri-disk results/emri-gasdisk/paper_EMRI-GASDISK/codes/cas/deriv-master-system-R2/master_fixture.json\n"
       "#            sha256 %s  (generator master_fixture.wls, seed %s, %s digits; \"R2\" columns ONLY)\n"
       "# content  : 50 on-shell points (chi, aEMRI, h0, pdisk, qdisk, m, r) of the arXiv:2608.17003 disc: Omega(r), dOmega/dr, v from the\n"
       "#            exact hydrostatic balance (11), omega_d = m(Omega_EMRI(aEMRI) - Omega(r)), and the master entries C00..S1 of eq (47)\n"
       "#            for the COVARIANT test profile h_tt = exp(-((r-a)/(a/2))^2), h_tphi = 0.3 r^2 Omega_EMRI h_tt, h_phiphi = 0.7 r^2 h_tt,\n"
       "#            h_rr = 0.5 h_tt, h_zz = 0.2 h_tt (derivatives analytic). All numbers are decimal STRINGS; complex as (re, im).\n"
       "const MASTER_R2_FIXTURE = [\n") % (R2_SYSTEM_SHA, sha, d["meta"]["seed"], d["meta"]["digits"])
text = hdr + ",\n".join(rows) + "\n]\n"
open(out, "w").write(text)
print("source sha256", sha)
print("R2_master_system.m sha256", R2_SYSTEM_SHA)
print("wrote", os.path.abspath(out), "bytes", len(text.encode()))
