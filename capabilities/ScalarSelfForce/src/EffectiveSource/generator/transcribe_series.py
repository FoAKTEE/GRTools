#!/usr/bin/env python3
"""transcribe_series.py — MECHANICAL C -> Julia transcription of the near-worldline
power-series effective source (ESSF node `effective-source-series`, schedule W8-B).

Source of truth (READ-ONLY mirror):
    ref-code/EffectiveSource/kerr-circular-comoving.c, function
    `effsource_calc_m_series` (lines 643-677 at the pinned SHA
    1cb9834f8a8f666aaf7b1f924b3197743efb4120).

Output (DO NOT EDIT BY HAND):
    src/EffectiveSource/series_coeffs.jl -- `_series_coeffs(M, a, r0, Om, zc, mm)`
    returning the ELEVEN named coefficients of arXiv:2606.04998 main.tex:493-504.

The C writes the series as `srcO0 + srcO1 + srcO2` with the lambda structure inlined.
This script DIVIDES OUT that structure, term by term, so each named coefficient is a
function of (M, a, r0, Omega, z_c, m) only:

    srcO0                              -> S0
    srcO1            /  dr             -> S1_1
    srcO2 term 1     /  rho^2          -> S2_0rat     (the non-harmonic part of S^[2]_0)
    srcO2 term 2     /  dr^2           -> S2_2rat     (the non-harmonic part of S^[2]_2)
    srcO2 terms 3-7  /  dr^{2k}/rho^{2k-2}, k=2..6    -> S2_4, S2_6, S2_8, S2_10, S2_12
    H_{-m-1/2} block /  {rho^2, dr^2}  -> S2_0L, S2_2L

and ASSERTS that the C's `log(rho^2/zc^2)` block is the H block with both denominators
doubled -- i.e. the C's combination is exactly

    (S2_0L*rho^2 + S2_2L*dr^2) * (H_{-m-1/2} + log(rho/zc)) ,

which is the paper's `(S^[2]_{0L} rho^2 + S^[2]_{2L} dr^2) log rho` with the
`H_{-m-1/2} - log z_c` remainder absorbed into S^[2]_0 / S^[2]_2 (main.tex:493-504).

Mechanical rewrite rules (the COMPLETE set; term order and parenthesization are otherwise
preserved verbatim so the Float64 evaluation tracks the C oracle):
    1. `pow(X, n)` -> `X^n` (bare identifier base) or `(X)^n`, innermost first.
    2. C trailing-dot double literals (`24.`, `48.`, `96.`, `3.`, `32.`, `64.`) -> the same
       exact integer literal (keeps the Julia code type-generic).
    3. C symbol renames: `rp` -> `r0` (the C's `rp` is the PARTICLE radius, not the horizon),
       `Omega` -> `Om`, `m` -> `mm`.

Usage:  python3 transcribe_series.py <out.jl> [<kerr-circular-comoving.c>]
"""
import hashlib
import io
import re
import sys
from pathlib import Path

DEFAULT_CSRC = Path("/data/haiyangw/claude/emri-disk/ref-code/EffectiveSource/"
                    "kerr-circular-comoving.c")
# sha256 of kerr-circular-comoving.c at the pinned upstream commit
# 1cb9834f8a8f666aaf7b1f924b3197743efb4120. The generator REFUSES to emit from any other
# file: the surgery table below is a set of exact substring matches against those very
# bytes, so a drifted mirror must fail loudly rather than half-match.
CSRC_SHA256 = "5b4f5647e0c673a29c0be9637804aba103d8153194617caa28951598ccf75d1f"
LO, HI = 643, 677                      # 1-based inclusive: the whole function
IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
GUARD = "(dr==0&&dth==0)? 0 : "

# name -> (source line offset within the function, sign, the single substitution that
# divides out the lambda structure).  Offsets are 0-based from the `void` line.
SURGERY = [
    ("S0",      15, "+", None),
    ("S1_1",    16, "+", ("*zc*dr)", "*zc)")),
    ("S2_0rat", 17, "+", ("pow(rho,2)*", "")),
    ("S2_2rat", 18, "+", ("*zc*pow(dr,2))", "*zc)")),
    ("S2_4",    19, "+", ("*zc*pow(dr,4))/(96.*pow(rho,2)*", "*zc)/(96.*")),
    ("S2_6",    20, "-", ("*zc*pow(dr,6))/(48.*pow(rho,4)*", "*zc)/(48.*")),
    ("S2_8",    21, "-", ("*zc*pow(dr,8))/(3.*pow(rho,6)*", "*zc)/(3.*")),
    ("S2_10",   22, "+", ("*zc*pow(dr,10))/(pow(rho,8)*", "*zc)/(")),
    ("S2_12",   23, "-", ("*zc*pow(dr,12))/(pow(rho,10)*", "*zc)/(")),
]


def sub1(s, old, new):
    n = s.count(old)
    assert n == 1, "expected exactly 1 occurrence of %r, found %d" % (old, n)
    return s.replace(old, new)


def split_top(s, sep=" + "):
    out, depth, last, i = [], 0, 0, 0
    while i < len(s):
        c = s[i]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif depth == 0 and s.startswith(sep, i):
            out.append(s[last:i])
            last = i + len(sep)
            i += len(sep) - 1
        i += 1
    out.append(s[last:])
    return out


def convert_pow(expr):
    while True:
        i = expr.find("pow(")
        if i < 0:
            return expr
        depth, comma = 0, None
        for j in range(i + 3, len(expr)):
            c = expr[j]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            elif c == "," and depth == 1:
                comma = j
        assert depth == 0 and comma is not None, "unbalanced pow at %d" % i
        base = convert_pow(expr[i + 4:comma].strip())
        exp = expr[comma + 1:j].strip()
        assert re.fullmatch(r"\d+", exp), "non-integer exponent %r" % exp
        rep = "%s^%s" % (base, exp) if IDENT.match(base) else "(%s)^%s" % (base, exp)
        expr = expr[:i] + rep + expr[j + 1:]


def to_julia(expr):
    expr = convert_pow(expr)
    expr = re.sub(r"(?<![\w.])(\d+)\.(?!\d)", r"\1", expr)
    expr = re.sub(r"\brp\b", "r0", expr)
    expr = re.sub(r"\bOmega\b", "Om", expr)
    expr = re.sub(r"\bm\b", "mm", expr)
    assert "pow" not in expr and expr.count("(") == expr.count(")")
    assert not re.search(r"\b(rho|dr|dth|th|alpha)\b", expr), expr[:200]
    return expr.strip()


def main(outpath, csrc=DEFAULT_CSRC):
    raw = Path(csrc).read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != CSRC_SHA256:
        raise SystemExit(
            "kerr-circular-comoving.c sha256 mismatch: expected %s, got %s (%s).\n"
            "The pinned upstream is EffectiveSource @ 1cb9834f8a8f666aaf7b1f924b3197743efb4120."
            % (CSRC_SHA256, got, csrc))
    self_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    lines = raw.decode().split("\n")[LO - 1:HI]
    assert lines[0].startswith("void effsource_calc_m_series("), lines[0]

    def body(idx, prefix):
        s = lines[idx].strip()
        assert s.startswith(prefix), (idx, s[:60])
        s = s[len(prefix):].strip()
        if s.endswith(";"):
            s = s[:-1].strip()
        if s.startswith(GUARD):
            s = s[len(GUARD):].strip()
        return s

    prefixes = {15: "const double srcO0 =", 16: "const double srcO1 =",
                17: "const double srcO2 ="}
    raw = {}
    for name, idx, sign, sub in SURGERY:
        e = body(idx, prefixes.get(idx, sign + " "))
        raw[name] = (("-" if sign == "-" else "") +
                     (sub1(e, *sub) if sub else e))

    # ---- the H_{-m-1/2} block and its log twin -----------------------------------
    hp = "+ HarmonicNumber(-0.5 - m) * "
    XH = body(24, hp)
    assert XH.startswith("(") and XH.endswith(")")
    hparts = split_top(XH[1:-1])
    assert len(hparts) == 2
    raw["S2_0L"] = sub1(hparts[0], "pow(rho,2)*", "")
    raw["S2_2L"] = sub1(hparts[1], "*zc*pow(dr,2))", "*zc)")

    lp = "+ log(pow(rho,2)/pow(zc,2)) * "
    XL = body(25, lp)
    assert XL.startswith("(") and XL.endswith(")")
    lparts = split_top(XL[1:-1])
    assert len(lparts) == 2
    assert sub1(lparts[0], "pow(rho,2)*", "") == \
        sub1(raw["S2_0L"], ")/(64.*", ")/(128.*"), "log block != H block/2 (rho^2 part)"
    assert sub1(lparts[1], "*zc*pow(dr,2))", "*zc)") == \
        sub1(raw["S2_2L"], ")/(32.*", ")/(64.*"), "log block != H block/2 (dr^2 part)"

    names = ["S0", "S1_1", "S2_0rat", "S2_2rat", "S2_0L", "S2_2L",
             "S2_4", "S2_6", "S2_8", "S2_10", "S2_12"]
    jl = {n: to_julia(raw[n]) for n in names}

    # ---- the paper's m-dependence claims, asserted on the transcription ----------
    for n in ("S0", "S1_1"):
        occ = re.findall(r"\bmm\b", jl[n])
        assert len(occ) == 2 and jl[n].count("(-1 + 2*mm)*(1 + 2*mm)") == 1, \
            "%s must depend on m ONLY through 1/((2m-1)(2m+1))" % n
    for n in ("S2_0rat", "S2_2rat"):
        assert jl[n].count("(-1 + 2*mm)*(1 + 2*mm)") == 1 and "mm^2" in jl[n]
    for n in ("S2_0L", "S2_2L", "S2_4", "S2_6", "S2_8", "S2_10", "S2_12"):
        assert not re.search(r"\bmm\b", jl[n]), "%s must be m-INDEPENDENT" % n

    out = io.StringIO()
    out.write(('''\
# series_coeffs.jl — GENERATED FILE, DO NOT EDIT BY HAND.
# Regenerate (byte-identical) with:
#   python3 src/EffectiveSource/generator/transcribe_series.py \\
#           src/EffectiveSource/series_coeffs.jl <path-to>/kerr-circular-comoving.c
#
# PROVENANCE (both hashes verified by the generator at emit time):
#   input  kerr-circular-comoving.c sha256 %s
#   generator transcribe_series.py  sha256 %s
#
# Mechanical transcription of `effsource_calc_m_series`,
# ref-code/EffectiveSource/kerr-circular-comoving.c:643-677 (READ-ONLY mirror, pinned SHA
# 1cb9834f8a8f666aaf7b1f924b3197743efb4120; https://github.com/barrywardell/EffectiveSource,
# MIT license), with the λ-structure DIVIDED OUT so each name is one coefficient of the
# near-worldline series of arXiv:2606.04998 main.tex:493-504. The C's `rp` is the PARTICLE
# radius (renamed `r0` here; the package's `orbit.rp` is the OUTER HORIZON), `zc` = z_c of
# eq:s2, `Om` = Ω, `mm` = |m|. Rewrite rules and the surgery table: the generator header.
#
# The paper's m-dependence (main.tex:503-504) is ASSERTED by the generator on the
# transcription itself: S^[0] and S^[1]_1 carry m ONLY through 1/((2m−1)(2m+1)); S^[2]_0 and
# S^[2]_2 are rational in m² over the same factor (their harmonic-number parts are the
# SEPARATE S2_0L/S2_2L blocks, added in series.jl); every other coefficient is m-independent.

"""
    _series_coeffs(M, a, r0, Om, zc, mm) -> NamedTuple

The eleven near-worldline series coefficients (main.tex:493-504) for a circular equatorial
Kerr geodesic: `M`/`a` the background mass and spin, `r0` the orbital radius, `Om` = Ω the
orbital frequency, `zc` = z_c (`sqrt(orbit.zc2)`, eq:s2), `mm` = |m| ≥ 0.

`S2_0rat`/`S2_2rat` are the NON-harmonic parts of S^[2]_0/S^[2]_2; `S2_0L`/`S2_2L` multiply
`log ϱ` and (as `S2_*L·(H_{−m−1/2} − log z_c)`) complete S^[2]_0/S^[2]_2 —
[`effective_source_series_coefficients`](@ref) does that assembly. Type-generic (integer
literals only), so BigFloat arguments give full-precision BigFloat coefficients.
"""
function _series_coeffs(M, a, r0, Om, zc, mm)
''') % (CSRC_SHA256, self_sha))
    for n in names:
        out.write("    %s = %s\n" % (n, jl[n]))
    out.write("    return (; " + ", ".join(names) + ")\nend\n")
    Path(outpath).write_text(out.getvalue())
    print("wrote %s (%d coefficients, %d bytes)"
          % (outpath, len(names), len(out.getvalue())))


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
