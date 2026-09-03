#!/usr/bin/env python3
"""capture_series_rational.py — EXACT-RATIONAL reference values for the near-worldline
series coefficients of `effsource_calc_m_series`
(ref-code/EffectiveSource/kerr-circular-comoving.c:643-677, pinned SHA 1cb9834f).

Baked output:  test/fixtures/effsource/series_coeffs_rational.jl

WHY IT EXISTS (the "double transcription" gate of ESSF node `effective-source-series`).
The shipped `src/EffectiveSource/series_coeffs.jl` is Julia CODE GENERATED from the C text
by `src/EffectiveSource/generator/transcribe_series.py` (pow(X,n) -> X^n, `24.` -> `24`,
`rp`->`r0`, `Omega`->`Om`, `m`->`mm`) after a TEXTUAL surgery that divides out the
λ-structure. This script reaches the SAME eleven numbers by two INDEPENDENT routes:

  * ARITHMETIC: the C expression text is evaluated DIRECTLY as Python (the C body is valid
    Python once `24.`-style literals are made integers and the `?:` guard is stripped) over
    `fractions.Fraction`, so no pow/precedence/literal rewrite is shared with the Julia side.
  * EXTRACTION: instead of the generator's textual surgery, the coefficients are obtained by
    SUBSTITUTION AND EXACT LINEAR SOLVE. With `zc = 1` and `rho = 1` the C's
    `log(rho^2/zc^2)` block vanishes identically and

        srcO2(dr) = [A0 + A2 dr^2 + A4 dr^4 + ... + A12 dr^12] + H·[X0 + X2 dr^2] ,

    which is read off from 7 values of `dr` by an exact Vandermonde solve; `H` is carried
    as a formal symbol by a two-component (rational, H-coefficient) number type, so nothing
    irrational ever enters. `srcO0` is S0; `srcO1(dr=1)` is S1_1.

Evaluation points: rational (a, r0) with r0 a perfect square, so the geodesic frequency
Omega = M/(a M + sqrt(M r0^3)) is EXACTLY rational and every coefficient is a rational
number (each C term carries exactly ONE linear factor of z_c, hence zc = 1 gives
coefficient / z_c). These points are ALGEBRAIC probes, not physical orbits — (0, 4) is
inside the Schwarzschild ISCO; only the transcription is under test.

Usage:  python3 capture_series_rational.py <out.jl> [<kerr-circular-comoving.c>]
"""
import re
import sys
from fractions import Fraction as F
from pathlib import Path

DEFAULT_CSRC = Path("/data/haiyangw/claude/emri-disk/ref-code/EffectiveSource/"
                    "kerr-circular-comoving.c")
LO, HI = 643, 677
GUARD = "(dr==0&&dth==0)? 0 : "

# (a, r0) with r0 = k^2 so that Omega = 1/(a + r0^(3/2)) is exactly rational (M = 1).
POINTS = [(F(0), F(4)), (F(1, 2), F(9)), (F(-1, 2), F(16))]
MS = [0, 1, 2, 5, 10, 20]
NAMES = ["S0", "S1_1", "S2_0rat", "S2_2rat", "S2_0L", "S2_2L",
         "S2_4", "S2_6", "S2_8", "S2_10", "S2_12"]


class HNum:
    """x = a + b*H with a, b exact rationals and H = HarmonicNumber(-1/2-m) a symbol.
    Only linear-in-H arithmetic occurs in the C expression, which is asserted."""

    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a = F(a)
        self.b = F(b)

    @staticmethod
    def _c(o):
        return o if isinstance(o, HNum) else HNum(o)

    def __add__(self, o):
        o = HNum._c(o)
        return HNum(self.a + o.a, self.b + o.b)
    __radd__ = __add__

    def __neg__(self):
        return HNum(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-HNum._c(o))

    def __rsub__(self, o):
        return HNum._c(o) + (-self)

    def __mul__(self, o):
        o = HNum._c(o)
        assert self.b == 0 or o.b == 0, "H appears nonlinearly"
        return HNum(self.a * o.a, self.a * o.b + self.b * o.a)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = HNum._c(o)
        assert o.b == 0, "division by an H-carrying value"
        return HNum(self.a / o.a, self.b / o.a)

    def __rtruediv__(self, o):
        assert self.b == 0, "division by an H-carrying value"
        return HNum(F(o) / self.a)

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        r = HNum(1)
        for _ in range(n):
            r = r * self
        return r


def c_body(csrc):
    lines = Path(csrc).read_text().split("\n")[LO - 1:HI]
    assert lines[0].startswith("void effsource_calc_m_series(")
    def get(i, prefix):
        s = lines[i].strip()
        assert s.startswith(prefix), (i, s[:60])
        s = s[len(prefix):].strip()
        return s[:-1].strip() if s.endswith(";") else s
    e0 = get(15, "const double srcO0 =")
    e1 = get(16, "const double srcO1 =")
    e2 = get(17, "const double srcO2 =")
    assert e1.startswith(GUARD) and e2.startswith(GUARD)
    e1 = e1[len(GUARD):]
    e2 = "".join([e2[len(GUARD):]] + [l.strip() + " " for l in lines[18:26]])
    assert e2.rstrip().endswith(";")
    e2 = e2.rstrip()[:-1]
    # the ONLY rewrite: C trailing-dot double literals -> exact integers (Python `24.` is a
    # float and would poison the Fraction arithmetic). `pow`, `log` and the operator
    # precedence of the C text are used AS IS by the Python evaluator.
    return [re.sub(r"(?<![\w.])(\d+)\.(?!\d)", r"\1", e) for e in (e0, e1, e2)]


def make_env(a, r0, m, dr):
    """Evaluation environment: zc = 1, rho = 1 (so log(rho^2/zc^2) == 0), dth = 0."""
    Om = F(1) / (a + F(int(round(float(r0) ** 1.5))))
    assert Om == F(1) / (a + F(int(round(float(r0) ** 1.5)))) and \
        int(round(float(r0) ** 1.5)) ** 2 == int(r0) ** 3, "r0^(3/2) must be exact"
    return dict(M=HNum(1), a=HNum(a), rp=HNum(r0), Omega=HNum(Om), m=m,
                zc=HNum(1), rho=HNum(1), dr=HNum(dr), dth=HNum(0),
                pow=lambda x, n: HNum._c(x) ** n,
                log=lambda x: (_assert_zero(x)),
                HarmonicNumber=lambda x: _harm(x, m))


def _assert_zero(x):
    assert isinstance(x, HNum) and x.a == 1 and x.b == 0, "log() argument must be 1 here"
    return HNum(0)


def _harm(x, m):
    assert x == F(-1, 2) - m, "unexpected HarmonicNumber argument"
    return HNum(0, 1)


def solve_vandermonde(pts, vals):
    """Exact solve of sum_k c_k * x^k == v over Fractions (Gaussian elimination)."""
    n = len(pts)
    A = [[p ** k for k in range(n)] + [v] for p, v in zip(pts, vals)]
    for i in range(n):
        piv = next(r for r in range(i, n) if A[r][i] != 0)
        A[i], A[piv] = A[piv], A[i]
        A[i] = [x / A[i][i] for x in A[i]]
        for r in range(n):
            if r != i and A[r][i] != 0:
                f = A[r][i]
                A[r] = [x - f * y for x, y in zip(A[r], A[i])]
    return [A[i][n] for i in range(n)]


def coeffs(exprs, a, r0, m):
    e0, e1, e2 = exprs
    S0 = eval(e0, {}, make_env(a, r0, m, 0))
    assert S0.b == 0
    S1 = eval(e1, {}, make_env(a, r0, m, 1))
    assert S1.b == 0
    # srcO2(dr) at rho = zc = 1 is  P(dr^2) + H * Q(dr^2), deg P = 6, deg Q = 1 in dr^2
    xs, pv, qv = [], [], []
    for d in range(7):
        v = eval(e2, {}, make_env(a, r0, m, d))
        xs.append(F(d) ** 2)
        pv.append(v.a)
        qv.append(v.b)
    A = solve_vandermonde(xs, pv)                    # A[k] multiplies dr^(2k)
    X = solve_vandermonde(xs[:2], qv[:2])            # X[k] multiplies dr^(2k)
    assert all(qv[i] == X[0] + X[1] * xs[i] for i in range(7)), "H block is not linear in dr^2"
    return dict(S0=S0.a, S1_1=S1.a, S2_0rat=A[0], S2_2rat=A[1], S2_0L=X[0], S2_2L=X[1],
                S2_4=A[2], S2_6=A[3], S2_8=A[4], S2_10=A[5], S2_12=A[6])


def main(outpath, csrc=DEFAULT_CSRC):
    exprs = c_body(csrc)
    rows = []
    for ip, (a, r0) in enumerate(POINTS, start=1):
        for m in MS:
            c = coeffs(exprs, a, r0, m)
            for n in NAMES:
                v = c[n]
                rows.append('    (%d, %d, :%s, "%d", "%d"),' %
                            (ip, m, n, v.numerator, v.denominator))
    pts = ",\n".join(
        '    (a = %d // %d, r0 = %d // %d, Omega = %d // %d)'
        % (a.numerator, a.denominator, r0.numerator, r0.denominator,
           (F(1) / (a + F(int(round(float(r0) ** 1.5))))).numerator,
           (F(1) / (a + F(int(round(float(r0) ** 1.5))))).denominator)
        for a, r0 in POINTS)
    Path(outpath).write_text('''\
# series_coeffs_rational.jl — BAKED EXACT-RATIONAL reference for the near-worldline series
# coefficients (ESSF node `effective-source-series`, W8-B). DO NOT EDIT BY HAND.
# Regenerate (byte-identical, no third-party Python packages needed):
#   python3 test/fixtures/effsource/generator/capture_series_rational.py \\
#           test/fixtures/effsource/series_coeffs_rational.jl \\
#           <path-to>/kerr-circular-comoving.c
#
# CITATION: https://github.com/barrywardell/EffectiveSource
#   pinned SHA 1cb9834f8a8f666aaf7b1f924b3197743efb4120 (2026-04-11, MIT license)
#   read-only mirror: ref-code/EffectiveSource (PROVENANCE.elliptic-ssf.md)
# ORACLE: `effsource_calc_m_series` (kerr-circular-comoving.c:643-677) evaluated EXACTLY in
#   Python `fractions.Fraction` — the C body is used as Python source verbatim except that
#   trailing-dot double literals (`24.`) are made integers; `pow`, `log` and the C operator
#   precedence are consumed as is. The eleven coefficients are extracted NOT by the shipped
#   generator's textual surgery but by SUBSTITUTION + exact Vandermonde solve at zc = rho = 1
#   (where the C's `log(rho^2/zc^2)` block vanishes) with `HarmonicNumber(-1/2-m)` carried as
#   a formal symbol. This is therefore an INDEPENDENT second transcription of the same C
#   text: it shares no rewrite rule and no extraction step with
#   src/EffectiveSource/generator/transcribe_series.py.
#
# CONVENTION: values are the coefficients DIVIDED BY z_c (every C term carries exactly one
#   linear z_c factor), i.e. exactly `_series_coeffs(1, a, r0, Omega, 1, m)`.
#
# GRID: rational (a, r0) with r0 a perfect square so Omega = 1/(a + r0^(3/2)) is exactly
#   rational (M = 1). These are ALGEBRAIC PROBES, not physical orbits — (a, r0) = (0, 4) is
#   inside the Schwarzschild ISCO; only the transcription is under test. m in %s.
#
# COLUMNS: (ipoint, m, coefficient name, numerator, denominator) — decimal strings, parsed
#   into `Rational{BigInt}` by the test.
const EFFSOURCE_SERIES_RATIONAL_POINTS = [
%s,
]

const EFFSOURCE_SERIES_RATIONAL = [
%s
]
''' % (MS, pts, "\n".join(rows)))
    print("wrote %s (%d rows)" % (outpath, len(rows)))


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
