#!/usr/bin/env python3
"""transcribe_field4d.py — MECHANICAL C -> Julia transcription of the 4D Detweiler-Whiting
puncture A_ijk coefficient block (ESSF node `puncture-field-4d`, schedule W3-D).

Source of truth (READ-ONLY mirror):
    ref-code/EffectiveSource/kerr-circular.c, function effsource_set_particle,
    the `/* Compute A coefficients */` block (lines 792-841 at the pinned SHA
    1cb9834f8a8f666aaf7b1f924b3197743efb4120 -- see fixture headers under
    test/fixtures/effsource/ for the full citation).

Output (DO NOT EDIT BY HAND; regenerate with this script):
    src/Puncture/field4d_coeffs.jl  -- `_puncture_Aijk(M, a, r, E, L)` returning the 50
    named coefficients A006..A900 (the C static names, kerr-circular.c:18) as a NamedTuple.

Mechanical rules (the COMPLETE rewrite set -- nothing else is changed, term order and
parenthesization are preserved verbatim so the FP evaluation matches the C oracle
ulp-for-ulp up to pow(x,n) vs Julia ^):
    1. `pow(X, n)` -> `X^n` (bare-identifier base) or `(X)^n` (composite base), innermost
       first via balanced-parenthesis parsing; all exponents in the block are integers.
    2. C trailing-dot double literals `2.` `3.` `6.` `12.` `24.` `48.` -> the same integer
       (exact; keeps the Julia code type-generic so BigFloat precision is inherited).
    3. Leading tab and trailing `;` stripped.

Run from the package root:
    python3 src/Puncture/generator/transcribe_field4d.py
It verifies structure (exactly 50 assignments; the name set == {A_ijk : j, k even,
6 <= i+j+k <= 9}; no `pow` left; balanced parens) and rewrites field4d_coeffs.jl.
"""
import io
import itertools
import re
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[3]          # .../grtools/capabilities/ScalarSelfForce
REPO = PKG.parents[2]                              # .../emri-disk
CSRC = REPO / "ref-code" / "EffectiveSource" / "kerr-circular.c"
OUT = PKG / "src" / "Puncture" / "field4d_coeffs.jl"

IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def extract_block(text):
    """Return the assignment lines of the `/* Compute A coefficients */ { ... }` block."""
    tail = text[text.index("/* Compute A coefficients */"):]
    body = tail[tail.index("{") + 1 : tail.index("}")]
    lines = [ln.strip() for ln in body.splitlines() if ln.strip()]
    assert all(ln.endswith(";") for ln in lines), "every A line must end with ';'"
    return lines


def convert_pow(expr):
    """Rewrite every pow(X, n) as X^n (innermost first, balanced-paren parsing)."""
    while True:
        i = expr.find("pow(")
        if i < 0:
            return expr
        # find the matching close paren and the top-level comma
        depth, j, comma = 0, i + 3, None
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
        assert depth == 0 and comma is not None, f"unbalanced pow at {i}"
        base = convert_pow(expr[i + 4 : comma].strip())
        exp = expr[comma + 1 : j].strip()
        assert re.fullmatch(r"\d+", exp), f"non-integer exponent {exp!r}"
        rep = f"{base}^{exp}" if IDENT.match(base) else f"({base})^{exp}"
        expr = expr[:i] + rep + expr[j + 1 :]


def convert_line(line):
    lhs, rhs = line.rstrip(";").split("=", 1)
    rhs = convert_pow(rhs.strip())
    rhs = re.sub(r"(?<![\w.])(\d+)\.(?!\d)", r"\1", rhs)   # 3.* -> 3*, /24. -> /24 (exact)
    assert "pow" not in rhs and rhs.count("(") == rhs.count(")")
    return lhs.strip(), rhs


def expected_names():
    return sorted(
        f"A{i}{j}{k}"
        for i, j, k in itertools.product(range(10), range(0, 9, 2), range(0, 9, 2))
        if 6 <= i + j + k <= 9
    )


def main():
    pairs = [convert_line(ln) for ln in extract_block(CSRC.read_text())]
    names = [n for n, _ in pairs]
    assert names == sorted(names) == expected_names(), "A_ijk name-set mismatch"
    assert len(pairs) == 50

    out = io.StringIO()
    out.write(f'''\
# field4d_coeffs.jl — GENERATED FILE, DO NOT EDIT BY HAND.
# Regenerate (byte-identical) with:  python3 src/Puncture/generator/transcribe_field4d.py
#
# Mechanical transcription of the `/* Compute A coefficients */` block of
# effsource_set_particle, ref-code/EffectiveSource/kerr-circular.c (READ-ONLY mirror,
# pinned SHA 1cb9834f8a8f666aaf7b1f924b3197743efb4120; https://github.com/barrywardell/EffectiveSource,
# MIT license). Rewrite rules (complete list): pow(X,n) -> X^n, C double literals
# 2./3./6./12./24./48. -> exact integer literals, layout otherwise verbatim -- see the
# generator script header. These are the A_ijk of eq:puncture (arxiv-2606.04998
# main.tex:440-457): Phi^P = R^-7 sum A_ijk dr^i dtheta^j Q^k, 6 <= i+j+k <= 9, and
# A_ijk = 0 for odd j or k on circular equatorial orbits -- hence exactly these 50
# even-j, even-k names (the C statics of kerr-circular.c:18).

"""
    _puncture_Aijk(M, a, r, E, L) -> NamedTuple (A006, A008, ..., A900)

The 50 nonzero 4th-order Detweiler-Whiting puncture coefficients A_ijk for a circular
equatorial Kerr geodesic: `M`/`a` background mass and spin, `r` the orbital radius, `E`/`L`
the specific energy and axial angular momentum (kerr-circular.c effsource_set_particle).
Type-generic: all arithmetic promotes from the argument types (integer literals only), so
BigFloat inputs yield full-precision BigFloat coefficients.
"""
function _puncture_Aijk(M, a, r, E, L)
''')
    for name, rhs in pairs:
        out.write(f"    {name} = {rhs}\n")
    out.write("    return (; " + ", ".join(names) + ")\nend\n")
    OUT.write_text(out.getvalue())
    print(f"wrote {OUT} ({len(pairs)} coefficients)")


if __name__ == "__main__":
    sys.exit(main())
