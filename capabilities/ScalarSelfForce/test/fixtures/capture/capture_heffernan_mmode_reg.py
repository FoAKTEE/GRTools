#!/usr/bin/env python3
"""Capture script for test/fixtures/heffernan_mmode_reg.jl (ESSF node
`heffernan-regularization`, GLOBAL_DAG row 45).

ORACLE (independent of the tex transcription shipped in
src/Extraction/regularization.jl): the `Fmr[4]` definition of the ancillary
Mathematica notebook

    <mirror>/arxiv-1610.09319/src/anc/m-ModeKerrRegularizationParameters.nb

"Kerr Scalar Field m-mode Regularization parameters", B. Wardell and
J. Thornburg (supplement to arXiv:1610.09319) -- the "electronic form" that
Heffernan, Ottewill & Wardell 2014 (arXiv:1211.6446) point to at
KerrRegularization.tex:1685-1687 for the expressions they do not print.  It
gives F^m_{r[4]} for a GENERIC equatorial Kerr geodesic in terms of
(E, L, r, a, M, m) INCLUDING the 24/(pi (2m-3)(2m-1)(2m+1)(2m+3)) m-factor;
this script specialises it to circular orbits by substituting the circular
(E, L) and evaluates it in 60-digit arithmetic.

The notebook is a READ-ONLY mirror living outside this repository, so the
fixture is baked; re-run with

    python3 capture_heffernan_mmode_reg.py [--mirror <ref-paper dir>] \
        > ../heffernan_mmode_reg.jl

Requires: mpmath.  Deliberately NOT a live test dependency (dep-budget gate
test/hygiene.jl: no Python/GSL/ccall channel in the Julia test env).
"""
import argparse, hashlib, json, os, re, sys

DEFAULT_MIRROR = "/data/haiyangw/claude/emri-disk/ref-paper"
NB_REL = "arxiv-1610.09319/src/anc/m-ModeKerrRegularizationParameters.nb"
CELL_FIRST, CELL_LAST = 2806, 9600          # 1-based: the `Fmr[4]` Input cell

# --------------------------------------------------------------- notebook box parser
TOKEN_RE = re.compile(r'''
      (?P<str>"(?:[^"\\]|\\.)*")
    | (?P<sym>\\\[[A-Za-z]+\]|[A-Za-z$][A-Za-z0-9$]*)
    | (?P<num>[0-9]+\.?[0-9]*(?:`[0-9.]*)?(?:\*\^-?[0-9]+)?)
    | (?P<op>->|\[|\]|\{|\}|,|;)
    | (?P<ws>\s+)
''', re.X)

def tokenize(s):
    out, i, n = [], 0, len(s)
    while i < n:
        m = TOKEN_RE.match(s, i)
        if not m:
            raise ValueError("tokenizer stuck at %d: %r" % (i, s[i:i + 60]))
        i = m.end()
        if m.lastgroup != "ws":
            out.append((m.lastgroup, m.group()))
    return out

class Parser:
    def __init__(self, toks): self.t, self.i = toks, 0
    def peek(self): return self.t[self.i] if self.i < len(self.t) else (None, None)
    def next(self): tok = self.t[self.i]; self.i += 1; return tok
    def parse(self):
        k, s = self.next()
        if k == "str":
            node = ("str", json.loads(s.replace("\\[", "\\\\[")))
        elif k == "num":
            node = ("num", s)
        elif k == "sym":
            if self.peek()[1] == "[":
                self.next()
                node = ("call", s, self.seq("]"))
            else:
                node = ("sym", s)
        elif s == "{":
            node = ("list", self.seq("}"))
        else:
            raise ValueError("unexpected token %r" % (s,))
        while self.peek()[1] == "->":
            self.next(); self.parse()          # swallow option values
        return node
    def seq(self, close):
        args = []
        if self.peek()[1] == close:
            self.next(); return args
        while True:
            args.append(self.parse())
            _, s = self.next()
            if s == close: return args
            if s != ",": raise ValueError("bad separator %r" % (s,))

def render(node):
    kind = node[0]
    if kind in ("str", "num", "sym"):
        return node[1]
    if kind == "list":
        return "".join(render(a) for a in node[1])
    head, args = node[1], node[2]
    if head == "RowBox":       return "".join(render(a) for a in args)
    if head == "FractionBox":  return "((%s)/(%s))" % (render(args[0]), render(args[1]))
    if head == "SuperscriptBox": return "((%s)**(%s))" % (render(args[0]), render(args[1]))
    if head == "SqrtBox":      return "(sqrt(%s))" % (render(args[0]),)
    if head in ("StyleBox", "TagBox", "FormBox", "InterpretationBox", "BoxData"):
        return render(args[0])
    raise ValueError("unhandled box head %r" % (head,))

_L = set("0123456789)]_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
_R = set("0123456789([_\\abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")

def fix_spaces(s):
    """TraditionalForm juxtaposition -> explicit `*`."""
    out, i, n = [], 0, len(s)
    while i < n:
        if s[i] == " ":
            j = i
            while j < n and s[j] == " ":
                j += 1
            if (out[-1] if out else "") in _L and (s[j] if j < n else "") in _R:
                out.append("*")
            i = j
            continue
        out.append(s[i]); i += 1
    return "".join(out)

def extract_fmr4(nb_path):
    with open(nb_path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    blob = "".join(lines[CELL_FIRST - 1:CELL_LAST]).strip()
    assert blob.startswith("Cell[BoxData["), blob[:40]
    assert blob.endswith('";"}]], "Input",'), blob[-40:]
    node = Parser(tokenize(blob[len("Cell["): -len(', "Input",')])).parse()
    assert node[0] == "call" and node[1] == "BoxData"
    txt = fix_spaces(render(node[2][0]))
    for k, v in ((r"\[ScriptCapitalE]", "En"), (r"\[ScriptCapitalL]", "Lz"), (r"\[Pi]", "PI")):
        txt = txt.replace(k, v)
    assert txt.startswith("Fmr[4]="), txt[:40]
    expr = txt[len("Fmr[4]="):].rstrip().rstrip(";")
    left = re.findall(r"\\\[[A-Za-z]+\]", expr)
    assert not left, left
    idents = sorted(set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expr)))
    assert idents == ["En", "Lz", "M", "PI", "a", "m", "r"], idents
    return expr.replace(")(", ")*(").replace("^", "**")

# ------------------------------------------------------------------------ grid + main
A_GRID  = [(0, 1), (1, 2), (-1, 2), (3, 5), (9, 10), (-9, 10), (499, 500), (-499, 500)]
R0_GRID = [6, 7, 10, 20, 50]
M_GRID  = [0, 1, 2, 5, 20]
DIGITS  = 50

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mirror", default=DEFAULT_MIRROR)
    args = ap.parse_args()
    nb = os.path.join(args.mirror, NB_REL)
    if not os.path.isfile(nb):
        sys.exit("read-only mirror not present: %s" % nb)
    sha = hashlib.sha256(open(nb, "rb").read()).hexdigest()
    code = compile(extract_fmr4(nb), "<Fmr4>", "eval")

    from mpmath import mp, mpf, sqrt, pi, nstr
    mp.dps = 60

    def circ_EL(a, r, M):
        v = sqrt(M * r)
        den = r * sqrt(r**2 - 3*M*r + 2*a*v)
        return (r**2 - 2*M*r + a*v) / den, v * (r**2 - 2*a*v + a**2) / den

    rows = []
    for (an, ad) in A_GRID:
        a = mpf(an) / ad
        for r0 in R0_GRID:
            r = mpf(r0)
            En, Lz = circ_EL(a, r, mpf(1))
            for m in M_GRID:
                val = eval(code, {"__builtins__": {}},
                           {"a": a, "r": r, "m": mpf(m), "M": mpf(1),
                            "PI": pi, "En": En, "Lz": Lz})
                rows.append(((an, ad), r0, m, nstr(val, DIGITS, strip_zeros=False)))

    w = sys.stdout.write
    w("""# heffernan_mmode_reg.jl — BAKED m-mode radial-self-force regularization-parameter
# ORACLE for circular equatorial Kerr geodesics (ESSF node `heffernan-regularization`,
# GLOBAL_DAG row 45).  GENERATED FILE — do not edit by hand; regenerate with
# test/fixtures/capture/capture_heffernan_mmode_reg.py (see its docstring).
#
# CITATION / ORACLE (INDEPENDENT of the tex transcription in
#   src/Extraction/regularization.jl):
#   B. Wardell and J. Thornburg, "Kerr Scalar Field m-mode Regularization parameters",
#   ancillary Mathematica notebook `m-ModeKerrRegularizationParameters.nb` of
#   J. Thornburg and B. Wardell, "Scalar self-force for highly eccentric equatorial
#   orbits in Kerr spacetime", Phys. Rev. D 95, 084043 (2017); arXiv:1610.09319.
#   Read-only mirror: ref-paper/arxiv-1610.09319/src/anc/m-ModeKerrRegularizationParameters.nb
#   sha256 = HEFF_MMODE_NB_SHA256 below.  Cell `Fmr[4]` (notebook lines %d-%d) gives
#   F^m_{r[4]} for a GENERIC equatorial Kerr geodesic in terms of (E, L, r, a, M, m)
#   INCLUDING the 24/(pi (2m-3)(2m-1)(2m+1)(2m+3)) m-factor; it is the "electronic form"
#   cited by Heffernan, Ottewill & Wardell 2014 (arXiv:1211.6446) at
#   KerrRegularization.tex:1685-1687 in place of the generic printed expressions.
#
# SPECIALISATION: circular equatorial orbits, M = 1, u^r = 0 (u^r does not appear in
#   `Fmr[4]` at all), with the Bardeen-Press-Teukolsky circular constants of motion
#     E = (r^2 - 2Mr + a sqrt(Mr)) / (r sqrt(r^2 - 3Mr + 2a sqrt(Mr))),
#     L = sqrt(Mr) (r^2 - 2a sqrt(Mr) + a^2) / (r sqrt(r^2 - 3Mr + 2a sqrt(Mr))),
#   signed a (a < 0 = retrograde, motion always +phi) — the package-wide convention.
#
# VALUES: %d significant decimal digits from 60-digit mpmath arithmetic, stored as EXACT
#   DECIMAL STRINGS so the consumer chooses the precision (parse(Float64, .) or
#   parse(BigFloat, .) at any setprecision).  Rows are (a::Rational{Int}, r0::Int,
#   m::Int, dec::String) with a and r0 in units of M = 1.
#
# CONSUMER: test/test_regularization.jl — cross-SOURCE gate on
#   `mmode_reg_Fr(4, CircularOrbit(a, r0), m)`.

const HEFF_MMODE_NB_SHA256 = "%s"
const HEFF_MMODE_NB_PATH = "ref-paper/%s"
const HEFF_MMODE_NB_CELL_LINES = %d:%d
const HEFF_MMODE_ORACLE_DIGITS = %d

const HEFF_MMODE_ORACLE_A  = (%s)
const HEFF_MMODE_ORACLE_R0 = (%s)
const HEFF_MMODE_ORACLE_MS = (%s)

\"\"\"
`HEFF_MMODE_FMR4[(a, r0, m)]` — the baked decimal string for `F^m_{r[4]}` (Wardell-Thornburg
notebook oracle, see this file's header).  `a::Rational{Int}` and `r0::Int` in units of
`M = 1`, `m::Int >= 0`; the value is the FULL m-mode parameter, m-factor included.
\"\"\"
const HEFF_MMODE_FMR4 = Dict{Tuple{Rational{Int},Int,Int},String}(
""" % (CELL_FIRST, CELL_LAST, DIGITS, sha, NB_REL, CELL_FIRST, CELL_LAST, DIGITS,
       ", ".join("%d//%d" % ad for ad in A_GRID),
       ", ".join(str(r) for r in R0_GRID),
       ", ".join(str(m) for m in M_GRID)))
    for ((an, ad), r0, m, dec) in rows:
        w('    (%d//%d, %d, %d) => "%s",\n' % (an, ad, r0, m, dec))
    w(")\n")

main()
