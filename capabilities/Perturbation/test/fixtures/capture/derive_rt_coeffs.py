# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/derive_rt_coeffs.py (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = 985de8be15981b87a17e4935b3a1eb730dea8d7e2d9bbeee7ff719c71045b251
# fixture(s)  : rt_coeffs_reference.jl
# oracle      : INDEPENDENT sympy re-derivation (sympy 1.12) of arXiv:1301.1591 Eq. (15); no external oracle file
# run         : python3 test/fixtures/capture/derive_rt_coeffs.py
#
# The script body is VERBATIM from the original one-shot capture EXCEPT for the output
# path (single deviation, GRR::g2): the original resolved it relative to its scratchpad
# location as `<script dir>/../grtools/capabilities/Perturbation/test/fixtures/
# rt_coeffs_reference.jl`, which no longer exists from the tracked location.  The tracked
# copy resolves the SAME fixture through GRTOOLS_FIXTURE_OUT_DIR (default: the sibling
# `test/fixtures/` directory).  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run
# WITHOUT overwriting the committed fixture, then diff.
# ==========================================================================================================

#!/usr/bin/env python3
# ONE-SHOT capture (DG-P4 transcription gate): INDEPENDENT sympy re-derivation of the
# RT-chart m-mode Teukolsky coefficients of arXiv:1301.1591 Eq. (15) (appendix list),
# baked as numeric reference values into a Perturbation test fixture.
#
# Derivation (no appendix input — only the master equation + chart definitions):
#   master TE in BL coords [arXiv:2307.01349, Eq. (1)] (≡ arXiv:1301.1591 Eq. (10) modulo
#     its (r²+a²)Δ⁻¹ → (r²+a²)²Δ⁻¹ typo),
#   ingoing-Kerr:  τ = t − r + r*_K,  φ̃ = φ + a∫dr/Δ      [1301.1591, Eqs. (3),(4)]
#   RT:            τ = T − 4M log(1−R²) + (1+R²)/(1−R²),  r = 2R/(1−R²)   [Eqs. (5),(6)]
#   rescaling+mode: Ψ = Δ⁻ˢ r⁻¹ e^{imφ̃} ψ(T(t,r), R(r), θ)                [Eqs. (12),(15)]
# so ∂t = ∂T, ∂φ = im (φ̃-mode), ∂r|_t acts via k'(r)∂T + R'(r)∂R + prefactor-log-derivative
# L = −sΔ'/Δ − 1/r + im·a/Δ, with T = t + k(r),
#   k'(r) = (r²+a²)/Δ − 1 + [−8MR/(1−R²) − 4R/(1−R²)²]·R'(r),  R'(r) = (1−R²)²/(2(1+R²)).
# Collect the jet coefficients of {ψ, ψ_T, ψ_R, ψ_θ, ψ_TT, ψ_TR, ψ_RR, ψ_θθ}, normalize to
# the paper's convention by scaling coeff(ψ_θθ) → −1/R, and evaluate numerically.
#
# Usage: python3 scratchpad/derive_rt_coeffs.py
import sympy as sp

M, a, R, r, th = sp.symbols('M a R r theta', positive=True)
s, m = sp.symbols('s m', integer=True)

Delta = r**2 - 2*M*r + a**2

# chart maps
r_of_R = 2*R/(1 - R**2)
Rp_of_r = (sp.sqrt(1 + r**2) - 1)/r          # inverse of r = 2R/(1-R^2), R in (0,1)
dR_dr = sp.simplify(sp.diff(Rp_of_r, r)).subs(r, r_of_R)  # R'(r) as function of R
dR_dr = sp.simplify(dR_dr)

# k(r): T = t + k,  k = r*_K - r + 4M log(1-R^2) - (1+R^2)/(1-R^2);  only k', k'' needed
kp_R = (r_of_R**2 + a**2)/Delta.subs(r, r_of_R) - 1 \
    + (-8*M*R/(1 - R**2) - 4*R/(1 - R**2)**2)*dR_dr
kp_r = kp_R.subs(R, Rp_of_r)                  # k'(r) as function of r

# prefactor log-derivative: P = Delta^-s r^-1 e^{im*Phi(r)}, Phi' = a/Delta
L = -s*sp.diff(Delta, r)/Delta - 1/r + sp.I*m*a/Delta

# jet symbols
p0, pT, pR, pth, pTT, pTR, pRR, pthth = sp.symbols(
    'psi0 psiT psiR psith psiTT psiTR psiRR psithth')

Rr = Rp_of_r  # R(r)


def Dt(expr):
    return expr.subs({pR: pTR, pT: pTT, p0: pT}, simultaneous=True)


def Dr(expr):
    """d/dr of (prefactor * expr(jet)) divided by prefactor: chain rule + L*expr."""
    # coefficient derivative (explicit r dependence)
    out = sp.diff(expr, r) + L*expr
    # jet shifts: d/dr psi = k' psi_T + R' psi_R  etc.
    kpr = kp_r
    Rpr = sp.diff(Rp_of_r, r)
    out += sp.diff(expr, p0)*0  # p0 handled below
    jetmap = {
        p0: kpr*pT + Rpr*pR,
        pT: kpr*pTT + Rpr*pTR,
        pR: kpr*pTR + Rpr*pRR,
    }
    for k_, v_ in jetmap.items():
        out += sp.diff(expr, k_)*(v_ - 0)  # replaces the "jet derivative" of k_
        out -= sp.diff(expr, k_)*sp.diff(k_, r)  # (zero; jets are symbols)
    return out


# ---- master equation [arXiv:2307.01349, Eq. (1)], vacuum, applied to P*psi ----------------------
# -[(r²+a²)²/Δ - a²sin²θ]∂tt Ψ - (4Mar/Δ)∂tφ Ψ - 2s[r - M(r²-a²)/Δ + ia cosθ]∂t Ψ
# + Δ^{-s}∂r(Δ^{s+1}∂r Ψ) + (1/sinθ)∂θ(sinθ ∂θ Ψ) + [1/sin²θ - a²/Δ]∂φφ Ψ
# + 2s[a(r-M)/Δ + i cosθ/sin²θ]∂φ Ψ - (s²cot²θ - s)Ψ = 0
# with Ψ = P·ψ(T,R,θ), ∂t→ψ_T-shift, ∂φ→im, ∂r via Dr (reduced by P).

psi = p0
term_tt = -((r**2 + a**2)**2/Delta - a**2*sp.sin(th)**2) * Dt(Dt(psi))
term_tphi = -(4*M*a*r/Delta) * (sp.I*m) * Dt(psi)
term_t = -2*s*(r - M*(r**2 - a**2)/Delta + sp.I*a*sp.cos(th)) * Dt(psi)
# Δ^{-s}∂r(Δ^{s+1}∂r Ψ)/P = Δ^{-s}[ ∂r(Δ^{s+1})·(DrΨ/P) + Δ^{s+1}·(Dr² ...) ]
inner = Dr(psi)
term_rr = (s + 1)*sp.diff(Delta, r)*inner + Delta*Dr(inner)
term_thth = pthth + sp.cos(th)/sp.sin(th)*pth
term_phiphi = (1/sp.sin(th)**2 - a**2/Delta) * (sp.I*m)**2 * psi
term_phi = 2*s*(a*(r - M)/Delta + sp.I*sp.cos(th)/sp.sin(th)**2) * (sp.I*m) * psi
term_0 = -(s**2*(sp.cos(th)/sp.sin(th))**2 - s) * psi

eq = term_tt + term_tphi + term_t + term_rr + term_thth + term_phiphi + term_phi + term_0

jet = [p0, pT, pR, pth, pTT, pTR, pRR, pthth]
coeffs = {j: sp.expand(eq).coeff(j) for j in jet}

# normalize to the paper convention: coeff(psithth) -> -1/R  (C_θθ = -1/R)
norm = coeffs[pthth]  # should be a function of r only


def bake(a_val, s_val, m_val, R_val, th_val, prec=80):
    subs0 = {M: 1, a: sp.nsimplify(a_val), s: s_val, m: m_val}
    r_val = (2*R_val/(1 - R_val**2)).evalf(prec) if R_val != 1 else None
    out = {}
    nrm = None
    for j in jet:
        c = coeffs[j].subs(subs0)
        c = c.subs({th: th_val})
        # substitute r as a high-precision number
        cval = c.subs(r, r_val).evalf(prec)
        out[j] = cval
    nrm = out[pthth]
    scale = (-1/sp.Float(R_val, prec))/nrm
    return {j: sp.simplify(out[j]*scale).evalf(20) for j in jet}


def main():
    samples = []
    # (a, s, m): coverage of a=0/0.9 × s∈{-2,-1,0,1,2} × m∈{0,±1,2,3}
    cases = [(0, -2, 0), (0, 2, 1), (sp.Rational(9, 10), -2, 0), (sp.Rational(9, 10), -1, 1),
             (sp.Rational(9, 10), 1, 0), (sp.Rational(9, 10), 2, 2), (sp.Rational(9, 10), 0, 3),
             (sp.Rational(1, 2), -1, -2)]
    for (av, sv, mv) in cases:
        # R samples: strictly inside + very near both ends (R=1 exactly is singular in the
        # intermediate r-substitution; probe 1-1e-6 for the scri limit; R_+ + 1e-6 likewise)
        rp = sp.nsimplify(1 + sp.sqrt(1 - av**2), rational=False)
        Rp_val = (sp.sqrt(1 + rp**2) - 1)/rp
        Rlist = [(Rp_val + sp.Rational(1, 10**4)).evalf(80),
                 ((Rp_val + 1)/2).evalf(80), sp.Rational(4, 5).evalf(80),
                 (1 - sp.Rational(1, 10**4)).evalf(80)]
        for Rv in Rlist:
            for thv in [sp.pi/2, sp.pi/3]:
                vals = bake(av, sv, mv, Rv, thv)
                samples.append(((av, sv, mv, Rv, thv), vals))
    _os = __import__('os')
    _fixdir = _os.environ.get('GRTOOLS_FIXTURE_OUT_DIR') or _os.path.join(
        _os.path.dirname(_os.path.abspath(__file__)), '..')
    out = _os.path.join(_fixdir, 'rt_coeffs_reference.jl')
    with open(out, 'w') as f:
        f.write("# AUTO-GENERATED one-shot fixture — DO NOT EDIT.\n")
        f.write("# INDEPENDENT sympy re-derivation (scratchpad/derive_rt_coeffs.py, sympy %s)\n"
                % sp.__version__)
        f.write("# of the RT-chart m-mode Teukolsky coefficients [arXiv:1301.1591, Eq. (15)\n")
        f.write("# appendix list] from the master equation [arXiv:2307.01349, Eq. (1)] +\n")
        f.write("# chart maps [arXiv:1301.1591, Eqs. (3)-(6),(12)] ONLY (no appendix input).\n")
        f.write("# Normalization: coeff(psi_thth) = -1/R (paper convention). M = 1.\n")
        f.write("# Tuple layout: (a, s, m, R, theta, C0, CT, CR, Cth, CTT, CTR, CRR)\n")
        f.write("# (each C complex, order matches the jet {psi, psi_T, psi_R, psi_th, psi_TT,\n")
        f.write("#  psi_TR, psi_RR}; C_thth == -1/R by construction, not stored).\n")
        f.write("const RT_COEFFS_REF = [\n")
        for (key, vals) in samples:
            av, sv, mv, Rv, thv = key
            def jl(z):
                zz = sp.N(z, 17)
                re, im = float(sp.re(zz)), float(sp.im(zz))
                return "complex(%r, %r)" % (re, im)
            f.write("    (%r, %d, %d, %r, %r, %s, %s, %s, %s, %s, %s, %s),\n" % (
                float(av), sv, mv, float(Rv), float(thv),
                jl(vals[p0]), jl(vals[pT]), jl(vals[pR]), jl(vals[pth]),
                jl(vals[pTT]), jl(vals[pTR]), jl(vals[pRR])))
        f.write("]\n")
    print("wrote", out, "samples:", len(samples))


if __name__ == '__main__':
    main()
