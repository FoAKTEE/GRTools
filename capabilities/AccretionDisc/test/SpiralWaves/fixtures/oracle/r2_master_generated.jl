# ============================================================================================
# ORACLE FIXTURE -- NOT PART OF THE BUILD.  Route R2 (OUR independent Mathematica derivation).
#
# The BUILD source of the master entries is `src/SpiralWaves/master_derived.jl`, emitted by route
# R3 (the Julia-native `Symbolics.jl` derivation of `capabilities/AccretionDisc/derivation/`).
# This file is loaded ONLY by `test/SpiralWaves/test_derived_equivalence.jl` (through
# `fixtures/oracle/oracle_module.jl`), where it gates route R3 FOREVER at arbitrary ON-SHELL
# points.  Nothing under `src/` may include it.
#
# It replaces the notebook-derived oracle (route R1, `master_generated.jl`) that was REMOVED from
# this repository on 2026-08-22 by owner decision (EGD::jd-notebook-purge).  R1 == R2 on-shell was
# certified by EGD::R1-vs-R2-gate, so the equivalence gate keeps its strength and only changes
# lineage; the R1 lineage and the full R1-vs-R2-vs-R3 evidence are archived in the mission
# repository (emri-disk results/emri-gasdisk paper_EMRI-GASDISK, evidence/{R1-vs-R2-gate,
# jd-equivalence-gate,jd-notebook-purge}).
#
# The routes agree ON-SHELL only -- where the background satisfies W = 1/sqrt(1-v^2) (8) AND the
# hydrostatic balance (11).  Off-shell each route carries its own continuation.
# ============================================================================================
# r2_master_generated.jl -- GENERATED FILE, DO NOT EDIT BY HAND (AccretionDisc.SpiralWaves test oracle, arXiv:2608.17003 eq 47).
# generator : emri-disk results/emri-gasdisk/paper_EMRI-GASDISK/codes/cas/R2_master_to_julia.wls
#             (shared tree walker codes/cas/cas_to_julia.wl; rerun: wolframscript -file R2_master_to_julia.wls)
# input     : codes/cas/deriv-master-system-R2/R2_master_system.m  sha256 67825ad049e45ed995f1a087cde6ff9f0d64101ea1ac52f457286b645de26bf6
#             (route R2 = the INDEPENDENT CAS derivation of EGD::deriv-master-system-R2)
# tool      : 14.3.0 for Linux x86 (64-bit) (July 8, 2025)
# content   : dy/dr = C.y + S with y = (deltaU, urbyi) = (delta p-hat/(eps+p), -i delta u-hat^r);
#             C = [C00 C01; C10 C11], S = [S0, S1].  ON-SHELL forms: nu' has been eliminated with the
#             hydrostatic balance (11), which is also what routes R1 and R3 emit.
# signature : (inp::MasterInputs, omt, vp) -- omt = omega-tilde(r) and vp = dv/dr are the two extra
#             background scalars of the R2 alphabet (MasterInputs carries Omega and Omega' instead,
#             related by Omega = omt + v e^{nu-psi} (10)); `r2_context` of master.jl assembles them.
#             This is the SAME convention as the _r2_* diagnostics kernels of r2_generated.jl.
# symbols   : inp::MasterInputs of master.jl (m, nu, psi, mu, psi_r, mu_r, omt_r, v, eps, eps_r, p, p_r,
#             Gamma1, omega_d, htt, htt_r, htphi, htphi_r, hphiphi, hphiphi_r, hrr, hzz -- COVARIANT
#             Fourier modes, ERRATUM E-M2).  W and Omega, Omega' are NOT read (W is on-shell-eliminated).
#             Rationals are emitted as n//d and exp(a nu + b psi + c mu) as enu^a*epsi^b*emu^c (hoisted);
#             the only other non-arithmetic call is sqrt (the W = 1/sqrt(1-v^2) factors of the R2
#             alphabet), and no non-integer `^` survives -- so every float type T (Float64, BigFloat,
#             ...) inherits its own precision.

"""
    _r2_master_entries(inp::MasterInputs, omt, vp) -> (C00, C01, C10, C11, S0, S1)

Entries of the relativistic master system dy/dr = C y + S (arXiv:2608.17003 eq 47) as DERIVED by
route R2 (our independent Mathematica derivation).  ORACLE ONLY -- the build source is the route-R3
`_master_entries` of `src/SpiralWaves/master_derived.jl`.  Valid ON-SHELL (see the file header).
Pure, allocation-free, type-generic.
"""
@inline function _r2_master_entries(inp::MasterInputs, omt, vp)
    m = inp.m
    nu = inp.nu
    psi = inp.psi
    mu = inp.mu
    psip = inp.psi_r
    mup = inp.mu_r
    omtp = inp.omt_r
    v = inp.v
    eps = inp.eps
    epsp = inp.eps_r
    p = inp.p
    pp = inp.p_r
    Gam1 = inp.Gamma1
    omd = inp.omega_d
    htt = inp.htt
    http = inp.htt_r
    htp = inp.htphi
    htpp = inp.htphi_r
    hpp = inp.hphiphi
    hppp = inp.hphiphi_r
    hrr = inp.hrr
    hzz = inp.hzz
    enu = exp(nu)
    epsi = exp(psi)
    emu = exp(mu)
    C00 = ((Gam1)^(-1) * (omd)^(-1) * (p)^(-1) * ((eps + p))^(-1) * (((-1) * epsp * Gam1 * omd * p) + ((enu^(-1) * epsi^1) * Gam1 * omd * omtp * p * (eps + p) * v) + ((-2) * (enu^1 * epsi^(-1)) * Gam1 * m * p * (pp + ((eps + p) * psip)) * v * ((-1) + (v)^2)) + (p * ((omd * pp) + ((-2) * Gam1 * omd * (pp + (p * psip)) * (v)^2) + (Gam1 * m * omtp * p * ((-1) + (v)^2)))) + (eps * ((omd * pp) + (Gam1 * p * (((-2) * omd * psip * (v)^2) + (m * omtp * ((-1) + (v)^2))))))))
    C01 = ((enu^(-1)) * (Gam1)^(-1) * (omd)^(-1) * (p)^(-1) * ((eps + p))^(-2) * (sqrt((1 + ((-1) * (v)^2))))^(-1) * (((-1) * (emu^2) * Gam1 * (omd)^2 * p * ((eps + p))^2) + ((enu^1 * epsi^1) * Gam1 * omtp * p * (eps + p) * ((pp * v * ((-1) + (v)^2)) + ((eps + p) * ((psip * v * ((-1) + (v)^2)) + ((-1) * vp))))) + ((enu^2) * (((-1) * epsp * Gam1 * p * pp * ((-1) + (v)^2)) + (2 * (eps)^2 * Gam1 * p * psip * v * ((psip * (v + ((-1) * (v)^3))) + vp)) + (p * (((-1) * (pp)^2 * ((-1) + (v)^2) * ((-1) + (2 * Gam1 * (v)^2))) + (2 * Gam1 * p * pp * v * ((2 * psip * v) + ((-2) * psip * (v)^3) + vp)) + (2 * Gam1 * (p)^2 * psip * v * ((psip * v) + ((-1) * psip * (v)^3) + vp)))) + (eps * (((pp)^2 * ((-1) + (v)^2)) + (2 * Gam1 * p * pp * v * ((2 * psip * v) + ((-2) * psip * (v)^3) + vp)) + (4 * Gam1 * (p)^2 * psip * v * ((psip * v) + ((-1) * psip * (v)^3) + vp))))))))
    C10 = ((enu^(-1) * epsi^(-2)) * (Gam1)^(-1) * (omd)^(-1) * (p)^(-1) * (sqrt((1 + ((-1) * (v)^2))))^(-1) * (((epsi^2) * eps * (omd)^2) + ((-1) * p * (((-1) * (epsi^2) * (omd)^2) + (Gam1 * ((((epsi^1) * omd * v) + ((enu^1) * m * ((-1) + (v)^2))))^2)))))
    C11 = (((enu^(-1) * epsi^1) * omtp * v) + ((enu^1 * epsi^(-1)) * m * (omd)^(-1) * ((eps + p))^(-1) * ((pp * v * ((-1) + (v)^2)) + ((eps + p) * ((psip * v * ((-1) + (v)^2)) + ((-1) * vp))))) + ((-1) * (Gam1)^(-1) * (p)^(-1) * ((eps + p))^(-1) * (((-1) + (v)^2))^(-1) * ((eps * pp * ((-1) + (v)^2)) + (p * pp * ((-1) + (v)^2)) + (eps * Gam1 * p * ((mup * ((-1) + (v)^2)) + (psip * ((-1) + (v)^2)) + (v * vp))) + (Gam1 * p * (pp + ((-1) * pp * (v)^2) + (mup * p * ((-1) + (v)^2)) + (p * psip * ((-1) + (v)^2)) + (p * v * vp))))))
    S0 = ((1//2) * (enu^(-3) * epsi^(-3)) * (omd)^(-1) * (((-1) + (v)^2))^(-1) * (((epsi^4) * omd * (htt + (omt * ((2 * htp) + (hpp * omt)))) * omtp * v) + ((-2) * (enu^3) * hpp * m * ((eps + p))^(-1) * (v)^3 * (((-1) * (epsi^1) * omtp * (eps + p) * v) + ((enu^1) * (pp + ((eps + p) * psip)) * ((-1) + (v)^2)))) + ((-1) * (enu^2 * epsi^2) * v * ((2 * htpp * omd) + (2 * hppp * omd * omt) + (2 * htp * m * omtp) + (2 * hpp * omd * omtp) + (2 * hpp * m * omt * omtp) + ((-4) * htp * omd * psip) + ((-4) * htp * m * omt * psip) + ((-4) * hpp * omd * omt * psip) + ((-2) * hpp * m * (omt)^2 * psip) + (2 * htp * m * omtp * (v)^2) + (hpp * omd * omtp * (v)^2) + (2 * hpp * m * omt * omtp * (v)^2) + (4 * htp * m * omt * ((v * (((-1) * (enu^(-1) * epsi^1) * omtp) + (psip * v))) + (((eps + p))^(-1) * pp * ((-1) + (v)^2)))) + (2 * hpp * m * (omt)^2 * ((v * (((-1) * (enu^(-1) * epsi^1) * omtp) + (psip * v))) + (((eps + p))^(-1) * pp * ((-1) + (v)^2)))) + (2 * htt * m * (((-1) * (enu^(-1) * epsi^1) * omtp * v) + (((eps + p))^(-1) * pp * ((-1) + (v)^2)) + (psip * ((-1) + (v)^2)))))) + ((-1) * (enu^3 * epsi^1) * (v)^2 * ((hppp * omd) + (4 * htp * m * (((-1) * (enu^(-1) * epsi^1) * omtp * v) + (((eps + p))^(-1) * pp * ((-1) + (v)^2)) + (psip * ((-1) + (v)^2)))) + (hpp * (((-4) * omd * psip) + (m * (omtp + ((-4) * omt * psip) + (omtp * (v)^2))) + (2 * (omd + (2 * m * omt)) * ((v * (((-1) * (enu^(-1) * epsi^1) * omtp) + (psip * v))) + (((eps + p))^(-1) * pp * ((-1) + (v)^2)))))))) + ((-1) * (enu^1 * epsi^3) * ((http * omd) + (2 * htpp * omd * omt) + (hppp * omd * (omt)^2) + (2 * htp * omd * omtp) + (2 * htp * m * omt * omtp) + (2 * hpp * omd * omt * omtp) + (hpp * m * (omt)^2 * omtp) + (2 * htp * m * omt * omtp * (v)^2) + (hpp * m * (omt)^2 * omtp * (v)^2) + ((-4) * htp * omd * omt * ((v * (((-1) * (enu^(-1) * epsi^1) * omtp) + (psip * v))) + (((eps + p))^(-1) * pp * ((-1) + (v)^2)))) + ((-2) * hpp * omd * (omt)^2 * ((v * (((-1) * (enu^(-1) * epsi^1) * omtp) + (psip * v))) + (((eps + p))^(-1) * pp * ((-1) + (v)^2)))) + (htt * ((m * omtp * (1 + (v)^2)) + (2 * omd * ((v * (((enu^(-1) * epsi^1) * omtp) + ((-1) * psip * v))) + (((eps + p))^(-1) * (pp + ((-1) * pp * (v)^2)))))))))))
    S1 = ((1//2) * (enu^(-1) * epsi^(-4) * emu^(-2)) * (omd)^(-1) * (sqrt((1 + ((-1) * (v)^2))))^(-1) * (((epsi^4) * hrr * (omd)^2) + ((epsi^4 * emu^2) * hzz * (omd)^2) + ((-1) * (epsi^2 * emu^2) * ((htt * (m)^2) + ((omd + (m * omt)) * ((2 * htp * m) + (hpp * (omd + (m * omt)))))) * ((-1) + (v)^2)) + ((-2) * (enu^1 * epsi^1 * emu^2) * m * ((htp * m) + (hpp * (omd + (m * omt)))) * v * ((-1) + (v)^2)) + ((-1) * (enu^2 * emu^2) * hpp * (m)^2 * (v)^2 * ((-1) + (v)^2))))
    return (C00, C01, C10, C11, S0, S1)
end
