"""
Option A pressure law reduced by the SAME reduce.jl as the adiabatic default.
Only the residual selector differs. The fluid equations and alphabet are shared
with R3Derivation; no adiabatic diagnostic tensors are used as LI diagnostics.
"""
module LocallyIsothermal

using Symbolics
using ..R3Derivation: dp, dp_r, dur, dur_r, dut, duphi, deps,
    p, p_r, eps, eps_r, Gamma1, XI_R,
    y0, y1, y0_r, y1_r, nu_r, psi_r, epsi, v, omt_r, enu, W,
    DEULER_UP, DENERGY, DNORM, DUPHI_LOW, IT, IR, IPHI, IZ,
    iszero_entry, to_omega_d, no_second_derivatives, ELIMINATED_DERIVS,
    htt, htr, htphi, hrr, hrphi, hphiphi, hzz,
    htt_r, htr_r, htphi_r, hrr_r, hrphi_r, hphiphi_r, hzz_r,
    Certificate, certify_symbolic, certify_exact, certify_bool, certify_pole_order,
    lie_u_xi_r, onshell
import ..R3Derivation as R3

const PRESSURE_CLOSURE = Val(:locally_isothermal)
include("reduce.jl")
# Reuse backsub, MASTER_NAMES and structural inspection from the existing stage.
# Its adiabatic-labelled suite is not the LI gate; closure_certificates is.
include("stage2.jl")
include("li_certificates.jl")

end
