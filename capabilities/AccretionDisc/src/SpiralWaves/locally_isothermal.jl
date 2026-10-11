# OPTION A: prescribed Eulerian p=Theta(r)*epsilon and conserved stress energy.
# These are runtime adapters to the certified emitted arithmetic, not a second
# closure derivation. There is no adiabatic U/V/W wave-flux substitution here.
include("generated/eccentric_li.jl")
include("generated/eccentric_li_stress.jl")

function _context_operator(::Val{LocallyIsothermalOptionA}, inp)
    [_eccentric_li_C00(inp) _eccentric_li_C01(inp);
     _eccentric_li_C10(inp) _eccentric_li_C11(inp)]
end
_context_source(::Val{LocallyIsothermalOptionA}, inp) = (_eccentric_li_S0(inp),_eccentric_li_S1(inp))
function _context_reconstruction(::Val{LocallyIsothermalOptionA}, inp)
    (; dut=_eccentric_li_dut(inp),dur=im*inp.y1,
       duphi_up=_eccentric_li_duphi_up(inp),deps=_eccentric_li_deps(inp),
       du_phi=_eccentric_li_duphi_low(inp),dp=(inp.eps+inp.p)*inp.y0)
end

"""
    linear_stress(context, r, y; hmodes, landau=true)

OPTION A Eulerian contravariant stress perturbation dT^{mu nu} (coordinate order
t,r,phi,z) and mixed coordinate currents dJt^mu=dT^mu_t,
dJphi^mu=dT^mu_phi (t,r,phi components). Includes the perturbed metric in lowering
the current index. Complex amplitudes are retained, with no reality projection.

These are linear stress jets, not the connected quadratic wave flux or its
mean-transport subtraction; those belong to the LI diagnostics certificate.
"""
function linear_stress(c::ModeContext, r, y; hmodes, landau::Bool=true)
    c.closure == LocallyIsothermalOptionA || throw(ArgumentError("linear_stress requires OPTION A context"))
    x = master_inputs(c,r;hmodes,landau,y)
    tt,tr,tp = _eccentric_li_stress_dT_tt(x),_eccentric_li_stress_dT_tr(x),_eccentric_li_stress_dT_tphi(x)
    rr,rp = _eccentric_li_stress_dT_rr(x),_eccentric_li_stress_dT_rphi(x)
    pp,zz = _eccentric_li_stress_dT_phiphi(x),_eccentric_li_stress_dT_zz(x)
    z = zero(tt)
    dT = [tt tr tp z; tr rr rp z; tp rp pp z; z z z zz]
    dJt = (_eccentric_li_stress_dJt_t(x),_eccentric_li_stress_dJt_r(x),_eccentric_li_stress_dJt_phi(x))
    dJphi = (_eccentric_li_stress_dJphi_t(x),_eccentric_li_stress_dJphi_r(x),_eccentric_li_stress_dJphi_phi(x))
    (;dT,dJt,dJphi)
end
