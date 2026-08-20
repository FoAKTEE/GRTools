# ONE-SHOT capture script for `mr19_newtonian.jl` — the characterization oracle of the MR19/MR20
# Newtonian planet-forced disc-mode solver (`AccretionDisc.Modes`, arXiv:1811.09628).
#
# The oracle is the LEGACY `DiscModes` package plus the figure-local homogeneous solver of the
# legacy MR19 reproduction, both in the emri-disk mission repository. Neither is a dependency of
# GrTools: this script exists so the fixture can be REGENERATED and audited, and it is run by hand
# against a checkout of the legacy tree, never by the test suite.
#
#   MR19_LEGACY_ROOT=/path/to/emri-disk \
#     julia --project=$MR19_LEGACY_ROOT/packages/DiscModes -t1 capture_mr19_newtonian.jl [outfile]
#
# Wall time ≈ 25 min at -t1 (the m = 10 forced solve alone is 733 s).
using DiscModes, Printf
const DM = DiscModes
const LEGACY = get(ENV, "MR19_LEGACY_ROOT", joinpath(@__DIR__, "..", "..", "..", "..", "..", "emri-disk"))
include(joinpath(LEGACY, "results", "emri-disk", "figures", "1811.09628", "homog_solver.jl"))
const OUT = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "mr19_newtonian.jl")

d = NewtonianDisc(h0=0.1, q=1.0, p=1.0, gamma=7/5, beta=80.0, soft=0.6*0.1, indirect=0, pole_disp=1e-6)
dc = NewtonianDisc(h0=0.1, q=1.0, p=1.0, gamma=7/5, beta=0.1,  soft=0.6*0.1, indirect=0, pole_disp=1e-6)
dh = NewtonianDisc(h0=0.1, q=1.0, p=1.0, gamma=7/5, beta=0.0,  soft=0.6*0.1, indirect=0, pole_disp=1e-6)

io = open(OUT, "w")
println(io, "# ONE-SHOT capture of the legacy `DiscModes` Newtonian planet-forced disc-mode solver")
println(io, "# (emri-disk packages/DiscModes: newtonian_disc.jl, master_equation_cool_newt.jl, solver_newt.jl,")
println(io, "#  integrator.jl) + the figure-local homogeneous solver")
println(io, "# (results/emri-disk/figures/1811.09628/homog_solver.jl), captured for the GrTools port.")
println(io, "# Disc: MR19 fiducial h0=0.1, q=1, p=1, gamma=7/5, soft=0.06, pole_disp=1e-6; beta as noted.")
println(io, "const MR19 = (")

RS = [0.1, 0.3, 0.7, 1.0, 1.5, 3.0, 4.5]
println(io, "  rs = ", RS, ",")
for (nm, f) in (("Sigma", DM.Sigma), ("dSigmadr", DM.dSigmadr), ("OmegaK", DM.OmegaK),
                ("Omega", DM.Omega), ("dOmegadr", DM.dOmegadr), ("kappa", DM.kappa),
                ("dkappadr", DM.dkappadr), ("csiso", DM.csiso), ("dcsisodr", DM.dcsisodr),
                ("csadi", DM.csadi), ("Nr2", DM.Nr2), ("dNr2dr", DM.dNr2dr),
                ("LSinv", DM.LSinv), ("dLSinvdr", DM.dLSinvdr), ("LTinv", DM.LTinv),
                ("dLTinvdr", DM.dLTinvdr), ("tc", DM.tc), ("dtcdr", DM.dtcdr))
    println(io, "  ", nm, " = ", repr([f(d, r) for r in RS]), ",")
end
println(io, "  Dfreq_m5 = ", repr([DM.Dfreq(d, r, 5) for r in RS]), ",")
println(io, "  dDdr_m5 = ", repr([DM.dDdr(d, r, 5) for r in RS]), ",")
println(io, "  Phi_m5 = ", repr([DM.Phi_m(d, r, 5) for r in RS]), ",")
println(io, "  dPhidr_m5 = ", repr([DM.dPhidr_m(d, r, 5) for r in RS]), ",")
println(io, "  d2Phidr2_m5 = ", repr([DM.d2Phidr2_m(d, r, 5) for r in RS]), ",")
println(io, "  kwave_m5 = ", repr([DM.kwave(d, r, 5) for r in RS]), ",")
println(io, "  wkbC_m5 = ", repr([DM.wkbC(d, r, 5) for r in RS]), ",")

# master-equation coefficients at three (r, m)
meq = Any[]
for (r, m) in ((0.5, 2), (1.3, 5), (3.0, 10))
    Om=DM.Omega(d,r); dOm=DM.dOmegadr(d,r); Sig=DM.Sigma(d,r); dSig=DM.dSigmadr(d,r)
    cs=DM.csadi(d,r); ot=DM.omega_pat(d,m)-m*Om
    Dv=DM.Dfreq(d,r,m); dD=DM.dDdr(d,r,m)
    LSi=DM.LSinv(d,r); dLSi=DM.dLSinvdr(d,r); LTi=DM.LTinv(d,r); dLTi=DM.dLTinvdr(d,r)
    N2=DM.Nr2(d,r); dN2=DM.dNr2dr(d,r); t=DM.tc(d,r); dt=DM.dtcdr(d,r)
    Ph=DM.Phi_m(d,r,m); dPh=DM.dPhidr_m(d,r,m); d2Ph=DM.d2Phidr2_m(d,r,m)
    push!(meq, collect(DM.master_equation_cool(r, d.gamma, cs, Om, dOm, Sig, dSig, LTi, dLTi,
        LSi, dLSi, ot, t, dt, Dv, dD, N2, dN2, d.pole_disp, Ph, dPh, d2Ph, m)))
end
println(io, "  meq_cases = ", repr([(0.5,2),(1.3,5),(3.0,10)]), ",")
println(io, "  meq = ", repr(meq), ",")

# forced modes on the MR19 fig-1 grid
x = DM.loggrid(d, 0.05, 5.0, 2000)
idx = round.(Int, range(1, 2000; length=21))
println(io, "  grid_n = 2000, grid_rin = 0.05, grid_rout = 5.0,")
println(io, "  grid_idx = ", repr(collect(idx)), ",")
println(io, "  grid_r = ", repr(x[idx]), ",")
for m in (2, 5, 10)
    t0=time(); fm = DM.solve_m_mode(d, m, x); el=time()-t0
    _, _, sig = DM.vars_from_enthalpy(d, fm, x, m)
    frac = [sig[i]/DM.Sigma(d, x[i]) for i in eachindex(x)]
    println(io, "  dh_m$(m) = ", repr(fm[idx,1]), ",")
    println(io, "  dhp_m$(m) = ", repr(fm[idx,2]), ",")
    println(io, "  frac_m$(m) = ", repr(frac[idx]), ",")
    @printf("forced m=%d solved in %.1f s\n", m, el); flush(stdout)
end
# cooling regime (beta = 0.1, the MR20 committed value) — one mode, to gate the beta path
fm = DM.solve_m_mode(dc, 5, x)
println(io, "  dh_cool_m5 = ", repr(fm[idx,1]), ",")

# homogeneous (passive-wake) mode, MR19 §passive / figs 10-12
dho = NewtonianDisc(h0=0.05, q=1.0, p=1.0, gamma=7/5, beta=0.0, soft=0.06, indirect=0, pole_disp=1e-6)
xh = DM.loggrid(dho, 0.05, 1.0, 2000)
idxh = round.(Int, range(1, 2000; length=21))
println(io, "  homog_h0 = 0.05, homog_rin = 0.05, homog_rout = 1.0, homog_n = 2000,")
println(io, "  homog_idx = ", repr(collect(idxh)), ",")
println(io, "  homog_r = ", repr(xh[idxh]), ",")
for m in (5, 10)
    fh = homog_mode(dho, m, xh)
    sg = sigma_from_enthalpy0(dho, fh, xh, m)
    println(io, "  homog_dh_m$(m) = ", repr(fh[idxh,1]), ",")
    println(io, "  homog_sig_m$(m) = ", repr(sg[idxh]), ",")
    println("homog m=$m done"); flush(stdout)
end
println(io, ")")
close(io)
println("FIXTURE WRITTEN")
