# Policy sources (text only): ref-paper/arXiv-2407.20398/src/Sections/
# 3_circular_benchmark.tex:13; 4_eccentric_torques.tex:19;
# 2_linear_framework.tex:192,204-213. Unprinted choices are explicit below.
using ..ForcingModes: NewtonianMode, InstantaneousScaleHeight

"""
    newtonian_mode_set(orbit; mmax, dlmax)

Ordered positive-frequency Newtonian harmonics, m >= 1, l >= 1,
|l-m| <= dlmax. This is the landed FR22 mask, generalized only in its cutoffs.
The excluded l <= 0 sector is a policy, not a claim of vanishing forcing.
"""
function newtonian_mode_set(orbit::EccentricOrbit; mmax::Integer, dlmax::Integer)
    mmax >= 0 && dlmax >= 0 || throw(ArgumentError("mode truncations must be nonnegative"))
    [NewtonianMode(m,l,orbit.n) for m in 1:mmax for l in max(1,m-dlmax):m+dlmax]
end

"""
    fr24_mode_set(orbit)

FR24 v1 mask: circular m=l=1:150; eccentric mmax=150 for e<=0.06,
175 above it, |l-m|<=40. Retains m=1. Adopts FR22's l>=1 policy because
FR24 refers to that method but does not specify l<=0; static/retrograde modes
have no prograde corotation shooting origin and are outside this anchor path.
No amplitude pruning, exterior-resonance exclusion or manual m=1 filtering.
"""
fr24_mode_set(orbit::EccentricOrbit) = newtonian_mode_set(orbit;
    mmax=orbit.e <= oftype(orbit.e,3)/50 ? 150 : 175, dlmax=iszero(orbit.e) ? 0 : 40)

# Adopted resolution, not a printed FR24 value. Keep output-grid, cache and
# angular convergence separate. The solver consumes FR22_SOLVER_SETTINGS.
const FR24_GRID_SETTINGS = (; radii=4000, cache_nodes=3997, Nt=1024, Nphi=1024)

"""
    fr24_mode_domain(d, mode)

Per-mode domain plan in FR24 code units a=n=1. Report on [0.05,5]. Shoot
from the grid point nearest the mode's exact pressure-supported corotation.

The measured FR22 near-edge modes (6:9,1) use r_out=10, checked against 12.
For other modes with a nominal resonance within 2H of the outer boundary,
or with corotation outside it, use r_out=max(10,2r_OLR), checked against
1.2r_out. These latter extensions are a declared numerical policy requiring
consumer convergence; they are not certified by the four-mode FR22 sweep.
No positive-frequency mode is dropped. A missing prograde corotation (l<=0)
raises ArgumentError, consistent with the mode-set exclusion. The return
keeps solve/check/report domains distinct and flags the measured subset.
Softening is b H(r_p), b=0.3 (FR24 2_linear_framework.tex:202).
"""
function fr24_mode_domain(d::NewtonianDisc, mode::HarmonicMode)
    mode.m > 0 && mode.m+mode.k > 0 && mode.omega == mode.m+mode.k ||
        throw(ArgumentError("FR24 domain policy requires m,l>=1 and n=1"))
    iszero(d.beta) || throw(ArgumentError("FR24 requires locally isothermal beta=0"))
    T=typeof(float(d.h0+mode.omega)); a=one(T)
    report=(a/20,5a)
    orb=EccentricOrbit(a,zero(T),a)
    roots=lindblad_locations(orb,mode.m,mode.m+mode.k)
    rc_out=mode.m*Omega(d,report[2]) > mode.omega
    near=any(r !== nothing && r+2csiso(d,r)/OmegaK(d,r) >= report[2] for r in roots)
    measured=6 <= mode.m <= 9 && mode.m+mode.k == 1
    rout=measured ? 10a : (near || rc_out ? max(10a,2roots[2]) : report[2])
    solve_domain=(report[1],rout)
    check_domain=(report[1],rout > report[2] ? 6rout/5 : report[2])
    residual(r)=mode.m*Omega(d,r)-mode.omega
    residual(solve_domain[1]) > 0 > residual(solve_domain[2]) ||
        throw(ArgumentError("mode needs a different corotation bracket"))
    rc=brent(residual,solve_domain...;atol=eps(T),rtol=8eps(T))
    (; report_domain=report, solve_domain, check_domain, corotation=rc,
       anchor=:corotation, corotation_in_report=report[1] < rc < report[2],
       edge_policy=measured ? :fr22_edge_extension : (rout>report[2] ? :requires_boundary_check : :published),
       grid=FR24_GRID_SETTINGS, softening=InstantaneousScaleHeight(T(3)/10))
end
