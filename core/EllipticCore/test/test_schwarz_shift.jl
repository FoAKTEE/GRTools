# test_schwarz_shift.jl — grounded gates for src/schwarz_shift.jl (DAG node:
# schwarz-complex-shift, elliptic-ssf GLOBAL_DAG row 18 + TOPOLOGICAL_SCHEDULE W7-A).
# Written BEFORE the implementation (TDD red run captured in
# evidence/schwarz-complex-shift/red_targeted.log).
#
# Target: the CRP main.tex:680 complex-shift POLICY — the operator-level shift value
# s = (mΩ)² with its slicing policy (t-slicing: essential; hyperboloidal/vtu: disabled
# — the shipped vtu production yaml sets ComplexShift: 0.0, anc yaml:194,203), its
# MASS-WEIGHTED application route on the massive Gauss operator (smoother built on the
# damped β − i·s coefficient ⇒ subdomain matrices A_ss − i·s·M_ss), and the
# Schwarz-preconditioned-GMRES glue. The literal −i·s·I MECHANISM (spectre
# ExplicitInverse convention) is the sealed W6-B schwarz-smoother scope, gated bitwise
# in test_schwarz.jl at s = 0.7; gates (c)/(e) below CITE those and extend them to the
# POLICY value.
#
# Model problems (coefficient callables only — no SSF dependency; the real m-mode
# operator integration is the solver-stack-ssf node's scope):
#   * t-sliced CHARACTER model: 2D complex Helmholtz on [0,1]² (2×2 elements, N = 14,
#     ndofs = 784), A = −∇·∇ + β with indefinite real β = −k²·(1 + 0.3x) and the
#     oscillatory manufactured solution e^{ik(x + y/2)}·sin(πy) at the operator's own
#     wavenumber k = mΩ (m = 4-like scale: Ω = k/4 ≈ 3.15). The solver-killing regime
#     is k AT A SUBDOMAIN RESONANCE (calibration rounds 2-3: away from resonances the
#     unshifted smoother is fine; on one, the patch matrix is numerically singular —
#     σ_min ~ 8e-14 measured — and its solves amplify the resonant mode: the
#     main.tex:680 "noise ... excites oscillations" pathology). The test computes the
#     resonance from ITS OWN patch matrices (generalized eigenproblem, public
#     subdomain_matrix API), so the knife-edge tracks the platform's arithmetic.
#   * vtu CHARACTER model: the null-slicing phase absorption cancels the −k² bulk term
#     and leaves complex first-order advection γ¹ = −2ik·(1 + 0.3x) (piecewise-null
#     damped, non-oscillatory tails), same mesh, same k.
#
# Oracle gates (worker node contract; constants calibrated PRE-RED in
# evidence/schwarz-complex-shift/calibrate_shift{1..4}.jl + calibration_shift{1..4}.log
# — THROWAWAY prototypes composing sealed W6 APIs, NOT the implementation):
#   (a) at k = k_res: unshifted Schwarz-pGMRES (overlap 2, 2 passes, reltol 1e-10)
#       FAILS its full 150-iteration budget with the residual plateaued ≥ 1e6× above
#       tolerance (round 4, measured: res 1206.7 ≈ 2.7·‖b‖ after the full 150 its,
#       res/tol = 2.7e10); the mass-weighted policy shift s = (mΩ)² = k² converges
#       (measured 47 its, gate CAL4_A_GATE); 3× shifted its fit inside the unshifted
#       budget (the row-18 >3x CONTRAST); 1% off resonance the unshifted run converges
#       again (measured 26 its) — the failure is the resonance, not the model.
#   (b) vtu model, same mesh/config, s = helmholtz_shift(m, Ω; slicing = :vtu) = 0:
#       converges in ≤ CAL4_B_GATE its (measured 15).
#   (c) literal seam at the policy value: SchwarzSmoother(op; shift = s) patch ==
#       unshifted − i·s·I BITWISE (extends the s = 0.7 gate of test_schwarz.jl);
#       mass-weighted route: smoother-on-(β − i·s) patch == unshifted − i·s·M_ss at
#       rel ≤ 1e-12, M_ss extracted as the patch difference of (β+1) vs β operators.
#   (d) the shift does not change the CONVERGED solution (right preconditioning):
#       u(s) vs u(2s) at reltol 1e-12 agree ≤ CAL4_D_GATE rel (measured 9.66e-11), both
#       true-residual-verified, and both match the manufactured discrete solution
#       (b := A·u_ex) to ≤ CAL4_D_UEX (measured 9.69e-11).
#   (e) cache rekey at the policy value (extends the test_schwarz.jl builds/hits
#       counter gate): 0 → s rebuilds all, repeat s is all hits.
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting siblings).

using Test
using LinearAlgebra
using EllipticCore

const SHF_C = 1.5                                    # CRP main.tex:595 penalty parameter
const SHF_RTOL = 1e-10                               # solver-behavior gate tolerance

# Oracle constants, calibrated PRE-RED in round 4 (evidence/schwarz-complex-shift/
# calibrate_shift4.jl + calibration_shift4.log — the round runs EXACTLY the configurations
# gated below, with this file's helpers copied verbatim). Each is the measured value plus a
# small margin; the measurement is quoted so a drift is diagnosable, not just red.
const CAL4_A_GATE = 49        # (a) policy-shift solve at k_res: measured 47 its at reltol 1e-10.
                              #     Capped at 49, NOT higher: the row-18 >3x contrast gate
                              #     `3*its < 150` (the failed unshifted budget) admits at most 49.
const CAL4_OFF_GATE = 30      # (a) unshifted 1% OFF resonance, maxiter 100: measured 26 its.
const CAL4_B_GATE = 19        # (b) vtu-character model at s = 0, maxiter 100: measured 15 its.
const CAL4_D_GATE = 1e-9      # (d) |u(s) − u(2s)| / |u(s)| at reltol 1e-12: measured 9.66e-11.
const CAL4_D_UEX = 1e-9       # (d) |u(s) − u_ex| / |u_ex| at reltol 1e-12: measured 9.69e-11.

shf_one2(x, y) = 1.0
shf_zero2(x, y) = 0.0
shf_profile(x, y) = 1 + 0.3x                         # smooth O(1) radial-like profile

# 2×2 level-1 mesh on the unit block, N×N Gauss points per element
shf_mesh(N) = Mesh2D(Block(81, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0)),
                     vec([ElementId(81, SegmentId(1, ix), SegmentId(1, iy))
                          for ix in 0:1, iy in 0:1]), (N, N))

function shf_zerobcs(mesh)
    bcs = Dict{Any,Any}()
    for (el, d, sidx) in boundary_faces(mesh)
        nf = mesh.orders[mesh.index[el]][3-d]
        bcs[(el, d, sidx)] = DirichletGhost(zeros(ComplexF64, nf))
    end
    return bcs
end

# t-sliced-character operator with an OPERATOR-level complex shift s (0 = original):
# β = −k²·profile − i·s ⇒ subdomain matrices A_ss − i·s·M_ss (the mass-weighted route)
shf_top(mesh, k; s = 0.0, βoff = 0.0) =
    MeshOperator(ComplexF64, mesh, (shf_one2, shf_one2),
                 (x, y) -> -k^2 * shf_profile(x, y) - im * s + βoff,
                 (shf_zero2, shf_zero2); C = SHF_C, bcs = shf_zerobcs(mesh))
# vtu-character operator: null-phase absorption cancels the −k² bulk term, complex
# first-order advection remains (non-oscillatory damped tails)
shf_vop(mesh, k) =
    MeshOperator(ComplexF64, mesh, (shf_one2, shf_one2), shf_zero2,
                 ((x, y) -> -2im * k * shf_profile(x, y), shf_zero2);
                 C = SHF_C, bcs = shf_zerobcs(mesh))

# nodal values of f(x, y) on the mesh (physical tensor grid per element)
function shf_fieldvec(mesh, f)
    u = zeros(ComplexF64, mesh.ndofs)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        blk = mesh.blocks[mesh.blockindex[el.block]]
        Nx, Ny = mesh.orders[e]
        grid = ntuple(d -> begin
            lo, hi = logical_interval(Float64, el, d)
            [map_forward(blockmap(blk, d), (lo * (1 - ξ) + hi * (1 + ξ)) / 2)
             for ξ in EllipticCore.refelement(Float64, d == 1 ? Nx : Ny).nodes]
        end, 2)
        for j in 1:Ny, i in 1:Nx
            u[mesh.offsets[e] + tensor_lin(i, j, Nx)] = f(grid[1][i], grid[2][j])
        end
    end
    return u
end

# manufactured RHS: b := A·u_ex (linearized apply) makes u_ex the EXACT discrete
# solution — converged GMRES must reproduce it (oscillatory for the t-sliced model)
shf_uex_t(mesh, k) = shf_fieldvec(mesh, (x, y) -> cis(k * (x + 0.5y)) * sinpi(y))
shf_uex_vtu(mesh) = shf_fieldvec(mesh, (x, y) -> exp(-2.0x) * (1 + 0.5x) * sinpi(y))
shf_rhs(op, uex) = apply!(zeros(ComplexF64, op.mesh.ndofs), op, uex)

# the patch-1 resonant wavenumber nearest 14: β enters the patch matrix as the exactly
# diagonal −k²·(profile-mass), so A_ss(k) = A_ss(0) − k²·B with B = A_ss(0) − A_ss(1)
# and the resonances are the generalized eigenvalues k² of (A_ss(0), B)
function shf_kres(mesh)
    A0 = subdomain_matrix(SchwarzSmoother(shf_top(mesh, 0.0); overlap = 2), 1)
    B = A0 .- subdomain_matrix(SchwarzSmoother(shf_top(mesh, 1.0); overlap = 2), 1)
    ks = [sqrt(real(λ)) for λ in eigvals(Matrix(A0), Matrix(B))
          if abs(imag(λ)) < 1e-8 * abs(λ) && real(λ) > 0 && 10 ≤ sqrt(real(λ)) ≤ 20]
    return ks[argmin(abs.(ks .- 14.0))]              # ≈ 12.588 on the calibration host
end

# Schwarz-preconditioned GMRES on op, smoother built on smop (== op ⇒ unshifted)
function shf_solve(op, smop, b; maxiter = 100, reltol = SHF_RTOL)
    sm = SchwarzSmoother(smop; overlap = 2)
    return gmres(op, b; abstol = 0.0, reltol = reltol, maxiter = maxiter,
                 Minv! = schwarz_preconditioner(sm; passes = 2))
end

@testset "schwarz shift policy (schwarz-complex-shift)" begin

@testset "helmholtz_shift policy values" begin
    Ω = 0.3
    @test helmholtz_shift(4, Ω) == abs2(4 * Ω)               # t-slicing is the default
    @test helmholtz_shift(4, Ω; slicing = :t) == abs2(4 * Ω)
    # the CRP claim (main.tex:680) + shipped vtu yaml ComplexShift: 0.0: disabled
    # under hyperboloidal/null slicing
    @test helmholtz_shift(4, Ω; slicing = :vtu) == 0
    @test helmholtz_shift(4, Ω; slicing = :hyperboloidal) == 0
    @test helmholtz_shift(0, Ω) == 0                          # m = 0: nothing oscillates
    @test helmholtz_shift(4, -Ω) == helmholtz_shift(4, Ω)     # retrograde: even in mΩ
    @test helmholtz_shift(4, big"0.25") isa BigFloat          # precision-inheriting
    @test helmholtz_shift(4, big"0.25") == 1
    @test_throws ArgumentError helmholtz_shift(4, Ω; slicing = :radial)
end

@testset "glue: schwarz_preconditioner == passes× smooth! from zero (bitwise)" begin
    op = shf_top(shf_mesh(8), 12.0)
    n = op.mesh.ndofs
    sm = SchwarzSmoother(op; overlap = 2)
    v = shf_rhs(op, shf_uex_t(op.mesh, 12.0))
    for p in (1, 2, 3)
        pc = schwarz_preconditioner(sm; passes = p)
        z = fill(complex(7.0, -7.0), n)                       # stale data must be reset
        @test pc(z, v) === z
        @test z == smooth!(zeros(ComplexF64, n), sm, v; passes = p)   # bitwise
    end
    @test_throws ArgumentError schwarz_preconditioner(sm; passes = 0)
end

@testset "(a) t-sliced resonance: unshifted pGMRES fails, policy shift converges" begin
    mesh = shf_mesh(14)
    kres = shf_kres(mesh)
    @test 10 < kres < 20
    op = shf_top(mesh, kres)
    # precondition: ON the resonance the unshifted patch is numerically singular —
    # the main.tex:680 spurious-oscillation amplifier (calibrated σ_min = 8.2e-14)
    σmin = minimum(svdvals(subdomain_matrix(SchwarzSmoother(op; overlap = 2), 1)))
    @test σmin < 1e-8
    uex = shf_uex_t(mesh, kres)
    b = shf_rhs(op, uex)
    tol = SHF_RTOL * norm(b)
    # NODE GATE: unshifted burns the FULL 150-iteration budget, plateaued ≥ 1e6× tol
    r0 = shf_solve(op, op, b; maxiter = 150)
    @test !r0.converged
    @test r0.iterations == 150
    @test r0.residual > 1e6 * tol
    # NODE GATE: the mass-weighted policy shift converges well inside the budget
    s = helmholtz_shift(4, kres / 4)                          # m = 4-like: Ω = k/4
    @test s == kres^2                                         # exact (×4 is a pow2)
    rs = shf_solve(op, shf_top(mesh, kres; s = s), b; maxiter = 150)
    @test rs.converged
    @test rs.residual ≤ tol                                   # true residual, recomputed
    @test rs.iterations ≤ CAL4_A_GATE
    # row-18 CONTRAST > 3x: three shifted solves fit inside the failed unshifted budget
    @test 3 * rs.iterations < 150
    # falsifier of misattribution: 1% OFF the resonance the unshifted run converges
    # (measured 26 its) — the failure above is the resonance, not the model
    koff = kres * 1.01
    opoff = shf_top(mesh, koff)
    roff = shf_solve(opoff, opoff, shf_rhs(opoff, shf_uex_t(mesh, koff)))
    @test roff.converged
    @test roff.iterations ≤ CAL4_OFF_GATE
end

@testset "(b) vtu-character model: converges fast with the policy s = 0" begin
    mesh = shf_mesh(14)
    kres = shf_kres(mesh)
    op = shf_vop(mesh, kres)
    b = shf_rhs(op, shf_uex_vtu(mesh))
    s = helmholtz_shift(4, kres / 4; slicing = :vtu)
    @test s == 0                                              # the disabled-shift policy
    rv = shf_solve(op, op, b)                                 # s = 0 ⇒ smoother on op
    @test rv.converged
    @test rv.residual ≤ SHF_RTOL * norm(b)
    @test rv.iterations ≤ CAL4_B_GATE
end

@testset "(c) shifted patches at the policy value: literal −i·s·I and −i·s·M_ss" begin
    op = shf_top(shf_mesh(8), 14.0)
    s = helmholtz_shift(4, 3.5)                               # (mΩ)² = 196
    # literal seam (sealed W6-B mechanism, extends the s = 0.7 bitwise gate of
    # test_schwarz.jl to the policy value)
    cache = SubdomainCache()
    sm0 = SchwarzSmoother(op; overlap = 2, cache = cache)
    sml = SchwarzSmoother(op; overlap = 2, shift = s, cache = cache)
    for sd in eachindex(sml.subdomains)
        A0 = subdomain_matrix(sm0, sd)
        @test subdomain_matrix(sml, sd) == A0 .- (im * s) .* Matrix(I, size(A0)...)
    end
    # mass-weighted route: the β − i·s smoother operator's patches equal
    # A_ss − i·s·M_ss with M_ss the (exactly diagonal) patch mass, extracted as the
    # patch difference of the (β+1) and β operators
    smw = SchwarzSmoother(shf_top(shf_mesh(8), 14.0; s = s); overlap = 2)
    smm = SchwarzSmoother(shf_top(shf_mesh(8), 14.0; βoff = 1.0); overlap = 2)
    for sd in eachindex(smw.subdomains)
        A0 = subdomain_matrix(sm0, sd)
        Mss = subdomain_matrix(smm, sd) .- A0                 # = patch mass diagonal
        @test Mss == Diagonal(Mss)                            # exactly diagonal
        @test all(>(0), real.(diag(Mss))) && all(==(0), imag.(diag(Mss)))
        As = subdomain_matrix(smw, sd)
        @test maximum(abs, As .- (A0 .- (im * s) .* Mss)) ≤ 1e-12 * maximum(abs, A0)
    end
end

@testset "(d) the shift does not change the CONVERGED solution" begin
    mesh = shf_mesh(14)
    kres = shf_kres(mesh)
    op = shf_top(mesh, kres)
    uex = shf_uex_t(mesh, kres)
    b = shf_rhs(op, uex)
    s = helmholtz_shift(4, kres / 4)
    rs = shf_solve(op, shf_top(mesh, kres; s = s), b; reltol = 1e-12, maxiter = 150)
    r2 = shf_solve(op, shf_top(mesh, kres; s = 2s), b; reltol = 1e-12, maxiter = 150)
    # both are residual-verified converged solutions of the ORIGINAL (undamped) system
    tol12 = 1e-12 * norm(b)
    @test rs.converged && r2.converged
    @test rs.residual ≤ tol12 && r2.residual ≤ tol12
    # NODE GATE: same answer — the shift lives ONLY in the preconditioner
    @test norm(rs.x .- r2.x) ≤ CAL4_D_GATE * norm(rs.x)
    # and it is the manufactured exact discrete solution (b := A·u_ex)
    @test norm(rs.x .- uex) ≤ CAL4_D_UEX * norm(uex)
end

@testset "(e) cache rekey at the policy value (builds/hits counters)" begin
    # extends the test_schwarz.jl counter gate to the policy value on the model mesh
    op = shf_top(shf_mesh(8), 14.0)
    nel = length(op.mesh.elements)
    s = helmholtz_shift(4, 3.5)
    cache = SubdomainCache()
    sm0 = SchwarzSmoother(op; overlap = 2, cache = cache)
    @test cache.builds == nel && cache.hits == 0
    sms = SchwarzSmoother(op; overlap = 2, shift = s, cache = cache)
    @test cache.builds == 2nel && cache.hits == 0             # s is a NEW cache key
    smr = SchwarzSmoother(op; overlap = 2, shift = s, cache = cache)
    @test cache.builds == 2nel && cache.hits == nel           # same s: all hits
end

end # testset schwarz shift policy
