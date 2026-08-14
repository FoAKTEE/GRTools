# Grounded tests for SelfForce.LorenzField boundary-condition series (homogeneous in/out asymptotic
# seeds for the Lorenz-gauge radial ODE).
#
# GROUNDING — characterization vs legacy + analytic identities (no heavy live dep):
#
# (1) CHARACTERIZATION: the baked R / R_rs_deriv literals below were captured ONE-SHOT from the legacy
#     H1Lorenz boundary functions (a HEAVY HDF5 package — never a live test-dep) via
#     scratchpad/capture_bc.jl (Pkg.activate packages/H1Lorenz), across every parity sector and both
#     boundaries, including the eigenproblem (ωp)-scaling selector branches. The port is a line-for-
#     line transcription, so it reproduces them to ~1e-12 (near-zero entries via atol).
# (2) ANALYTIC ASYMPTOTIC IDENTITY: at the far/near boundary the leading behaviour is the pure
#     e^{±iωr*} mode, so R_rs_deriv ≈ ±iω·R (outgoing +, ingoing −) up to the tiny 1/r corrections.
# (3) ODE-RESIDUAL via the CONSUMED RHS: each series, substituted into
#     MetricReconstruction.Lorenz.LorenzRadial.{even_ode!,odd_ode!}, satisfies the field equation — the
#     analytic R'' from the RHS matches a finite-difference second derivative of the series across rs.
# (4) NULL-BASIS Wronskian/independence: the consumed `_nullbasis` returns a linearly-independent set
#     spanning the kernel of a singular matrix (A·v ≈ 0), and the even/odd leading-power eigenvalues
#     are the expected integers.

using Test
using LinearAlgebra: norm, dot, det
using SelfForce.LorenzField
const LF = SelfForce.LorenzField
import MetricReconstruction
const LR = MetricReconstruction.Lorenz.LorenzRadial

# helper: compare a returned complex vector to baked (re,im) pairs (atol covers ~1e-12 roundoff zeros)
function _cmp(got, pairs; rtol=1e-11, atol=1e-9)
    @test length(got) == length(pairs)
    for (g, (re, im)) in zip(got, pairs)
        @test real(g) ≈ re rtol=rtol atol=atol
        @test imag(g) ≈ im rtol=rtol atol=atol
    end
end

@testset "even (l=2,ω=0.1,p=10) out/in vs baked legacy H1Lorenz" begin
    R,dR,rsb = LF.even_out_bcs([1.0,0,0,0,0], 2, 0.1, 10.0)
    @test rsb == 110000.0
    _cmp(R, [(-0.089778729509831887,-0.30556840251996431),(0.089777091570309484,0.30556284728940913),
             (0.35910837082820063,1.2222513878202017),(0.089777091570309525,0.3055628472894093),
             (-0.17955418162496509,-0.61112569502463088)])
    _cmp(dR, [(0.030556837740701229,-0.0089778721580057782),(-0.030556282202750362,0.0089777084145640836),
              (-0.12222512867734069,0.035910834112976908),(-0.030556282202750376,0.0089777084145640888),
              (0.061112564450054453,-0.017955416677554684)])
    R2,_,_ = LF.even_out_bcs([0,1.0,0,0,0], 2, 0.1, 10.0)
    _cmp(R2, [(-0.12285050427026908,-0.41714873667750707),(0.13819891709969237,0.46926550247730203),
              (0.0041682065359989143,0.014153503611262363),(-0.26521512458147734,-0.90055920606335227),
              (-0.016672852456894375,-0.056614006696765057)])
    Ri,dRi,rsi = LF.even_in_bcs([1.0,0,0,0,0], 2, 0.1)
    @test rsi == -50.0
    _cmp(Ri, [(0.28366218547791983,-0.95892427468175689),(0,0),(0,0),(0,0),(0,0)])
    _cmp(dRi, [(-0.095892427460828919,-0.02836621855710119),(0,0),(0,0),(0,0),(0,0)])
end

@testset "monopole (l=0,ω=0.05) out/in vs baked legacy" begin
    R,dR,rsb = LF.monopole_out_bcs([1.0,0,0], 0, 0.05)
    _cmp(R, [(-0.59899202314035582,0.80075497072747048),(0.00014562071660863602,0.00010892926118399734),
             (0.00014562071660863602,0.00010892926118399734)])
    Ri,dRi,_ = LF.monopole_in_bcs([0,1.0,0], 0, 0.05)
    _cmp(Ri, [(0,0),(-0.80114361554241009,0.5984721441018036),(0,0)])
end

@testset "odd (l=2,ω=0.1,p=10) out/in vs baked legacy" begin
    R,dR,_ = LF.odd_out_bcs([1.0,0], 2, 0.1, 10.0)
    _cmp(R, [(-0.23274324542951563,-0.79109456794870825),(0.23274324150545544,0.79109456910421716)])
    R2,_,_ = LF.odd_out_bcs([0,1.0], 2, 0.1, 10.0)
    _cmp(R2, [(-0.079953625295565711,-0.27130687139013593),(-0.31981450252814436,-1.0852274851640038)])
    Ri,dRi,_ = LF.odd_in_bcs([1.0,0], 2, 0.1)
    _cmp(Ri, [(0.28366218546910371,-0.95892427467058594),(0,0)])
end

@testset "odd-dipole (l=1,ω=0.05,p=10) out/in vs baked legacy" begin
    R,dR,_ = LF.odd_dipole_out_bcs([1.0], 1, 0.05, 10.0)
    _cmp(R, [(-0.59928325993845666,0.80053711567063779)])
    Ri,_,_ = LF.odd_dipole_in_bcs([1.0], 1, 0.05)
    _cmp(Ri, [(-0.80114361553788649,0.59847214409965077)])
end

@testset "even-dipole (l=1) out (both ω branches) / in vs baked legacy" begin
    Rs,_,_ = LF.even_dipole_out_bcs([1.0,0,0,0], 1, 0.005, 10.0)   # small-ω fixed-seed branch
    _cmp(Rs, [(0.97523464824100503,0.22126242153350945),(-0.97521691350412165,-0.22125839785846166),
              (-1.9504338669302892,-0.44251661826175681),(3.9008677822146494,0.88503302706776688)])
    Rl,_,_ = LF.even_dipole_out_bcs([1.0,0,0,0], 1, 0.5, 10.0)     # else branch
    _cmp(Rl, [(-0.99056861453040013,-0.1370175873598902),(-2.4920509384837591e-6,1.8013868046278868e-5),
              (-4.9841018650527343e-6,3.6027736094205594e-5),(-2.4920509384837591e-6,1.8013868046278868e-5)])
    Ri,_,_ = LF.even_dipole_in_bcs([1.0,0,0,0], 1, 0.05)
    _cmp(Ri, [(-0.8011436155514573,0.59847214410610927),(0,0),(0,0),(0,0)])
end

@testset "even-static (ω=0,l=2) out/in vs baked legacy" begin
    R,dR,rsb = LF.even_static_out_bcs([1.0,0,0], 2, 0.0)
    @test rsb == 110000.0
    @test all(imag.(R) .== 0)                       # static ⇒ purely real
    _cmp(R, [(8.2678191626204922e-11,0),(8.267969524791826e-11,0),(1.5036696520307697e-15,0)]; atol=1e-13)
    Ri,_,rsi = LF.even_static_in_bcs([1.0,0,0], 2, 0.0)
    @test rsi == -30.0
    _cmp(Ri, [(3.7992491912117358e-14,0),(0.99999999999997469,0),(7.5984961971727324e-14,0)]; atol=1e-10)
end

@testset "asymptotic identity: R_rs_deriv ≈ ±iω·R at the boundary" begin
    R,dR,_ = LF.even_out_bcs([1.0,0,0,0,0], 2, 0.1, 10.0)
    @test dR ≈ im*0.1 .* R rtol=1e-4              # outgoing leading e^{+iωr*}
    Ri,dRi,_ = LF.even_in_bcs([1.0,0,0,0,0], 2, 0.1)
    @test dRi[1] ≈ -im*0.1*Ri[1] rtol=1e-4        # ingoing leading e^{−iωr*}
    Ro,dRo,_ = LF.odd_out_bcs([1.0,0], 2, 0.1, 10.0)
    @test dRo ≈ im*0.1 .* Ro rtol=1e-4
end

# ODE-residual: feed the series (R, R') at a moderate boundary into the CONSUMED LorenzRadial RHS and
# check the analytic R'' it returns against a central finite difference of the series across rs. The
# match radius for ω=0.1 is ≈250 M, so r≈400 M (rs≈411) is safely inside the convergent regime.
function _ode_residual(bcfn, nf, odefn, l, ω; rs0=411.0, h=0.05)
    R0,dR0,b0 = bcfn(rs0); Rp,_,_ = bcfn(rs0+h); Rm,_,_ = bcfn(rs0-h)
    @test b0 ≈ rs0                                   # no boundary-push fallback (series converged)
    y = zeros(Float64, 4nf)
    for k in 0:(nf-1)
        y[1+4k]=real(R0[k+1]); y[2+4k]=real(dR0[k+1]); y[3+4k]=imag(R0[k+1]); y[4+4k]=imag(dR0[k+1])
    end
    g = similar(y); odefn(g, rs0, y, l, ω)
    for k in 0:(nf-1)
        Rpp_ode = complex(g[2+4k], g[4+4k])
        Rpp_fd  = (Rp[k+1] - 2R0[k+1] + Rm[k+1]) / h^2   # O(h²) central 2nd diff of the series
        @test Rpp_ode ≈ Rpp_fd rtol=1e-4 atol=1e-8       # tol = FD truncation floor (residual → 0)
    end
end

@testset "ODE-residual vs consumed LorenzRadial RHS (even & odd, l=2, ω=0.1)" begin
    _ode_residual((rs)->LF.even_out_bcs([1.0,0,0,0,0],2,0.1,10.0; rs_out=rs), 5, LR.even_ode!, 2, 0.1)
    _ode_residual((rs)->LF.odd_out_bcs([1.0,0],2,0.1,10.0; rs_out=rs), 2, LR.odd_ode!, 2, 0.1)
end

@testset "null-basis Wronskian/independence + leading-power spectrum" begin
    # singular 3×3 with a 1-dim kernel
    A = Float64[1 2 3; 2 4 6; 1 1 1]
    vecs = LF._nullbasis(A)
    @test length(vecs) == 1
    @test maximum(abs, A*vecs[1]) < 1e-10            # A·v ≈ 0
    @test abs(norm(vecs[1]) - 1) < 1e-12             # normalized
    # 2-dim kernel: independence (Gram matrix full rank)
    B = Float64[1 0 0; 0 0 0; 0 0 0]
    v2 = LF._nullbasis(B)
    @test length(v2) == 2
    G = [dot(v2[i], v2[j]) for i in 1:2, j in 1:2]
    @test abs(det(G)) > 1e-6                          # linearly independent
    # even 5×5 leading-power eigenvalues are integers (the asymptotic powers)
    l = 2; L = l*(l+1); lam = (l+2)*(l-1)
    a_data = Float64[ -(2+L) 2 2 2 0; 2 -(2+L) -2 -2 0; 2L -2L -(4+L) -2L 2; 2 -2 -2 -(2+L) 0; 0 0 2lam 0 -(L-2) ]
    λs,_,_ = LF._eig_sorted(permutedims(a_data))
    @test all(λ -> abs(real(λ) - round(real(λ))) < 1e-8, λs)
end

@testset "boundary BigFloat precision path (type-generic {T})" begin
    setprecision(BigFloat, 200) do
        R,dR,rsb = LF.odd_out_bcs(BigFloat[1,0], 2, BigFloat(1)/10, BigFloat(10))
        @test eltype(R) == Complex{BigFloat}
        @test all(isfinite, R) && all(isfinite, dR)
        # The eigenproblem seed is a HOMOGENEOUS solution defined only up to scale (the BigFloat
        # _eig_sorted recomputes eigenvectors via nullspace, so its overall sign/scale may differ from
        # LAPACK's). Physical invariant: BigFloat and Float64 span the SAME ray — the per-component
        # ratio Rᵇⁱᵍ/Rᶠˡᵒᵃᵗ is one common (real) constant. (And dR ≈ iω R asymptotically, as above.)
        Rf,dRf,_ = LF.odd_out_bcs([1.0,0], 2, 0.1, 10.0)
        ratio = ComplexF64(R[1]) / Rf[1]
        @test ComplexF64(R[2]) / Rf[2] ≈ ratio rtol=1e-9
        @test abs(imag(ratio)) < 1e-9                 # real scale factor (same homogeneous solution)
        @test dR ≈ im*(BigFloat(1)/10) .* R rtol=1e-4
    end
end
