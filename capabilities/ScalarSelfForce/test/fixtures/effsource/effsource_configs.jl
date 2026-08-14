# effsource_configs.jl — BAKED EffectiveSource-C oracle fixture (ESSF node effsource-oracle-fixtures, W2-D).
# DO NOT EDIT BY HAND. Regenerate (one command, byte-identical on this toolchain):
#   make -C test/fixtures/effsource/generator fixtures
#
# CITATION: https://github.com/barrywardell/EffectiveSource
#   pinned SHA 1cb9834f8a8f666aaf7b1f924b3197743efb4120 (2026-04-11, MIT license)
#   read-only mirror: ref-code/EffectiveSource (PROVENANCE.elliptic-ssf.md)
# VARIANT: elliptic = kerr-circular.c (complete-elliptic-integral E/K m-mode path)
# TOOLCHAIN: gcc 13.3.0 | GSL 2.7.1
# COMPILE: cc -O2 -std=c99 $(pkg-config --cflags gsl) [-DEFFSOURCE_VARIANT_COMOVING]
#   -I$BUILD -o capture_elliptic capture.c kerr-circular.c $(pkg-config --libs gsl) -lm
# RUN: ./capture_elliptic <outdir>   (driver source: generator/capture.c)
#
# FROZEN GRID: 6 configs (a, r0) = (0,10), (0.5,10), (0.9,isco+0.1), (0.998,isco),
#   (-0.5,10), (-0.998,isco) [M=1; r0/E/Lz = solid Spacetime values, 17 digits];
#   particle at (t=0, r0, pi/2, phi_p=0.31), u^r=0; per-config points npts =
#   [40,40,38,37,40,40] in blocks r-ray-minus | r-ray-plus | theta-ray | diagonal
#   (|dr| down to 1e-6 M, out to +10 M; minus block clipped above the horizon;
#   both sides of the comoving z=0.995 switch covered in dr AND dtheta);
#   m in {0,1,2,5,10,20} (comoving PhiS_m adds {-2,-1}); dphi in {0, 0.5, 2.0}.
#
# PHASE CONVENTIONS (the phase bridge, measured by test_effsource_fixtures.jl):
#   * PhiS_m / calc_m (BOTH variants): complex quantities carry exp(-i m phi_p) —
#     pairs are (v cos(m phi_p), -v sin(m phi_p)) for real-type v [PhiS, d/dr,
#     d/dtheta, d2/dt2, d2/dtdphi, d2/dr2, d2/dtheta2, d2/dphi2, src] and
#     (v sin(m phi_p), v cos(m phi_p)) for pure-imaginary-type v [d/dt, d/dphi].
#   * calc_m_series (comoving only) uses the OPPOSITE convention exp(+i m phi_p):
#     src pairs are (s cos(m phi_p), +s sin(m phi_p)) — kerr-circular-comoving.c:673-674
#     vs :634-635. Measured bridge numbers: see effsource_configs.jl header.
#
# NaN SENTINELS: calc_m mixed second derivatives d2/dtdr, d2/dtdtheta, d2/drdtheta,
#   d2/drdphi, d2/dthetadphi are NOT provided by the C code (kerr-circular.c:709-713,
#   kerr-circular-comoving.c:570-574) and are recorded as literal NaN.
#
# COLUMNS (fixture row k = matrix column M[:, k]):
#   EFFSOURCE_CONFIGS: (idx, a, r0, E, Lz, risco, rp, alpha20, alpha02, beta, npts);
#   EFFSOURCE_POINTS: 1 cfg | 2 ipt | 3 dr | 4 dtheta.  alpha20 = r0^2/Delta,
#   alpha02 = r0^2, beta = z_c^2 = 4(Lz^2 + r0^2 + a^2(r0+2M)/r0) (kerr-circular.c:844-847;
#   main.tex:445-455). MEASURED PHASE BRIDGE + conditioning (this committed capture,
#   gated in test_effsource_fixtures.jl):
#   * PhiS_m elliptic (E/K) vs comoving (3P Legendre): SAME exp(-i m phi_p) convention,
#     worst relative difference 2.8e-10 over the overlap window C1 = alpha/beta in
#     [1e-8, 1e-2) (test gate 3e-9). Outside the window the two paths are COMPLEMENTARY
#     oracles: the elliptic E/K path loses digits near the worldline (C1 < 1e-8;
#     k' = sqrt(C1) cancellation entering Kcomp — up to ~1.2e-3 relative at C1 ~ 1e-15,
#     values stay finite), and in the far field (C1 >~ 0.1) it is the ELLIPTIC
#     E/K-polynomial route (kerr-circular.c) that cancels catastrophically for
#     m >= 2 — BigFloat(256)-quadrature adjudication (W4-D): elliptic-vs-truth up
#     to ~5.8e91 relative at C1 ~ 120, m = 20, while the comoving/3P route tracks
#     truth to ~1e-11 over the whole frozen grid. An earlier revision of this
#     header mis-attributed the far-field cancellation to the comoving route.
#     Evidence: results/elliptic-ssf/paper_ELLIPTIC-SSF/evidence/
#     mmode-puncture-integrals/verification_summary.json.
#   * calc_m DIRECT src vs calc_m_series (comoving): the direct Box[Phi^P] evaluation
#     is DESTROYED by cancellation near the worldline (e.g. |src_direct| ~ 1.4e2 vs
#     true ~ 7.0e-5 at dr = 1e-6, cfg 1 m 1) — calc_m_series is the near-worldline
#     oracle; the direct src is trustworthy only for rho^2 >~ 1e-5. Handoff agreement:
#     for EVERY (config, m) the minimum-over-points relative |src| difference is
#     <= 6.0e-4 (worst at m = 20; test gate 5e-3).
#   * calc_m_series phase is EXACTLY opposite: series = s exp(+i m phi_p) vs
#     calc_m = s' exp(-i m phi_p) — the +sin(m phi_p) bridge; all phase-rotation
#     residuals <= 1.3e-16 (test gates 1e-12).
#   * conjugate symmetry (comoving): |Psi_{-m} - conj(Psi_m)|/|Psi_m| <= 4.1e-15 for
#     m = 1, 2 (test gate 1e-12).
const EFFSOURCE_PHIP = 0.31

const EFFSOURCE_CONFIGS = [
    (idx = 1, a = 0, r0 = 10, E = 0.9561828874675149, Lz = 3.7796447300922735, risco = 6, rp = 2, alpha20 = 1.25, alpha02 = 100, beta = 457.14285714285717, npts = 40),
    (idx = 2, a = 0.5, r0 = 10, E = 0.95377548362550246, Lz = 3.5893904590904584, risco = 4.2330025295308262, rp = 1.8660254037844386, alpha20 = 1.2461059190031152, alpha02 = 100, beta = 452.73489547123842, npts = 40),
    (idx = 3, a = 0.90000000000000002, r0 = 2.4208830417618872, E = 0.8449777316180217, Lz = 2.1031359739931039, risco = 2.3208830417618871, rp = 1.4358898943540672, alpha20 = 3.2044655720003883, alpha02 = 5.8606747018902876, beta = 47.052131745150042, npts = 38),
    (idx = 4, a = 0.998, r0 = 1.2369706551751847, E = 0.67900583438380036, Lz = 1.3918131634600921, risco = 1.2369706551751847, rp = 1.0632139225171164, alpha20 = 29.335181274826571, alpha02 = 1.5300964017645255, beta = 24.294546275416153, npts = 37),
    (idx = 5, a = -0.5, r0 = 10, E = 0.95920119269209247, Lz = 4.0000118938067279, risco = 7.5545847145123579, rp = 1.8660254037844386, alpha20 = 1.2461059190031152, alpha02 = 100, beta = 465.20038060238113, npts = 40),
    (idx = 6, a = -0.998, r0 = 8.9943744548035696, E = 0.96222637470790651, Lz = 4.2326025907146425, risco = 8.9943744548035696, rp = 1.0632139225171164, alpha20 = 1.2659020710193396, alpha02 = 80.898771833223009, beta = 400.12469272078317, npts = 40),
]

const EFFSOURCE_POINTS = reshape(Float64[
    1, 1, -9.9999999999999995e-07, 0,
    1, 2, -1.0000000000000001e-05, 0,
    1, 3, -0.0001, 0,
    1, 4, -0.001, 0,
    1, 5, -0.01, 0,
    1, 6, -0.029999999999999999, 0,
    1, 7, -0.10000000000000001, 0,
    1, 8, -0.29999999999999999, 0,
    1, 9, -1, 0,
    1, 10, -3, 0,
    1, 11, 9.9999999999999995e-07, 0,
    1, 12, 1.0000000000000001e-05, 0,
    1, 13, 0.0001, 0,
    1, 14, 0.001, 0,
    1, 15, 0.01, 0,
    1, 16, 0.029999999999999999, 0,
    1, 17, 0.10000000000000001, 0,
    1, 18, 0.29999999999999999, 0,
    1, 19, 1, 0,
    1, 20, 3, 0,
    1, 21, 10, 0,
    1, 22, 0, 9.9999999999999995e-07,
    1, 23, 0, -9.9999999999999995e-07,
    1, 24, 0, 0.0001,
    1, 25, 0, -0.0001,
    1, 26, 0, 0.001,
    1, 27, 0, -0.001,
    1, 28, 0, 0.01,
    1, 29, 0, -0.01,
    1, 30, 0, 0.10000000000000001,
    1, 31, 0, -0.10000000000000001,
    1, 32, 0, 1,
    1, 33, 0, -1,
    1, 34, 0.001, 0.001,
    1, 35, -0.001, -0.001,
    1, 36, 0.01, 0.0050000000000000001,
    1, 37, -0.01, -0.0050000000000000001,
    1, 38, 0.10000000000000001, 0.050000000000000003,
    1, 39, -0.10000000000000001, -0.050000000000000003,
    1, 40, 0.5, 0.20000000000000001,
    2, 1, -9.9999999999999995e-07, 0,
    2, 2, -1.0000000000000001e-05, 0,
    2, 3, -0.0001, 0,
    2, 4, -0.001, 0,
    2, 5, -0.01, 0,
    2, 6, -0.029999999999999999, 0,
    2, 7, -0.10000000000000001, 0,
    2, 8, -0.29999999999999999, 0,
    2, 9, -1, 0,
    2, 10, -3, 0,
    2, 11, 9.9999999999999995e-07, 0,
    2, 12, 1.0000000000000001e-05, 0,
    2, 13, 0.0001, 0,
    2, 14, 0.001, 0,
    2, 15, 0.01, 0,
    2, 16, 0.029999999999999999, 0,
    2, 17, 0.10000000000000001, 0,
    2, 18, 0.29999999999999999, 0,
    2, 19, 1, 0,
    2, 20, 3, 0,
    2, 21, 10, 0,
    2, 22, 0, 9.9999999999999995e-07,
    2, 23, 0, -9.9999999999999995e-07,
    2, 24, 0, 0.0001,
    2, 25, 0, -0.0001,
    2, 26, 0, 0.001,
    2, 27, 0, -0.001,
    2, 28, 0, 0.01,
    2, 29, 0, -0.01,
    2, 30, 0, 0.10000000000000001,
    2, 31, 0, -0.10000000000000001,
    2, 32, 0, 1,
    2, 33, 0, -1,
    2, 34, 0.001, 0.001,
    2, 35, -0.001, -0.001,
    2, 36, 0.01, 0.0050000000000000001,
    2, 37, -0.01, -0.0050000000000000001,
    2, 38, 0.10000000000000001, 0.050000000000000003,
    2, 39, -0.10000000000000001, -0.050000000000000003,
    2, 40, 0.5, 0.20000000000000001,
    3, 1, -9.9999999999999995e-07, 0,
    3, 2, -1.0000000000000001e-05, 0,
    3, 3, -0.0001, 0,
    3, 4, -0.001, 0,
    3, 5, -0.01, 0,
    3, 6, -0.029999999999999999, 0,
    3, 7, -0.10000000000000001, 0,
    3, 8, -0.29999999999999999, 0,
    3, 9, 9.9999999999999995e-07, 0,
    3, 10, 1.0000000000000001e-05, 0,
    3, 11, 0.0001, 0,
    3, 12, 0.001, 0,
    3, 13, 0.01, 0,
    3, 14, 0.029999999999999999, 0,
    3, 15, 0.10000000000000001, 0,
    3, 16, 0.29999999999999999, 0,
    3, 17, 1, 0,
    3, 18, 3, 0,
    3, 19, 10, 0,
    3, 20, 0, 9.9999999999999995e-07,
    3, 21, 0, -9.9999999999999995e-07,
    3, 22, 0, 0.0001,
    3, 23, 0, -0.0001,
    3, 24, 0, 0.001,
    3, 25, 0, -0.001,
    3, 26, 0, 0.01,
    3, 27, 0, -0.01,
    3, 28, 0, 0.10000000000000001,
    3, 29, 0, -0.10000000000000001,
    3, 30, 0, 1,
    3, 31, 0, -1,
    3, 32, 0.001, 0.001,
    3, 33, -0.001, -0.001,
    3, 34, 0.01, 0.0050000000000000001,
    3, 35, -0.01, -0.0050000000000000001,
    3, 36, 0.10000000000000001, 0.050000000000000003,
    3, 37, -0.10000000000000001, -0.050000000000000003,
    3, 38, 0.5, 0.20000000000000001,
    4, 1, -9.9999999999999995e-07, 0,
    4, 2, -1.0000000000000001e-05, 0,
    4, 3, -0.0001, 0,
    4, 4, -0.001, 0,
    4, 5, -0.01, 0,
    4, 6, -0.029999999999999999, 0,
    4, 7, -0.10000000000000001, 0,
    4, 8, 9.9999999999999995e-07, 0,
    4, 9, 1.0000000000000001e-05, 0,
    4, 10, 0.0001, 0,
    4, 11, 0.001, 0,
    4, 12, 0.01, 0,
    4, 13, 0.029999999999999999, 0,
    4, 14, 0.10000000000000001, 0,
    4, 15, 0.29999999999999999, 0,
    4, 16, 1, 0,
    4, 17, 3, 0,
    4, 18, 10, 0,
    4, 19, 0, 9.9999999999999995e-07,
    4, 20, 0, -9.9999999999999995e-07,
    4, 21, 0, 0.0001,
    4, 22, 0, -0.0001,
    4, 23, 0, 0.001,
    4, 24, 0, -0.001,
    4, 25, 0, 0.01,
    4, 26, 0, -0.01,
    4, 27, 0, 0.10000000000000001,
    4, 28, 0, -0.10000000000000001,
    4, 29, 0, 1,
    4, 30, 0, -1,
    4, 31, 0.001, 0.001,
    4, 32, -0.001, -0.001,
    4, 33, 0.01, 0.0050000000000000001,
    4, 34, -0.01, -0.0050000000000000001,
    4, 35, 0.10000000000000001, 0.050000000000000003,
    4, 36, -0.10000000000000001, -0.050000000000000003,
    4, 37, 0.5, 0.20000000000000001,
    5, 1, -9.9999999999999995e-07, 0,
    5, 2, -1.0000000000000001e-05, 0,
    5, 3, -0.0001, 0,
    5, 4, -0.001, 0,
    5, 5, -0.01, 0,
    5, 6, -0.029999999999999999, 0,
    5, 7, -0.10000000000000001, 0,
    5, 8, -0.29999999999999999, 0,
    5, 9, -1, 0,
    5, 10, -3, 0,
    5, 11, 9.9999999999999995e-07, 0,
    5, 12, 1.0000000000000001e-05, 0,
    5, 13, 0.0001, 0,
    5, 14, 0.001, 0,
    5, 15, 0.01, 0,
    5, 16, 0.029999999999999999, 0,
    5, 17, 0.10000000000000001, 0,
    5, 18, 0.29999999999999999, 0,
    5, 19, 1, 0,
    5, 20, 3, 0,
    5, 21, 10, 0,
    5, 22, 0, 9.9999999999999995e-07,
    5, 23, 0, -9.9999999999999995e-07,
    5, 24, 0, 0.0001,
    5, 25, 0, -0.0001,
    5, 26, 0, 0.001,
    5, 27, 0, -0.001,
    5, 28, 0, 0.01,
    5, 29, 0, -0.01,
    5, 30, 0, 0.10000000000000001,
    5, 31, 0, -0.10000000000000001,
    5, 32, 0, 1,
    5, 33, 0, -1,
    5, 34, 0.001, 0.001,
    5, 35, -0.001, -0.001,
    5, 36, 0.01, 0.0050000000000000001,
    5, 37, -0.01, -0.0050000000000000001,
    5, 38, 0.10000000000000001, 0.050000000000000003,
    5, 39, -0.10000000000000001, -0.050000000000000003,
    5, 40, 0.5, 0.20000000000000001,
    6, 1, -9.9999999999999995e-07, 0,
    6, 2, -1.0000000000000001e-05, 0,
    6, 3, -0.0001, 0,
    6, 4, -0.001, 0,
    6, 5, -0.01, 0,
    6, 6, -0.029999999999999999, 0,
    6, 7, -0.10000000000000001, 0,
    6, 8, -0.29999999999999999, 0,
    6, 9, -1, 0,
    6, 10, -3, 0,
    6, 11, 9.9999999999999995e-07, 0,
    6, 12, 1.0000000000000001e-05, 0,
    6, 13, 0.0001, 0,
    6, 14, 0.001, 0,
    6, 15, 0.01, 0,
    6, 16, 0.029999999999999999, 0,
    6, 17, 0.10000000000000001, 0,
    6, 18, 0.29999999999999999, 0,
    6, 19, 1, 0,
    6, 20, 3, 0,
    6, 21, 10, 0,
    6, 22, 0, 9.9999999999999995e-07,
    6, 23, 0, -9.9999999999999995e-07,
    6, 24, 0, 0.0001,
    6, 25, 0, -0.0001,
    6, 26, 0, 0.001,
    6, 27, 0, -0.001,
    6, 28, 0, 0.01,
    6, 29, 0, -0.01,
    6, 30, 0, 0.10000000000000001,
    6, 31, 0, -0.10000000000000001,
    6, 32, 0, 1,
    6, 33, 0, -1,
    6, 34, 0.001, 0.001,
    6, 35, -0.001, -0.001,
    6, 36, 0.01, 0.0050000000000000001,
    6, 37, -0.01, -0.0050000000000000001,
    6, 38, 0.10000000000000001, 0.050000000000000003,
    6, 39, -0.10000000000000001, -0.050000000000000003,
    6, 40, 0.5, 0.20000000000000001,
], 4, :)
