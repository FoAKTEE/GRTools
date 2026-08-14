/* capture.c — frozen-grid characterization capture for barrywardell/EffectiveSource.
 *
 * ESSF node `effsource-oracle-fixtures` (schedule W2-D): bakes the C effective-source
 * oracles as Julia-readable literal fixture files for the ScalarSelfForce test suite.
 * The C code itself is NEVER a live test dependency (hygiene §7.7) — this driver is the
 * one-command regeneration path (`make fixtures`, see Makefile in this directory).
 *
 * Compiled TWICE against the two upstream variants (identical public API, effsource.h):
 *   capture_elliptic — kerr-circular.c          (m-mode via complete elliptic E/K series)
 *   capture_comoving — kerr-circular-comoving.c (m-mode via half-integer Legendre P /
 *                       hypergeometric 2F1 with the z = 0.995 branch switch, line 27;
 *                       adds effsource_calc_m_series)   [-DEFFSOURCE_VARIANT_COMOVING]
 *
 * SOURCE: https://github.com/barrywardell/EffectiveSource
 *         pinned SHA 1cb9834f8a8f666aaf7b1f924b3197743efb4120 (2026-04-11, MIT)
 *         read-only mirror: ref-code/EffectiveSource (PROVENANCE.elliptic-ssf.md)
 *
 * Determinism: output is pure %.17g of double results of a fixed evaluation order —
 * recompiling + rerunning on the same toolchain reproduces every file byte-identically
 * (gated by test_effsource_fixtures.jl "re-run reproducibility").
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_version.h>
#include "effsource.h"

#ifndef M_PI_2 /* strict -std=c99 hides it */
#define M_PI_2 1.57079632679489661923
#endif

#ifdef EFFSOURCE_VARIANT_COMOVING
#define VARIANT "comoving"
#define VARIANT_UC "COMOVING"
#define VARIANT_SRC "kerr-circular-comoving.c"
#else
#define VARIANT "elliptic"
#define VARIANT_UC "ELLIPTIC"
#define VARIANT_SRC "kerr-circular.c"
#endif

/* ---------------- frozen grid (BINDING; mirrored in test_effsource_fixtures.jl) ------ */

/* Particle: circular equatorial geodesic at (t=0, r0, theta=pi/2, phi=PHIP), u^r = 0.
 * PHIP != 0 so the exp(-i m phi_p) phase factors are non-trivial and measurable.       */
static const double PHIP = 0.31;

/* r0 for the ISCO-anchored configs and every (E, Lz) below are the SOLID Spacetime
 * values, hard-coded to 17 significant digits (transcription is gated by the test):
 *   r0    : kerr_isco(a) (+0.1 for a=0.9)
 *   E, Lz : constants_of_motion(a, r0, 0.0, 1.0)  [equatorial x = +1, M = 1]
 *   rp    : horizons(a)[1]   risco : kerr_isco(a)
 * a < 0 encodes retrograde orbits (motion always +phi — Spacetime sign convention).   */
struct config { double a, r0, E, L, risco, rp; };
static const struct config CFG[6] = {
  { 0.0,   10.0,                0.9561828874675149,  3.7796447300922735,
    6.0,                2.0 },
  { 0.5,   10.0,                0.95377548362550246, 3.5893904590904584,
    4.2330025295308262, 1.8660254037844386 },
  { 0.9,   2.4208830417618872,  0.8449777316180217,  2.1031359739931039,
    2.3208830417618871, 1.4358898943540672 },
  { 0.998, 1.2369706551751847,  0.67900583438380036, 1.3918131634600921,
    1.2369706551751847, 1.0632139225171164 },
  {-0.5,   10.0,                0.95920119269209247, 4.0000118938067279,
    7.5545847145123579, 1.8660254037844386 },
  {-0.998, 8.9943744548035696,  0.96222637470790651, 4.2326025907146425,
    8.9943744548035696, 1.0632139225171164 },
};

/* Point offsets (dr, dtheta) off the worldline, in blocks:
 *   r-ray minus | r-ray plus | theta-ray | diagonal.
 * Near-worldline |dr| reaches 1e-6 M; far field reaches r0+10 M. The minus block is
 * CLIPPED per config so every field point stays safely outside the horizon:
 *   NMINUS = {10, 10, 8, 7, 10, 10}  (a=0.9: drop -1,-3; a=0.998: drop -0.3,-1,-3).
 * Points straddle the comoving z=0.995 hypergeometric switch in BOTH dr (dr* in
 * [1.1e-3, 2.4e-2] across configs) and dtheta (dth* in [2.7e-3, 5.0e-3]).           */
static const double DR_MINUS[10] = {-1e-6,-1e-5,-1e-4,-1e-3,-1e-2,-3e-2,-1e-1,-3e-1,-1.0,-3.0};
static const double DR_PLUS[11]  = { 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0};
static const int    NMINUS[6]    = {10, 10, 8, 7, 10, 10};
static const double DTH[12] = { 1e-6,-1e-6, 1e-4,-1e-4, 1e-3,-1e-3, 1e-2,-1e-2, 1e-1,-1e-1, 1.0,-1.0};
static const double DIAG[7][2] = {
  { 1e-3, 1e-3}, {-1e-3,-1e-3}, { 1e-2, 5e-3}, {-1e-2,-5e-3},
  { 1e-1, 5e-2}, {-1e-1,-5e-2}, { 5e-1, 2e-1} };
#define MAXPTS 40

/* m-mode lists. The elliptic variant supports 0 <= m <= 20 only (ReEI/ImEI tables);
 * the comoving Legendre path also evaluates m < 0, captured for the conjugate-symmetry
 * gate Psi_{-m} = conj(Psi_m).                                                        */
static const int MS[6] = {0, 1, 2, 5, 10, 20};
#ifdef EFFSOURCE_VARIANT_COMOVING
static const int MS_PHIS_M[8] = {-2, -1, 0, 1, 2, 5, 10, 20};
#define N_MS_PHIS_M 8
#else
#define MS_PHIS_M MS
#define N_MS_PHIS_M 6
#endif

/* 4D singular field is captured at three azimuthal offsets dphi = phi - phi_p.        */
static const double DPHI[3] = {0.0, 0.5, 2.0};

/* ---------------- emission helpers --------------------------------------------------- */

/* Julia-readable double literal; NaN/Inf sentinels spelled the Julia way. */
static void jd(FILE *f, double x)
{
  if (isnan(x))      fputs("NaN", f);
  else if (isinf(x)) fputs(x > 0 ? "Inf" : "-Inf", f);
  else               fprintf(f, "%.17g", x);
}

static void row_begin(FILE *f) { fputs("    ", f); }
static void cell(FILE *f, double x) { jd(f, x); fputs(", ", f); }
static void row_end(FILE *f, double x) { jd(f, x); fputs(",\n", f); }

/* Shared CITATION header. KEEP FREE OF TIMESTAMPS (byte-identical regeneration). */
static void header(FILE *f, const char *title, const char *cols)
{
  fprintf(f,
    "# %s — BAKED EffectiveSource-C oracle fixture (ESSF node effsource-oracle-fixtures, W2-D).\n"
    "# DO NOT EDIT BY HAND. Regenerate (one command, byte-identical on this toolchain):\n"
    "#   make -C test/fixtures/effsource/generator fixtures\n"
    "#\n"
    "# CITATION: https://github.com/barrywardell/EffectiveSource\n"
    "#   pinned SHA 1cb9834f8a8f666aaf7b1f924b3197743efb4120 (2026-04-11, MIT license)\n"
    "#   read-only mirror: ref-code/EffectiveSource (PROVENANCE.elliptic-ssf.md)\n"
    "# VARIANT: " VARIANT " = " VARIANT_SRC
#ifdef EFFSOURCE_VARIANT_COMOVING
    " (3P half-integer Legendre / 2F1 m-mode path;\n"
    "#   hypergeometric branch switch at z = 0.995, kerr-circular-comoving.c:27)\n"
#else
    " (complete-elliptic-integral E/K m-mode path)\n"
#endif
    "# TOOLCHAIN: gcc %s | GSL %s\n"
    "# COMPILE: cc -O2 -std=c99 $(pkg-config --cflags gsl) [-DEFFSOURCE_VARIANT_COMOVING]\n"
    "#   -I$BUILD -o capture_" VARIANT " capture.c " VARIANT_SRC " $(pkg-config --libs gsl) -lm\n"
    "# RUN: ./capture_" VARIANT " <outdir>   (driver source: generator/capture.c)\n"
    "#\n"
    "# FROZEN GRID: 6 configs (a, r0) = (0,10), (0.5,10), (0.9,isco+0.1), (0.998,isco),\n"
    "#   (-0.5,10), (-0.998,isco) [M=1; r0/E/Lz = solid Spacetime values, 17 digits];\n"
    "#   particle at (t=0, r0, pi/2, phi_p=0.31), u^r=0; per-config points npts =\n"
    "#   [40,40,38,37,40,40] in blocks r-ray-minus | r-ray-plus | theta-ray | diagonal\n"
    "#   (|dr| down to 1e-6 M, out to +10 M; minus block clipped above the horizon;\n"
    "#   both sides of the comoving z=0.995 switch covered in dr AND dtheta);\n"
    "#   m in {0,1,2,5,10,20} (comoving PhiS_m adds {-2,-1}); dphi in {0, 0.5, 2.0}.\n"
    "#\n"
    "# PHASE CONVENTIONS (the phase bridge, measured by test_effsource_fixtures.jl):\n"
    "#   * PhiS_m / calc_m (BOTH variants): complex quantities carry exp(-i m phi_p) —\n"
    "#     pairs are (v cos(m phi_p), -v sin(m phi_p)) for real-type v [PhiS, d/dr,\n"
    "#     d/dtheta, d2/dt2, d2/dtdphi, d2/dr2, d2/dtheta2, d2/dphi2, src] and\n"
    "#     (v sin(m phi_p), v cos(m phi_p)) for pure-imaginary-type v [d/dt, d/dphi].\n"
    "#   * calc_m_series (comoving only) uses the OPPOSITE convention exp(+i m phi_p):\n"
    "#     src pairs are (s cos(m phi_p), +s sin(m phi_p)) — kerr-circular-comoving.c:673-674\n"
    "#     vs :634-635. Measured bridge numbers: see effsource_configs.jl header.\n"
    "#\n"
    "# NaN SENTINELS: calc_m mixed second derivatives d2/dtdr, d2/dtdtheta, d2/drdtheta,\n"
    "#   d2/drdphi, d2/dthetadphi are NOT provided by the C code (kerr-circular.c:709-713,\n"
    "#   kerr-circular-comoving.c:570-574) and are recorded as literal NaN.\n"
    "#\n"
    "# COLUMNS (fixture row k = matrix column M[:, k]):\n"
    "#   %s\n",
    title, __VERSION__, GSL_VERSION, cols);
}

/* ---------------- grid construction -------------------------------------------------- */

static int build_points(int c, double pts[MAXPTS][2])
{
  int n = 0, i;
  for (i = 0; i < NMINUS[c]; i++) { pts[n][0] = DR_MINUS[i]; pts[n][1] = 0.0; n++; }
  for (i = 0; i < 11; i++)        { pts[n][0] = DR_PLUS[i];  pts[n][1] = 0.0; n++; }
  for (i = 0; i < 12; i++)        { pts[n][0] = 0.0;         pts[n][1] = DTH[i]; n++; }
  for (i = 0; i < 7; i++)         { pts[n][0] = DIAG[i][0];  pts[n][1] = DIAG[i][1]; n++; }
  return n;
}

static FILE *openf(const char *outdir, const char *name)
{
  char path[4096];
  snprintf(path, sizeof path, "%s/%s", outdir, name);
  FILE *f = fopen(path, "w");
  if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
  return f;
}

int main(int argc, char *argv[])
{
  if (argc != 2) { fprintf(stderr, "usage: %s <outdir>\n", argv[0]); return 1; }
  const char *outdir = argv[1];

  /* GSL default handler aborts on roundoff-level errors near the worldline; upstream
   * test/effsource.c disables it the same way. */
  gsl_set_error_handler_off();

  FILE *f4d = openf(outdir, "phis4d_" VARIANT ".jl");
  header(f4d, "phis4d_" VARIANT ".jl",
    "1 cfg | 2 ipt | 3 dr | 4 dtheta | 5 dphi | 6 PhiS  [effsource_PhiS, 4D real]");
  fputs("const EFFSOURCE_PHIS4D_" VARIANT_UC " = reshape(Float64[\n", f4d);

  FILE *fpm = openf(outdir, "phis_m_" VARIANT ".jl");
  header(fpm, "phis_m_" VARIANT ".jl",
    "1 cfg | 2 m | 3 ipt | 4 dr | 5 dtheta | 6 Re Psi^P_m | 7 Im Psi^P_m  [effsource_PhiS_m]");
  fputs("const EFFSOURCE_PHIS_M_" VARIANT_UC " = reshape(Float64[\n", fpm);

  FILE *fcm = openf(outdir, "calc_m_" VARIANT ".jl");
  header(fcm, "calc_m_" VARIANT ".jl",
    "1 cfg | 2 m | 3 ipt | 4 dr | 5 dtheta | 6-7 PhiS re,im | 8-15 dPhiS (t,r,theta,phi)x(re,im)\n"
    "#   | 16-35 d2PhiS (tt,tr,ttheta,tphi,rr,rtheta,rphi,thetatheta,thetaphi,phiphi)x(re,im)\n"
    "#   | 36-37 src re,im  [effsource_calc_m; raw C array order]");
  fputs("const EFFSOURCE_CALC_M_" VARIANT_UC " = reshape(Float64[\n", fcm);

#ifdef EFFSOURCE_VARIANT_COMOVING
  FILE *fse = openf(outdir, "calc_m_series_comoving.jl");
  header(fse, "calc_m_series_comoving.jl",
    "1 cfg | 2 m | 3 ipt | 4 dr | 5 dtheta | 6 Re src | 7 Im src  [effsource_calc_m_series;\n"
    "#   NOTE the +sin(m phi_p) convention — see PHASE CONVENTIONS above]");
  fputs("const EFFSOURCE_CALC_M_SERIES_COMOVING = reshape(Float64[\n", fse);
#else
  /* configs + point table are emitted once, by the elliptic build (identical content
   * would be produced by either variant — the constants are driver-side). */
  FILE *fcfg = openf(outdir, "effsource_configs.jl");
  header(fcfg, "effsource_configs.jl",
    "EFFSOURCE_CONFIGS: (idx, a, r0, E, Lz, risco, rp, alpha20, alpha02, beta, npts);\n"
    "#   EFFSOURCE_POINTS: 1 cfg | 2 ipt | 3 dr | 4 dtheta.  alpha20 = r0^2/Delta,\n"
    "#   alpha02 = r0^2, beta = z_c^2 = 4(Lz^2 + r0^2 + a^2(r0+2M)/r0) (kerr-circular.c:844-847;\n"
    "#   main.tex:445-455). MEASURED PHASE BRIDGE + conditioning (this committed capture,\n"
    "#   gated in test_effsource_fixtures.jl):\n"
    "#   * PhiS_m elliptic (E/K) vs comoving (3P Legendre): SAME exp(-i m phi_p) convention,\n"
    "#     worst relative difference 2.8e-10 over the overlap window C1 = alpha/beta in\n"
    "#     [1e-8, 1e-2) (test gate 3e-9). Outside the window the two paths are COMPLEMENTARY\n"
    "#     oracles: the elliptic E/K path loses digits near the worldline (C1 < 1e-8;\n"
    "#     k' = sqrt(C1) cancellation entering Kcomp — up to ~1.2e-3 relative at C1 ~ 1e-15,\n"
    "#     values stay finite), and in the far field (C1 >~ 0.1) it is the ELLIPTIC\n"
    "#     E/K-polynomial route (kerr-circular.c) that cancels catastrophically for\n"
    "#     m >= 2 — BigFloat(256)-quadrature adjudication (W4-D): elliptic-vs-truth up\n"
    "#     to ~5.8e91 relative at C1 ~ 120, m = 20, while the comoving/3P route tracks\n"
    "#     truth to ~1e-11 over the whole frozen grid. An earlier revision of this\n"
    "#     header mis-attributed the far-field cancellation to the comoving route.\n"
    "#     Evidence: results/elliptic-ssf/paper_ELLIPTIC-SSF/evidence/\n"
    "#     mmode-puncture-integrals/verification_summary.json.\n"
    "#   * calc_m DIRECT src vs calc_m_series (comoving): the direct Box[Phi^P] evaluation\n"
    "#     is DESTROYED by cancellation near the worldline (e.g. |src_direct| ~ 1.4e2 vs\n"
    "#     true ~ 7.0e-5 at dr = 1e-6, cfg 1 m 1) — calc_m_series is the near-worldline\n"
    "#     oracle; the direct src is trustworthy only for rho^2 >~ 1e-5. Handoff agreement:\n"
    "#     for EVERY (config, m) the minimum-over-points relative |src| difference is\n"
    "#     <= 6.0e-4 (worst at m = 20; test gate 5e-3).\n"
    "#   * calc_m_series phase is EXACTLY opposite: series = s exp(+i m phi_p) vs\n"
    "#     calc_m = s' exp(-i m phi_p) — the +sin(m phi_p) bridge; all phase-rotation\n"
    "#     residuals <= 1.3e-16 (test gates 1e-12).\n"
    "#   * conjugate symmetry (comoving): |Psi_{-m} - conj(Psi_m)|/|Psi_m| <= 4.1e-15 for\n"
    "#     m = 1, 2 (test gate 1e-12).");
  fputs("const EFFSOURCE_PHIP = 0.31\n\nconst EFFSOURCE_CONFIGS = [\n", fcfg);
  for (int c = 0; c < 6; c++) {
    double pts[MAXPTS][2];
    int npts = build_points(c, pts);
    double a = CFG[c].a, r0 = CFG[c].r0, L = CFG[c].L;
    double alpha20 = r0 * r0 / (a * a + r0 * (r0 - 2.0));
    double alpha02 = r0 * r0;
    double beta = 4.0 * (L * L + r0 * r0 + a * a * (r0 + 2.0) / r0);
    fprintf(fcfg, "    (idx = %d, a = ", c + 1); jd(fcfg, a);
    fputs(", r0 = ", fcfg); jd(fcfg, r0);
    fputs(", E = ", fcfg); jd(fcfg, CFG[c].E);
    fputs(", Lz = ", fcfg); jd(fcfg, L);
    fputs(", risco = ", fcfg); jd(fcfg, CFG[c].risco);
    fputs(", rp = ", fcfg); jd(fcfg, CFG[c].rp);
    fputs(", alpha20 = ", fcfg); jd(fcfg, alpha20);
    fputs(", alpha02 = ", fcfg); jd(fcfg, alpha02);
    fputs(", beta = ", fcfg); jd(fcfg, beta);
    fprintf(fcfg, ", npts = %d),\n", npts);
  }
  fputs("]\n\nconst EFFSOURCE_POINTS = reshape(Float64[\n", fcfg);
  for (int c = 0; c < 6; c++) {
    double pts[MAXPTS][2];
    int npts = build_points(c, pts);
    for (int p = 0; p < npts; p++) {
      row_begin(fcfg);
      cell(fcfg, c + 1); cell(fcfg, p + 1); cell(fcfg, pts[p][0]);
      row_end(fcfg, pts[p][1]);
    }
  }
  fputs("], 4, :)\n", fcfg);
  fclose(fcfg);
#endif

  for (int c = 0; c < 6; c++) {
    double pts[MAXPTS][2];
    int npts = build_points(c, pts);

    effsource_init(1.0, CFG[c].a);
    struct coordinate xp = {0.0, CFG[c].r0, M_PI_2, PHIP};
    effsource_set_particle(&xp, CFG[c].E, CFG[c].L, 0.0);

    for (int p = 0; p < npts; p++) {
      struct coordinate x = {0.0, CFG[c].r0 + pts[p][0], M_PI_2 + pts[p][1], PHIP};

      /* 4D singular field at three dphi offsets */
      for (int q = 0; q < 3; q++) {
        double phis;
        x.phi = PHIP + DPHI[q];
        effsource_PhiS(&x, &phis);
        row_begin(f4d);
        cell(f4d, c + 1); cell(f4d, p + 1); cell(f4d, pts[p][0]); cell(f4d, pts[p][1]);
        cell(f4d, DPHI[q]);
        row_end(f4d, phis);
      }
      x.phi = PHIP; /* not used by the m-mode entry points */

      /* standalone m-mode puncture */
      for (int q = 0; q < N_MS_PHIS_M; q++) {
        double P2[2];
        effsource_PhiS_m(MS_PHIS_M[q], &x, P2);
        row_begin(fpm);
        cell(fpm, c + 1); cell(fpm, MS_PHIS_M[q]); cell(fpm, p + 1);
        cell(fpm, pts[p][0]); cell(fpm, pts[p][1]); cell(fpm, P2[0]);
        row_end(fpm, P2[1]);
      }

      /* m-mode puncture + derivatives + effective source */
      for (int q = 0; q < 6; q++) {
        double P2[2], d1[8], d2[20], s2[2];
        effsource_calc_m(MS[q], &x, P2, d1, d2, s2);
        row_begin(fcm);
        cell(fcm, c + 1); cell(fcm, MS[q]); cell(fcm, p + 1);
        cell(fcm, pts[p][0]); cell(fcm, pts[p][1]);
        cell(fcm, P2[0]); cell(fcm, P2[1]);
        for (int j = 0; j < 8; j++) cell(fcm, d1[j]);
        for (int j = 0; j < 20; j++) cell(fcm, d2[j]);
        cell(fcm, s2[0]);
        row_end(fcm, s2[1]);
      }

#ifdef EFFSOURCE_VARIANT_COMOVING
      /* near-worldline series form of the effective source */
      for (int q = 0; q < 6; q++) {
        double s2[2];
        effsource_calc_m_series(MS[q], &x, s2);
        row_begin(fse);
        cell(fse, c + 1); cell(fse, MS[q]); cell(fse, p + 1);
        cell(fse, pts[p][0]); cell(fse, pts[p][1]); cell(fse, s2[0]);
        row_end(fse, s2[1]);
      }
#endif
    }
  }

  fputs("], 6, :)\n", f4d); fclose(f4d);
  fputs("], 7, :)\n", fpm); fclose(fpm);
  fputs("], 37, :)\n", fcm); fclose(fcm);
#ifdef EFFSOURCE_VARIANT_COMOVING
  fputs("], 7, :)\n", fse); fclose(fse);
#endif
  return 0;
}
