# BAKED MMA cross-oracle fixture for the Q2-H hyperboloidal QNM radial eigenfunction.
# Source: scratchpad/qqnm_phase3/mma_radial_fixture.txt (QNM-Mathematica HEAD a74aadd,
# BHPToolkit Teukolsky+SWSH, WL 14.3.0). Mode (s=-2,l=2,m=2,n=0,a=0.7,M=1). Hyperboloidal
# Ripley chart R=1/r in [0,1/r+]; ef normalized so horizon value = 1; I = ef(scri)/ef(hor).
# Chart map to the GrTool hphc operator: SAME hyperboloidal chart; the GrTool 1/2-form
# operator = 1/2 x the QNM.wl ScriptCapitalM (== R22/2), so the NULL EIGENFUNCTION is
# identical (global scale is a null-space invariant). Comparison is on the scri/horizon
# R-RATIO (normalization-invariant): GrTool R(rho)/R(scri) vs baked ef(rho)/ef(scri).
module QNMHphcMMARef
const OMEGA   = 0.5326002435337337 - 0.08079287315478063im   # 32-pt hyperboloidal solve
const OMEGA_LEAVER = 0.532600243551019 - 0.08079287315501084im  # Leaver crosscheck (|Δ|=1.73e-11)
const LAMBDA  = 1.547682710237845 + 0.367208599709584im     # SWSH[-2,2,2,a*omega] (BHPT)
const I_SCRI  = 0.29252750963763074 + 0.011718684202705693im # scri amplitude ef(scri)/ef(hor)
const RPLUS   = 1.7141428428542849
const RP_COMP = 0.5833819533586021   # 1/r+ = R_+ (compactified horizon)
# (rho, ef) rows: rho=1/r in [0,1/r+], ef horizon-normalized (ef(hor)=1). 32-pt-omega source.
const RHO = Float64[
  0.0,   # idx 1
  0.01759111834392735,   # idx 12
  0.07430406146901833,   # idx 24
  0.16216075273765398,   # idx 36
  0.26857447936195195,   # idx 48
  0.3782999709547902,   # idx 60
  0.4756174999889236,   # idx 72
  0.5465849549473855,   # idx 84
  0.5810352448556491,   # idx 96
  0.5833819533586021,   # idx 100
]
const EF = ComplexF64[
  0.29252750963763074 + 0.011718684202705679im,
  0.30711177320896027 + 0.011854046529727645im,
  0.35627531565761644 + 0.012615370923589865im,
  0.43951725461064645 + 0.014287562189162845im,
  0.5534840458882648 + 0.015913688948061545im,
  0.6885619869265923 + 0.015356917636294825im,
  0.8257774042684607 + 0.01116006648894179im,
  0.9376568712620967 + 0.004759992411751604im,
  0.9959311543167918 + 0.0003368665981787364im,
  1.0 + -5.0551961381894235e-18im,
]
end # module
