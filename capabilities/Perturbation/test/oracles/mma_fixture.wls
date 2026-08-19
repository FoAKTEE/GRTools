(* TRACKED ORACLE SCRIPT  (GRR::g14-oracle-scripts, repro-audit gap G14)
   
   origin : emri-disk/scratchpad/qqnm_phase3/mma_fixture.wls (untracked one-shot)
   sha256 : 8dbeba9cd820af29fa21ef833787da78ae849fda4f7e34b682554c953f3338c8
   backs  : test/fixtures/qnm_hphc_mma_ref.jl (via scratchpad/qqnm_phase3/mma_radial_fixture.txt)
   engine : wolframscript + QNM-Mathematica HEAD a74aadd, WL 14.3.0
   
   Tracked for provenance/regeneration only: nothing in the Julia test suite executes it
   (no Wolfram/py dependency is added).  Re-running it requires the engine above. *)

(* Q2H-PREP-B task 2: export the validated QNM.wl (s=-2,l=2,m=2,n=0,a=0.7) radial
   eigenfunction on its native hyperboloidal grid + eigenvalue, to a plain-text fixture.
   READ-ONLY: runs against a byte-for-byte scratch copy of ref-code/QNM-Mathematica. *)
PacletDirectoryLoad["/data/haiyangw/claude/emri-disk/scratchpad/qqnm_phase3/qnm_scratch/QNM-Mathematica"];
Needs["QNM`"];
Off[QNMFrequency::acc]; Off[QNMFrequency::noacc];

s = -2; l = 2; m = 2; n = 0; a = 0.7;

(* ---- frequency: default (Automatic = Ripley hyperboloidal, 32-pt) + cross-checks ---- *)
wAuto = QNMFrequency[s,l,m,n,N[a,30]];                 (* default method *)
wLeav = QNMFrequency[s,l,m,n,N[a,30], Method->"Leaver"];
wInt  = Quiet@QNMFrequency[s,l,m,n,a,   Method->"Interpolation"];
w = wAuto;
lam = SpinWeightedSpheroidalEigenvalue[s,l,m,a w];

(* ---- radial eigenfunction on its native grid: replicate QNMRadialHyperboloidal(:701) ---- *)
numpoints = 100;                                        (* QNMRadialHyperboloidal default *)
rp = QNM`Private`rp[a, 1];
rho = (1/rp) Reverse[1/2 (1 + Cos[Pi Subdivide[numpoints-1]])];   (* ScriptCapitalM grid, asc 0->1/rp *)
Mat = QNM`Private`\[ScriptCapitalM][s, m, a, w, lam, numpoints];
ns  = NullSpace[Mat];
ef  = If[Length[ns]==0, Eigensystem[Mat][[2, -1]], First[ns]];
efN = ef/ef[[-1]];                                      (* normalize: horizon (rho=1/rp) amplitude = 1 *)
ampScri = ef[[1]]/ef[[-1]];                             (* I = scri amplitude, matches QNMRadialFunction *)

(* ---- public cross-check via QNMRadial ---- *)
qr = QNMRadial[s,l,m,n,a, "Frequency"->w];
lamPub = qr["Eigenvalue"]; ampPub = qr["Amplitudes"];

(* round-trippable, human-readable complex formatter (machine precision, ~16-17 sig) *)
cs[z_] := With[{zz=N[z]},
  ToString[Re[zz],InputForm] <> (If[Im[zz]>=0,"+","-"]) <> ToString[Abs[Im[zz]],InputForm] <> "im"];
rs[x_] := ToString[N[x], InputForm];

f = OpenWrite["/data/haiyangw/claude/emri-disk/scratchpad/qqnm_phase3/mma_radial_fixture.txt"];
WriteLine[f, "# ==================================================================="];
WriteLine[f, "# MMA radial-eigenfunction fixture  (Q2H-PREP-B task 2)"];
WriteLine[f, "# Generated: " <> DateString[]];
WriteLine[f, "# Source pkg: ref-code/QNM-Mathematica  HEAD a74aadd  (scratch byte-copy)"];
WriteLine[f, "#   Kernel file used: " <> FindFile["QNM`"]];
WriteLine[f, "#   " <> $Version];
WriteLine[f, "# Deps: SpinWeightedSpheroidalHarmonics, Teukolsky (BHPToolkit)"];
WriteLine[f, "# Mode: s=-2, l=2, m=2, n=0, a=0.7, M=1  (fundamental gravitational QNM)"];
WriteLine[f, "# Method: QNMFrequency default = Ripley hyperboloidal eigenvalue solve, 32 pts;"];
WriteLine[f, "#         radial eigenfunction = QNMRadialHyperboloidal replica, numpoints=100."];
WriteLine[f, "# CHART: hyperboloidal Ripley (Ripley:2022ypi). Compactified R=rho=1/r in [0,1/rp];"];
WriteLine[f, "#   rho=0 = null infinity (scri), rho=1/rp = horizon. Operator ScriptCapitalM"];
WriteLine[f, "#   (QNM.wl:97-107): -rho^2*Delta*d2R + A*dR + (B+lam+2 a m w-(a w)^2)*R, Delta=1-2M rho+a^2 rho^2."];
WriteLine[f, "#   Eigenfunction normalized so horizon value ef[[-1]] = 1 (I=scri amplitude below)."];
WriteLine[f, "# --- eigenvalues ---"];
WriteLine[f, "# omega_default_32pt_hyperboloidal = " <> cs[wAuto]];
WriteLine[f, "# omega_Leaver_crosscheck          = " <> cs[wLeav]];
WriteLine[f, "# omega_h5Interp_crosscheck        = " <> cs[wInt]];
WriteLine[f, "# |omega_auto - omega_Leaver|      = " <> ToString[N@Abs[wAuto-wLeav], InputForm]];
WriteLine[f, "# angular_separation_const_lambda  = " <> cs[lam] <> "   (SpinWeightedSpheroidalEigenvalue[-2,2,2,a*omega])"];
WriteLine[f, "# QNMRadial_public_Eigenvalue      = " <> cs[lamPub]];
WriteLine[f, "# scri_amplitude_I = ef[[1]]/ef[[-1]] = " <> cs[ampScri]];
WriteLine[f, "# QNMRadial_public_Amplitudes: H=" <> cs[ampPub["\[ScriptCapitalH]"]] <> "  I=" <> cs[ampPub["\[ScriptCapitalI]"]]];
WriteLine[f, "# rp = r_+ = " <> rs[rp] <> " ;  1/rp = " <> rs[1/rp]];
WriteLine[f, "# ==================================================================="];
WriteLine[f, "# columns:  idx    rho(=1/r)    r(=1/rho)    Re(ef_norm)    Im(ef_norm)"];
WriteLine[f, "# row 1: rho=0 => r=Infinity (scri).  last row: rho=1/rp => horizon, ef_norm=1+0im"];
Do[
  With[{rr = If[rho[[i]]==0, "Infinity", rs[1/rho[[i]]]]},
    WriteLine[f, StringRiffle[{
        ToString[i], rs[rho[[i]]], rr, rs[Re[efN[[i]]]], rs[Im[efN[[i]]]] }, "    "]]],
  {i, 1, numpoints}
];
Close[f];
Print["omega_auto = ", cs[wAuto]];
Print["omega_leav = ", cs[wLeav]];
Print["|wAuto-wLeav| = ", ScientificForm[N@Abs[wAuto-wLeav]]];
Print["lambda = ", cs[lam], "   lamPub = ", cs[lamPub]];
Print["scri amp I = ", cs[ampScri], "   pub I = ", cs[ampPub["\[ScriptCapitalI]"]]];
Print["wrote fixture with ", numpoints, " grid rows"];
Print["done"];
