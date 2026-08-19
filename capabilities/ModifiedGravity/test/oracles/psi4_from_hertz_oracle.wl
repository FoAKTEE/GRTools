(* TRACKED ORACLE SCRIPT  (GRR::g14-oracle-scripts, repro-audit gap G14)
   
   origin : emri-disk/packages/ModifiedTeukolsky/oracle/psi4_from_hertz_oracle.wl @ 9be1f980
   sha256 : f1d9bceb2d485bfd841471152b3354f7ebe0ede60fdf3c344fc39a14f86ae35f
   backs  : test/test_weyl_from_hertz.jl (necessary/asymptotic branch)
   engine : wolframscript (Wolfram Language)
   
   Tracked for provenance/regeneration only: nothing in the Julia test suite executes it
   (no Wolfram/py dependency is added).  Re-running it requires the engine above. *)

(* Oracle + generator for psi4_from_hertz: Radpsi4fromHzpotmode (Sourceincoordinates.nb).
   The reconstructed psi4 closed form is built xAct-FREE from Psi4check (Doperdag^4, SC l.408)
   reduced via DerHzpotsubstf (radial-Teukolsky derivative elimination, SC l.356), times
   signcoef[2,mi] (= Psi4fromHzpotlmomega, proven equal at SC IN176). It is grounded by the
   reconstruction-vs-Teukolsky identity (Shah & van de Meent 2015): Radpsi4fromHzpotmode equals
   the point-particle psi4 amplitude / r^4 from the BHPToolkit Teukolsky paclet.

   Prints (a) the reduced symbolic kernel in InputForm (for the Julia port) and
   (b) numeric reconstructed psi4 vs toolkit psi4/r^4 at sample radii. *)

SetDirectory[DirectoryName[$InputFileName]];
ToExpression[Import["/tmp/claude-1012/-data-haiyangw-claude-emri-disk/a271c4fb-a0c6-4324-b1ef-aec2ea6ab28d/scratchpad/derhz.txt", "Text"]];  (* binds DerHzpotsubstf via side effect *)
Doperdag[func_] := D[func, r[]] - (r[]^2/(r[]^2 - 2*MB*r[]))*D[func, t[]];
Psi4check = (1/32)*(1/r[]^4)*(r[]^2 - 2*MB*r[])^2*(Doperdag[Doperdag[Doperdag[Doperdag[(r[]^2 - 2*MB*r[])^2*RadHZ\[Dagger][li, mi, omega, r[]]*E^(I*omega*t[])]]]]/E^(I*omega*t[]));
signcoef[s_, m_] := (-1)^(s + m);
psi4lm = Simplify[signcoef[2, mi]*(Psi4check /. DerHzpotsubstf)];
psi4lmrad = Simplify[psi4lm /. r[] -> rad];

Print["=== psi4lm InputForm (rhzb=RadHZdag, drhzb=RadHZdag') ==="];
out = psi4lmrad /. {RadHZ\[Dagger][li, mi, omega, rad] -> rhzb,
                    Derivative[0, 0, 0, 1][RadHZ\[Dagger]][li, mi, omega, rad] -> drhzb};
Print[ToString[InputForm[out]]];

(* ---- numeric grounding vs BHPToolkit Teukolsky ---- *)
Needs["Teukolsky`"];
Needs["KerrGeodesics`"];
RadHzpotmode[l_, m_, n_, k_, orbit_] := Module[{kloc, T, freq, amp, Radfunc},
  If[orbit["Inclination"] == 1, kloc = 0, kloc = k];
  T = TeukolskyPointParticleMode[-2, l, m, orbit];
  freq = T["\[Omega]"]; amp = T["Amplitudes"][[1]];
  amp = (-1)^(l + kloc + m)*(2/freq^4)*amp;
  Radfunc = TeukolskyRadial[2, l, m, 0, freq];
  {Function[rr, amp*Radfunc["Up"][rr]],
   Function[rr, (D[amp*Radfunc["Up"][rad], rad] /. rad -> rr)], freq}];

reconpsi4[l_, m_, n_, k_, orbit_, rr_] := Module[{rhz, drhz, freq, expr},
  {rhz, drhz, freq} = RadHzpotmode[l, m, n, k, orbit];
  expr = psi4lmrad /. {MB -> 1, omega -> freq, li -> l, mi -> m,
     RadHZ\[Dagger][li, mi, omega, rad] -> Conjugate[rhz[rad]],
     Derivative[0, 0, 0, 1][RadHZ\[Dagger]][li, mi, omega, rad] -> Conjugate[drhz[rad]]};
  expr /. rad -> rr];

aa = 0.; pp = 10.; ee = 0.; xx = 1; ll = 2; mm = 2; nn = 0; kk = 0;
orbit = KerrGeoOrbit[aa, pp, ee, xx];
psitk = TeukolskyPointParticleMode[-2, ll, mm, orbit];
Print["=== identity: Radpsi4fromHzpotmode[l,-m,-n,-k] vs psi4toolkit[l,m,n,k]/r^4 ==="];
Do[
  rec = reconpsi4[ll, -mm, -nn, -kk, orbit, rv];
  tk = psitk[rv]/rv^4;
  Print["r=", rv, " recon=", ToString[N[rec, 12], InputForm],
        " toolkit/r^4=", ToString[N[tk, 12], InputForm],
        " reldiff=", N[Abs[rec - tk]/Abs[tk]]];
  , {rv, {20., 50., 100.}}];
