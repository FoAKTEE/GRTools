(* Oracle for h44_from_hertz (Radh44fromHzpotmode, Sourceincoordinates.nb SC l.368/430/432).
   The reconstructed (4,4)=m̄m̄ metric mode is assembled from the recovered operator
     h44metric = -signcoef[2,m] SWSH[-2,l,-m] ( (3ε'+3ε'†+2ρ') þ'Ψ + Ψ B + þ'þ'Ψ )
   with Ψ = Radf†[2,l,m] = RadHZ†/ksim1^4 = RadHZ† r^4 (FreqsubstHz, ksim1=1/r), signcoef[2,m]=(-1)^m,
     B = 2ε'^2+4ε'ε'†+2ε'†^2+ε'ρ'+ε'†ρ'-3ρ'^2+2þ'ε'+2þ'ε'†+3þ'ρ',
   and þ'≡PDNP[{2,-NP}] the GHP thorn-prime, expanded on the equatorial Schwarzschild Kinnersley
   tetrad as þ'F = (1/2)∂_t F - (Δ/2r²)∂_r F + (p+q)ε' F (weight-sum p+q). Schwspincoeff: ε'=ε'†=
   -M/2r², ρ'=(r-2M)/2r². This script composes the operator SYMBOLICALLY (D[] chain/product rule,
   independent of the Julia jet-tower implementation) and evaluates it numerically with the BHPToolkit
   Teukolsky paclet RadHZ† (RadHzpotmode Rup radial), the SWSH divided out (radial mode). *)

Needs["Teukolsky`"]; Needs["KerrGeodesics`"];
MB = 1;
Δ = r[]^2 - 2 MB r[];
esp = -MB/(2 r[]^2); rsp = (r[] - 2 MB)/(2 r[]^2);
tp[F_, w_] := (1/2) D[F, t[]] - (Δ/(2 r[]^2)) D[F, r[]] + w esp F;    (* GHP þ' with weight-sum w *)
bpesp = tp[esp, 0];    (* þ'ε' : ε' charge-free (w=0) *)
bprsp = tp[rsp, -2];   (* þ'ρ' : ρ' type {-1,-1} (w=-2) *)
signcoef[s_, m_] := (-1)^(s + m);
Aco = 3 esp + 3 esp + 2 rsp;
Bco = 2 esp^2 + 4 esp esp + 2 esp^2 + esp rsp + esp rsp - 3 rsp^2 + 2 bpesp + 2 bpesp + 3 bprsp;
Psi = RadHZd[r[]] r[]^4 E^(I omega t[]);   (* Radf†[2,l,m]; ∂_t -> Iω *)
tpPsi = tp[Psi, 4];                        (* Ψ weight-sum 4 *)
tp2Psi = tp[tpPsi, 2];                     (* þ'Ψ weight-sum 2 *)
h44mode = -signcoef[2, mi] (Aco tpPsi + Psi Bco + tp2Psi)/E^(I omega t[]);  (* SWSH-free radial mode *)
h44mode = h44mode /. t[] -> 0;
(* expand derivatives of the field into RadHZd, RadHZd', RadHZd'' at r[] *)
h44mode = h44mode /. { Derivative[2][RadHZd][r[]] -> rh2, Derivative[1][RadHZd][r[]] -> rh1,
                       RadHZd[r[]] -> rh0 };

(* ---- numeric: paclet Hertz Rup radial + its r-derivatives, daggered ---- *)
RadHzpotmode[l_, m_, n_, k_, orbit_] := Module[{kloc, T, freq, amp, Rup},
  If[orbit["Inclination"] == 1, kloc = 0, kloc = k];
  T = TeukolskyPointParticleMode[-2, l, m, orbit];
  freq = T["\[Omega]"]; amp = T["Amplitudes"][[1]];
  amp = (-1)^(l + kloc + m)*(2/freq^4)*amp;
  Rup = TeukolskyRadial[2, l, m, 0, freq]["Up"];
  {amp, Rup, freq}];

aa = 0.; pp = 10.; ee = 0.; xx = 1; ll = 2; mm = -2; nn = 0; kk = 0;
orbit = KerrGeoOrbit[aa, pp, ee, xx];
{amp, Rup, freq} = RadHzpotmode[ll, mm, nn, kk, orbit];
RH[rr_?NumericQ] := amp*Rup[rr];
dRH[rr_?NumericQ] := (D[amp*Rup[x], x] /. x -> rr);
d2RH[rr_?NumericQ] := (D[amp*Rup[x], {x, 2}] /. x -> rr);
(* RadHZ†(r) = Conjugate[RH(r)] for real r; derivs conjugate too *)
h44val[rr_] := h44mode /. { rh0 -> Conjugate[RH[rr]], rh1 -> Conjugate[dRH[rr]],
   rh2 -> Conjugate[d2RH[rr]], omega -> freq, mi -> mm, r[] -> rr };

Print["ω = ", ToString[N[freq, 16], InputForm]];
Do[
  h = h44val[rv];
  psi4amp = TeukolskyPointParticleMode[-2, ll, -mm, orbit][rv]/rv^4;  (* for asymptotic self-check *)
  Print["r=", rv, "  h44=", ToString[N[h, 14], InputForm]];
  , {rv, {12., 20., 50., 100., 40000.}}];
