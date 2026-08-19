(* TRACKED ORACLE SCRIPT  (GRR::g14-oracle-scripts, repro-audit gap G14)
   
   origin : emri-disk/packages/ModifiedTeukolsky/oracle/ring_functions_oracle.wl @ 9be1f980
   sha256 : 39368fa491b5b5717ed9eceeebdcc6cb10d1a7b537738871eb29737832ca23e1
   backs  : test/test_ring_dynamics.jl (211-mode ring observables, density/continuity literals)
   engine : wolframscript (Wolfram Language)
   
   Tracked for provenance/regeneration only: nothing in the Julia test suite executes it
   (no Wolfram/py dependency is added).  Re-running it requires the engine above. *)

(* Oracle for ring_functions: ResultingfunctionsRing (Sourceincoordinates.nb SC l.1134),
   composed from Urhoamplitudes + BLtoNPruleU (closed forms transcribed VERBATIM from the
   notebook). Pure closed-form algebra (no xAct). Imports a shipped .dat source table, builds
   the (t,phi) reconstruction Sum e^{-i(omega t - m phi)} + c.c., and prints the deformed-ring
   observables at a sample (tc,phic). Args: <datfile> <radius> <tc> <phic>. *)

Dcommon = 12*MB^2 + MB*(-4 + mi^2)*R - 2*Sqrt[MB]*mi*omega*R^(5/2) + omega^2*R^4;
PG = {
 fourvelocitycoef[1] -> (Sqrt[-3*MB + R]*((8*I)*MB^3*S1 - (2*I)*omega^2*R^5*S2 + I*Sqrt[MB]*omega*R^(7/2)*(4*mi*S2 + Sqrt[2]*(S3 - S4)) + 2*MB^(3/2)*R^(3/2)*(2*mi*S2 + Sqrt[2]*(1 - I*omega*R)*(S3 - S4)) + 4*Sqrt[2]*MB^(5/2)*Sqrt[R]*(-S3 + S4) - (2*I)*MB^2*R*(4*S1 + 4*S2 + Sqrt[2]*mi*(-S3 + S4)) + I*MB*R^2*(2*S1 + (4 - 2*mi^2 + (4*I)*omega*R)*S2 + Sqrt[2]*mi*(-S3 + S4))))/(2*(Sqrt[MB]*mi - omega*R^(3/2))*Dcommon),
 fourvelocitycoef[2] -> (R*Sqrt[-3*MB + R]*((-4*I)*MB^2*S1 - I*omega^2*R^4*S1 + MB*R*((2*I - I*mi^2 + 2*omega*R)*S1 + I*(4*S2 + Sqrt[2]*mi*(S3 - S4))) - 2*MB^(3/2)*Sqrt[R]*(mi*S1 + Sqrt[2]*(-S3 + S4)) + I*Sqrt[MB]*omega*R^(5/2)*(2*mi*S1 + Sqrt[2]*(-S3 + S4))))/((Sqrt[MB]*mi - omega*R^(3/2))*Dcommon),
 fourvelocitycoef[3] -> (Sqrt[-3*MB + R]*(-4*MB^(5/2)*Sqrt[R]*S1 + 2*MB^(3/2)*R^(3/2)*(S1 + I*omega*R*S1 + 2*S2) + I*Sqrt[2]*omega^2*R^5*S4 - I*Sqrt[MB]*omega*R^(7/2)*(S1 - 2*S2 + 2*Sqrt[2]*mi*S4) + I*MB*R^2*(mi*(S1 - 2*S2) + Sqrt[2]*mi^2*S4 - 2*Sqrt[2]*(S3 + S4)) - (2*I)*MB^2*R*(mi*S1 - 2*Sqrt[2]*(S3 + 2*S4))))/(Sqrt[2]*(Sqrt[MB]*mi - omega*R^(3/2))*Dcommon),
 fourvelocitycoef[4] -> (Sqrt[-3*MB + R]*(4*Sqrt[2]*MB^(5/2)*Sqrt[R]*S1 + 2*Sqrt[2]*MB^(3/2)*R^(3/2)*((-1 - I*omega*R)*S1 - 2*S2) + (2*I)*omega^2*R^5*S3 + I*Sqrt[MB]*omega*R^(7/2)*(Sqrt[2]*S1 - 2*Sqrt[2]*S2 - 4*mi*S3) + (2*I)*MB^2*R*(Sqrt[2]*mi*S1 + 8*S3 + 4*S4) - I*MB*R^2*(Sqrt[2]*mi*(S1 - 2*S2) - 2*mi^2*S3 + 4*(S3 + S4))))/(2*(Sqrt[MB]*mi - omega*R^(3/2))*Dcommon)};
normAfter = ((2*I)*MB*Sqrt[R]*S1 - I*R^(3/2)*(S1 + 2*S2) + Sqrt[2]*Sqrt[MB]*R*(-S3 + S4))/(-(Sqrt[MB]*mi) + omega*R^(3/2));
densRule = {
 drhoA -> (-4^(-1))*((md*(-3*MB + R)^(3/2)*(2*R*((-4*I)*MB + R*(2*I + omega*R))*fourvelocitycoef[1] - I*(2*MB - R)*((4*MB + R*(-2 - I*omega*R))*fourvelocitycoef[2] + Sqrt[2]*mi*R*(-fourvelocitycoef[3] + fourvelocitycoef[4]))))/(Pi*R^3*(-2*MB + R)^2*(-(Sqrt[MB]*mi) + omega*R^(3/2)))),
 dthetaA -> ((-I)*Sqrt[-3*MB + R]*(fourvelocitycoef[3] + fourvelocitycoef[4]))/(Sqrt[2]*(Sqrt[MB]*mi - omega*R^(3/2))),
 drA -> ((-I)*Sqrt[-3*MB + R]*(2*R*fourvelocitycoef[1] + (2*MB - R)*fourvelocitycoef[2]))/(2*Sqrt[MB]*mi - 2*omega*R^(3/2))};

args = $ScriptCommandLine;
datfile = args[[2]];
radius = ToExpression[args[[3]]];
tc = ToExpression[args[[4]]];
phic = ToExpression[args[[5]]];
mass = 10^(-4);
dat = Import[datfile];

Urho[row_] := Module[{vals, u, uv, dens, nc},
  vals = {MB -> 1, R -> radius, omega -> row[[4]], mi -> row[[1]], md -> mass,
          S1 -> row[[5]] + I*row[[6]], S2 -> row[[7]] + I*row[[8]],
          S3 -> row[[9]] + I*row[[10]], S4 -> row[[11]] + I*row[[12]]};
  u = Table[fourvelocitycoef[i] /. PG /. vals, {i, 1, 4}];
  uv = Thread[{fourvelocitycoef[1], fourvelocitycoef[2], fourvelocitycoef[3], fourvelocitycoef[4]} -> u];
  dens = ({drhoA, dthetaA, drA} /. densRule) /. uv /. vals;
  nc = normAfter /. vals;
  <|"m" -> row[[1]], "omega" -> row[[4]], "U" -> u, "dens" -> dens, "normchk" -> nc, "normSrc" -> 0|>];

amps = Urho /@ dat;
ph[a_] := a["omega"]*tc - a["m"]*phic;
sumdir[get_] := Sum[get[a]*Exp[-I*ph[a]], {a, amps}];
sumcc[get_] := Sum[Conjugate[get[a]]*Exp[I*ph[a]], {a, amps}];
U1e = sumdir[#["U"][[1]] &] + sumcc[#["U"][[1]] &];
U2e = sumdir[#["U"][[2]] &] + sumcc[#["U"][[2]] &];
U3e = sumdir[#["U"][[3]] &] + sumcc[#["U"][[4]] &];   (* +c.c. uses U4 *)
U4e = sumdir[#["U"][[4]] &] + sumcc[#["U"][[3]] &];   (* +c.c. uses U3 *)
drhoe = sumdir[#["dens"][[1]] &] + sumcc[#["dens"][[1]] &];
dthetae = sumdir[#["dens"][[2]] &] + sumcc[#["dens"][[2]] &];
dre = sumdir[#["dens"][[3]] &] + sumcc[#["dens"][[3]] &];
Normdife = sumdir[(#["normSrc"] + #["normchk"]) &] + sumcc[(#["normSrc"] + #["normchk"]) &];
MBv = 1;
Ute = -((radius*U1e)/(2*MBv - radius)) + U2e/2;
Ure = U1e + (-1/2 + MBv/radius)*U2e;
Uthetae = (U3e + U4e)/(Sqrt[2]*radius);
Uphie = (I*(U3e - U4e))/(Sqrt[2]*radius);

pr[x_] := ToString[Re[N[x, 20]], InputForm] <> " , " <> ToString[Im[N[x, 20]], InputForm];
Print["U1 ", pr[U1e]];
Print["U3 ", pr[U3e]];
Print["drho ", pr[drhoe]];
Print["dtheta ", pr[dthetae]];
Print["dr ", pr[dre]];
Print["Ut ", pr[Ute]];
Print["Ur ", pr[Ure]];
Print["Utheta ", pr[Uthetae]];
Print["Uphi ", pr[Uphie]];
Print["Normdif ", pr[Normdife]];
