(* TRACKED ORACLE SCRIPT  (GRR::g14-oracle-scripts, repro-audit gap G14)
   
   origin : emri-disk/packages/ModifiedTeukolsky/oracle/ring_urho_oracle.wl @ 9be1f980
   sha256 : 65b79df91d527ee4aae8cd84d53c13fd14d27df34657b8f0fd09f1ac666a33b0
   backs  : test/test_ring_dynamics.jl 'four-velocity vs wolframscript Urhoamplitudes' (the test header calls this file scratchpad/ora_urho.wl — that scratchpad copy is LOST; this is the packages/ copy of the same oracle)
   engine : wolframscript (Wolfram Language)
   
   Tracked for provenance/regeneration only: nothing in the Julia test suite executes it
   (no Wolfram/py dependency is added).  Re-running it requires the engine above. *)

(* Oracle for ring_urho_amplitudes: PertGeodesiceqsolution + densitysubstrule +
   LHsnormaftersubst, transcribed VERBATIM from Sourceincoordinates.nb (SC.txt lines
   1074/1078/1100). Closed-form algebra only (no xAct). Evaluated on row 5 of
   i1.0_p12.0_e0_r20.0.dat. Prints u1..u4, drho, dtheta, dr, normcheck. *)

Dcommon = 12*MB^2 + MB*(-4 + mi^2)*R - 2*Sqrt[MB]*mi*omega*R^(5/2) + omega^2*R^4;
PertGeodesiceqsolution = {
 fourvelocitycoef[1] -> (Sqrt[-3*MB + R]*((8*I)*MB^3*S1 - (2*I)*omega^2*R^5*S2 + I*Sqrt[MB]*omega*R^(7/2)*(4*mi*S2 + Sqrt[2]*(S3 - S4)) + 2*MB^(3/2)*R^(3/2)*(2*mi*S2 + Sqrt[2]*(1 - I*omega*R)*(S3 - S4)) + 4*Sqrt[2]*MB^(5/2)*Sqrt[R]*(-S3 + S4) - (2*I)*MB^2*R*(4*S1 + 4*S2 + Sqrt[2]*mi*(-S3 + S4)) + I*MB*R^2*(2*S1 + (4 - 2*mi^2 + (4*I)*omega*R)*S2 + Sqrt[2]*mi*(-S3 + S4))))/(2*(Sqrt[MB]*mi - omega*R^(3/2))*Dcommon),
 fourvelocitycoef[2] -> (R*Sqrt[-3*MB + R]*((-4*I)*MB^2*S1 - I*omega^2*R^4*S1 + MB*R*((2*I - I*mi^2 + 2*omega*R)*S1 + I*(4*S2 + Sqrt[2]*mi*(S3 - S4))) - 2*MB^(3/2)*Sqrt[R]*(mi*S1 + Sqrt[2]*(-S3 + S4)) + I*Sqrt[MB]*omega*R^(5/2)*(2*mi*S1 + Sqrt[2]*(-S3 + S4))))/((Sqrt[MB]*mi - omega*R^(3/2))*Dcommon),
 fourvelocitycoef[3] -> (Sqrt[-3*MB + R]*(-4*MB^(5/2)*Sqrt[R]*S1 + 2*MB^(3/2)*R^(3/2)*(S1 + I*omega*R*S1 + 2*S2) + I*Sqrt[2]*omega^2*R^5*S4 - I*Sqrt[MB]*omega*R^(7/2)*(S1 - 2*S2 + 2*Sqrt[2]*mi*S4) + I*MB*R^2*(mi*(S1 - 2*S2) + Sqrt[2]*mi^2*S4 - 2*Sqrt[2]*(S3 + S4)) - (2*I)*MB^2*R*(mi*S1 - 2*Sqrt[2]*(S3 + 2*S4))))/(Sqrt[2]*(Sqrt[MB]*mi - omega*R^(3/2))*Dcommon),
 fourvelocitycoef[4] -> (Sqrt[-3*MB + R]*(4*Sqrt[2]*MB^(5/2)*Sqrt[R]*S1 + 2*Sqrt[2]*MB^(3/2)*R^(3/2)*((-1 - I*omega*R)*S1 - 2*S2) + (2*I)*omega^2*R^5*S3 + I*Sqrt[MB]*omega*R^(7/2)*(Sqrt[2]*S1 - 2*Sqrt[2]*S2 - 4*mi*S3) + (2*I)*MB^2*R*(Sqrt[2]*mi*S1 + 8*S3 + 4*S4) - I*MB*R^2*(Sqrt[2]*mi*(S1 - 2*S2) - 2*mi^2*S3 + 4*(S3 + S4))))/(2*(Sqrt[MB]*mi - omega*R^(3/2))*Dcommon)};

LHsnormaftersubst = ((2*I)*MB*Sqrt[R]*S1 - I*R^(3/2)*(S1 + 2*S2) + Sqrt[2]*Sqrt[MB]*R*(-S3 + S4))/(-(Sqrt[MB]*mi) + omega*R^(3/2));

densitysubstrule = {
 drho -> (-4^(-1))*((md*(-3*MB + R)^(3/2)*(2*R*((-4*I)*MB + R*(2*I + omega*R))*fourvelocitycoef[1] - I*(2*MB - R)*((4*MB + R*(-2 - I*omega*R))*fourvelocitycoef[2] + Sqrt[2]*mi*R*(-fourvelocitycoef[3] + fourvelocitycoef[4]))))/(Pi*R^3*(-2*MB + R)^2*(-(Sqrt[MB]*mi) + omega*R^(3/2)))),
 dtheta -> ((-I)*Sqrt[-3*MB + R]*(fourvelocitycoef[3] + fourvelocitycoef[4]))/(Sqrt[2]*(Sqrt[MB]*mi - omega*R^(3/2))),
 dr -> ((-I)*Sqrt[-3*MB + R]*(2*R*fourvelocitycoef[1] + (2*MB - R)*fourvelocitycoef[2]))/(2*Sqrt[MB]*mi - 2*omega*R^(3/2))};

(* Row 5 of i1.0_p12.0_e0_r20.0.dat: m=5,n=0,k=0 *)
dat = Import["/data/haiyangw/claude/emri-disk/packages/ModifiedTeukolsky/test/data/i1.0_p12.0_e0_r20.0.dat"];
row = dat[[5]];
vals = {MB -> 1, R -> 20, omega -> row[[4]], mi -> row[[1]], md -> 10^(-4),
        S1 -> row[[5]] + I*row[[6]], S2 -> row[[7]] + I*row[[8]],
        S3 -> row[[9]] + I*row[[10]], S4 -> row[[11]] + I*row[[12]]};
u = Table[fourvelocitycoef[i] /. PertGeodesiceqsolution /. vals, {i, 1, 4}];
uvals = Thread[{fourvelocitycoef[1], fourvelocitycoef[2], fourvelocitycoef[3], fourvelocitycoef[4]} -> u];
dens = ({drho, dtheta, dr} /. densitysubstrule) /. uvals /. vals;
normchk = LHsnormaftersubst /. vals;
pr[x_] := ToString[Re[N[x,20]], InputForm] <> " , " <> ToString[Im[N[x,20]], InputForm];
Print["U1 ", pr[u[[1]]]];
Print["U2 ", pr[u[[2]]]];
Print["U3 ", pr[u[[3]]]];
Print["U4 ", pr[u[[4]]]];
Print["drho ", pr[dens[[1]]]];
Print["dtheta ", pr[dens[[2]]]];
Print["dr ", pr[dens[[3]]]];
Print["normcheck ", pr[normchk]];
