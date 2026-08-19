(* TRACKED ORACLE SCRIPT  (GRR::g14-oracle-scripts, repro-audit gap G14)
   
   origin : emri-disk/scratchpad/qqnm_phase3/run_qnm.wls (untracked one-shot)
   sha256 : 7e14865d6a233e55415bb7841cc17d7b51cc97bd11bc818f6b61eaea8a703a3a
   backs  : QNM-Mathematica cross-oracle runs behind test/fixtures/qnm_hphc_mma_ref.jl and the QNM frequency tolerances in test_qnm_frequency.jl
   engine : wolframscript + QNM-Mathematica (HEAD a74aadd) + BHPToolkit Teukolsky/SWSH, WL 14.3.0
   
   Tracked for provenance/regeneration only: nothing in the Julia test suite executes it
   (no Wolfram/py dependency is added).  Re-running it requires the engine above. *)

(* Load the scratch QNM paclet + BHPT deps, run all compute methods *)
PacletDirectoryLoad["/data/haiyangw/claude/emri-disk/scratchpad/qqnm_phase3/qnm_scratch/QNM-Mathematica"];
Needs["QNM`"];
Print["Loaded QNM. Kernel file used: ", FindFile["QNM`"]];

fmt[z_] := If[Head[z]===Complex || NumericQ[z],
  ToString[NumberForm[Re[z], 17], InputForm] <> (If[Im[z]>=0,"+","-"]) <> ToString[NumberForm[Abs[Im[z]], 17], InputForm] <> "im",
  ToString[z]];

anchors = <|0.0 -> (0.37367168441804177 - 0.08896231568893546 I)|>;

avals = {0.0, 0.3, 0.6, 0.9};
methods = {"Interpolation", Automatic, "Leaver", "IncidenceAmplitude"};

Do[
  Print["================ a = ", a, " ================"];
  Do[
    Module[{res, t},
      {t, res} = AbsoluteTiming[
        Quiet[Check[QNMFrequency[-2, 2, 2, 0, N[a, 30], Method -> mth], "ERR"]]
      ];
      Print["  Method ", ToString[mth], "  ->  ",
        If[NumericQ[res], fmt[res], ToString[res]],
        "    [", ToString[NumberForm[t,3]], "s]"];
    ],
    {mth, methods}
  ],
  {a, avals}
];

Print["done"];
