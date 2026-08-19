(* TRACKED ORACLE SCRIPT  (GRR::g14-oracle-scripts, repro-audit gap G14)
   
   origin : emri-disk/results/gw-flux/reference/gen/gen_pn_series_ops.wl @ 29b9fa7d (tracked in emri-disk, not in GRTools)
   sha256 : 68c545a1e465e6703ebe8af1567a659c67dc235a4458e0152b7d195b7e632119
   backs  : test/fixtures/pn_series_ops_oracle.json (vendored verbatim); consumed by test/test_series.jl
   engine : wolframscript (Wolfram Language); needs the BHPToolkit PN series bundle
   note   : This is the ONLY tracked wolframscript generator behind PostNewtonian's baked literals. test_catalog.jl / test_selfforce.jl / test_mmainputform.jl cite one-shot INTERACTIVE wolframscript evaluations against the bundled SeriesData *.m files (which ARE in-tree); no script for those was ever written [HOLE].
   
   Tracked for provenance/regeneration only: nothing in the Julia test suite executes it
   (no Wolfram/py dependency is added).  Re-running it requires the engine above. *)

(* ::Package:: *)
(* Oracle generator for TeukolskyPN slice-1 (PNSeries type + Tools.wl ops).
   Builds SeriesData test cases (symbolic coeffs in v) and applies each op.
   Dumps SeriesData slots (coeffs InputForm strings + nmin/nmax/den) and op results
   to results/gw-flux/reference/pn_series_ops_oracle.json.

   The Tools.wl ops are transcribed here directly from PN.wl lines 1020-1164 so the
   oracle is the REAL semantics (not a reimplementation that could drift). *)

(* ---- Tools.wl ops, transcribed verbatim from Teukolsky`PN`Private (PN.wl) ---- *)
ClearAll[SeriesMinOrder, SeriesMaxOrder, SeriesLength, SeriesCollect, SeriesTake,
         SeriesTerms, polyToSeries, ChangeSeriesParameter, PowerCounting, ExpandLog, ExpandGamma];

SeriesMinOrder[series_SeriesData] := series[[4]]/series[[6]];
SeriesMinOrder[1] = 0;
Attributes[SeriesMinOrder] = {Listable};

SeriesMaxOrder[series_SeriesData] := series[[5]]/series[[6]];
Attributes[SeriesMaxOrder] = {Listable};

SeriesLength[series_SeriesData] := SeriesMaxOrder[series] - SeriesMinOrder[series];

SeriesCollect[series_SeriesData, var__, func_:Identity] := Collect[#, var, func] & /@ series;
SeriesCollect[list_List, var__, func_:Identity] := SeriesCollect[#, var, func] & /@ list;
SeriesCollect[expr_, var__, func_:Identity] := Collect[#, var, func] & @ expr;

SeriesTake[series_SeriesData, order_Integer:1] := series (1 + O[series[[1]]]^(order));
SeriesTake[series_SeriesData, 0] := Block[{minOrder}, minOrder = series // SeriesMinOrder; O[series[[1]]]^(minOrder)];
SeriesTake[series_O, order_Integer:1] := series;
Attributes[SeriesTake] = {Listable};

SeriesTerms[expr_, {x_, x0_, termOrder_}] := Module[{minOrder},
  minOrder = Series[expr, x -> x0] // SeriesMinOrder;
  Series[expr, {x, x0, minOrder + termOrder - 1}]];
SeriesTerms[expr___] := Series[expr];

polyToSeries[poly_, x_:\[Eta], x0_:0] := Block[{maxPower},
  maxPower = poly // Exponent[#, x] & // Ceiling;
  poly + O[(x - x0)]^(maxPower + 1)];
polyToSeries[0, x_:\[Eta], x0_:0] := 0;

ChangeSeriesParameter[series_SeriesData, var_Symbol] := ReplacePart[series, 1 -> var];
ChangeSeriesParameter[series_SeriesData, var_/;MatchQ[var, Power[_Symbol, __]]] := Module[
  {sym, pow, powNum, powDen, oldMin, oldMax, oldDen, oldCoeffs, newCoeffs},
  {sym, pow} = var /. Power[b_Symbol, c__] :> {b, c};
  {powNum, powDen} = {pow // Numerator, pow // Denominator};
  {oldCoeffs, oldMin, oldMax, oldDen} = series[[#]] & /@ {3, 4, 5, 6};
  newCoeffs = {oldCoeffs}~Join~ConstantArray[ConstantArray[0, Length[oldCoeffs]], powNum-1] // Transpose // Flatten;
  ReplacePart[series, {1 -> sym, 3 -> newCoeffs, 4 -> oldMin powNum, 5 -> oldMax powNum, 6 -> oldDen powDen}]];
ChangeSeriesParameter[series_SeriesData, var_/;MatchQ[var, __ _Symbol]] := Module[
  {sym, fac, oldCoeffs, oldMin, oldMax, oldDen, exps, newCoeffs},
  {fac, sym} = var /. b__ c_Symbol :> {b, c};
  {oldCoeffs, oldMin, oldMax, oldDen} = series[[#]] & /@ {3, 4, 5, 6};
  exps = Range[oldMin/oldDen, oldMin/oldDen + (Length[oldCoeffs]-1)/oldDen, 1/oldDen];
  newCoeffs = oldCoeffs (fac^# & /@ exps);
  ReplacePart[series, {1 -> sym, 3 -> newCoeffs}]];
ChangeSeriesParameter[series_SeriesData, var_/;MatchQ[var, __ Power[_Symbol, __]]] := Module[
  {sym, fac, pow, powDen, powNum, oldCoeffs, oldMin, oldMax, oldDen, exps, newCoeffs},
  {fac, sym, pow} = var /. a__ Power[b_Symbol, c__] :> {a, b, c};
  {powNum, powDen} = {pow // Numerator, pow // Denominator};
  {oldCoeffs, oldMin, oldMax, oldDen} = series[[#]] & /@ {3, 4, 5, 6};
  exps = Range[oldMin/oldDen, oldMin/oldDen + (Length[oldCoeffs]-1)/oldDen, 1/oldDen];
  newCoeffs = oldCoeffs (fac^# & /@ exps);
  newCoeffs = {newCoeffs}~Join~ConstantArray[ConstantArray[0, Length[newCoeffs]], powNum-1] // Transpose // Flatten;
  ReplacePart[series, {1 -> sym, 3 -> newCoeffs, 4 -> oldMin powNum, 5 -> oldMax powNum, 6 -> oldDen powDen}]];
ChangeSeriesParameter[a_, b_] := a;

PowerCounting[series_SeriesData, var_] := Module[
  {sym, oldCoeffs, oldMin, oldMax, oldDen, exps, newCoeffs},
  {sym, oldCoeffs, oldMin, oldMax, oldDen} = (series[[#]] &) /@ {1, 3, 4, 5, 6};
  exps = Range[oldMin/oldDen, oldMin/oldDen + (Length[oldCoeffs]-1)/oldDen, 1/oldDen];
  newCoeffs = oldCoeffs (sym^# & /@ exps);
  ReplacePart[series, {1 -> var, 3 -> newCoeffs}]];

ExpandLog[expr_] := (expr /. Log[a_] :> PowerExpand[Log[a]]);

ExpandGamma[expr_Gamma, j_Integer:0] := Module[{arg, x, n},
  arg = expr /. Gamma[argg_] :> argg;
  {n, x} = If[MatchQ[arg, a_Integer + b__], arg /. nn_Integer + xx__ :> {nn, xx}, {0, arg}];
  (Product[(x + j + i), {i, 0, n - j - 1}]) (Product[(x + j + i)^(-1), {i, n - j, -1}]) Gamma[x + j]];
ExpandGamma[expr_, j_Integer:0] /; FreeQ[expr, Gamma] := expr;
ExpandGamma[expr_, j_Integer:0] := ExpandGamma[#, j] & /@ expr;

(* ---- helpers to serialize ---- *)
inp[e_] := ToString[e, InputForm];

(* Dump a SeriesData (or plain expr) to an assoc with structural slots + a numeric sample. *)
seriesToAssoc[s_SeriesData] := <|
  "head" -> "SeriesData",
  "var" -> inp[s[[1]]],
  "x0" -> inp[s[[2]]],
  "coeffs" -> (inp /@ (s[[3]])),
  "nmin" -> s[[4]], "nmax" -> s[[5]], "den" -> s[[6]],
  "normal" -> inp[Normal[s]]
|>;
seriesToAssoc[e_] := <| "head" -> "Expr", "value" -> inp[e], "normal" -> inp[e] |>;

(* numeric evaluation of Normal[result] at given (v, eta) -- both real and complex probes *)
numAt[s_, vval_, etaval_, etasym_] := Module[{nf},
  nf = Normal[s] /. {v -> vval, etasym -> etaval};
  inp[N[nf, 30]]];

(* ============================ BUILD TEST CASES ============================ *)
results = <||>;

(* Base series in eta with symbolic coeffs in v *)
A = Series[(1 + v \[Eta] + v^2 \[Eta]^2/2)/(1 - \[Eta]), {\[Eta], 0, 5}];   (* A: nmin 0 *)
B = Series[(v + \[Eta])/(1 - v \[Eta]), {\[Eta], 0, 5}];                    (* B *)
C0 = Series[(1 + \[Eta])^v, {\[Eta], 0, 4}];                               (* C: binomial-like, symbolic exponent *)
D0 = Series[Exp[v \[Eta]]/\[Eta], {\[Eta], 0, 4}];                          (* D: Laurent, nmin = -1 *)
F0 = Series[Sqrt[1 + v \[Eta]], {\[Eta], 0, 4}];                            (* F: half-integer? no, integer powers; for den test use G *)
G0 = Series[Sqrt[\[Eta]] (1 + v \[Eta]), {\[Eta], 0, 3}];                   (* G: den = 2 (half-integer powers) *)

results["inputs"] = <|
  "A" -> seriesToAssoc[A], "B" -> seriesToAssoc[B], "C" -> seriesToAssoc[C0],
  "D" -> seriesToAssoc[D0], "F" -> seriesToAssoc[F0], "G" -> seriesToAssoc[G0]
|>;

(* ---------- ARITHMETIC (the crux: O[] propagation) ---------- *)
results["arith"] = <|
  "A_plus_B" -> seriesToAssoc[A + B],
  "A_times_B" -> seriesToAssoc[A B],
  "inv_A"     -> seriesToAssoc[1/A],
  "A_pow_3"   -> seriesToAssoc[A^3],
  "A_pow_m2"  -> seriesToAssoc[A^(-2)],
  "A_div_B"   -> seriesToAssoc[A/B],
  "A_minus_B" -> seriesToAssoc[A - B],
  "scalar_plus_A" -> seriesToAssoc[3 + A],
  "num_times_A"   -> seriesToAssoc[v^2 A],
  "A_plus_num"    -> seriesToAssoc[A + v],
  "D_times_A"     -> seriesToAssoc[D0 A],   (* Laurent * regular *)
  "D_pow_2"       -> seriesToAssoc[D0^2],
  "C_times_C"     -> seriesToAssoc[C0 C0],
  "G_times_G"     -> seriesToAssoc[G0 G0],  (* half-integer den arithmetic *)
  "inv_G"         -> seriesToAssoc[1/G0]
|>;

(* ---------- SeriesTake ---------- *)
results["take"] = <|
  "A_take_1" -> seriesToAssoc[SeriesTake[A, 1]],
  "A_take_3" -> seriesToAssoc[SeriesTake[A, 3]],
  "A_take_0" -> seriesToAssoc[SeriesTake[A, 0]],
  "D_take_2" -> seriesToAssoc[SeriesTake[D0, 2]],
  "B_take_4" -> seriesToAssoc[SeriesTake[B, 4]],
  "G_take_2" -> seriesToAssoc[SeriesTake[G0, 2]]
|>;

(* ---------- SeriesCollect ---------- *)
results["collect"] = <|
  "A_collect_v"        -> seriesToAssoc[SeriesCollect[A, v]],
  "B_collect_v"        -> seriesToAssoc[SeriesCollect[B, v]],
  "A_collect_v_expand" -> seriesToAssoc[SeriesCollect[A B, v, Expand]]
|>;

(* ---------- SeriesTerms ---------- *)
results["terms"] = <|
  "expr1_3terms" -> seriesToAssoc[SeriesTerms[(1 + v \[Eta])/(1 - \[Eta]), {\[Eta], 0, 3}]],
  "expr2_2terms" -> seriesToAssoc[SeriesTerms[Exp[v \[Eta]]/\[Eta], {\[Eta], 0, 2}]],
  "expr3_4terms" -> seriesToAssoc[SeriesTerms[Sqrt[1 + v \[Eta]], {\[Eta], 0, 4}]]
|>;

(* ---------- polyToSeries ---------- *)
results["polyToSeries"] = <|
  "p1" -> seriesToAssoc[polyToSeries[1 + v \[Eta] + v^2 \[Eta]^3, \[Eta], 0]],
  "p2" -> seriesToAssoc[polyToSeries[v \[Eta]^2, \[Eta], 0]],
  "p0" -> seriesToAssoc[polyToSeries[0, \[Eta], 0]]
|>;

(* ---------- ChangeSeriesParameter (4 forms) ---------- *)
results["change"] = <|
  "A_rename"     -> seriesToAssoc[ChangeSeriesParameter[A, x]],            (* symbol *)
  "A_rescale"    -> seriesToAssoc[ChangeSeriesParameter[A, 2 \[Omega]]],   (* fac*sym *)
  "A_rescale_v"  -> seriesToAssoc[ChangeSeriesParameter[A, v x]],          (* fac*sym, symbolic fac *)
  "A_pow2"       -> seriesToAssoc[ChangeSeriesParameter[A, x^2]],          (* sym^pow *)
  "A_facpow"     -> seriesToAssoc[ChangeSeriesParameter[A, 3 x^2]],        (* fac*sym^pow *)
  "G_pow2"       -> seriesToAssoc[ChangeSeriesParameter[G0, x^2]],         (* sym^pow on half-int den *)
  "A_facpow_v"   -> seriesToAssoc[ChangeSeriesParameter[A, v x^2]]         (* fac*sym^pow, symbolic fac *)
|>;

(* ---------- PowerCounting ---------- *)
results["powercounting"] = <|
  "A_pc" -> seriesToAssoc[PowerCounting[A, x]],
  "D_pc" -> seriesToAssoc[PowerCounting[D0, x]],
  "G_pc" -> seriesToAssoc[PowerCounting[G0, x]]
|>;

(* ---------- min/max/length ---------- *)
results["orders"] = <|
  "A" -> <| "min" -> inp[SeriesMinOrder[A]], "max" -> inp[SeriesMaxOrder[A]], "len" -> inp[SeriesLength[A]] |>,
  "D" -> <| "min" -> inp[SeriesMinOrder[D0]], "max" -> inp[SeriesMaxOrder[D0]], "len" -> inp[SeriesLength[D0]] |>,
  "G" -> <| "min" -> inp[SeriesMinOrder[G0]], "max" -> inp[SeriesMaxOrder[G0]], "len" -> inp[SeriesLength[G0]] |>
|>;

(* ---------- Normal ---------- *)
results["normal"] = <|
  "A" -> inp[Normal[A]], "D" -> inp[Normal[D0]], "G" -> inp[Normal[G0]]
|>;

(* ---------- ExpandLog ---------- *)
results["expandlog"] = <|
  "l1" -> inp[ExpandLog[Log[a b c]]],
  "l2" -> inp[ExpandLog[Log[x^n]]],
  "l3" -> inp[ExpandLog[Log[a^2 b/c]]],
  "l4" -> inp[ExpandLog[2 Log[v \[Eta]] + Log[(1 - \[Eta])^3]]]
|>;

(* ---------- ExpandGamma ---------- *)
results["expandgamma"] = <|
  "g1" -> inp[ExpandGamma[Gamma[3 + x]]],
  "g2" -> inp[ExpandGamma[Gamma[x]]],
  "g3" -> inp[ExpandGamma[Gamma[2 + x] + Gamma[1 + y]]],
  "g4" -> inp[ExpandGamma[Gamma[n + 1]]],     (* n symbolic, no integer head -> unchanged-ish *)
  "g5" -> inp[ExpandGamma[v Gamma[4 + x]]],
  "g6shift" -> inp[ExpandGamma[Gamma[3 + x], 1]]
|>;

(* numeric samples of arithmetic results for tolerance-based test (avoids InputForm parsing drift) *)
probesV  = {3/10, 11/10, -7/10};
probesEt = {1/100, 1/20, -3/100};
numProbe[s_] := Table[
  inp[N[Normal[s] /. {v -> probesV[[i]], \[Eta] -> probesEt[[i]]}, 30]],
  {i, 1, Length[probesV]}];

results["arith_numeric"] = <|
  "A_plus_B" -> numProbe[A + B], "A_times_B" -> numProbe[A B],
  "inv_A" -> numProbe[1/A], "A_pow_3" -> numProbe[A^3],
  "A_pow_m2" -> numProbe[A^(-2)], "A_div_B" -> numProbe[A/B],
  "A_minus_B" -> numProbe[A - B], "scalar_plus_A" -> numProbe[3 + A],
  "num_times_A" -> numProbe[v^2 A], "A_plus_num" -> numProbe[A + v],
  "D_times_A" -> numProbe[D0 A], "D_pow_2" -> numProbe[D0^2],
  "C_times_C" -> numProbe[C0 C0], "G_times_G" -> numProbe[G0 G0],
  "inv_G" -> numProbe[1/G0]
|>;

Export["/data/haiyangw/claude/gw-flux/results/gw-flux/reference/pn_series_ops_oracle.json", results, "JSON"];
Print["WROTE pn_series_ops_oracle.json"];
Print["A slots: ", {A[[4]], A[[5]], A[[6]]}];
Print["A*B nmax/den: ", {(A B)[[5]], (A B)[[6]]}];
Print["1/A nmax: ", (1/A)[[5]]];
Print["A^3 nmax: ", (A^3)[[5]]];
Print["G*G slots: ", {(G0 G0)[[4]], (G0 G0)[[5]], (G0 G0)[[6]]}];
Print["A_pow2 (change x^2) slots: ", With[{r=ChangeSeriesParameter[A,x^2]}, {r[[4]],r[[5]],r[[6]], r[[3]]}]];
