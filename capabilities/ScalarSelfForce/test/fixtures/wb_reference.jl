# wb_reference.jl — BAKED Warburton–Barack printed reference-SSF tables (ESSF node
# `wb-reference-fixtures`, schedule W5-F, GLOBAL_DAG row 6). DO NOT EDIT BY HAND: every
# literal below is asserted STRING-EQUAL to the tex mirror by test_wb_fixtures.jl
# (independent re-extraction parser), so any hand edit that deviates from the print fails CI.
#
# CITATION: N. Warburton and L. Barack, "Self force on a scalar charge in Kerr spacetime:
#   circular equatorial orbits", Phys. Rev. D 81, 084039 (2010); arXiv:1003.1860.
#   Read-only tex mirror: ref-paper/arxiv-1003.1860/src/Warburton.tex (PROVENANCE.md;
#   sha256 = WB_TEX_SHA256 below, re-pinned by the test when the mirror is present).
#   * WB_FT — Table II, \label{table:temporal-sf-results}: (M^2/q^2) F_t.
#       tex: table env lines 567-592, tabular 569-588, DATA ROWS 575-586
#       (= WB_FT_TEX_DATA_LINES); caption lines 589-591.
#   * WB_FR — Table III, \label{table:radial-sf-results}: (M^2/q^2) F_r.
#       tex: table env lines 620-645, tabular 622-641, DATA ROWS 628-639
#       (= WB_FR_TEX_DATA_LINES); caption lines 642-644.
#   * Table I (\label{table:fluxes}, tex lines 513-550) prints ENERGY FLUXES Edot, not SSF
#     components — deliberately OUT OF SCOPE here. It is the only table with |a|=0.998 rows.
#
# GENERATION RECIPE (one-shot; the test supersedes it as the reproducibility guarantee):
#   split each data row of the two tabulars on '&', strip math-mode '$', map the printed
#   exponent macro "e{-}k" -> "e-k" (macro defined in the Table I caption, tex line 546),
#   bake '-' cells as WB_BELOW_ISCO. No other transform: significant figures are EXACTLY
#   as printed ("All figures presented are significant", captions tex lines 589/642).
#
# CONVENTIONS: G = c = 1, M = black-hole mass, q = scalar charge; a < 0 encodes retrograde
#   orbits (motion always +phi); tables store the dimensionless printed products
#   (M^2/q^2) F_t and (M^2/q^2) F_r on the grid r0 in WB_R0_OVER_M (rows) x a in
#   WB_A_OVER_M (columns), both in the printed order.
#   Empty printed cells ('-') are orbits below the ISCO (captions, tex lines 589/642) and
#   are stored as the typed sentinel WB_BELOW_ISCO :: WBBelowISCO.
#
# PRINTED PARITY DOMAIN [recorded decision Q2]: the printed F_t/F_r grids cover
#   |a| <= 0.9 M (= WB_MAX_ABS_A) and r0 >= 4 M (= WB_MIN_PRINTED_R0; r0 = 4M printed only
#   for a = 0.7M, 0.9M). The CRP's |a| = 0.998/near-ISCO WB comparison used PRIVATE WB data
#   that is NOT in the paper — parity claims against this fixture must stay inside the
#   printed domain; wb_cell throws ArgumentError outside it.
#
# STANDALONE: pure Base — no package or stdlib deps; test_wb_fixtures.jl loads this file
#   into a fresh module to enforce that.

# A printed table cell: `dec` is the EXACT printed decimal string (sole source of truth,
# normalized "e{-}k" -> "e-k" only); `value` is derived from it at load time.
struct WBRefCell
    dec::String
    value::Float64
    WBRefCell(dec::AbstractString) = new(String(dec), parse(Float64, dec))
end

# Typed sentinel for the printed '-' cells: circular orbit below the ISCO, no SSF printed.
struct WBBelowISCO end
const WB_BELOW_ISCO = WBBelowISCO()

# sha256 of the tex mirror the literals were extracted from (drift pin, test gate A).
const WB_TEX_SHA256 = "334cb564a7c8763d701508cc92a671dbaa0af09ca8078baef82e9700ba8ae408"

# 1-based tex line ranges of the data rows (re-checked by the extraction-equality test).
const WB_FT_TEX_DATA_LINES = 575:586
const WB_FR_TEX_DATA_LINES = 628:639

# Printed grid axes, in printed order (columns / rows of the matrices below).
const WB_A_OVER_M  = [-0.9, -0.7, -0.5, 0.0, 0.5, 0.7, 0.9]
const WB_R0_OVER_M = [4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 14.0, 20.0, 30.0, 50.0, 70.0, 100.0]

# Printed coverage bounds [recorded decision Q2].
const WB_MAX_ABS_A = 0.9
const WB_MIN_PRINTED_R0 = 4.0

# Table-literal helpers (BI is shorthand for the sentinel INSIDE the literals below only).
const BI = WB_BELOW_ISCO
function _wb_bake(tbl::AbstractMatrix)::Matrix{Union{WBRefCell, WBBelowISCO}}
    return map(x -> x isa WBBelowISCO ? x : WBRefCell(x), tbl)
end

# Table II — (M^2/q^2) F_t, tex data rows 575-586 (67 printed cells, 17 dashes).
const WB_FT = _wb_bake(Any[
    # a=-0.9M              a=-0.7M              a=-0.5M              a=0                  a=0.5M               a=0.7M               a=0.9M
    BI                   BI                   BI                   BI                   BI                   "1.35921815e-3"      "1.14204820e-3"   # r0 = 4M
    BI                   BI                   BI                   BI                   "6.07684087e-4"      "5.35768561e-4"      "4.79634985e-4"   # r0 = 5M
    BI                   BI                   BI                   "3.60907254e-4"      "2.78394798e-4"      "2.55161013e-4"      "2.35733853e-4"   # r0 = 6M
    BI                   BI                   BI                   "1.76732019e-4"      "1.46366447e-4"      "1.37103703e-4"      "1.29046747e-4"   # r0 = 7M
    BI                   BI                   "1.15781360e-4"      "9.77204485e-5"      "8.44876316e-5"      "8.02407393e-5"      "7.64519160e-5"   # r0 = 8M
    "4.60475173e-5"      "4.38590519e-5"      "4.18429073e-5"      "3.75022727e-5"      "3.40410532e-5"      "3.28611197e-5"      "3.17760168e-5"   # r0 = 10M
    "1.03173965e-5"      "1.00539090e-5"      "9.80387438e-6"      "9.23672660e-6"      "8.74728207e-6"      "8.57077224e-6"      "8.40373578e-6"   # r0 = 14M
    "2.28457108e-6"      "2.25311511e-6"      "2.22274047e-6"      "2.15159216e-6"      "2.08709237e-6"      "2.06300902e-6"      "2.03980574e-6"   # r0 = 20M
    "4.30761267e-7"      "4.27729235e-7"      "4.24767592e-7"      "4.17678576e-7"      "4.11035602e-7"      "4.08496912e-7"      "4.06021007e-7"   # r0 = 30M
    "5.43419839e-8"      "5.41729302e-8"      "5.40064364e-8"      "5.36016621e-8"      "5.32132722e-8"      "5.30623647e-8"      "5.29138807e-8"   # r0 = 50M
    "1.40256823e-8"      "1.39999178e-8"      "1.39744575e-8"      "1.39121644e-8"      "1.38518165e-8"      "1.38282103e-8"      "1.38048982e-8"   # r0 = 70M
    "3.35072295e-9"      "3.34717963e-9"      "3.34366914e-9"      "3.33503895e-9"      "3.32661812e-9"      "3.32330755e-9"      "3.32002917e-9"   # r0 = 100M
])

# Table III — (M^2/q^2) F_r, tex data rows 628-639 (67 printed cells, 17 dashes).
# Sign structure: repulsive (F_r > 0) for a <= 0 at every printed radius; attractive
# (F_r < 0) at every printed prograde cell (the critical radius r_c(a) < 5M for a >= 0.5M).
const WB_FR = _wb_bake(Any[
    # a=-0.9M              a=-0.7M              a=-0.5M              a=0                  a=0.5M               a=0.7M               a=0.9M
    BI                   BI                   BI                   BI                   BI                   "-5.24194e-4"        "-9.5941e-4"      # r0 = 4M
    BI                   BI                   BI                   BI                   "-4.160235e-5"       "-2.044174e-4"       "-3.63448e-4"     # r0 = 5M
    BI                   BI                   BI                   "1.677283e-4"        "-2.421685e-5"       "-9.528095e-5"       "-1.645525e-4"    # r0 = 6M
    BI                   BI                   BI                   "7.850679e-5"        "-1.467677e-5"       "-4.980678e-5"       "-8.410331e-5"    # r0 = 7M
    BI                   BI                   "9.642777e-5"        "4.082502e-5"        "-9.21907e-6"        "-2.829488e-5"       "-4.696081e-5"    # r0 = 8M
    "4.939995e-5"        "4.100712e-5"        "3.28942e-5"         "1.378448e-5"        "-4.03517e-6"        "-1.091819e-5"       "-1.768232e-5"    # r0 = 10M
    "9.968208e-6"        "8.303689e-6"        "6.67043e-6"         "2.720083e-6"        "-1.07573e-6"        "-2.561183e-6"       "-4.02935e-6"     # r0 = 14M
    "1.878548e-6"        "1.565128e-6"        "1.2550019e-6"       "4.93790e-7"         "-2.50260e-7"        "-5.43942e-7"        "-8.35474e-7"     # r0 = 20M
    "2.873310e-7"        "2.389538e-7"        "1.90843e-7"         "7.1719e-8"          "-4.595209e-8"       "-9.26682e-8"        "-1.391883e-7"    # r0 = 30M
    "2.74358e-8"         "2.272902e-8"        "1.803392e-8"        "6.3467e-9"          "-5.27419e-9"        "-9.90589e-9"        "-1.452810e-8"    # r0 = 50M
    "5.87543e-9"         "4.8525e-9"          "3.8312e-9"          "1.2845e-9"          "-1.25352e-9"        "-2.26649e-9"        "-3.27820e-9"     # r0 = 70M
    "1.1508e-9"          "9.4715e-10"         "7.4364e-10"         "2.356e-10"          "-2.7134e-10"        "-4.7388e-10"        "-6.7625e-10"     # r0 = 100M
])

# wb_cell(table, a, r0): printed cell of `table` (WB_FT or WB_FR) at spin a/M and orbital
# radius r0/M. Returns a WBRefCell (printed value) or WB_BELOW_ISCO (printed '-'). Throws
# ArgumentError for any (a, r0) off the printed grids — in particular the |a| = 0.998 /
# near-ISCO domain of the CRP comparison (private WB data) [recorded decision Q2].
function wb_cell(table::AbstractMatrix, a::Real, r0::Real)
    j = findfirst(==(Float64(a)), WB_A_OVER_M)
    i = findfirst(==(Float64(r0)), WB_R0_OVER_M)
    (i === nothing || j === nothing) && throw(ArgumentError(
        "(a=$a, r0=$r0) is not on the printed WB grid: a/M in $(WB_A_OVER_M), " *
        "r0/M in $(WB_R0_OVER_M) [see PRINTED PARITY DOMAIN note, decision Q2]"))
    return table[i, j]
end
