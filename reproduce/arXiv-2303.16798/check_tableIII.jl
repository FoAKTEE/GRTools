# arXiv:2303.16798 Table III (tab:generic, paper.tex L878-903): the linear-in-spin amplitude parts
# C^±_{S,lmnk} at infinity (C^+) and at the horizon (C^-) for 24 published (m,l,n,k) modes of the
# generic orbit a = 0.9 M, p = 12, e = 0.2, I = 30° (x = cos 30°).
#
# THIN driver: every number comes from Perturbation.SpinFluxes
# (`KerrSpinOrbitCorrection` + `TeukolskySpinModeCorrection`); the script only arranges them.
# Sign convention: C_GrTools = −C_paper (the flux |C|² is convention-independent).
#
#   julia --project=reproduce -t4 reproduce/arXiv-2303.16798/check_tableIII.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Perturbation: KerrSpinOrbitCorrection, TeukolskySpinModeCorrection

const A, P, E, X = 0.9, 12.0, 0.2, cos(pi / 6)

# published (m, l, n, k, Re C^+, Im C^+, Re C^-, Im C^-)
const TABLE_III = [
    (1, 2, 0, 1,  4.8962e-6, -1.6020e-6, -5.2716e-6, -2.7823e-7),
    (1, 2, 1, 1,  9.9514e-6, -2.7846e-6,  3.8592e-6,  8.7391e-7),
    (1, 2, 2, 1,  7.3027e-6, -2.1468e-6,  3.7407e-6,  7.8065e-7),
    (1, 2, 3, 1,  3.6008e-6, -1.1232e-6,  1.7742e-6,  3.8842e-7),
    (1, 3, 0, 2, -9.4587e-8,  1.3025e-7, -5.8770e-7, -4.1971e-7),
    (1, 3, 1, 2, -8.8801e-7, -1.7932e-6, -3.2098e-7, -2.1754e-7),
    (1, 3, 2, 2, -9.8569e-7, -1.9468e-6, -1.1060e-7, -6.4897e-8),
    (1, 3, 3, 2, -6.6574e-7, -1.2388e-6, -3.7354e-8, -1.4854e-8),
    (2, 2, 0, 0, -1.9890e-5,  5.8986e-6, -7.5636e-6, -3.6977e-6),
    (2, 2, 1, 0, -3.6535e-5,  1.0473e-5, -2.8727e-5, -9.5114e-6),
    (2, 2, 2, 0, -2.8239e-5,  8.6430e-6, -2.1354e-5, -7.5839e-6),
    (2, 2, 3, 0, -1.5302e-5,  4.9730e-6, -1.0601e-5, -4.1408e-6),
    (2, 3, 0, 1,  8.2420e-7,  9.2288e-7,  6.8449e-7,  1.3381e-6),
    (2, 3, 1, 1,  3.8727e-6,  7.7467e-6, -1.2893e-7, -2.4943e-7),
    (2, 3, 2, 1,  4.3094e-6,  8.2104e-6, -4.2520e-7, -8.2810e-7),
    (2, 3, 3, 1,  3.0636e-6,  5.4471e-6, -3.3219e-7, -6.4545e-7),
    (3, 3, 0, 0, -2.2312e-6, -3.3582e-6, -3.0967e-7, -1.2954e-6),
    (3, 3, 1, 0, -8.9746e-6, -1.7781e-5,  5.0351e-7,  3.8648e-6),
    (3, 3, 2, 0, -1.0099e-5, -1.8845e-5,  6.2534e-7,  4.8041e-6),
    (3, 3, 3, 0, -7.3942e-6, -1.2830e-5,  4.0686e-7,  3.3295e-6),
    (3, 4, 0, 1,  1.1671e-6, -4.1185e-7,  1.9807e-7, -2.0808e-7),
    (3, 4, 1, 1, -3.6720e-6,  2.6844e-6,  4.2934e-8, -4.2764e-8),
    (3, 4, 2, 1, -5.6428e-6,  4.2331e-6, -8.4756e-8,  9.4155e-8),
    (3, 4, 3, 1, -4.6912e-6,  3.7655e-6, -1.0502e-7,  1.1655e-7),
]

orbit = KerrSpinOrbitCorrection(A, P, E, X, 8, 8)

nrows = length(TABLE_III)
Cplus  = Vector{ComplexF64}(undef, nrows)
Cminus = Vector{ComplexF64}(undef, nrows)
Threads.@threads for i in 1:nrows
    m, l, n, k = TABLE_III[i][1], TABLE_III[i][2], TABLE_III[i][3], TABLE_III[i][4]
    md = TeukolskySpinModeCorrection(l, m, n, k, orbit)
    Cplus[i]  = -ComplexF64(md.AmplitudesCorrection.I)     # C_GrTools = −C_paper
    Cminus[i] = -ComplexF64(md.AmplitudesCorrection.H)
end

pubP = [complex(r[5], r[6]) for r in TABLE_III]
pubM = [complex(r[7], r[8]) for r in TABLE_III]
relP = abs.(Cplus  .- pubP) ./ abs.(pubP)
relM = abs.(Cminus .- pubM) ./ abs.(pubM)

println("[table III] worst |ΔC^+|/|C^+| = ", maximum(relP), " at row ", argmax(relP))
println("[table III] worst |ΔC^-|/|C^-| = ", maximum(relM), " at row ", argmax(relM))

@repro_result "modes"            [string(r[1], ",", r[2], ",", r[3], ",", r[4]) for r in TABLE_III]
@repro_result "Cplus_computed"   Cplus
@repro_result "Cminus_computed"  Cminus
@repro_result "Cplus_published"  pubP
@repro_result "Cminus_published" pubM
@repro_result "max_rel_Cplus"    maximum(relP)
@repro_result "max_rel_Cminus"   maximum(relM)
@repro_result "all_rows_within_0p5pct" (maximum(relP) < 5e-3) && (maximum(relM) < 5e-3)

exit(ReproHarness.conclude!(opts; paper = "2303.16798", id = "tableIII", script = @__FILE__,
                            rtol = 1e-8,
                            tolerances = Dict("max_rel_Cplus"  => (rtol = 1e-6, atol = 1e-12),
                                              "max_rel_Cminus" => (rtol = 1e-6, atol = 1e-12))))
