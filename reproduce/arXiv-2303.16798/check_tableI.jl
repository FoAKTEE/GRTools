# arXiv:2303.16798 Table I (tab:error_generic, paper.tex L613-619): the linear-in-spin total energy
# flux F^E_{S,m} from GENERIC (eccentric + inclined) orbits at a = 0.9 M, summed over l, n and k,
# for the six published (p, e, I, m) rows.
#
# THIN driver: Perturbation.SpinFluxes supplies the orbit (`KerrSpinOrbitCorrection`) and every mode
# (`TeukolskySpinModeCorrection`); the script only sums and compares.
#
# SCOPE: this FAST driver reproduces ROW 1 only (p = 10, e = 0.1, I = 15°, m = 2). The sum runs
# l = 2..3, |n| ≤ 2, |k| ≤ 2 — the truncation the legacy lock used, which is convergent at e = 0.1
# and lands 0.30 % from the published −2.8259e-6 (the fully converged Julia value is −2.819e-6,
# 0.24 %, the MST-analytic generic-mode level).
#
# MEASURED FINDING (rows 2–6): at this truncation the remaining published rows are NOT converged —
# |n| ≤ 2 misses the radial harmonics that carry the flux at larger eccentricity. Measured
# rel diffs at l ≤ 3, |n| ≤ 2, |k| ≤ 2 (4 min at -t4): row2 2.7 %, row3 3.8 %, row4 17 %,
# row5 6.3 %, row6 (e = 0.5) 96 %. Those rows therefore need the e-scaled truncation of the
# HEAVY driver `check_tableI_full.jl` (registered, not yet run — see MANIFEST notes). Reporting
# the under-truncated numbers as a reproduction would be dishonest, so they are not registered here.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2303.16798/check_tableI.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Perturbation: KerrSpinOrbitCorrection, TeukolskySpinModeCorrection

const A = 0.9
const NS = -2:2
const KS = -2:2
const LMAX = 3

# published (p, e, I°, m, F^E_{S,m}, ΔF^E vs the time-domain Teukode)
const TABLE_I = [
    (10.0, 0.1, 15.0, 2, -2.8259e-6, 1e-3),      # row 1 — the only row convergent at |n| ≤ 2
]

# one orbit per distinct (p, e, I): rows 2–4 share it
orbits = Dict{NTuple{3,Float64},Any}()
for (p, e, Ideg, _, _, _) in TABLE_I
    get!(orbits, (p, e, Ideg)) do
        KerrSpinOrbitCorrection(A, p, e, cos(Ideg * pi / 180), 8, 8)
    end
end

FE = zeros(length(TABLE_I))
for (row, (p, e, Ideg, m, _, _)) in enumerate(TABLE_I)
    oc = orbits[(p, e, Ideg)]
    modes = [(l, n, k) for l in max(2, m):LMAX for n in NS for k in KS]
    term = Vector{ComplexF64}(undef, length(modes))
    Threads.@threads for i in eachindex(modes)
        l, n, k = modes[i]
        md = TeukolskySpinModeCorrection(l, m, n, k, oc)
        term[i] = ComplexF64(md.FluxesCorrection.Energy.I + md.FluxesCorrection.Energy.H)
    end
    FE[row] = real(sum(term))     # fixed-order reduction ⇒ thread-count invariant
    println("  p=", p, " e=", e, " I=", Ideg, "° m=", m, "  F^E=", FE[row],
            "  (published ", TABLE_I[row][5], ")")
    flush(stdout)
end

pub = [r[5] for r in TABLE_I]
rel = abs.(FE .- pub) ./ abs.(pub)
println("[table I] worst |ΔF^E|/|F^E| = ", maximum(rel), " at row ", argmax(rel))

@repro_result "rows" [string("p=", r[1], ",e=", r[2], ",I=", Int(r[3]), ",m=", r[4]) for r in TABLE_I]
@repro_result "FE_computed"          FE
@repro_result "FE_published"         pub
@repro_result "rel_diff_vs_published" rel
@repro_result "max_rel"              maximum(rel)
@repro_result "FE_all_negative"      all(FE .< 0)
@repro_result "row1_within_0p5pct"   rel[1] < 5e-3
@repro_result "truncation"          "l<=3, |n|<=2, |k|<=2"

exit(ReproHarness.conclude!(opts; paper = "2303.16798", id = "tableI", script = @__FILE__,
                            rtol = 1e-8,
                            tolerances = Dict("max_rel" => (rtol = 1e-6, atol = 1e-12),
                                              "rel_diff_vs_published" => (rtol = 1e-6, atol = 1e-12))))
