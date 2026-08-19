# arXiv:2303.16798 Table II (tab:spherical, paper.tex L850-865): the linear-in-spin total energy and
# angular-momentum fluxes F^E_{S,m}, F^{Jz}_{S,m} from nearly SPHERICAL orbits (e = 0) at a = 0.9 M,
# p = 10, for I ∈ {30°, 60°, 120°, 150°} and m ∈ {1, 2, 3}, summed over l = max(2,m)..5 and k = -4..4.
#
# THIN driver: Perturbation.SpinFluxes supplies the spherical orbit correction, its r/x derivatives
# and every mode (`KerrSpinOrbitCorrectionSpherical`, `_orbit_derivatives_spherical`,
# `TeukolskySpinModeSphericalCorrection`); the script only sums and compares.
#
# PAPER TYPO (carried from the legacy lock, not silently fixed): the I = 150°, m = 2 row prints
# F^{Jz} = 1.658e-5 but the value is 1.658e-4 — that row's published entry is corrected here and
# both the printed and the corrected numbers are registered in the golden.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2303.16798/check_tableII.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Perturbation: KerrSpinOrbitCorrectionSpherical, _orbit_derivatives_spherical,
                            TeukolskySpinModeSphericalCorrection

const A, P = 0.9, 10.0
const INCLINATIONS = (30.0, 60.0, 120.0, 150.0)
const MS = (1, 2, 3)
const LSMAX = 5          # l runs max(2,m)..LSMAX (|m| ≤ l)
const KS = -4:4

# published (I°, m, F^E, F^Jz-as-printed, F^Jz-corrected)
const TABLE_II = [
    ( 30.0, 1, -2.642e-7, -2.446e-6, -2.446e-6),
    ( 30.0, 2, -2.702e-6, -6.431e-5, -6.431e-5),
    ( 30.0, 3, -3.921e-7, -1.016e-5, -1.016e-5),
    ( 60.0, 1, -1.533e-6, -1.891e-5, -1.891e-5),
    ( 60.0, 2, -2.177e-6, -5.110e-5, -5.110e-5),
    ( 60.0, 3, -2.223e-7, -5.463e-6, -5.463e-6),
    (120.0, 1, -4.175e-6,  3.021e-5,  3.021e-5),
    (120.0, 2, -1.796e-6,  3.597e-5,  3.597e-5),
    (120.0, 3, -1.730e-7,  4.020e-6,  4.020e-6),
    (150.0, 1, -2.859e-6,  3.280e-5,  3.280e-5),
    (150.0, 2, -6.930e-6,  1.658e-5,  1.658e-4),   # printed 1.658e-5 — PAPER TYPO, is 1.658e-4
    (150.0, 3, -1.069e-6,  2.723e-5,  2.723e-5),
]

FE = zeros(length(TABLE_II))
FJ = zeros(length(TABLE_II))

for (iI, Ideg) in enumerate(INCLINATIONS)
    x   = cos(Ideg * pi / 180)
    oc  = KerrSpinOrbitCorrectionSpherical(A, P, 0.0, x, 8)
    odv = _orbit_derivatives_spherical(A, P, 0.0, x, 8; h = 1e-5)
    for (im, m) in enumerate(MS)
        row = (iI - 1) * length(MS) + im
        # per-mode into distinct slots, then a fixed-order reduction ⇒ thread-count invariant
        modes = [(l, k) for l in max(2, m):LSMAX for k in KS]
        eterm = Vector{ComplexF64}(undef, length(modes))
        jterm = Vector{ComplexF64}(undef, length(modes))
        Threads.@threads for i in eachindex(modes)
            l, k = modes[i]
            md = TeukolskySpinModeSphericalCorrection(l, m, k, oc, odv)
            eterm[i] = ComplexF64(md.FluxesCorrection.Energy.I + md.FluxesCorrection.Energy.H)
            jterm[i] = ComplexF64(md.FluxesCorrection.AngularMomentum.I +
                                  md.FluxesCorrection.AngularMomentum.H)
        end
        FE[row] = real(sum(eterm))
        FJ[row] = real(sum(jterm))
        println("  I=", Ideg, "° m=", m, "  F^E=", FE[row], "  F^Jz=", FJ[row])
        flush(stdout)
    end
end

pubE = [r[3] for r in TABLE_II]
pubJ = [r[5] for r in TABLE_II]          # typo-corrected column
relE = abs.(FE .- pubE) ./ abs.(pubE)
relJ = abs.(FJ .- pubJ) ./ abs.(pubJ)

println("[table II] worst |ΔF^E|/|F^E|   = ", maximum(relE), " at row ", argmax(relE))
println("[table II] worst |ΔF^Jz|/|F^Jz| = ", maximum(relJ), " at row ", argmax(relJ))

@repro_result "rows"             [string(Int(r[1]), "deg,m=", r[2]) for r in TABLE_II]
@repro_result "FE_computed"      FE
@repro_result "FJz_computed"     FJ
@repro_result "FE_published"     pubE
@repro_result "FJz_published_corrected" pubJ
@repro_result "FJz_published_as_printed" [r[4] for r in TABLE_II]
@repro_result "max_rel_FE"       maximum(relE)
@repro_result "max_rel_FJz"      maximum(relJ)
# the paper's qualitative claims: F^E < 0 always, and F^Jz flips sign prograde → retrograde
@repro_result "FE_all_negative"  all(FE .< 0)
@repro_result "FJz_sign_flip_prograde_to_retrograde" all(FJ[1:6] .< 0) && all(FJ[7:12] .> 0)

exit(ReproHarness.conclude!(opts; paper = "2303.16798", id = "tableII", script = @__FILE__,
                            rtol = 1e-8,
                            tolerances = Dict("max_rel_FE"  => (rtol = 1e-6, atol = 1e-12),
                                              "max_rel_FJz" => (rtol = 1e-6, atol = 1e-12))))
