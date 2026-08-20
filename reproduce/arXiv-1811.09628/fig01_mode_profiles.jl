# arXiv:1811.09628 (Miranda & Rafikov 2019) Fig. 1 (`mode_profiles`) — the fractional surface-density
# perturbation profiles δΣ_m/Σ of the planet-forced linear modes m = 2, 5, 10, in the figure's unit
# (M_p/M_th)^{-1} δΣ_m/Σ = h₀³ × (the solver's per-unit-(M_p/M_*) value).
#
# THIN driver: `AccretionDisc.Modes` supplies the disc (`NewtonianDisc`), the grid (`loggrid`), the
# mode solve (`solve_modes`) and the reconstruction (`vars_from_enthalpy`); this script only scales,
# reduces and registers.
#
# PHYSICS REGIME (a measured finding of the legacy reproduction, carried here): MR19 is the ADIABATIC,
# NO-COOLING linear theory, i.e. the limit β = Ωt_c → ∞. The cooling extension's shipped
# `parameters.py` value β = 0.1 damps the outward-propagating waves and destroys the profile ordering
# of this figure; β = 80 already realises the adiabatic limit (β ≳ 10³ is converged).
#
# RESOLUTION: 20 000 radial points (the published figure used ~10⁵ — a denser red fill at high m,
# cosmetic). The golden stores the profiles decimated to 400 points per mode plus the envelope
# amplitudes the figure is read off.
#
# tier HEAVY: three forced solves on a 20 000-point grid, ≈ 8 min at -t4.
#
#   julia --project=reproduce -t4 reproduce/arXiv-1811.09628/fig01_mode_profiles.jl --check --tier heavy

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

using GrTools
using GrTools.AccretionDisc: NewtonianDisc, loggrid, solve_modes, vars_from_enthalpy
const ADM = GrTools.AccretionDisc.Modes
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))

const H0 = 0.1
const BETA = 80.0                 # adiabatic limit (MR19 has no cooling)
const NPTS = 20_000
const MS = (2, 5, 10)
const SCALE = H0^3                # M_th/M_* with M_* = 1: the figure's (M_p/M_th)^{-1} unit

d = NewtonianDisc(h0 = H0, q = 1.0, p = 1.0, gamma = 7 / 5, beta = BETA, soft = 0.6 * H0,
                  indirect = 0, pole_disp = 1e-6)
x = loggrid(d, 0.05, 5.0, NPTS)
println("[repro] grid: $(length(x)) pts, r ∈ [$(x[1]), $(x[end])], beta=$(BETA)"); flush(stdout)

t0 = time()
modes = solve_modes(d, collect(MS), x)
println("[repro] solve_modes wall=$(round(time()-t0, digits=1)) s"); flush(stdout)

@repro_result "npts" NPTS
@repro_result "beta" BETA
@repro_result "scale_Mth" SCALE
idx = 1:div(NPTS, 400):NPTS
@repro_result "r" collect(x[idx])

fracs = Dict{Int,Vector{ComplexF64}}()
for m in MS
    _, _, sig = vars_from_enthalpy(d, modes[m], x, m)
    frac = ComplexF64[SCALE * sig[i] / ADM.Sigma(d, x[i]) for i in eachindex(x)]
    fracs[m] = frac
    @repro_result "re_m$(m)" real.(frac[idx])
    @repro_result "im_m$(m)" imag.(frac[idx])
    amp = abs.(frac)
    env(lo, hi) = maximum(amp[i] for i in eachindex(x) if lo < x[i] < hi; init = 0.0)
    @repro_result "env_interior_m$(m)" env(0.3, 0.9)
    @repro_result "env_exterior_m$(m)" env(1.2, 2.5)
    @repro_result "env_far_m$(m)" env(3.0, 4.7)
    @repro_result "amp_max_m$(m)" maximum(amp)
    @repro_result "r_at_amp_max_m$(m)" x[argmax(amp)]
end

# ---- legacy-campaign anchors (emri-disk results/emri-disk/figures/1811.09628/repro_fig1.spec.md,
# "Rule 3 — VLM COMPARISON", tracked): the same quantities measured by the legacy DiscModes run on
# the same 20 000-point grid, in the SAME unscaled code units the spec quotes (× SCALE here).
const LEGACY_ENV_INTERIOR = Dict(2 => 36.3, 5 => 56.7, 10 => 37.8)
const LEGACY_AMP_MAX      = Dict(2 => 269.5, 5 => 274.8, 10 => 181.6)
for m in MS
    @repro_result "legacy_env_interior_m$(m)" SCALE * LEGACY_ENV_INTERIOR[m]
    @repro_result "legacy_amp_max_m$(m)" SCALE * LEGACY_AMP_MAX[m]
end

# ---- figure: the published 3-panel stack, Re black / Im red ---------------------------------------
if opts.mode === :regen
    emri_defaults!()
    plt = plot(layout = (3, 1), size = (560, 640), link = :x, legend = false,
               bottom_margin = 0Plots.mm, left_margin = 5Plots.mm)
    for (k, m) in enumerate(MS)
        plot!(plt[k], x, real.(fracs[m]), c = :black, lw = 0.6, xscale = :log10,
              xlim = (0.05, 5.0), ylim = (-0.4, 0.4), yticks = [-0.4, -0.2, 0.0, 0.2, 0.4])
        plot!(plt[k], x, imag.(fracs[m]), c = :red, lw = 0.6)
        vline!(plt[k], [1.0], c = :black, lw = 0.8)
        annotate!(plt[k], 0.06, 0.32, text("($("abc"[k])) m = $(m)", :left, 11))
        if k == 1                                   # manual Re/Im legend, as in the published panel (a)
            plot!(plt[k], [0.115, 0.16], [0.08, 0.08], c = :black, lw = 1.2)
            plot!(plt[k], [0.115, 0.16], [-0.02, -0.02], c = :red, lw = 1.2)
            annotate!(plt[k], 0.17, 0.08, text("Re", :left, 9))
            annotate!(plt[k], 0.17, -0.02, text("Im", :left, 9))
        end
        if k < length(MS)
            plot!(plt[k], xformatter = _ -> "")
        else
            plot!(plt[k], xticks = ([0.1, 0.3, 1.0, 3.0], ["0.1", "0.3", "1.0", "3.0"]),
                  xlabel = "r / r_p")
        end
    end
    plot!(plt[2], ylabel = "(M_p/M_th)^-1 dSigma_m / Sigma")
    png = joinpath(@__DIR__, "fig01_mode_profiles.png")
    savefig(plt, png)
    println("[repro] wrote ", png, " (", filesize(png), " bytes)")
end

exit(ReproHarness.conclude!(opts; paper = "1811.09628", id = "fig01_mode_profiles",
                            script = @__FILE__, rtol = 1e-7))
