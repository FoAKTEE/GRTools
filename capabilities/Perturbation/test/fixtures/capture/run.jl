# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/npghp_moveD/TeukEvolution_run/minrun/run.jl (untracked one-shot, mtime 2026-07-04)
#               sha256(original) = ec1cec1e58674ab45a0374d9df0d946908267790c6fde7520095f00c997137e9
# fixture(s)  : td_hphc_upstream.jl (via build_fixture.jl)
# oracle      : ref-code/TeukEvolution.jl @ 3c8e47f — minimal linear psi4 (s=-2) Kerr a=0.7 Gaussian-ID run (nx=48, ny=24, nt=600)
# run         : julia --project=<TeukEvolution env> test/fixtures/capture/run.jl   # writes out_min/ CSVs
#
# Paths rewritten to the ENV-overridable GRTOOLS_TEUKEVOLUTION_ROOT / GRTOOLS_REFCODE_ROOT
# (defaults reproduce the capture host).  Body otherwise VERBATIM from the original one-shot.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const TEUK_ROOT = get(ENV, "GRTOOLS_TEUKEVOLUTION_ROOT",
                      joinpath(get(ENV, "GRTOOLS_REFCODE_ROOT", "/data/haiyangw/claude/emri-disk/ref-code"),
                               "TeukEvolution.jl"))
# ----------------------------------------------------------------------------------------------

# Minimal Move-D characterization run: linear psi4 (spin -2), Kerr a=0.7, Gaussian ID.
# Self-contained: cd(TEUK_ROOT) so launch's `cp("./run.jl",...)` and relative outdir resolve here.
cd(TEUK_ROOT)

include(joinpath(TEUK_ROOT, "src") * "/TeukEvolution.jl")
import .TeukEvolution as TE

params = Dict(
    "outdir" => "out_min",
    "nx" => 48,      # radial grid points (>=6 for FD boundary stencils)
    "ny" => 24,      # angular collocation points (>20 required by Sphere.num_l)
    "nt" => 600,     # number of time steps
    "ts" => 200,     # save every ts steps -> snapshots at tc=200,400,600
    "psi_spin" => -2,
    "id_kind" => "gaussian",
    "runtype" => "linear_field",
    "m_vals" => [2],
    "id_l_ang" => [2],
    "id_ru" => [4.0],
    "id_rl" => [2.5],
    "id_width" => [1.5],
    "ingoing" => true,
    "id_amp" => [[0.1, 0.0]],   # [re, im] per m
    "cl" => 1.0,
    "cfl" => 0.5,
    "bhs" => 0.7,               # a = 0.7
    "bhm" => 1,
    "save_l" => false,
    "l_to_save" => 2,
    "precision" => Float64,
)

@time TE.launch(params)
