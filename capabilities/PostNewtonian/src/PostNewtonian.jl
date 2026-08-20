"""
    PostNewtonian

PN-series engine for Teukolsky-mode post-Newtonian expansions and PN self-force — the GrTools
**Symbolics quarantine** (design §4.6). Fuses the two legacy PN leaf packages (TeukolskyPN +
PostNewtonianSelfForce) into one capability. `Symbolics` is a HARD dep here and ONLY here, because
the Laurent/Taylor coefficient vector and the resummed self-force scalar are `Symbolics.Num` struct
field types; the capability is a top-of-DAG leaf so the Symbolics precompile floor never leaks onto
the numeric/GPU waveform path.

Submodules:
- `Series`      — the `PNSeries` Num-algebra engine: Laurent/Taylor arithmetic, `Normal`
                  materialization, change-of-parameter, power-counting, `expand_log`/`expand_gamma`.
- `Catalog`     — PN `.m` discovery / substring browse (`SERIES_DATA_DIR`, `pn_series_catalog`,
                  `pn_browse`).
- `MMAInputForm`— the unified `read_mma_inputform_symbolic` MMA-InputForm→`Num` symbolic walker (Association/
                  SeriesData/Module text → `PNSeries`/`Dict`, no eval/world-age).
- `SelfForce`   — `PostNewtonianData` accessor, `postnewtonian_coefficient`, `ResummedSeriesData`,
                  `postnewtonian_expansion`.
- `MST`        — symbolic Mano-Suzuki-Takasugi coefficient engine: truncated-`ε` Laurent algebra,
                  the renormalized-angular-momentum `ν`-solve, and the `aₙ` `ε`-series lifted to
                  `PNSeries` (`η`-form and `γ`-FreqRep form).
- `RadialMode` — η-graded symbolic radial Coulomb/MST solution (`TeukolskyRadialPN` In/Up pair) and
                  the asymptotic mode amplitudes, from `MST`'s `ε`-series.
- `SourceMode` — point-particle source δ-coefficient series, sourced mode assembly, and the
                  asymptotic energy flux (`TeukolskyPointParticleModePN`).
"""
module PostNewtonian

# Submodule includes (built tests-first, grounded by characterization vs the legacy byte-exact
# .m oracle / analytic PN identities).
include("Series.jl")
include("Catalog.jl")
include("MMAInputForm.jl")
include("SelfForce.jl")
include("MST.jl")
include("RadialMode.jl")
include("SourceMode.jl")

using .Series
using .Catalog
using .MMAInputForm
using .SelfForce
using .MST
using .RadialMode
using .SourceMode

end # module PostNewtonian
