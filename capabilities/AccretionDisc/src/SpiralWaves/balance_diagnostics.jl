# ECC::gr-balance-diagnostics — the ONE closure-tagged entry point for the work/torque/AMF diagnostics of
# a solved two-index mode, for both primary closures. A thin aggregate: every number comes from its
# certified owner, diagnostics.jl (ECC::gr-adiabatic-balance-diagnostics) or li_diagnostics.jl
# (ECC::gr-li-balance-diagnostics), selected by the solution's closure through the same `Val` dispatch
# those files use. Nothing is derived, weighted or integrated here.

# Header of a generated kernel file, as `emit_eccentric` (derivation/src/eccentric_emission.jl) writes it:
# the emitted closure tag, the expression artifact and its sha256, and the R3 derivation sha256. Read
# when this file is compiled, so the record describes the compiled kernels (the generated files are
# `include`d, so editing one recompiles this file as well).
function _emission_header(file)
    head = open(io -> [readline(io) for _ in 1:4], joinpath(@__DIR__, "generated", file))
    tag = match(r"^# GENERATED: closure = (\w+); do not edit by hand\.$", head[1])
    input = match(r"^# input: (\S+) sha256 ([0-9a-f]{64})$", head[3])
    derivation = match(r"^# derivation: sha256 ([0-9a-f]{64})$", head[4])
    any(isnothing, (tag, input, derivation)) && error("generated/$file: not an emit_eccentric header")
    (; file, tag = String(tag[1]), artifact = String(input[1]), artifact_sha256 = String(input[2]),
       derivation_sha256 = String(derivation[1]))
end

# The generated kernels that each closure's solution and diagnostics run: the adiabatic operator, source
# and reconstruction (master.jl); the OPTION A operator, source and reconstruction, and its stress jets
# (locally_isothermal.jl).
const _ADIABATIC_EMISSION = (_emission_header("eccentric_adiabatic.jl"),)
const _LI_EMISSION = (_emission_header("eccentric_li.jl"), _emission_header("eccentric_li_stress.jl"))
_emission(::Val{Adiabatic}) = _ADIABATIC_EMISSION
_emission(::Val{LocallyIsothermalOptionA}) = _LI_EMISSION

"""
    balance_diagnostics(sol::EccentricModeSolution) -> NamedTuple

Work, torque and angular-momentum-flux diagnostics of a solved two-index mode, tagged with the closure
that produced it, for the adiabatic and the locally isothermal OPTION A closures alike:

  * `closure`, `mode`: the solution's closure tag (taken from its validated context) and harmonic;
  * `emitted`: one `(file, tag, artifact, artifact_sha256, derivation_sha256)` record per generated
    kernel file that closure runs, from the file headers at compile time;
  * `torque_density`, `power_density`: direct gravitational work at every node
    ([`gravitational_work`](@ref)), positive into the disc;
  * `torque`, `power`: their integrals over the grid ([`total_torque`](@ref), [`total_power`](@ref));
  * `amf`: the closure's wave angular-momentum flux at every node ([`amf`](@ref)), the adiabatic (72a)
    flux or the OPTION A connected Reynolds wave flux (`li_flux_decomposition(...).wave`).

For OPTION A the edge difference of `amf` is not a Lindblad torque; the wave/mean-transport split is in
[`li_flux_decomposition`](@ref). Precision follows the solution.
"""
function balance_diagnostics(sol::EccentricModeSolution)
    work = gravitational_work(sol)
    torque_density = getproperty.(work, :torque)
    power_density = getproperty.(work, :power)
    (; closure = sol.closure, mode = sol.context.mode, emitted = _emission(Val(sol.closure)),
       torque_density, power_density,
       torque = total_torque(sol; density = torque_density),
       power = total_power(sol; density = power_density),
       amf = amf(sol))
end
