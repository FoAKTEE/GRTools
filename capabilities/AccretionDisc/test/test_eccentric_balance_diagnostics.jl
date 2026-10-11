module BalanceDiagnosticsTests
# ECC::gr-balance-diagnostics: the one closure-tagged entry point over the certified adiabatic and
# OPTION A work/torque/AMF owners. This file checks aggregation and tagging only; the physics of every
# returned number is certified in test_eccentric_diagnostics.jl.
using Test
using AccretionDisc, Spacetime
using AccretionDisc.ForcingModes
const SW = AccretionDisc.SpiralWaves
include("eccentric_test_support.jl")

const GENERATED = joinpath(@__DIR__, "../src/SpiralWaves/generated")
const KERNELS = Dict(Adiabatic => ["eccentric_adiabatic.jl"],
                     LocallyIsothermalOptionA => ["eccentric_li.jl", "eccentric_li_stress.jl"])
const TAGS = Dict(Adiabatic => ["adiabatic"], LocallyIsothermalOptionA => ["li", "li_stress"])
const PROBE = Dict("adiabatic" => "_C00", "li" => "_C00", "li_stress" => "_dT_tt")
kernel(tag, name) = getfield(SW, Symbol("_eccentric_", tag, name))

# The header fields, read with plain string operations (independently of the source's regexes).
function header(file)
    l = readlines(joinpath(GENERATED, file))[1:4]
    tag = split(split(l[1], "closure = ")[2], ";")[1]
    _, _, artifact, _, artifact_sha256 = split(l[3])
    (; file, tag, artifact, artifact_sha256, derivation_sha256 = split(l[4])[end])
end

# Production path: one forced (2,1) mode through the eccentric source bank, Float64.
function forced(closure)
    d = SW.DiscModel(; χ=0.9, r_ref=30.0, h_ref=0.08, pSigma=1.5, qTemp=0.75,
                     closure, evolution=:fixed, εp=1e-3)
    pert = SW.EccentricPerturber(EquatorialOrbit(0.9, 30.0, 0.2); mass=1e-5,
                                 softening=FixedSoftening(0.72))
    c = SW.ModeContext(d, pert, 2, 1; closure, domain=(12.0, 75.0))
    r = collect(range(12.0, 75.0; length=129))
    grid = SW.EccentricGrid(r, c.r_c, 1:0, SW.mode_cache_key(c), (; balance_aggregate=true))
    SW.solve_mode(c, pert; grid, bank_options=(; Nq=32))
end

# Two-cell BigFloat160 unit with corotation inside the domain and the Landau shift off, so a
# precision cast or a dropped `landau` setting changes the bits (no production BigFloat run).
function unit(closure)
    d = SW.DiscModel(; χ=big"0.4", r_ref=big"30", h_ref=big"0.08", pSigma=big"1.5", qTemp=big"0.75",
                     closure, evolution=:fixed)
    omega = 2 / (big"29.5"^(3 // 2) + big"0.4")
    c = SW.ModeContext(d, HarmonicMode(2, 1, omega); closure, domain=(big"28", big"32"))
    grid = SW.EccentricGrid(collect(range(c.domain...; length=3)), c.r_c, 1:0, SW.mode_cache_key(c),
                            (; balance_unit=true))
    SW.solve_mode(c; grid, source=eccentric_test_source, landau=false,
                  boundary_rows=SW.mode_boundary_rows(c))
end

function owner_checks(sol, b, closure)
    @test b.closure === closure === sol.closure === sol.context.closure
    @test b.mode === sol.context.mode
    # every number is its certified owner's, bit for bit
    @test b.torque_density == SW.torque_density(sol)
    @test b.power_density == SW.power_density(sol)
    @test b.torque == SW.total_torque(sol) && b.power == SW.total_power(sol)
    @test b.amf == SW.amf(sol)
end

@testset "ECC::gr-balance-diagnostics closure-tagged aggregate" begin
    @test isdefined(SW, :balance_diagnostics)   # registered in the module, not only in this file
    derivations = String[]
    for closure in (Adiabatic, LocallyIsothermalOptionA)
        sol = forced(closure)
        b = SW.balance_diagnostics(sol)
        owner_checks(sol, b, closure)
        # the emitted kernels the closure runs, as recorded in their headers
        @test [e.file for e in b.emitted] == KERNELS[closure]
        @test [e.tag for e in b.emitted] == TAGS[closure]
        for e in b.emitted
            @test e == header(e.file)
            @test e.artifact == "ECC_$(e.tag).jls.txt"
            @test endswith(Base.functionloc(only(methods(kernel(e.tag, PROBE[e.tag]))))[1], "generated/" * e.file)
            push!(derivations, e.derivation_sha256)
        end
        # ... and they are what the closure's runtime dispatch calls
        i = cld(length(sol.r), 2); r = sol.r[i]; h = sol.source(r); y = view(sol.y, :, i)
        inp = SW.master_inputs(sol.context, r; hmodes=h, landau=true, y)
        tag = first(b.emitted).tag
        @test SW._context_operator(Val(closure), inp) == [kernel(tag, "_C$p$q")(inp) for p in 0:1, q in 0:1]
        @test SW.reconstruct_fluid(sol.context, r, y; hmodes=h).deps == kernel(tag, "_deps")(inp)  # closure-dependent
        if closure == LocallyIsothermalOptionA   # the work diagnostics' stress is the emitted li_stress jet
            @test SW.gravitational_stress(sol.context, r, y; hmodes=h).eulerian[1, 1] ==
                  kernel("li_stress", "_dT_tt")(inp)
        end
        println("PASS balance aggregate ", closure, ": torque ", b.torque, ", power ", b.power,
                ", emitted ", join((e.tag for e in b.emitted), "+"))
        setprecision(BigFloat, 160) do
            u = unit(closure)
            @test u.context.r_c !== nothing && !u.params.landau
            bu = SW.balance_diagnostics(u)
            owner_checks(u, bu, closure)
            @test eltype(bu.torque_density) == eltype(bu.amf) == BigFloat
            @test bu.torque isa BigFloat && precision(bu.torque) == precision(bu.power) == 160
            @test bu.emitted === b.emitted
        end
        println("PASS balance aggregate ", closure, ": BigFloat160 two-cell unit, landau=false, corotation inside")
    end
    # both closures are emitted from one R3 derivation state
    @test length(derivations) == 3 && allequal(derivations)
end

end # module
