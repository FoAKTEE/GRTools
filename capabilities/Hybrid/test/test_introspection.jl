# Characterization of model introspection (Hybrid.Introspection, §4.9).
#
# GROUNDING (no WaSABI live dep — literal asset metadata + grounded forcing/kinematics + filled
# forcing-driven ODE assembly self-consistent against the grounded HybridC/Models primitives):
#   * cite_model reads the verbatim "Citation" bibtex from the shipped MetadataModels/<m>.m asset
#     (escapes resolved) — checked against literal substrings of 0PASchwarz.m.
#   * list_metadata_models == the 7 asset FileBaseNames.
#   * initial_conditions == the literal _ICS_FORMATS lists; error on un-ported.
#   * forcing_terms: per-model ODE structure + rhs evaluator == omega_dot / kepler_Ω (grounded).
#   * list_modes == sorted amp keys.
#   * get_inspiral_equations / inspiral_evaluate (filled H-WASABI-STUBS): the WaSABI-C 9-var
#     forcing-driven rhs at M=1 == the grounded wc_forcing channels (oracle16b "FCH" literals);
#     inspiral_evaluate(solution,t) == wc_state_at; rhs == finite-difference dω/dt self-consistency.

using Test
using Hybrid.Introspection
using Hybrid.Models: omega_dot, Forcing0PASchwarz
using Hybrid.HybridC: int_inspiral, duration, wc_state_at, WaSABICModel
using Spacetime: kepler_Ω

# Baked fixture copy of the upstream WaSABI data (provenance: test/fixtures/wasabi/README.md).
const _IN_REF = joinpath(@__DIR__, "fixtures", "wasabi")
const _IN_INSP = joinpath(_IN_REF, "InspiralModels", "sf_data")

@testset "Introspection" begin
    @testset "cite_model literal asset metadata" begin
        c = cite_model("0PASchwarz")
        @test startswith(c, "@article{Wardell:2021fyy,")
        @test occursin("eprint = \"2112.12265\"", c)      # escapes resolved to real quotes
        @test occursin("year = \"2023\"", c)
        @test !occursin("\\\"", c)                         # no leftover \" escapes
        @test_throws ErrorException cite_model("nonsense")
        @test list_metadata_models() ==
              sort(["0PAKerr", "0PASchwarz", "1PAT1", "1PAT1R", "1PAT1ea", "WaSABI-C", "WaSABI-C_v0.9"])
    end

    @testset "initial_conditions" begin
        @test initial_conditions("0PASchwarz") == ["M", "r0", "ν", "φ"]
        @test initial_conditions("1PAT1") == ["M", "r0", "ν", "φ"]
        @test initial_conditions("0PAKerr") == ["ω", "φ", "ν", "M", "χ1"]
        @test_throws ErrorException initial_conditions("WaSABI-C")
        @test_throws ErrorException initial_conditions("nonsense")
    end

    @testset "forcing_terms symbolic structure + rhs evaluator" begin
        rs = collect(range(6.0, 32.0; length = 16))
        fI(r) = 0.02 + 1.0e-3 * r - 3.0e-5 * r^2; fH(r) = 0.005 - 2.0e-4 * r
        sf = mkpath(joinpath(mktempdir(), "1SF_Flux", "Schwarz_Circ"))
        _wl(p, ys) = write(p, "{" * join(("{$(r), $(y)}" for (r, y) in zip(rs, ys)), ", ") * "}")
        _wl(joinpath(sf, "1SFCircShwarzDotEInf.m"), fI.(rs))
        _wl(joinpath(sf, "1SFCircShwarzDotEHorizon.m"), fH.(rs))
        sfroot = dirname(dirname(sf))

        ftm = forcing_terms("0PASchwarz", sfroot)
        @test ftm isa ForcingTermsModel
        @test ftm.odelhs == ["dΩ/dt", "dφ/dt", "dν/dt", "dM/dt"]
        @test ftm.aelhs == ["Ω"]
        @test ftm.parameters == ["r0", "φ", "ν", "M"]
        @test ftm.forcing isa Forcing0PASchwarz

        r0, ν, M = 10.0, 1.0e-3, 1.0
        out = ftm(Dict("r0" => r0, "φ" => 0.0, "ν" => ν, "M" => M))
        @test out["dΩ/dt"] == omega_dot(ftm.forcing, r0, ν, M)
        @test out["dφ/dt"] == kepler_Ω(r0) / M
        @test out["dν/dt"] == 0.0 && out["dM/dt"] == 0.0

        @test_throws ErrorException forcing_terms("nonsense", sfroot)   # unknown model
    end

    @testset "list_modes == sorted amp keys" begin
        rs = collect(range(6.0, 32.0; length = 16))
        amp = mkpath(joinpath(mktempdir(), "1SF", "Schwarz_Circ"))
        zf(l, m, r) = complex(0.01 + 0.001r, -0.002 + 5e-4r)
        for l in 2:5, mm in 1:l
            vals = [zf(l, mm, r) for r in rs]
            write(joinpath(amp, "1SFTeukampSchwarzCirc$(l)$(mm).m"),
                  "{" * join(("{$(r), $(real(z)) + $(imag(z))*I}" for (r, z) in zip(rs, vals)), ", ") * "}")
        end
        @test list_modes("0PASchwarz", dirname(dirname(amp))) == [(l, m) for l in 2:5 for m in 1:l]
    end

    @testset "get_inspiral_equations / inspiral_evaluate (filled H-WASABI-STUBS)" begin
        # The forcing-driven ODE system == the forcing-terms structure (AHA).
        eqs = get_inspiral_equations("WaSABI-C", _IN_INSP)
        @test eqs isa ForcingTermsModel
        @test eqs isa InspiralEquations
        @test eqs.forcing isa WaSABICModel
        @test eqs.odelhs == ["dω/dt", "dφ/dt", "dν/dt", "dM/dt", "dχ1/dt", "dχ2/dt",
                             "dδm/dt", "dδν/dt", "dδχ/dt"]
        @test eqs.parameters == ["ω", "φ", "ν", "M", "χ1", "χ2", "δm", "δν", "δχ"]

        # rhs at the M=1 state == the grounded wc_forcing channels (oracle16b "FCH" literals; the
        # du = F/M rates equal the raw forcing channels at M=1, except dM/dt=Fm which has no 1/M).
        st = Dict("ω" => 0.02, "φ" => 0.0, "ν" => 0.1, "M" => 1.0, "χ1" => 0.5, "χ2" => 0.3,
                  "δm" => 0.1, "δν" => 0.05, "δχ" => 0.2)
        du = eqs(st)
        @test isapprox(du["dω/dt"], 1.1270202541171997e-6; rtol = 1e-9)
        @test isapprox(du["dφ/dt"], 0.02; rtol = 1e-12)
        @test isapprox(du["dν/dt"], 2.681105260055963e-12; rtol = 1e-9, atol = 1e-24)
        @test isapprox(du["dM/dt"], -3.096924505082267e-11; rtol = 1e-9, atol = 1e-24)
        @test isapprox(du["dχ1/dt"], -1.942598793426137e-9; rtol = 1e-9, atol = 1e-24)
        @test isapprox(du["dχ2/dt"], -1.1385680120862304e-10; rtol = 1e-9, atol = 1e-24)
        @test isapprox(du["dδm/dt"], -3.1237355576828266e-10; rtol = 1e-9, atol = 1e-24)
        @test isapprox(du["dδν/dt"], 2.654294207455403e-10; rtol = 1e-9, atol = 1e-24)
        @test isapprox(du["dδχ/dt"], -1.9431350144781482e-8; rtol = 1e-9, atol = 1e-24)
        # inspiral_evaluate(model-eqs, state) == the rhs dict.
        @test inspiral_evaluate(eqs, st) == du

        # inspiral_evaluate(solution, t) == wc_state_at; rhs == finite-difference dω/dt along the
        # integrated trajectory (physical self-consistency vs the grounded HybridC int_inspiral).
        s = int_inspiral(eqs.forcing, (; ω = 0.02, ν = 0.1, χ1 = 0.5, χ2 = 0.3);
                         reltol = 1e-11, abstol = 1e-11)
        t0 = 0.3 * duration(s)
        @test inspiral_evaluate(s, t0) == wc_state_at(s, t0)
        u = wc_state_at(s, t0)
        stt = Dict("ω" => u.ω, "φ" => u.φ, "ν" => u.ν, "M" => u.M, "χ1" => u.χ1, "χ2" => u.χ2,
                   "δm" => u.δm, "δν" => u.δν, "δχ" => u.δχ)
        dt = 1e-3 * duration(s)
        fd = (wc_state_at(s, t0 + dt).ω - wc_state_at(s, t0 - dt).ω) / (2dt)
        @test isapprox(eqs(stt)["dω/dt"], fd; rtol = 1e-5)

        # MethodError for an unsupported solution type (typed dispatch, never a silent stub).
        @test_throws MethodError inspiral_evaluate(nothing, 0.0)
    end
end
