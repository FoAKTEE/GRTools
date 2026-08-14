# Characterization of the model registry (Hybrid.Registry, §4.9).
#
# GROUNDING (no WaSABI live dep — on-disk plugin discovery + the "not-yet-ported" contract +
# a synthetic-fixture ctor round-trip):
#   * _MODEL_CATALOG == the sorted 7 == FileBaseName of the baked fixture copy of upstream
#     WaSABI/{InspiralModels,AmplitudeModels}/*.m (the upstream FileNames["*.m"] discovery;
#     provenance: test/fixtures/wasabi/README.md).
#   * list_inspiral_models/list_amplitude_models/list_models all == the catalogue.
#   * get_inspiral/get_amplitudes error on an unknown name and on the catalogued-but-not-ported
#     WaSABI-C / WaSABI-C_v0.9 (the explicit porting contract).
#   * a wired ctor (0PASchwarz) builds the right model type from a synthetic datadir round-trip.

using Test
using Hybrid.Registry
using Hybrid.Models: Forcing0PASchwarz, Amp0PASchwarz, flux_1SF, Z_I
using Hybrid.HybridC: WaSABICModel, int_inspiral, duration, wc_state_at

# FileBaseName of *.m in a directory (the upstream FileNames["*.m", dir] discovery).
_basenames(dir) = sort([first(splitext(f)) for f in readdir(dir) if endswith(f, ".m")])
# Baked fixture copy of the upstream WaSABI tree (provenance: test/fixtures/wasabi/README.md).
_refdir(sub) = joinpath(@__DIR__, "fixtures", "wasabi", sub)

@testset "Registry" begin
    catalog = sort(["0PAKerr", "0PASchwarz", "1PAT1", "1PAT1R", "1PAT1ea",
                    "WaSABI-C", "WaSABI-C_v0.9"])

    @testset "catalogue == on-disk plugin discovery" begin
        @test list_inspiral_models() == catalog
        @test list_amplitude_models() == catalog
        @test list_models() == catalog
        @test length(catalog) == 7
        @test _basenames(_refdir("InspiralModels")) == catalog
        @test _basenames(_refdir("AmplitudeModels")) == catalog
    end

    @testset "error contract (unknown + amplitude-side not-yet-ported)" begin
        @test_throws ErrorException get_inspiral("nonsense", "/tmp")
        @test_throws ErrorException get_amplitudes("nonsense", "/tmp")
        # WaSABI-C / WaSABI-C_v0.9 inspiral side is now wired to HybridC.wasabic_model; only the
        # amplitude side remains the upstream-undefined contract.
        for c in ("WaSABI-C", "WaSABI-C_v0.9")
            @test_throws ErrorException get_amplitudes(c, "/tmp")
        end
    end

    @testset "WaSABI-C inspiral wired to HybridC.wasabic_model" begin
        # get_inspiral("WaSABI-C", sf) builds a WaSABICModel from the inspiral root + sibling
        # AmplitudeModels/sf_amp_data; its int_inspiral reproduces the module-3 oracle16b trajectory.
        insp = _refdir(joinpath("InspiralModels", "sf_data"))
        for name in ("WaSABI-C", "WaSABI-C_v0.9")
            m = get_inspiral(name, insp)
            @test m isa WaSABICModel
            s = int_inspiral(m, (; ω = 0.02, ν = 0.1, χ1 = 0.5, χ2 = 0.3);
                             reltol = 1e-11, abstol = 1e-11)
            @test isapprox(duration(s), 6173.477487710508; rtol = 1e-6)
            u = wc_state_at(s, 2057.825829236836)
            @test isapprox(u.ω, 0.022960044450056023; rtol = 1e-6)
            @test isapprox(u.φ, 43.94133361266952; rtol = 1e-6)
        end
    end

    @testset "wired ctor round-trip (synthetic 0PASchwarz fixtures)" begin
        fI(r) = 0.02 + 1.0e-3 * r - 5.0e-5 * r^2
        fH(r) = 0.005 - 2.0e-4 * r
        rs = collect(range(6.2, 32.0; length = 14))
        _wl(p, ys) = write(p, "{" * join(("{$(r), $(y)}" for (r, y) in zip(rs, ys)), ", ") * "}")

        # Inspiral: sf_datadir/1SF_Flux/Schwarz_Circ/...
        sf = mkpath(joinpath(mktempdir(), "1SF_Flux", "Schwarz_Circ"))
        _wl(joinpath(sf, "1SFCircShwarzDotEInf.m"), fI.(rs))
        _wl(joinpath(sf, "1SFCircShwarzDotEHorizon.m"), fH.(rs))
        m = get_inspiral("0PASchwarz", dirname(dirname(sf)))
        @test m isa Forcing0PASchwarz
        @test flux_1SF(m, 15.0) ≈ fI(15.0) + fH(15.0) rtol = 1e-10

        # Amplitudes: amp_datadir/1SF/Schwarz_Circ/...
        amp = mkpath(joinpath(mktempdir(), "1SF", "Schwarz_Circ"))
        zf(l, m, r) = complex(0.01 + 0.001r + 1e-4(l + m), -0.002 + 5e-4r)
        for l in 2:5, mm in 1:l
            vals = [zf(l, mm, r) for r in rs]
            write(joinpath(amp, "1SFTeukampSchwarzCirc$(l)$(mm).m"),
                  "{" * join(("{$(r), $(real(z)) + $(imag(z))*I}" for (r, z) in zip(rs, vals)), ", ") * "}")
        end
        a = get_amplitudes("0PASchwarz", dirname(dirname(amp)))
        @test a isa Amp0PASchwarz
        @test Z_I(a, 2, 2, 15.0) ≈ zf(2, 2, 15.0) rtol = 1e-9
    end
end
