# Module 1 — Modes/Background.jl. Kerr equatorial thin-disc background (Dyson & D'Orazio
# arXiv:2601.19123, M=1 BL). GROUNDING:
#   • disc_OmegaK & kappar2 — CHARACTERIZED vs KerrGeodesics.KerrGeoFrequencies
#     (KerrGeoFrequencies(a,r,0,1).Omegaph == disc_OmegaK rel<1e-12 ; .Omegar^2 == kappar2
#     rel<1e-10), verbatim from DiscModes/test/verify_rel_disc.jl over a∈{0,0.6,0.9}, r∈{8,20,50}.
#     The KerrGeodesics oracle (Omegaph, Omegar) is BAKED as Float64 fixtures (kg_freq below),
#     captured once via `julia --project=packages/KerrGeodesics -e
#     'using KerrGeodesics; KerrGeoFrequencies(a,r,0.0,1.0)'` — legacy no longer needed to test.
#   • risco_disc — Bardeen 1972 prograde ISCO: a=0 ⇒ 6 exactly; a=0.6/0.9 vs the DiscModes
#     oracle captured fixture (== Bardeen closed form).
#   • csound / eden0 / delta_r — captured one-shot fixtures from the DiscModes oracle
#     (packages/DiscModes, RelativisticDisc.csound/eden0/Delta_r) to full Float64 precision.

using Test
using AccretionDisc.Modes

# Baked KerrGeodesics.KerrGeoFrequencies(a, r, 0.0, 1.0) oracle (Omegaph, Omegar), captured once
# from `julia --project=packages/KerrGeodesics -e 'using KerrGeodesics; KerrGeoFrequencies(...)'`.
# (a, r) => (Omegaph, Omegar). Full Float64 precision preserved (repr round-trip).
const kg_freq = Dict(
    (0.0, 8.0)  => (0.044194173824159216, 0.022097086912079608),
    (0.0, 20.0) => (0.011180339887498949, 0.009354143466934852),
    (0.0, 50.0) => (0.0028284271247461896, 0.0026532998322843196),
    (0.6, 8.0)  => (0.04305257016255476, 0.028727939522846697),
    (0.6, 20.0) => (0.01110583965028321, 0.009624128865677045),
    (0.6, 50.0) => (0.0028236352568157294, 0.002668513768756883),
    (0.9, 8.0)  => (0.04250360335290111, 0.03094978326353972),
    (0.9, 20.0) => (0.011068960619761698, 0.009740832807923332),
    (0.9, 50.0) => (0.0028212454064164235, 0.0026755651969439477),
)

@testset "Modes.Background" begin

    @testset "disc_OmegaK & kappar2 vs KerrGeodesics (baked KerrGeoFrequencies fixture)" begin
        maxOm = 0.0; maxka = 0.0
        for a in (0.0, 0.6, 0.9), r in (8.0, 20.0, 50.0)
            d = RelativisticDisc(a = a)
            Om_kg, Omr = kg_freq[(a, r)]
            Om_dm = disc_OmegaK(d, r)
            ka_dm = kappar2(d, r); ka_kg = Omr^2
            maxOm = max(maxOm, abs(Om_dm - Om_kg) / abs(Om_kg))
            maxka = max(maxka, abs(ka_dm - ka_kg) / abs(ka_kg))
        end
        @test maxOm < 1e-12
        @test maxka < 1e-10
    end

    @testset "risco_disc (Bardeen 1972 prograde ISCO)" begin
        @test risco_disc(RelativisticDisc(a = 0.0)) ≈ 6.0 atol = 1e-13
        # captured DiscModes _risco fixtures (== Bardeen closed form)
        @test risco_disc(RelativisticDisc(a = 0.6)) ≈ 3.829069418813151 rtol = 1e-13
        @test risco_disc(RelativisticDisc(a = 0.9)) ≈ 2.320883041761887 rtol = 1e-13
    end

    @testset "csound / eden0 / delta_r captured DiscModes fixtures" begin
        # (a, h0, r) => (csound, eden0, delta_r), r_I=10, eden_I=2 ; captured from packages/DiscModes
        fix = [
            (0.0, 0.1, 8.0, 0.04476614810358452, 3.125, 48.0),
            (0.0, 0.1, 20.0, 0.024260699053001704, 0.5, 360.0),
            (0.0, 0.1, 50.0, 0.014588051152682587, 0.08000000000000002, 2400.0),
            (0.6, 0.1, 8.0, 0.04101173734110984, 3.125, 48.36),
            (0.6, 0.1, 20.0, 0.02377896564057417, 0.5, 360.36),
            (0.6, 0.1, 50.0, 0.014515417193603745, 0.08000000000000002, 2400.36),
            (0.9, 0.05, 8.0, 0.01974768394831828, 3.125, 48.81),
            (0.9, 0.05, 20.0, 0.011779959757203301, 0.5, 360.81),
            (0.9, 0.05, 50.0, 0.00724019504033, 0.08000000000000002, 2400.81),
        ]
        for (a, h0, r, cs, ed, dl) in fix
            d = RelativisticDisc(a = a, h0 = h0, r_I = 10.0, eden_I = 2.0)
            @test csound(d, r) ≈ cs rtol = 1e-14
            @test eden0(d, r) ≈ ed rtol = 1e-14
            @test delta_r(d, r) ≈ dl rtol = 1e-14
        end
    end

    @testset "BigFloat inheritance (type-generic)" begin
        d = RelativisticDisc(a = big"0.6", h0 = big"0.1", r_I = big"10.0", eden_I = big"2.0")
        @test disc_OmegaK(d, big"20.0") isa BigFloat
        @test kappar2(d, big"20.0") isa BigFloat
        @test csound(d, big"20.0") isa BigFloat
    end

end
