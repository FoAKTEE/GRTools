# Characterization of the hybrid binary-waveform assembly (Hybrid.BinaryWaveform, §4.9).
#
# GROUNDING (no WaSABI live dep — composition identity over already-grounded pieces + symmetry):
#   With the driver trajectory (Drivers) and the amplitude model (Models) each already grounded,
#   waveform(bim,l,m,t) MUST equal _mode_amplitude(amp,l,m,Ω(t),ics)·exp(-imφ(t)) — here
#   h0PA_ampΩ(amp,l,m,Ω(t),ν)·exp(-imφ(t)); plus the equatorial symmetry
#   h_{l,-m}=(-1)^l conj(h_{l,m}) (BinaryInspiral.m:172); mode_list == sorted amp keys;
#   trajectory accessors == omega_at/r0_at/phi_at; _outside_domain error for t<0 / t>duration;
#   binary_inspiral rejects mismatched ICS keys.

using Test
using Hybrid.BinaryWaveform
using Hybrid.Drivers: omega_at, phi_at, r0_at
using Hybrid.Models: h0PA_ampΩ

# Build a synthetic 0PASchwarz sf + amp data tree (degree-≤2 flux / degree-1 complex amps).
function _make_0pas_dirs()
    rs = collect(range(6.0, 32.0; length = 16))
    fI(r) = 0.02 + 1.0e-3 * r - 3.0e-5 * r^2
    fH(r) = 0.005 - 2.0e-4 * r
    sf = mkpath(joinpath(mktempdir(), "1SF_Flux", "Schwarz_Circ"))
    _wl(p, ys) = write(p, "{" * join(("{$(r), $(y)}" for (r, y) in zip(rs, ys)), ", ") * "}")
    _wl(joinpath(sf, "1SFCircShwarzDotEInf.m"), fI.(rs))
    _wl(joinpath(sf, "1SFCircShwarzDotEHorizon.m"), fH.(rs))
    amp = mkpath(joinpath(mktempdir(), "1SF", "Schwarz_Circ"))
    zf(l, m, r) = complex(0.01 + 0.001r + 1e-4(l + m), -0.002 + 5e-4r - 3e-4(l - m))
    _cx(z) = imag(z) >= 0 ? "$(real(z)) + $(imag(z))*I" : "$(real(z)) - $(-imag(z))*I"
    for l in 2:5, mm in 1:l
        vals = [zf(l, mm, r) for r in rs]
        write(joinpath(amp, "1SFTeukampSchwarzCirc$(l)$(mm).m"),
              "{" * join(("{$(r), $(_cx(z))}" for (r, z) in zip(rs, vals)), ", ") * "}")
    end
    return dirname(dirname(sf)), dirname(dirname(amp))
end

@testset "BinaryWaveform" begin
    sfroot, amproot = _make_0pas_dirs()
    ics = Dict("M" => 1.0, "r0" => 6.5, "ν" => 1.0e-2, "φ" => 0.0)
    bim = binary_inspiral(ics; model = "0PASchwarz", sf_datadir = sfroot, amp_datadir = amproot)

    @testset "mode_list == sorted amp keys" begin
        @test mode_list(bim) == [(l, m) for l in 2:5 for m in 1:l]
    end

    @testset "waveform composition identity h_lm = A_lm(Ω(t)) e^{-imφ(t)}" begin
        td = duration(bim)
        for t in (0.1td, 0.5td, 0.9td), (l, m) in ((2, 2), (3, 1), (4, 3))
            Ω = omega_at(bim.inspiral, t); φ = phi_at(bim.inspiral, t)
            expected = h0PA_ampΩ(bim.amp, l, m, Ω, ics["ν"]) * exp(-im * m * φ)
            @test waveform(bim, l, m, t) ≈ expected rtol = 0
        end
    end

    @testset "equatorial negative-m symmetry" begin
        td = duration(bim)
        for t in (0.3td, 0.7td), (l, m) in ((2, 2), (3, 2), (5, 4))
            @test waveform(bim, l, -m, t) == (-1)^l * conj(waveform(bim, l, m, t))
        end
    end

    @testset "trajectory accessors" begin
        td = duration(bim)
        for t in (0.2td, 0.6td)
            @test trajectory(bim, "Ω", t) == omega_at(bim.inspiral, t)
            @test trajectory(bim, "r0", t) == r0_at(bim.inspiral, t)
            @test trajectory(bim, "φ", t) == phi_at(bim.inspiral, t)
            @test trajectory(bim, "ν", t) == ics["ν"]
            @test trajectory(bim, "M", t) == ics["M"]
        end
    end

    @testset "domain / mode / ICS errors" begin
        td = duration(bim)
        @test_throws ErrorException waveform(bim, 2, 2, -1.0)
        @test_throws ErrorException waveform(bim, 2, 2, td * 1.0001 + 1.0)
        @test_throws ErrorException waveform(bim, 9, 1, 0.5td)        # mode not available
        @test_throws ErrorException trajectory(bim, "Ω", -1.0)
        bad = Dict("M" => 1.0, "r0" => 6.5, "ν" => 1.0e-2)           # missing φ
        @test_throws ErrorException binary_inspiral(bad; model = "0PASchwarz",
                                                    sf_datadir = sfroot, amp_datadir = amproot)
    end
end
