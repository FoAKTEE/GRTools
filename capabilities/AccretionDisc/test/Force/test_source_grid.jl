# Force.Source — source_F_grid, the WING-PRECISION array pass (high-order Fornberg + Chebyshev
# denoise layers over the per-point source_F). GROUNDING (fixtures only; legacy cov_kerr_grid
# assertions carried at the legacy frozen tolerances):
#   • BYTE-IDENTITY inside the protect bands: the array pass overwrites the resonance/corotation
#     spike cores with the per-point source_F ⇒ Fg[k] == source_F(k) exactly (capture verified
#     byte@OLR=true on these fixtures).
#   • smooth far wing: high-order-FD grid value ≈ per-point 2-pt-FD value to <1e-2 (SAME physics,
#     only the derivative stencil differs; capture measured 3.7e-3 at r≈14).
#   • spectral denoise path runs finite at full grid length.
#   • f5_source(highorder=true) peaks at the SAME OLR as the per-point path (legacy Δ<0.2 + OLR ==).
#   • _highorder_drH_grid ≈ central-FD(hmunu) at an interior point to <1e-3 (legacy).
using Test
using AccretionDisc
using AccretionDisc.Force
const FSG = AccretionDisc.Force

const FIXDIR_SG = joinpath(@__DIR__, "..", "fixtures")

@testset "Force.source_F_grid (wing-precision array pass)" begin
    d = disc_a0(); m = 2; r_p = 8.0
    ω0 = m * disc_OmegaK(d, r_p)
    modes = load_modes(2:6, m; dir = FIXDIR_SG)
    g = modes[m].r_right; N = length(g)

    Fg = source_F_grid(d, modes, m, r_p, ω0; lmax = 6, npts = 9, spectral = false)
    @test length(Fg) == N && all(isfinite, Fg)

    rI, rO = lindblad_locations(d, r_p, m)
    @testset "protect band: BYTE-IDENTICAL to per-point source_F at the OLR spike" begin
        kprot = argmin(abs.(g .- rO))
        @test Fg[kprot] == source_F(d, modes, kprot, m, r_p, ω0; lmax = 6)
    end

    @testset "smooth far wing: high-order ≈ per-point (same physics, different stencil)" begin
        ksm = argmin(abs.(g .- (r_p + 6.0)))
        Fpt = source_F(d, modes, ksm, m, r_p, ω0; lmax = 6)
        @test abs(Fg[ksm] - Fpt) / (abs(Fpt) + eps()) < 1e-2          # legacy gate (capture 3.7e-3)
    end

    @testset "spectral denoise path (finite, full length)" begin
        Fg_spec = source_F_grid(d, modes, m, r_p, ω0; lmax = 6, npts = 9, spectral = true,
                                smooth_lo = r_p + 3, smooth_hi = g[end] - 0.2, cheb_deg = 6)
        @test length(Fg_spec) == N && all(isfinite, Fg_spec)
    end

    @testset "f5_source highorder branch: same OLR peak as per-point" begin
        res_pt = f5_source(d, m, r_p, modes, 6; rg = 8.8:0.1:11.5, verbose = false)
        res_hi = f5_source(d, m, r_p, modes, 6; rg = 8.8:0.1:11.5, verbose = false,
                           highorder = true, npts = 9)
        @test abs(Float64(res_pt.r_peak) - Float64(res_hi.r_peak)) < 0.2   # legacy gate
        @test res_pt.OLR == res_hi.OLR
    end

    @testset "_highorder_drH_grid vs central-FD(hmunu)" begin
        dH = FSG._highorder_drH_grid(modes[2], 2, m, π/2, 0.0; npts = 9)
        @test length(dH) == N
        kmid = argmin(abs.(g .- 11.0))
        fd = (hmunu(modes[2], kmid+1, 2, m, π/2) .- hmunu(modes[2], kmid-1, 2, m, π/2)) ./
             (g[kmid+1] - g[kmid-1])
        @test maximum(abs.(dH[kmid] .- fd)) / maximum(abs.(fd)) < 1e-3    # legacy gate
    end
end
