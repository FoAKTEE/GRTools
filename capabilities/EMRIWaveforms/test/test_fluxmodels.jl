# FluxModels characterization — grid-interpolated flux ODE RHS (KB Schwarz-ecc / FEW Kerr-ecc-eq forks).
#
# Grounding (NO heavy live dep — FastEMRIWaveforms/KerrBackwards are never test-deps):
#  (A) ANALYTIC baked literals — the PN closed forms (EdotPN/LdotPN/kerr_PN_alt/kerr_pdot_PN/
#      kerr_edot_PN) and the Schwarz pdot/edot algebra are pure transcriptions of the verbatim
#      formulas in KerrBackwards/schwarzeccflux.jl + FEW/trajectory/kerrflux.jl; the expected values
#      are computed INDEPENDENTLY in python (the arithmetic oracle), baked here as literals.
#  (B) CAPTURED FIXTURE at real grid nodes — the checksum-verified on-disk FEW grids
#      FluxNewMinusPNScaled_fixed_y_order.dat (Schwarz) read via Data.read_dat and
#      KerrEccEqFluxData.h5 (Kerr) read via Data.read_fluxgrid(:A/:B): the interpolant surface
#      reproduces the raw grid node value exactly (GslBicubic / TricubicNotAKnot are node-exact),
#      and the loader composes read_fluxgrid + apex + PN-division correctly.
#  (C) separatrix-guard zeroing checked at p = 6 + 2e (KB fork upstream guard).

using Test
using EMRIWaveforms.FluxModels
using NumericsCore: GslBicubic, TricubicNotAKnot
using EMRIWaveforms.Mappings: apex_of_uwyz, apex_of_UWYZ
using Spacetime: get_separatrix
import Data
using HDF5    # fires DataHDF5Ext — read_fluxgrid body activates
using TestKit: oracle_root, oracle_or_skip

const KB = EMRIWaveforms.FluxModels.KerrBackwards
const FEW = EMRIWaveforms.FluxModels.FEW

const _FEWDATA = oracle_root("GRTOOLS_FEWDATA_DIR",
    "/data/haiyangw/claude/emri-disk/ref-code/FastEMRIWaveforms/.fewdata")
const _SCHWARZ_DAT = joinpath(_FEWDATA, "FluxNewMinusPNScaled_fixed_y_order.dat")
const _KERR_H5 = joinpath(_FEWDATA, "KerrEccEqFluxData.h5")

@testset "FluxModels" begin

    # =========================================================================
    # KB fork — SchwarzEccFlux
    # =========================================================================
    @testset "KB.SchwarzEccFlux" begin
        # ---- (A) PN closed-form pieces vs python oracle ----
        @testset "EdotPN/LdotPN analytic literals" begin
            # python: (96+292e²+37e⁴)/(15(1-e²)^3.5)·yPN^5 and 4(8+7e²)/(5(-1+e²)²)·yPN^3.5
            @test KB._EdotPN(0.3, 0.1) ≈ 0.00011367956994091843 rtol = 1e-14
            @test KB._LdotPN(0.3, 0.1) ≈ 0.002636440643135188  rtol = 1e-14
        end

        # ---- (A) pdot/edot algebra vs python oracle ----
        @testset "pdot/edot algebra literal" begin
            pdot, edot = KB._pdot_edot(-1e-4, -1e-3, 10.0, 0.3)
            @test pdot ≈ -0.0036213470912744563 rtol = 1e-13
            @test edot ≈ -0.003749153448430986  rtol = 1e-13
            # e = 0 ⇒ edot = 0 (upstream branch)
            @test KB._pdot_edot(-1e-4, -1e-3, 10.0, 0.0)[2] == 0.0
        end

        # ---- traits ----
        @testset "ODEModel traits" begin
            if oracle_or_skip(isfile(_SCHWARZ_DAT), "KB.SchwarzEccFlux ODEModel traits"; env = "GRTOOLS_FEWDATA_DIR")
                m = KB.SchwarzEccFlux(; file = _SCHWARZ_DAT)
                @test KB.convert_Y(m) == false
                @test KB.background(m) == :Schwarzschild
                @test KB.num_add_args(m) == 0
            end
        end

        # ---- (B) GslBicubic surface == raw .dat grid node value ----
        @testset "flux surface == raw grid node (Data.read_dat)" begin
            if oracle_or_skip(isfile(_SCHWARZ_DAT), "Schwarz flux surface == raw .dat grid node"; env = "GRTOOLS_FEWDATA_DIR")
                raw = Data.read_dat(_SCHWARZ_DAT)
                yc = Float64.(raw[:, 1]); ec = Float64.(raw[:, 2])
                Ec = Float64.(raw[:, 3]); Lc = Float64.(raw[:, 4])
                uy = sort(unique(yc)); ue = sort(unique(ec))
                nx = length(uy); ny = length(ue)
                ZE = Matrix{Float64}(undef, nx, ny)
                ZL = Matrix{Float64}(undef, nx, ny)
                for j in 1:ny, i in 1:nx
                    k = (j - 1) * nx + i
                    ZE[i, j] = Ec[k]; ZL[i, j] = Lc[k]
                end
                m = KB.SchwarzEccFlux(; file = _SCHWARZ_DAT)
                for (i, j) in ((1, 1), (1, ny), (nx, 1), (nx, ny),
                               (nx ÷ 2, ny ÷ 2), (2, ny - 1))
                    @test m.Edot(uy[i], ue[j]) ≈ ZE[i, j] rtol = 1e-12
                    @test m.Ldot(uy[i], ue[j]) ≈ ZL[i, j] rtol = 1e-12
                end
            end
        end

        # ---- (C) separatrix-guard zeroing at/below p = 6 + 2e ----
        @testset "separatrix guard zeroes at p <= 6 + 2e" begin
            if oracle_or_skip(isfile(_SCHWARZ_DAT), "Schwarz separatrix-guard zeroing"; env = "GRTOOLS_FEWDATA_DIR")
                m = KB.SchwarzEccFlux(; file = _SCHWARZ_DAT)
                e = 0.3
                psep = 6.0 + 2.0 * e
                d = KB.schwarz_ecc_flux_deriv(m, 1e-5, 0.0, psep - 0.5, e, 1.0)
                @test all(==(0.0), d)
                @test all(==(0.0), KB.schwarz_ecc_flux_deriv(m, 1e-5, 0.0, psep, e, 1.0))
                # ode_deriv alias agrees
                @test KB.ode_deriv(m, 1e-5, 0.0, psep - 0.5, e, 1.0) == d
            end
        end

        # ---- full deriv composition: finite, physical (inspiral ⇒ pdot < 0) ----
        @testset "deriv physical sanity above separatrix" begin
            if oracle_or_skip(isfile(_SCHWARZ_DAT), "Schwarz deriv physical sanity"; env = "GRTOOLS_FEWDATA_DIR")
                m = KB.SchwarzEccFlux(; file = _SCHWARZ_DAT)
                pdot, edot, xdot, Oph, Oth, Or = KB.schwarz_ecc_flux_deriv(m, 1e-5, 0.0, 10.0, 0.3, 1.0)
                @test all(isfinite, (pdot, edot, xdot, Oph, Oth, Or))
                @test pdot < 0          # energy loss shrinks p
                @test xdot == 0.0       # equatorial
                @test Oth == Oph        # Schwarzschild Ω_θ = Ω_φ
                @test Oph > 0 && Or > 0
            end
        end
    end

    # =========================================================================
    # FEW fork — KerrEccEqFlux
    # =========================================================================
    @testset "FEW.KerrEccEqFlux" begin
        # ---- (A) PN normalizations vs python oracle ----
        @testset "kerr PN normalizations analytic literals" begin
            Edot, Ldot = FEW.kerr_PN_alt(10.0, 0.3)
            @test Edot ≈ 7.09397058842768e-05 rtol = 1e-14
            @test Ldot ≈ 0.0018952341408474051 rtol = 1e-14
            @test FEW.kerr_pdot_PN(10.0, 0.3, 6.0, 6.6) ≈ 0.07664010977481024 rtol = 1e-13
            @test FEW.kerr_edot_PN(10.0, 0.3, 6.0, 6.6) ≈ 0.011651798072127264 rtol = 1e-13
        end

        # ---- (B) loader composition + TricubicNotAKnot node-exactness ----
        @testset "flux surface == raw grid node (Data.read_fluxgrid)" begin
            if oracle_or_skip(isfile(_KERR_H5), "Kerr flux surface == raw fluxgrid node"; env = "GRTOOLS_FEWDATA_DIR")
                # ≥4 nodes/axis (5×5×5) — not-a-knot tricubic is singular at 3 nodes; fast.
                ds = ((32, 16, 16), (16, 8, 8))
                flux = FEW.KerrEccEqFlux(; convention = "ELQ", downsample = ds, file = _KERR_H5)
                @test flux.convention == "ELQ"

                for (tag, region, dsr, is_A) in ((:A, :A, ds[1], true), (:B, :B, ds[2], false))
                    schema, ch = Data.read_fluxgrid(_KERR_H5, region; downsample = dsr)
                    u, w, z = schema.u, schema.w, schema.z
                    nu, nw, nz = length(u), length(w), length(z)
                    Eg = flux.grids[Symbol(:Edot, tag)]
                    Lg = flux.grids[Symbol(:Ldot, tag)]
                    @test size(Eg) == (nu, nw, nz)
                    for (i, j, k) in ((1, 1, 1), (2, 2, 2), (nu, nw, nz), (2, nw, 2), (nu, 2, nz))
                        # loader composition: stored grid = raw / PN(apex). Region A uses the close
                        # apex map; region B the far map (apex_of_UWYZ) — same split the loader uses.
                        ai, pi, ei, xi = is_A ? apex_of_uwyz(u[i], w[j], 1.0, z[k]) :
                                         apex_of_UWYZ(u[i], w[j], 1.0, z[k], true)
                        EPN, LPN = FEW.kerr_PN_alt(pi, ei)
                        @test Eg[i, j, k] ≈ ch[:Edot][i, j, k] / EPN rtol = 1e-12
                        @test Lg[i, j, k] ≈ ch[:Ldot][i, j, k] / LPN rtol = 1e-12
                        # tricubic node-exactness
                        @test flux.interp[Symbol(:Edot, tag)](u[i], w[j], z[k]) ≈ Eg[i, j, k] rtol = 1e-10
                        @test flux.interp[Symbol(:Ldot, tag)](u[i], w[j], z[k]) ≈ Lg[i, j, k] rtol = 1e-10
                    end
                end
            end
        end

        # ---- traits ----
        @testset "ODEModel traits" begin
            if oracle_or_skip(isfile(_KERR_H5), "FEW.KerrEccEqFlux ODEModel traits"; env = "GRTOOLS_FEWDATA_DIR")
                flux = FEW.KerrEccEqFlux(; convention = "ELQ",
                                         downsample = ((32, 16, 16), (16, 8, 8)), file = _KERR_H5)
                @test FEW.convert_Y(flux) == false
                @test FEW.background(flux) == :Kerr
                @test FEW.num_add_args(flux) == 0
            end
        end

        # ---- retrograde separatrix (FEW flux.py:640-641) ----
        # The default pLSO is the apex separatrix get_separatrix(a, e, x), not one taken at the grid's
        # signed spin a_in = -a; under the signed-spin chart (-a, x = -1) is the PROGRADE orbit.
        @testset "retrograde default pLSO is the apex separatrix" begin
            if oracle_or_skip(isfile(_KERR_H5), "FEW.KerrEccEqFlux retrograde pLSO"; env = "GRTOOLS_FEWDATA_DIR")
                for convention in ("ELQ", "pex")
                    flux = FEW.KerrEccEqFlux(; convention,
                                             downsample = ((32, 16, 16), (16, 8, 8)), file = _KERR_H5)
                    for (a, p, e) in ((0.9, 12.0, 0.3), (0.5, 14.0, 0.2))
                        psep = get_separatrix(a, e, -1.0)
                        @test psep > get_separatrix(a, e, 1.0)
                        @test FEW.interpolate_flux_grids(flux, p, e, -1.0, a) ===
                              FEW.interpolate_flux_grids(flux, p, e, -1.0, a; pLSO = psep)
                    end
                end
            end
        end
    end
end
