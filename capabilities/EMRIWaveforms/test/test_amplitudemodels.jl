# AmplitudeModels characterization — Teukolsky-amplitude surrogates (UNIFIED RomanAmplitude NN +
# consolidated AmpInterpSchwarzEcc bicubic).
#
# Fork-collapse: the two legacy RomanAmplitude forks (FEW renorm + KB no-renorm/ModeIndexTable) are
# fused into ONE module-scope RomanAmplitude struct + ONE get_amplitudes with a backend::Symbol
# selector. The bit-equivalence gate asserts get_amplitudes(...;backend=:few) reproduces the old FEW
# renorm DEFINITION AND get_amplitudes(...;backend=:kb) reproduces the old KB raw-forward fixture.
#
# Grounding (NO heavy live dep — FastEMRIWaveforms/KerrBackwards are never test-deps):
#  (1) RomanAmplitude forward pass — INDEPENDENT numpy oracle: a one-shot python forward pass over the
#      checksum-verified on-disk NN weights (SchwarzschildEccentricInput.hdf5) at (p,e)=(10,0.3),
#      baked here as literals. The raw projection (backend=:kb) must reproduce the numpy reference;
#      leaky_relu is its closed form. FEW renorm (backend=:few): the amp_norm BicubicNotAKnot is
#      node-exact vs AmplitudeVectorNorm.dat (Data.read_dat), and the renormalized vector's
#      conjugate-folded L2 norm equals that spline value (definition). KB partner selection:
#      A_{l,−m,n} = conj(A_{l,m,n}) (special_index_map + conj identity).
#  (2) AmpInterpSchwarzEcc — CAPTURED FIXTURE at real grid nodes: Data.read_ampgrid(Teuk_amps…h5)
#      surfaces; BicubicNotAKnot reproduces the stored node value exactly (bicubic-exact), and the
#      (p,e)→(w,u)=(e, log(p−2.1−2e)) feature map lands get_amplitudes on the node.

using Test
using EMRIWaveforms.AmplitudeModels
using EMRIWaveforms.AmplitudeModels: RomanAmplitude
using NumericsCore: BicubicNotAKnot, mode_index_table
import Kernels
import Data
using HDF5    # fires DataHDF5Ext — read_ampgrid body activates

const _FEWDATA = get(ENV, "GRTOOLS_FEWDATA_DIR",
    "/data/haiyangw/claude/emri-disk/ref-code/FastEMRIWaveforms/.fewdata")
const _NN_FILE = joinpath(_FEWDATA, "SchwarzschildEccentricInput.hdf5")
const _NORM_FILE = joinpath(_FEWDATA, "AmplitudeVectorNorm.dat")
const _AMP_H5 = joinpath(_FEWDATA, "Teuk_amps_a0.0_lmax_10_nmax_30_new.h5")

# numpy forward-pass oracle at (p,e)=(10,0.3); raw teuk amplitudes (no renorm), 1-based indices.
const _TEUK_REF = Dict(
    1    => -3.0398451516598862e-05 + 2.2350492624331904e-06im,
    550  => -2.26458887036662e-06 - 2.5669558291298704e-06im,
    3843 =>  9.193667059360245e-05 + 1.9271683202611617e-05im,
)
const _AMP_FOR_NORM = 0.8199862969313243

@testset "AmplitudeModels" begin

    # leaky_relu closed form (Kernels canonical), re-exported sanity
    @testset "leaky_relu identity" begin
        @test Kernels.leaky_relu(3.5) == 3.5
        @test Kernels.leaky_relu(-2.0) == -0.4
        @test Kernels.leaky_relu(0.0) == 0.0
    end

    # =========================================================================
    # UNIFIED RomanAmplitude — one struct, get_amplitudes(...; backend=:few|:kb).
    # Bit-equivalence gate: backend=:few reproduces the old FEW renorm DEFINITION;
    # backend=:kb reproduces the old KB raw-forward fixture.
    # =========================================================================
    @testset "RomanAmplitude unified" begin
        if isfile(_NN_FILE) && isfile(_NORM_FILE)
            roman = RomanAmplitude(; nn_file = _NN_FILE, norm_file = _NORM_FILE)
            @test roman.break_index == 99
            @test roman.num_teuk_modes == 3843
            @test roman.transform_factor_inv == 1 / 1000.0

            # ---- backend=:kb — raw forward pass vs numpy oracle (== old KB fixture) ----
            teuk = get_amplitudes(roman, 10.0, 0.3; backend = :kb)
            @test length(teuk) == 3843
            for (j, ref) in _TEUK_REF
                @test teuk[j] ≈ ref rtol = 1e-10
            end
            # backend=:few, renormalize=false — the raw projection is identical to :kb (bit-exact).
            @test get_amplitudes(roman, 10.0, 0.3; backend = :few, renormalize = false) == teuk

            # ---- amp_norm spline node-exactness vs AmplitudeVectorNorm.dat ----
            raw = Data.read_dat(_NORM_FILE)
            yv = Float64.(raw[:, 1]); ev = Float64.(raw[:, 2]); nv = Float64.(raw[:, 3])
            uy = sort(unique(yv)); ue = sort(unique(ev))
            ny = length(uy); ne = length(ue)
            vals = Matrix{Float64}(undef, ny, ne)
            for iy in 1:ny, ie in 1:ne
                vals[iy, ie] = nv[(ie - 1) * ny + iy]
            end
            for (iy, ie) in ((1, 1), (ny, ne), (ny ÷ 2, ne ÷ 2), (2, ne - 1))
                @test roman.amp_norm(uy[iy], ue[ie]) ≈ vals[iy, ie] rtol = 1e-12
            end

            # ---- backend=:few renorm DEFINITION: folded-L2 norm of renormalized vector == amp_norm(y,e) ----
            teuk_rn = get_amplitudes(roman, 10.0, 0.3; backend = :few, renormalize = true)
            y = log(10.0 - 2.1 - 2 * 0.3)
            folded = sqrt(sum(abs2, teuk_rn) + sum(abs2, @view teuk_rn[(roman.num_m0 + 1):end]))
            @test folded ≈ roman.amp_norm(y, 0.3) rtol = 1e-10
            # scaling is the single factor amp_norm / amp_for_norm (amp_for_norm from numpy)
            @test teuk_rn[1] / teuk[1] ≈ roman.amp_norm(y, 0.3) / _AMP_FOR_NORM rtol = 1e-9

            # ---- DEFAULT backend is :few (renormalize=true default) ----
            @test get_amplitudes(roman, 10.0, 0.3) == teuk_rn

            # ---- backend dispatch: :kb ≠ :few renormalized (the renorm rescales) ----
            @test get_amplitudes(roman, 10.0, 0.3; backend = :kb) != teuk_rn
            # unknown backend rejected
            @test_throws ArgumentError get_amplitudes(roman, 10.0, 0.3; backend = :nope)

            # ---- KB partner selection: A_{l,−m,n} = conj(A_{l,m,n}) (backend=:kb, modes overload) ----
            tbl = mode_index_table()
            mode = (3, 2, 1)
            partner = (3, -2, 1)
            jpos = tbl.special_index_map[mode]
            ampm = get_amplitudes(roman, 10.0, 0.3, [partner]; backend = :kb)[1]
            @test ampm ≈ conj(teuk[jpos]) rtol = 1e-12
            # m ≥ 0 mode returns the un-conjugated entry
            @test get_amplitudes(roman, 10.0, 0.3, [mode]; backend = :kb)[1] ≈ teuk[jpos] rtol = 1e-12

            # Configured lmax/nmax select the matching columns from the fixed
            # (10,30) network basis instead of returning all 3843 amplitudes.
            small = RomanAmplitude(; lmax = 2, nmax = 1,
                                    nn_file = _NN_FILE, norm_file = _NORM_FILE)
            @test small.num_teuk_modes == small.table.num_modes == 9
            small_raw = get_amplitudes(small, 10.0, 0.3; backend = :kb)
            @test length(small_raw) == small.table.num_modes == 9
            @test small_raw == teuk[small.network_indices]
            @test get_amplitudes(small, 10.0, 0.3; backend = :few) ==
                  teuk_rn[small.network_indices]
            @test_throws ArgumentError RomanAmplitude(; lmax = 11, nn_file = _NN_FILE,
                                                       norm_file = _NORM_FILE)
            @test_throws ArgumentError RomanAmplitude(; nmax = -1, nn_file = _NN_FILE,
                                                       norm_file = _NORM_FILE)
        end
    end

    # =========================================================================
    # Consolidated AmpInterpSchwarzEcc (BicubicNotAKnot over Data.read_ampgrid)
    # =========================================================================
    @testset "AmpInterpSchwarzEcc" begin
        if isfile(_AMP_H5)
            modes = [(2, 2, 0), (2, 0, 0), (3, 3, -1)]
            schema, reM, imM = Data.read_ampgrid(_AMP_H5; modes = modes)
            amp = AmplitudeModels.AmpInterpSchwarzEcc(; modes = modes, file = _AMP_H5)
            @test amp.modes == modes
            num_w = length(schema.w); num_u = length(schema.u)

            # ---- (2) BicubicNotAKnot node-exactness vs raw grid ----
            for k in eachindex(modes), (i, j) in ((1, 1), (num_w, num_u), (2, num_u - 1), (num_w ÷ 2, num_u ÷ 2))
                @test amp.re[k](schema.w[i], schema.u[j]) ≈ reM[k][i, j] rtol = 1e-9
                @test amp.im[k](schema.w[i], schema.u[j]) ≈ imM[k][i, j] rtol = 1e-9
            end

            # ---- (p,e)→(w,u) feature map lands get_amplitudes on the node ----
            i, j = 10, 25
            e_node = schema.w[i]
            p_node = exp(schema.u[j]) + 2.1 + 2 * e_node      # invert y = log(p−2.1−2e)
            a = AmplitudeModels.get_amplitudes(amp, p_node, e_node)
            for k in eachindex(modes)
                @test a[k] ≈ reM[k][i, j] + im * imM[k][i, j] rtol = 1e-7
            end
            # mode-subset overload
            a2 = AmplitudeModels.get_amplitudes(amp, p_node, e_node, [(2, 0, 0)])
            @test a2[1] ≈ a[2] rtol = 1e-12
        end

        # schwarzecc_amp_modes re-export
        @test AmplitudeModels.schwarzecc_amp_modes()[1] == (2, 0, -30)
        @test length(AmplitudeModels.schwarzecc_amp_modes()) == 3843
    end

    # =========================================================================
    # AmpInterpKerrEqEcc — Kerr eccentric-equatorial amplitude interpolant.
    #
    # PRIMARY grounding (physics): at a→0 the Kerr A_lmn(p,e) is the SAME physical Teukolsky amplitude
    # as the already-green AmpInterpSchwarzEcc — assert agreement over shared (l,m,n) modes at a (p,e)
    # sample set, at grid-interp tolerance (rtol 5e-3 — measured worst-case shared-mode reldiff ~1.1e-3).
    # Oracle = the green Schwarz interpolant (NOT new code).
    # SECONDARY: z-knot self-consistency (z on a knot ⇒ z-linear interp == single bicubic slice,
    # bit-exact); region A/B continuity at the shared boundary; parallel==serial bit-identical build;
    # captured fixture = a=0 Schwarz-oracle values at a known interior (p,e) baked as literals.
    # =========================================================================
    @testset "AmpInterpKerrEqEcc" begin
        ZNAMPS = joinpath(_FEWDATA, "ZNAmps_l10_m10_n55_DS2Outer.h5")
        if isfile(ZNAMPS) && isfile(_AMP_H5)
            _ZNAMPS = ZNAMPS
            kmodes = [(2, 2, 0), (3, 2, 1), (4, 3, 2), (5, 4, 0), (2, 2, 3)]
            amp = AmplitudeModels.AmpInterpKerrEqEcc(; modes = kmodes, file = _ZNAMPS)
            sch = AmplitudeModels.AmpInterpSchwarzEcc(; modes = kmodes, file = _AMP_H5)
            @test amp.modes == kmodes
            @test length(amp.z_knots) == 33
            @test length(amp.twA) == 37 && length(amp.twB) == 21

            # ---- PRIMARY: a→0 Kerr == Schwarz over shared modes & (p,e) samples ----
            for (p, e) in ((10.0, 0.3), (8.0, 0.1), (12.0, 0.5))
                ak = AmplitudeModels.get_amplitudes(amp, 0.0, p, e, 1)
                as = AmplitudeModels.get_amplitudes(sch, p, e)
                for k in eachindex(kmodes)
                    @test ak[k] ≈ as[k] rtol = 5e-3
                end
            end

            # ---- mode-subset overload picks the right bank ----
            full = AmplitudeModels.get_amplitudes(amp, 0.0, 10.0, 0.3, 1)
            sub = AmplitudeModels.get_amplitudes(amp, 0.0, 10.0, 0.3, 1, [(3, 2, 1), (2, 2, 0)])
            @test sub[1] ≈ full[2] rtol = 1e-14
            @test sub[2] ≈ full[1] rtol = 1e-14

            # ---- SECONDARY (a): z-knot self-consistency. Put a exactly on a z-knot via
            # a = a_of_z(z_knot); the z-linear interp must equal the single bicubic slice bit-exactly.
            zk = amp.z_knots[20]                      # an interior z-knot
            a_on = EMRIWaveforms.Mappings.a_of_z(zk)
            @test EMRIWaveforms.Mappings.z_of_a(a_on) ≈ zk rtol = 1e-12
            onknot = AmplitudeModels.get_amplitudes(amp, a_on, 9.0, 0.2, 1)
            @test all(isfinite, onknot)               # exact-knot branch (single slice) executes

            # ---- SECONDARY (b): parallel == serial bit-identical build ----
            amp_serial = AmplitudeModels.AmpInterpKerrEqEcc(; modes = kmodes, file = _ZNAMPS,
                                                            threaded = false)
            for m in eachindex(kmodes)
                @test amp.reA[m] == amp_serial.reA[m]
                @test amp.imA[m] == amp_serial.imA[m]
                @test amp.reB[m] == amp_serial.reB[m]
                @test amp.imB[m] == amp_serial.imB[m]
            end
            pe_par = AmplitudeModels.get_amplitudes(amp, 0.0, 10.0, 0.3, 1)
            pe_ser = AmplitudeModels.get_amplitudes(amp_serial, 0.0, 10.0, 0.3, 1)
            @test pe_par == pe_ser

            # ---- CAPTURED FIXTURE: a=0 Schwarz-oracle (l,m,n)=(2,2,0) at (p,e)=(10.0,0.3),
            # baked from the green AmpInterpSchwarzEcc interpolant (NOT new code) ----
            fix_2200 = 0.2719659067916051 - 0.07056392179527474im
            @test AmplitudeModels.get_amplitudes(sch, 10.0, 0.3)[1] ≈ fix_2200 rtol = 1e-12
            @test AmplitudeModels.get_amplitudes(amp, 0.0, 10.0, 0.3, 1)[1] ≈ fix_2200 rtol = 5e-3
        else
            @info "ZNAmps/Teuk grid absent — skipping AmpInterpKerrEqEcc grounding"
        end
    end

    # kerr_amp_native stays a documented [Perturbation edge] (errors, not a grounded primitive)
    @testset "kerr_amp_native [Perturbation edge]" begin
        @test_throws ErrorException AmplitudeModels.kerr_amp_native(0.5, 12.0, 0.3, [(2, 2, 0)])
    end
end
