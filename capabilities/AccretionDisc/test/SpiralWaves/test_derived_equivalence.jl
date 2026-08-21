# test_derived_equivalence.jl — EGD::jd-equivalence-gate: the admission bar of route R3, the
# JULIA-NATIVE derivation of the arXiv:2608.17003 master system (capabilities/AccretionDisc/derivation).
#
# Route R3 emits `master_derived.jl` (_master_entries) and `r2_derived.jl` (_r2_delta_u_phi, _r2_amf,
# _r2_torque_density) from `Symbolics.jl` alone.  Route R2 — OUR OWN independent Mathematica
# derivation — is the ORACLE: `fixtures/oracle/r2_master_generated.jl` (the master entries C, S) and
# `fixtures/oracle/r2_generated.jl` (the diagnostics kernels) gate R3 here, they are not in the build
# path.  This file loads the oracle explicitly by path, into its own module, so it does not care
# which sources `SpiralWaves` includes.
#
# EGD::jd-notebook-purge (2026-08-22, owner decision): the notebook-derived oracle of route R1 was
# REMOVED from this repository, and every comparison that used it is now made against route R2 at
# UNCHANGED tolerance — legitimate because EGD::R1-vs-R2-gate certified R1 == R2 on-shell
# (`Simplify[R1|(11) − R2] = 0`; all 300 entries of the 50-point fixture byte-identical) and because
# the purge re-measured the margins: over the SAME 200 random on-shell points and the same seed, the
# new oracle reproduces the old one to 3.0e-12 (Float64) / 5.7e-75 (BigFloat 256), and the gate
# deviation moved from 3.0e-12 / 1.9e-73 (R3 vs R1) to 1.7e-12 / 1.9e-73 (R3 vs R2) [SOLID, measured].
# The R1 lineage and the full R1-vs-R2-vs-R3 evidence are archived in the mission repository
# (emri-disk results/emri-gasdisk paper_EMRI-GASDISK, evidence/{R1-vs-R2-gate,jd-equivalence-gate,
# jd-notebook-purge}).
#
# WHAT MAY AND MAY NOT BE COMPARED (the central fact of this gate, established below in (e)):
# the routes derive the SAME master system, but each carries its own OFF-SHELL continuation —
# they agree only where the background satisfies the two algebraic constraints of the alphabet,
#   (i)  W = 1/sqrt(1-v^2)                                                              (8)
#   (ii) nu' = v^2 psi' - e^{psi-nu} v omt' - p'/((eps+p) W^2)   (hydrostatic balance)  (11)
# and differ by multiples of those residuals elsewhere (this is exactly what the R1-vs-R2 gate
# recorded: "R1 == R2 symbolically ON-SHELL; off-shell the two routes differ by multiples of the (11)
# residual").  Every comparison below is therefore made ON-SHELL:
#   (a) 200 random on-shell `MasterInputs`, derived vs the R2 oracle CODE (Float64 + BigFloat 256);
#   (b) the 50 PHYSICAL points of `fixtures/master_r2_fixture.jl` (Mathematica route R2, 30 digits);
#   (c) the 50 points of `fixtures/r2_diag_fixture.jl` (Mathematica, 30 digits) for the diagnostics;
#   (d) the Newtonian limit (50b-g) with the deviation shrinking as 1/a•;
#   (e) the FINDING, re-expressed on the R2 lineage: route R2 eliminated W and nu' with (8) and (11),
#       so its entries do not read those two fields at all, while route R3 keeps them — the two
#       therefore agree to roundoff ON the constraint surface and diverge by O(1) off it, which is
#       why an off-shell fixture cannot gate two on-shell-equivalent routes;
#   (f) allocation-free after warmup, type-generic, `@inferred`.
using Test
using AccretionDisc
using AccretionDisc.SpiralWaves
using AccretionDisc.SpiralWaves: MasterInputs, master_inputs
using Random

# route R3 — the Julia-native derivation (the BUILD SOURCE after EGD::jd-swap)
module DerivedR3
    using AccretionDisc.SpiralWaves: MasterInputs
    include(joinpath(@__DIR__, "..", "..", "src", "SpiralWaves", "master_derived.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "SpiralWaves", "r2_derived.jl"))
end

# route R2 — our independent Mathematica derivation, the external-CAS ORACLE, kept forever as a gate
isdefined(@__MODULE__, :SpiralWavesOracle) ||
    include(joinpath(@__DIR__, "fixtures", "oracle", "oracle_module.jl"))
const OracleR2 = SpiralWavesOracle

include(joinpath(@__DIR__, "fixtures", "master_r2_fixture.jl"))
include(joinpath(@__DIR__, "fixtures", "r2_diag_fixture.jl"))

_dfx(::Type{T}, s::AbstractString) where {T<:Real} = parse(T, s)
_dfx(::Type{T}, z::Tuple{<:AbstractString,<:AbstractString}) where {T<:Real} =
    Complex(parse(T, z[1]), parse(T, z[2]))
_dreldev(a, b) = abs(a - b) / max(abs(b), eps(real(typeof(float(abs(b))))))

const _ENTRY_NAMES = ("C00", "C01", "C10", "C11", "S0", "S1")

"""
    _onshell_inputs(T, rng) -> (inp, omt, vp)

A random `MasterInputs` that satisfies BOTH algebraic constraints exactly: `W` is computed from `v`
(8) and `nu_r` from the hydrostatic balance (11); `omt` and `v_r` (the two extra background scalars
of the route-R2 signature) are recovered from Ω, Ω′ by (10) and its r-derivative, so the point is
simultaneously a legal input of every route.
"""
function _onshell_inputs(::Type{T}, rng::AbstractRNG) where {T<:Real}
    u(a, b) = T(a) + (T(b) - T(a)) * T(rand(rng))
    nu = u(-0.4, 0.2); psi = u(0.2, 1.5); mu = u(-0.3, 0.4)
    psi_r = u(-1, 1); mu_r = u(-1, 1); omt_r = u(-1, 1)
    v = u(-0.8, 0.8); W = 1 / sqrt(1 - v^2)
    p = u(0.1, 2); p_r = u(-1, 1); eps = u(0.2, 3); eps_r = u(-1, 1); Gamma1 = u(1.05, 2)
    enu = exp(nu); epsi = exp(psi)
    nu_r = v^2 * psi_r - epsi * v * omt_r / enu - p_r / ((eps + p) * W^2)          # (11)
    Omega = u(0.1, 1.2)
    omt = Omega - v * enu / epsi                                                   # (10)
    v_r = u(-1, 1)
    Omega_r = omt_r + enu * (v_r + v * (nu_r - psi_r)) / epsi                      # d/dr of (10)
    ωd = Complex(u(-1, 1) + (rand(rng, Bool) ? 1 : -1) * T(2 // 10), u(-0.3, 0.3))
    inp = MasterInputs{T,Complex{T}}(rand(rng, 1:12), nu, nu_r, psi, psi_r, mu, mu_r, omt_r,
        W, v, Omega, Omega_r, ωd, p, p_r, eps, eps_r, Gamma1,
        u(-1, 1), u(-1, 1), u(-1, 1), u(-1, 1), u(-1, 1), u(-1, 1), u(-1, 1), u(-1, 1))
    (inp, omt, v_r)
end

"a copy of `inp` pushed OFF the constraint surface: W and nu' are scaled away from (8) and (11)."
function _offshell_copy(inp::MasterInputs{T,C}, fW, fnu) where {T,C}
    MasterInputs{T,C}(inp.m, inp.nu, T(fnu) * inp.nu_r, inp.psi, inp.psi_r, inp.mu, inp.mu_r,
        inp.omt_r, T(fW) * inp.W, inp.v, inp.Omega, inp.Omega_r, inp.omega_d,
        inp.p, inp.p_r, inp.eps, inp.eps_r, inp.Gamma1,
        inp.htt, inp.htt_r, inp.htphi, inp.htphi_r, inp.hphiphi, inp.hphiphi_r, inp.hrr, inp.hzz)
end

# the covariant analytic test profile of the on-shell fixture (identical to test_master.jl (f))
function _hmodes_gate(r, a, ΩE)
    htt = exp(-((r - a) / (a / 2))^2)
    dhtt = htt * (-2 * (r - a) / (a / 2)^2)
    (htt = htt, htt_r = dhtt,
     htphi = 3 // 10 * r^2 * ΩE * htt, htphi_r = 3 // 10 * ΩE * (2r * htt + r^2 * dhtt),
     hphiphi = 7 // 10 * r^2 * htt, hphiphi_r = 7 // 10 * (2r * htt + r^2 * dhtt),
     hrr = htt / 2, hzz = htt / 5)
end

# MasterInputs from the route-R2 fixture's own background columns (identical to test_r2_generated.jl)
function _inputs_r2_gate(::Type{T}, pt) where {T<:Real}
    g(k) = parse(T, getproperty(pt, k))
    MasterInputs{T,Complex{T}}(pt.m, g(:nu), g(:nup), g(:psi), g(:psip), g(:mu), g(:mup), g(:omtp),
        1 / sqrt(1 - g(:v)^2), g(:v), g(:Omega), g(:Omp), Complex{T}(g(:omd), zero(T)),
        g(:p), g(:pp), g(:eps), g(:epsp), g(:Gam1),
        g(:htt), g(:http), g(:htp), g(:htpp), g(:hpp), g(:hppp), g(:hrr), g(:hzz))
end

@testset "SpiralWaves R3-derived == route-R2 oracle (arXiv:2608.17003)" begin

    @testset "(0) the LIVE build source IS the derived one (EGD::jd-swap)" begin
        # the kernels `SpiralWaves` actually uses must be bitwise the emitted route-R3 ones; this is
        # what makes every other testset of the suite a gate on the Julia derivation
        rng = MersenneTwister(7)
        for _ in 1:20
            inp, omt, vp = _onshell_inputs(Float64, rng)
            y0 = 0.3 - 0.1im; y1 = -0.2 + 0.4im
            @test AccretionDisc.SpiralWaves._master_entries(inp) === DerivedR3._master_entries(inp)
            @test AccretionDisc.SpiralWaves._r2_delta_u_phi(inp, omt, vp, y0, y1) ===
                  DerivedR3._r2_delta_u_phi(inp, omt, vp, y0, y1)
            @test AccretionDisc.SpiralWaves._r2_amf(inp, omt, vp, y0, y1) ===
                  DerivedR3._r2_amf(inp, omt, vp, y0, y1)
            @test AccretionDisc.SpiralWaves._r2_torque_density(inp, omt, vp, y0, y1) ===
                  DerivedR3._r2_torque_density(inp, omt, vp, y0, y1)
        end
        @test !isfile(joinpath(@__DIR__, "..", "..", "src", "SpiralWaves", "master_generated.jl"))
        @test !isfile(joinpath(@__DIR__, "..", "..", "src", "SpiralWaves", "r2_generated.jl"))
        # EGD::jd-notebook-purge: no notebook-derived artefact may EVER come back into this package —
        # neither the R1 oracle code nor either of the two fixtures whose numbers were read from the
        # R1 columns. They are archived in emri-disk (paper_EMRI-GASDISK, evidence/R1-vs-R2-gate).
        for gone in ("fixtures/oracle/master_generated.jl", "fixtures/master_r1_fixture.jl",
                     "fixtures/master_onshell_fixture.jl",
                     "fixtures/capture/bake_master_onshell_fixture.py")
            @test !isfile(joinpath(@__DIR__, split(gone, "/")...))
        end
    end

    @testset "(a) 200 random ON-SHELL MasterInputs: derived vs the R2 oracle code" begin
        for (T, N, tolC, tolD) in ((Float64, 200, 1e-10, 1e-10), (BigFloat, 20, 1e-60, 1e-60))
            setprecision(BigFloat, 256) do
                rng = MersenneTwister(4242)
                worst = zeros(real(T), 6)
                wdup = zero(real(T)); wamf = zero(real(T)); wtq = zero(real(T))
                for _ in 1:N
                    inp, omt, vp = _onshell_inputs(T, rng)
                    got = DerivedR3._master_entries(inp)
                    ref = OracleR2._r2_master_entries(inp, omt, vp)
                    for j in 1:6
                        worst[j] = max(worst[j], _dreldev(got[j], ref[j]))
                    end
                    y0 = Complex{T}(T(3 // 10), T(-1 // 10)); y1 = Complex{T}(T(-2 // 10), T(4 // 10))
                    wdup = max(wdup, _dreldev(DerivedR3._r2_delta_u_phi(inp, omt, vp, y0, y1),
                                              OracleR2._r2_delta_u_phi(inp, omt, vp, y0, y1)))
                    wamf = max(wamf, _dreldev(DerivedR3._r2_amf(inp, omt, vp, y0, y1),
                                              OracleR2._r2_amf(inp, omt, vp, y0, y1)))
                    wtq = max(wtq, _dreldev(DerivedR3._r2_torque_density(inp, omt, vp, y0, y1),
                                            OracleR2._r2_torque_density(inp, omt, vp, y0, y1)))
                end
                @info "R3 vs R2 ($T, $N on-shell points): C,S $(Float64.(worst)); " *
                      "δû_φ $(Float64(wdup)), AMF $(Float64(wamf)), torque $(Float64(wtq))"
                for j in 1:6
                    @test worst[j] ≤ tolC
                end
                @test wdup ≤ tolD
                @test wamf ≤ tolD
                @test wtq ≤ tolD
            end
        end
    end

    @testset "(b) the 50 PHYSICAL on-shell points: derived vs the route-R2 fixture" begin
        @test length(MASTER_R2_FIXTURE) == 50
        worstC = 0.0; worstS = 0.0
        for pt in MASTER_R2_FIXTURE
            bg = DiscBackground(χ = parse(Float64, pt.chi), aEMRI = parse(Float64, pt.aEMRI),
                                h0 = parse(Float64, pt.h0), 𝔭 = parse(Float64, pt.pdisk),
                                𝔮 = parse(Float64, pt.qdisk))
            r = parse(Float64, pt.r)
            inp = master_inputs(bg, pt.m, r; hmodes = _hmodes_gate(r, bg.aEMRI, Omega_emri(bg)),
                                landau = false)
            got = DerivedR3._master_entries(inp)
            ref = (_dfx(Float64, pt.C00), _dfx(Float64, pt.C01), _dfx(Float64, pt.C10),
                   _dfx(Float64, pt.C11), _dfx(Float64, pt.S0), _dfx(Float64, pt.S1))
            for j in 1:4
                worstC = max(worstC, _dreldev(got[j], ref[j]))
            end
            for j in 5:6
                worstS = max(worstS, _dreldev(got[j], ref[j]))
            end
        end
        # the DERIVED form is better conditioned in Float64 than the CAS form of the same entries:
        # over these 50 points the R3 build measures 2.6e-12 where the route-R1 build measured
        # 7.0e-10 (EGD::jd-swap, archived in emri-disk) [SOLID, measured]
        @info "on-shell fixture (Float64, 50 points): derived worst rel dev C $worstC, S $worstS"
        @test worstC ≤ 1e-10
        @test worstS ≤ 1e-12
        setprecision(BigFloat, 256) do
            worst = big"0"
            for pt in MASTER_R2_FIXTURE[1:10:50]
                bg = DiscBackground{BigFloat}(χ = parse(BigFloat, pt.chi), aEMRI = parse(BigFloat, pt.aEMRI),
                                              h0 = parse(BigFloat, pt.h0), 𝔭 = parse(BigFloat, pt.pdisk),
                                              𝔮 = parse(BigFloat, pt.qdisk))
                r = parse(BigFloat, pt.r)
                inp = master_inputs(bg, pt.m, r; hmodes = _hmodes_gate(r, bg.aEMRI, Omega_emri(bg)),
                                    landau = false)
                got = DerivedR3._master_entries(inp)
                ref = (_dfx(BigFloat, pt.C00), _dfx(BigFloat, pt.C01), _dfx(BigFloat, pt.C10),
                       _dfx(BigFloat, pt.C11), _dfx(BigFloat, pt.S0), _dfx(BigFloat, pt.S1))
                for j in 1:6
                    worst = max(worst, _dreldev(got[j], ref[j]))
                end
            end
            @info "on-shell fixture (BigFloat 256, 5 points): derived worst rel dev = $(Float64(worst))"
            @test worst ≤ big"1e-27"      # the fixture itself carries 30 significant digits
        end
    end

    @testset "(c) the route-R2 diagnostics fixture: derived δû_φ, AMF (72a), torque density (72b)" begin
        for (T, n, tol) in ((Float64, 50, 1e-12), (BigFloat, 8, 1e-28))
            setprecision(BigFloat, 256) do
                wd = zero(T); wa = zero(T); wt = zero(T)
                for pt in R2_DIAG_FIXTURE[1:n]
                    inp = _inputs_r2_gate(T, pt)
                    omt = parse(T, pt.omt); vp = parse(T, pt.vp)
                    y0 = _dfx(T, pt.sample_y0); y1 = _dfx(T, pt.sample_y1)
                    c0 = _dfx(T, pt.deltaUphi_coeffs[1]); c1 = _dfx(T, pt.deltaUphi_coeffs[2])
                    ch = _dfx(T, pt.deltaUphi_coeffs[3])
                    wd = max(wd, _dreldev(DerivedR3._r2_delta_u_phi(inp, omt, vp, y0, y1),
                                          c0 * y0 + c1 * y1 + ch))
                    wa = max(wa, _dreldev(DerivedR3._r2_amf(inp, omt, vp, y0, y1),
                                          _dfx(T, pt.R2_AMF_at_sample)))
                    wt = max(wt, _dreldev(DerivedR3._r2_torque_density(inp, omt, vp, y0, y1),
                                          _dfx(T, pt.R2_TorqueDensity_at_sample)))
                end
                @info "r2 fixture ($T, $n points): derived δû_φ $(Float64(wd)), AMF $(Float64(wa)), torque $(Float64(wt))"
                @test wd ≤ tol
                @test wa ≤ tol
                @test wt ≤ tol
            end
        end
    end

    @testset "(d) Newtonian limit of the DERIVED entries: (50b-g), deviation ∝ 1/a•" begin
        m = 3; Φ = -2.0e-3; dΦ = 3.0e-7
        hm = (htt = -2Φ, htt_r = -2dΦ, htphi = 0.0, htphi_r = 0.0,
              hphiphi = 0.0, hphiphi_r = 0.0, hrr = 0.0, hzz = 0.0)
        devs = Dict{Float64,Float64}()
        for aEMRI in (1e4, 1e5)
            bg = DiscBackground(χ = 0.0, aEMRI = aEMRI, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0)
            worst = 0.0
            for x in (0.7, 1.1, 1.6)
                r = x * aEMRI
                got = DerivedR3._master_entries(master_inputs(bg, m, r; hmodes = hm, landau = false))
                Cn, Sn = master_matrices_newtonian(bg, m, r, Φ, dΦ; landau = false)
                ref = (Cn[1, 1], Cn[1, 2], Cn[2, 1], Cn[2, 2], Sn[1], Sn[2])
                for j in 1:6
                    worst = max(worst, _dreldev(got[j], ref[j]))
                end
                @test worst ≤ 25 / aEMRI        # the O(M/a•) bound measured for the CAS routes
            end
            devs[aEMRI] = worst
        end
        ratio = devs[1e4] / devs[1e5]
        @info "derived Newtonian limit: dev(a•=1e4) $(devs[1e4]), dev(a•=1e5) $(devs[1e5]), ratio $ratio"
        @test 7 ≤ ratio ≤ 13
    end

    @testset "(e) FINDING: the routes share only their ON-SHELL restriction" begin
        # Route R2 eliminated W with (8) and nu' with the hydrostatic balance (11) BEFORE emission, so
        # `_r2_master_entries` never reads `inp.W` or `inp.nu_r`; route R3 emits the on-shell form but
        # keeps both fields in its alphabet. On the constraint surface the two agree to roundoff; push
        # the SAME point off the surface by rescaling exactly those two fields and R2 cannot move while
        # R3 does — which is why an OFF-SHELL fixture (independent random inputs, e.g. the 40-point R1
        # fixture archived in emri-disk) can never gate two on-shell-equivalent routes. [FINDING]
        setprecision(BigFloat, 256) do
            rng = MersenneTwister(20260822)
            worst_on = big"0"; worst_off = big"0"; worst_r2_moved = big"0"
            for _ in 1:40
                inp, omt, vp = _onshell_inputs(BigFloat, rng)
                off = _offshell_copy(inp, 13 // 10, 7 // 10)      # W → 1.3 W, nu' → 0.7 nu'
                on3 = DerivedR3._master_entries(inp);  on2 = OracleR2._r2_master_entries(inp, omt, vp)
                of3 = DerivedR3._master_entries(off);  of2 = OracleR2._r2_master_entries(off, omt, vp)
                for j in 1:6
                    worst_on = max(worst_on, _dreldev(on3[j], on2[j]))
                    worst_off = max(worst_off, _dreldev(of3[j], of2[j]))
                    worst_r2_moved = max(worst_r2_moved, _dreldev(of2[j], on2[j]))
                end
            end
            @info "constraint surface, BigFloat(256): R3 vs R2 ON-SHELL $(Float64(worst_on)), " *
                  "the same points pushed OFF-SHELL $(Float64(worst_off)); R2's own response to " *
                  "(W, ν′) off-shell = $(Float64(worst_r2_moved))"
            @test worst_on ≤ big"1e-60"        # pure BigFloat roundoff — the routes ARE the same system
            @test worst_off > big"1e-3"        # the off-shell continuations genuinely differ [FINDING]
            @test worst_r2_moved == 0          # R2 is EXACTLY blind to W and nu' (they were eliminated)
        end
    end

    @testset "(f) allocation-free after warmup, inferred, BigFloat-generic" begin
        bg = DiscBackground(χ = 0.6, aEMRI = 15.0, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0)
        hm = (htt = -1e-3, htt_r = 1e-5, htphi = 2e-4, htphi_r = -3e-5,
              hphiphi = 1e-4, hphiphi_r = 2e-5, hrr = 5e-4, hzz = 1e-4)
        inp = master_inputs(bg, 4, 17.5; hmodes = hm, landau = true)
        omt = 0.01; vp = -1e-3; y0 = 0.3 + 0.1im; y1 = -0.2 + 0.05im
        # the kernels are `@inline`, so an `@allocated` on the bare call would compile the inlined
        # body inside the measurement closure; measure a warmed wrapper instead
        mw(i) = sum(abs, DerivedR3._master_entries(i))
        aw(i) = abs(DerivedR3._r2_amf(i, omt, vp, y0, y1))
        tw(i) = abs(DerivedR3._r2_torque_density(i, omt, vp, y0, y1))
        dw(i) = abs(DerivedR3._r2_delta_u_phi(i, omt, vp, y0, y1))
        ow(i) = sum(abs, OracleR2._r2_master_entries(i, omt, vp))
        mw(inp); aw(inp); tw(inp); dw(inp); ow(inp)
        @test (@allocated mw(inp)) == 0
        @test (@allocated aw(inp)) == 0
        @test (@allocated tw(inp)) == 0
        @test (@allocated dw(inp)) == 0
        @test (@allocated ow(inp)) == 0
        @test (@inferred DerivedR3._master_entries(inp)) isa NTuple{6,ComplexF64}
        @test (@inferred DerivedR3._r2_amf(inp, omt, vp, y0, y1)) isa ComplexF64
        @test (@inferred DerivedR3._r2_torque_density(inp, omt, vp, y0, y1)) isa ComplexF64
        @test (@inferred DerivedR3._r2_delta_u_phi(inp, omt, vp, y0, y1)) isa ComplexF64
        @test (@inferred OracleR2._r2_master_entries(inp, omt, vp)) isa NTuple{6,ComplexF64}
        setprecision(BigFloat, 256) do
            bgb = DiscBackground{BigFloat}(χ = 0.6, aEMRI = 15, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0)
            hmb = (htt = big"-1e-3", htt_r = big"1e-5", htphi = big"2e-4", htphi_r = big"-3e-5",
                   hphiphi = big"1e-4", hphiphi_r = big"2e-5", hrr = big"5e-4", hzz = big"1e-4")
            inb = master_inputs(bgb, 4, big"17.5"; hmodes = hmb, landau = true)
            gb = DerivedR3._master_entries(inb)
            @test all(x -> x isa Complex{BigFloat}, gb)
            @test all(isfinite, gb)
            ob = OracleR2._r2_master_entries(inb, big"0.01", big"-1e-3")
            @test all(x -> x isa Complex{BigFloat}, ob)
            @test precision(real(ob[1])) == 256          # the oracle inherits BigFloat precision too
            @test DerivedR3._r2_amf(inb, big"0.01", big"-1e-3",
                                    Complex{BigFloat}(big"0.3", big"0.1"),
                                    Complex{BigFloat}(big"-0.2", big"0.05")) isa Complex{BigFloat}
        end
    end
end
