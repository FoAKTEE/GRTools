# test_derived_equivalence.jl — EGD::jd-equivalence-gate: the admission bar of route R3, the
# JULIA-NATIVE derivation of the arXiv:2608.17003 master system (capabilities/AccretionDisc/derivation).
#
# Route R3 emits `master_derived.jl` (_master_entries) and `r2_derived.jl` (_r2_delta_u_phi, _r2_amf,
# _r2_torque_density) from `Symbolics.jl` alone.  The authors' supplementary notebook (route R1,
# `master_generated.jl`) and the independent Mathematica derivation (route R2, `r2_generated.jl`) are
# ORACLES: they gate R3 here, they are not in the build path.  This file loads BOTH generated sets
# explicitly by path, each into its own module, so it does not care which one `SpiralWaves` includes.
#
# WHAT MAY AND MAY NOT BE COMPARED (the central fact of this gate, established below in (e)):
# the three routes derive the SAME master system, but each carries its own OFF-SHELL continuation —
# they agree only where the background satisfies the two algebraic constraints of the alphabet,
#   (i)  W = 1/sqrt(1-v^2)                                                              (8)
#   (ii) nu' = v^2 psi' - e^{psi-nu} v omt' - p'/((eps+p) W^2)   (hydrostatic balance)  (11)
# and differ by multiples of those residuals elsewhere (this is exactly what the R1-vs-R2 gate
# recorded: "R1 == R2 symbolically ON-SHELL; off-shell the two routes differ by multiples of the (11)
# residual").  Every comparison below is therefore made ON-SHELL:
#   (a) 200 random on-shell `MasterInputs`, derived vs the R1/R2 oracle CODE (Float64 + BigFloat 256);
#   (b) the 50 PHYSICAL points of `fixtures/master_onshell_fixture.jl` (Mathematica, 30 digits);
#   (c) the 50 points of `fixtures/r2_diag_fixture.jl` (Mathematica, 30 digits) for the diagnostics;
#   (d) the Newtonian limit (50b-g) with the deviation shrinking as 1/a•;
#   (e) the FINDING: the 40-point `fixtures/master_r1_fixture.jl` is OFF-SHELL by construction (every
#       input an independent random rational) and sits O(1) away from BOTH constraints, so it cannot
#       compare two on-shell-equivalent routes; projected onto the constraint surface, R3 and R1 agree
#       to 3e-74 at all 40 points, which is the strongest single statement of this gate;
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

# routes R1 + R2 — the external-CAS ORACLES (notebook / Mathematica), kept forever as a gate
isdefined(@__MODULE__, :SpiralWavesOracle) ||
    include(joinpath(@__DIR__, "fixtures", "oracle", "oracle_module.jl"))
const OracleR1R2 = SpiralWavesOracle

include(joinpath(@__DIR__, "fixtures", "master_r1_fixture.jl"))
include(joinpath(@__DIR__, "fixtures", "master_onshell_fixture.jl"))
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

@testset "SpiralWaves R3-derived == notebook/R2 oracles (arXiv:2608.17003)" begin

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
    end

    @testset "(a) 200 random ON-SHELL MasterInputs: derived vs the R1/R2 oracle code" begin
        for (T, N, tolC, tolD) in ((Float64, 200, 1e-10, 1e-10), (BigFloat, 20, 1e-60, 1e-60))
            setprecision(BigFloat, 256) do
                rng = MersenneTwister(4242)
                worst = zeros(real(T), 6)
                wdup = zero(real(T)); wamf = zero(real(T)); wtq = zero(real(T))
                for _ in 1:N
                    inp, omt, vp = _onshell_inputs(T, rng)
                    got = DerivedR3._master_entries(inp)
                    ref = OracleR1R2._master_entries(inp)
                    for j in 1:6
                        worst[j] = max(worst[j], _dreldev(got[j], ref[j]))
                    end
                    y0 = Complex{T}(T(3 // 10), T(-1 // 10)); y1 = Complex{T}(T(-2 // 10), T(4 // 10))
                    wdup = max(wdup, _dreldev(DerivedR3._r2_delta_u_phi(inp, omt, vp, y0, y1),
                                              OracleR1R2._r2_delta_u_phi(inp, omt, vp, y0, y1)))
                    wamf = max(wamf, _dreldev(DerivedR3._r2_amf(inp, omt, vp, y0, y1),
                                              OracleR1R2._r2_amf(inp, omt, vp, y0, y1)))
                    wtq = max(wtq, _dreldev(DerivedR3._r2_torque_density(inp, omt, vp, y0, y1),
                                            OracleR1R2._r2_torque_density(inp, omt, vp, y0, y1)))
                end
                @info "R3 vs R1/R2 ($T, $N on-shell points): C,S $(Float64.(worst)); " *
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

    @testset "(b) the 50 PHYSICAL on-shell points: derived vs the Mathematica fixture" begin
        @test length(MASTER_ONSHELL_FIXTURE) == 50
        worstC = 0.0; worstS = 0.0
        for pt in MASTER_ONSHELL_FIXTURE
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
        # the DERIVED form is better conditioned in Float64 than the notebook form: the same fixture
        # measures 7.0e-10 for route R1 (test_master.jl (f)) and 2.6e-12 here [SOLID, measured]
        @info "on-shell fixture (Float64, 50 points): derived worst rel dev C $worstC, S $worstS"
        @test worstC ≤ 1e-10
        @test worstS ≤ 1e-12
        setprecision(BigFloat, 256) do
            worst = big"0"
            for pt in MASTER_ONSHELL_FIXTURE[1:10:50]
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
                @test worst ≤ 25 / aEMRI        # the O(M/a•) bound measured for route R1
            end
            devs[aEMRI] = worst
        end
        ratio = devs[1e4] / devs[1e5]
        @info "derived Newtonian limit: dev(a•=1e4) $(devs[1e4]), dev(a•=1e5) $(devs[1e5]), ratio $ratio"
        @test 7 ≤ ratio ≤ 13
    end

    @testset "(e) FINDING: the 40-point R1 fixture is OFF-SHELL and cannot gate route R3" begin
        # 1. the fixture sits O(1) away from BOTH algebraic constraints (every input is an
        #    independent random rational, by construction of R1_numeric_fixture.wls)
        wW = 0.0; wH = 0.0
        for pt in MASTER_R1_FIXTURE
            v = parse(Float64, pt.v); W = parse(Float64, pt.W)
            enu = exp(parse(Float64, pt.nu)); epsi = exp(parse(Float64, pt.psi))
            nu_on = v^2 * parse(Float64, pt.psi_r) - epsi * v * parse(Float64, pt.omt_r) / enu -
                    parse(Float64, pt.p_r) / ((parse(Float64, pt.eps) + parse(Float64, pt.p)) * W^2)
            wW = max(wW, abs(W - 1 / sqrt(1 - v^2)) * sqrt(1 - v^2))
            wH = max(wH, abs(parse(Float64, pt.nu_r) - nu_on) / max(abs(nu_on), 1e-300))
        end
        @info "R1 fixture distance from the constraints: |W√(1−v²) − 1| ≤ $wW, |ν′/ν′_(11) − 1| ≤ $wH"
        @test wW > 0.1                       # NOT on the W = 1/√(1−v²) branch
        @test wH > 0.1                       # NOT on the hydrostatic balance (11)

        # 2. projected onto the constraint surface (W and ν′ recomputed, everything else kept) the two
        #    routes agree to BigFloat(256) roundoff at all 40 points — the strongest statement here
        setprecision(BigFloat, 256) do
            worst_on = big"0"; worst_off = big"0"
            for pt in MASTER_R1_FIXTURE
                T = BigFloat
                v = _dfx(T, pt.v); nu = _dfx(T, pt.nu); psi = _dfx(T, pt.psi)
                enu = exp(nu); epsi = exp(psi)
                psi_r = _dfx(T, pt.psi_r); omt_r = _dfx(T, pt.omt_r)
                p = _dfx(T, pt.p); p_r = _dfx(T, pt.p_r); eps = _dfx(T, pt.eps)
                for onshell in (true, false)
                    W = onshell ? 1 / sqrt(1 - v^2) : _dfx(T, pt.W)
                    nu_r = onshell ? v^2 * psi_r - epsi * v * omt_r / enu - p_r / ((eps + p) * W^2) :
                                     _dfx(T, pt.nu_r)
                    inp = MasterInputs{T,Complex{T}}(pt.m, nu, nu_r, psi, psi_r, _dfx(T, pt.mu),
                        _dfx(T, pt.mu_r), omt_r, W, v, _dfx(T, pt.Omega), _dfx(T, pt.Omega_r),
                        _dfx(T, pt.omega_d), p, p_r, eps, _dfx(T, pt.eps_r), _dfx(T, pt.Gamma1),
                        _dfx(T, pt.htt), _dfx(T, pt.htt_r), _dfx(T, pt.htphi), _dfx(T, pt.htphi_r),
                        _dfx(T, pt.hphiphi), _dfx(T, pt.hphiphi_r), _dfx(T, pt.hrr), _dfx(T, pt.hzz))
                    got = DerivedR3._master_entries(inp); ref = OracleR1R2._master_entries(inp)
                    for j in 1:6
                        d = _dreldev(got[j], ref[j])
                        onshell ? (worst_on = max(worst_on, d)) : (worst_off = max(worst_off, d))
                    end
                end
            end
            @info "R1 fixture points, BigFloat(256): derived vs R1 code — ON-SHELL projection " *
                  "$(Float64(worst_on)), as printed (off-shell) $(Float64(worst_off))"
            @test worst_on ≤ big"1e-60"      # measured 3.5e-74 — pure BigFloat roundoff
            @test worst_off > big"1e-3"      # the off-shell continuations genuinely differ [FINDING]
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
        mw(inp); aw(inp); tw(inp); dw(inp)
        @test (@allocated mw(inp)) == 0
        @test (@allocated aw(inp)) == 0
        @test (@allocated tw(inp)) == 0
        @test (@allocated dw(inp)) == 0
        @test (@inferred DerivedR3._master_entries(inp)) isa NTuple{6,ComplexF64}
        @test (@inferred DerivedR3._r2_amf(inp, omt, vp, y0, y1)) isa ComplexF64
        @test (@inferred DerivedR3._r2_torque_density(inp, omt, vp, y0, y1)) isa ComplexF64
        @test (@inferred DerivedR3._r2_delta_u_phi(inp, omt, vp, y0, y1)) isa ComplexF64
        setprecision(BigFloat, 256) do
            bgb = DiscBackground{BigFloat}(χ = 0.6, aEMRI = 15, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0)
            hmb = (htt = big"-1e-3", htt_r = big"1e-5", htphi = big"2e-4", htphi_r = big"-3e-5",
                   hphiphi = big"1e-4", hphiphi_r = big"2e-5", hrr = big"5e-4", hzz = big"1e-4")
            inb = master_inputs(bgb, 4, big"17.5"; hmodes = hmb, landau = true)
            gb = DerivedR3._master_entries(inb)
            @test all(x -> x isa Complex{BigFloat}, gb)
            @test all(isfinite, gb)
            @test DerivedR3._r2_amf(inb, big"0.01", big"-1e-3",
                                    Complex{BigFloat}(big"0.3", big"0.1"),
                                    Complex{BigFloat}(big"-0.2", big"0.05")) isa Complex{BigFloat}
        end
    end
end
