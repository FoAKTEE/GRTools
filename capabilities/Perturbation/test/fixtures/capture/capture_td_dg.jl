# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/capture_td_dg.jl (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = ac3ed66c2449dde7bc887ea48fa336b471adb68fb6afedb793a4d87e66d9cc95
# fixture(s)  : td_mode1d.jl, td_flux1d.jl, td_selfforce.jl, td_kerrtail.jl, td_kerrflux.jl, td_substitution.jl
# oracle      : /data/haiyangw/claude/DG-Teulolsky @ 0a88967 (own env)
# run         : cd $DGTEUK_ROOT && julia --project=. <path>/capture_td_dg.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const DGTEUK_ROOT          = get(ENV, "GRTOOLS_DGTEUK_ROOT",          "/data/haiyangw/claude/DG-Teulolsky")
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

#!/usr/bin/env julia
# ============================================================================================
# ONE-SHOT characterization capture for the DG-P3 TimeDomain port (Perturbation.TimeDomain).
# Runs the ORIGINAL /data/haiyangw/claude/DG-Teulolsky (READ-ONLY; executed in ITS OWN env:
#   julia --project=/data/haiyangw/claude/DG-Teulolsky scratchpad/capture_td_dg.jl )
# and bakes fixtures into grtools/capabilities/Perturbation/test/fixtures/td_*.jl.
# Configs are shrunk from the repo's experiments/data_*.jl so each capture runs <~2-3 min;
# every fixture header records the EXACT config. NEVER regenerate from the ported code.
# ============================================================================================
using Printf
using ScalarTeukolskyDG
using ScalarTeukolskyDG.Hyperboloidal
using ScalarTeukolskyDG.Backgrounds
using ScalarTeukolskyDG.Mode1D
using ScalarTeukolskyDG.Flux1D
using ScalarTeukolskyDG.KerrModes
using ScalarTeukolskyDG.KerrFlux

const OUT = "$(FIXTURE_DIR)"
const REPO = "$(DGTEUK_ROOT)"
const COMMIT = strip(read(`git -C $REPO rev-parse --short HEAD`, String))

vec_str(v) = "[" * join(repr.(Float64.(v)), ", ") * "]"
ivec_str(v) = "[" * join(string.(v), ", ") * "]"

function header(io, what, cfg)
    println(io, "# AUTO-GENERATED one-shot fixture — DO NOT EDIT.")
    println(io, "# Baked ORIGINAL DG-Teulolsky ($REPO @ $COMMIT, its own env, julia $(VERSION))")
    println(io, "# by scratchpad/capture_td_dg.jl on $(Sys.MACHINE). $what")
    println(io, "# EXACT CONFIG: $cfg")
end

sample_idx(n, fracs) = unique(clamp.(round.(Int, fracs .* n), 2, n))

# ---------------------------------------------------------------------------------------------
# (e/substitution) Ylm_equator, orbital_params, tortoise/inverse, potential/coupling/friction,
#                  Hyperboloidal layer values — grounds the reuse substitutions.
# ---------------------------------------------------------------------------------------------
let t0 = time()
    lm = [(0,0),(1,0),(1,1),(2,0),(2,1),(2,2),(2,-2),(3,1),(3,3),(3,-1),(4,2),(5,3),(5,5),
          (6,6),(7,3),(9,7),(15,15),(10,-6)]
    orb = [(a, rp) for a in (0.0, 0.5, 0.9) for rp in (6.0, 8.0, 10.0, 16.0914363989845)]
    tor = [(a, r) for a in (0.0, 0.5, 0.9) for r in (2.5, 6.0, 10.0, 100.0, 1000.0)]
    rsv = [(a, rs) for a in (0.0, 0.5, 0.9) for rs in (-40.0, -5.0, 0.0, 12.0, 500.0)]
    bg  = [(a, r, l, m) for a in (0.0, 0.9) for r in (2.5, 6.0, 20.0) for (l, m) in ((0,0),(2,0),(2,2),(5,3))]
    L12 = Layer(150.0, 1200.0, 4); L35 = Layer(30.0, 50.0, 4)
    hyp = [(150.0,1200.0,4,ρ) for ρ in (10.0, 149.0, 150.0, 400.0, 900.0, 1199.999, 1200.0)]
    open(joinpath(OUT, "td_substitution.jl"), "w") do io
        header(io, "Reuse-substitution ground truth: repo Ylm_equator/orbital_params/tortoise/inverse_tortoise/potential_V/coupling_f/friction_mu/Hyperboloidal.", "M=1 everywhere (repo experiment convention)")
        println(io, "# TD_YLM_REF: (l, m, Ylm_equator(l,m))")
        println(io, "const TD_YLM_REF = [")
        for (l, m) in lm
            println(io, "    ($l, $m, $(repr(Ylm_equator(l, m)))),")
        end
        println(io, "]")
        println(io, "# TD_ORBIT_REF: (a, rp, E, Lz, ut, Omega)  [Kerr circular equatorial prograde]")
        println(io, "const TD_ORBIT_REF = [")
        for (a, rp) in orb
            op = orbital_params(rp; M=1.0, a=a)
            println(io, "    ($(repr(a)), $(repr(rp)), $(repr(op.E)), $(repr(op.Lz)), $(repr(op.ut)), $(repr(op.Ω))),")
        end
        println(io, "]")
        println(io, "# TD_TORTOISE_REF: (a, r, tortoise(r; M=1, a))")
        println(io, "const TD_TORTOISE_REF = [")
        for (a, r) in tor
            println(io, "    ($(repr(a)), $(repr(r)), $(repr(tortoise(r; M=1.0, a=a)))),")
        end
        println(io, "]")
        println(io, "# TD_INVTORT_REF: (a, rs, inverse_tortoise(rs; M=1, a))  [repo tol=1e-13]")
        println(io, "const TD_INVTORT_REF = [")
        for (a, rs) in rsv
            println(io, "    ($(repr(a)), $(repr(rs)), $(repr(inverse_tortoise(rs; M=1.0, a=a)))),")
        end
        println(io, "]")
        println(io, "# TD_BG_REF: (a, r, l, m, potential_V, coupling_f, friction_mu_imag)")
        println(io, "const TD_BG_REF = [")
        for (a, r, l, m) in bg
            V = potential_V(r, l, m; M=1.0, a=a); f = coupling_f(r; M=1.0, a=a)
            μ = friction_mu(r, m; M=1.0, a=a)
            println(io, "    ($(repr(a)), $(repr(r)), $l, $m, $(repr(V)), $(repr(f)), $(repr(imag(μ)))),")
        end
        println(io, "]")
        println(io, "# TD_HYP_REF: (R, s, P, rho, Omega, Omegaprime, denomOmega, H)  [Layer(R,s,P)]")
        println(io, "const TD_HYP_REF = [")
        for (R, s, P, ρ) in hyp
            L = Layer(R, s, P)
            println(io, "    ($(repr(R)), $(repr(s)), $P, $(repr(ρ)), $(repr(Ωfun(L, ρ))), $(repr(Ωprime(L, ρ))), $(repr(denomΩ(L, ρ))), $(repr(Hfun(L, ρ)))),")
        end
        println(io, "]")
        println(io, "# TD_HYP_TEST_REF: Layer(30,50,4) values from the repo's own unit testset")
        L = L35
        println(io, "const TD_HYP_TEST_REF = (H10 = $(repr(Hfun(L, 10.0))), Om10 = $(repr(Ωfun(L, 10.0))), " *
                    "H50 = $(repr(Hfun(L, 50.0))), Om50 = $(repr(Ωfun(L, 50.0))), " *
                    "dOm50 = $(repr(denomΩ(L, 50.0))), H30 = $(repr(Hfun(L, 30.0))))")
    end
    @printf("[substitution] done in %.1f s\n", time() - t0); flush(stdout)
end

# ---------------------------------------------------------------------------------------------
# (a) Mode1D Schwarzschild Price tail, ℓ=0 (repo's own runtests config, F1.3 class shrunk)
# ---------------------------------------------------------------------------------------------
let t0 = time()
    cfg = "Mode1DProblem(0; N=8, K=120, ρL=-200.0, R=150.0, s=1200.0, P=4, M=1, extract_rstar=[500.0, Inf]); solve_mode1d(tfinal=6000.0, dt=0.0866, record_every=100, momentum_ic=(amp=1/sqrt(2π*10^2), ρ0=30.0, w=10.0)); fit_tail window=(2500,6000)"
    p = Mode1DProblem(0; N=8, K=120, ρL=-200.0, R=150.0, s=1200.0, P=4,
                      extract_rstar=[500.0, Inf])
    τs, Ψ = solve_mode1d(p; tfinal=6000.0, dt=0.0866, record_every=100)
    n1 = fit_tail(τs, Ψ[1, :]; window=(2500.0, 6000.0))
    n2 = fit_tail(τs, Ψ[2, :]; window=(2500.0, 6000.0))
    js = sample_idx(length(τs), [0.02, 0.1, 0.25, 0.5, 0.75, 1.0])
    open(joinpath(OUT, "td_mode1d.jl"), "w") do io
        header(io, "(a) Mode1D Schwarzschild ℓ=0 Price tail: ψ probe samples at r*=500 and scri + fit_tail exponents.", cfg)
        println(io, "const TD_MODE1D = (")
        println(io, "    nsamp = $(length(τs)),")
        println(io, "    idx = $(ivec_str(js)),")
        println(io, "    tau = $(vec_str(τs[js])),")
        println(io, "    psi_r500 = $(vec_str(Ψ[1, js])),")
        println(io, "    psi_scri = $(vec_str(Ψ[2, js])),")
        println(io, "    n_r500 = $(repr(n1)),")
        println(io, "    n_scri = $(repr(n2)),")
        println(io, "    # problem-construction probes (grounds inverse_tortoise→Spacetime + verbatim C):")
        println(io, "    probe_rho = $(vec_str(p.probe_ρ)),")
        println(io, "    probe_k = $(ivec_str(p.probe_k)),")
        println(io, "    C_3_17 = $(repr(p.C[3, 17])),")
        println(io, "    C_5_60 = $(repr(p.C[5, 60])),")
        println(io, "    C_scri = $(repr(p.C[end, end])),")
        println(io, "    H_5_100 = $(repr(p.H[5, 100])),")
        println(io, ")")
    end
    @printf("[mode1d] n_r500=%.4f n_scri=%.4f  done in %.1f s\n", n1, n2, time() - t0); flush(stdout)
end

# ---------------------------------------------------------------------------------------------
# (b) Flux1D s=0 Schwarzschild circular-orbit flux (2,2) rp=10 (T2.1 class shrunk) + selfforce_ft
# ---------------------------------------------------------------------------------------------
let t0 = time()
    cfg = "Flux1DProblem(2, 2, 10.0; N=8, K=100, ρL=-150.0, R=100.0, s=200.0, P=4, M=1, q=1, ramp_τ=400.0, ramp_δ=0.00025); solve_flux(tfinal=1500.0, dt=0.01, record_every=200)"
    p = Flux1DProblem(2, 2, 10.0; N=8, K=100, ρL=-150.0, R=100.0, s=200.0, P=4,
                      ramp_τ=400.0, ramp_δ=0.00025)
    τs, Es, Eh = solve_flux(p; tfinal=1500.0, dt=0.01, record_every=200)
    js = sample_idx(length(τs), [0.3, 0.5, 0.7, 0.9, 1.0])
    op = orbital_params(10.0; M=1.0, a=0.0)
    open(joinpath(OUT, "td_flux1d.jl"), "w") do io
        header(io, "(b) Flux1D s=0 Schwarzschild circular-orbit (ℓ,m)=(2,2), rp=10: scri+horizon flux series and finals + selfforce_ft.", cfg)
        println(io, "const TD_FLUX1D = (")
        println(io, "    nsamp = $(length(τs)),")
        println(io, "    idx = $(ivec_str(js)),")
        println(io, "    tau = $(vec_str(τs[js])),")
        println(io, "    Escri = $(vec_str(Es[js])),")
        println(io, "    Ehor = $(vec_str(Eh[js])),")
        println(io, "    Escri_end = $(repr(Es[end])),")
        println(io, "    Ehor_end = $(repr(Eh[end])),")
        println(io, "    # construction probes (grounds tortoise/orbital_params/Ylm substitutions):")
        println(io, "    rsp = $(repr(p.rsp)),")
        println(io, "    kp = $(p.kp),")
        println(io, "    Gpre = $(repr(ComplexF64(p.Gpre))),")
        println(io, "    Omega = $(repr(p.Ω)),")
        println(io, "    ut = $(repr(op.ut)),")
        println(io, "    C_4_37 = $(repr(p.C[4, 37])),")
        println(io, "    C_scri = $(repr(p.C[end, end])),")
        println(io, ")")
    end
    @printf("[flux1d] Escri_end=%.8e Ehor_end=%.8e  done in %.1f s\n", Es[end], Eh[end], time() - t0); flush(stdout)

    t1 = time()
    cfg2 = "same Flux1DProblem; selfforce_ft(tfinal=1200.0, dt=0.01, twin=800.0)"
    ftL, ftR = selfforce_ft(p; tfinal=1200.0, dt=0.01, twin=800.0)
    open(joinpath(OUT, "td_selfforce.jl"), "w") do io
        header(io, "(b2) Flux1D dissipative self-force F_t^{22} at the particle (outflow/inflow nodes).", cfg2)
        println(io, "const TD_SELFFORCE = (ftL = $(repr(ftL)), ftR = $(repr(ftR)))")
    end
    @printf("[selfforce] ftL=%.10e ftR=%.10e  done in %.1f s\n", ftL, ftR, time() - t1); flush(stdout)
end

# ---------------------------------------------------------------------------------------------
# (c) KerrModes m=0 coupled tail, a=0.9, ℓs=[0,2,4], (2,0) seed (F1.4 class shrunk)
# ---------------------------------------------------------------------------------------------
let t0 = time()
    cfg = "KerrTailProblem([0,2,4], 0; N=8, K=120, ρL=-400.0, R=149.35, s=400.0, P=4, M=1, a=0.9, extract_rstar=[100.0]); solve_kerr_tail(init_mode=2, tfinal=3000.0, dt=0.04, record_every=100, ic=(amp=1/sqrt(2π*6^2), ρ0=25.0, w=6.0)); fit_tail_window (1200,3000)"
    p = KerrTailProblem([0, 2, 4], 0; N=8, K=120, ρL=-400.0, R=149.35, s=400.0, P=4,
                        a=0.9, extract_rstar=[100.0])
    τs, Ψ = solve_kerr_tail(p; init_mode=2, tfinal=3000.0, dt=0.04, record_every=100,
                            ic=(amp=1 / sqrt(2π * 6^2), ρ0=25.0, w=6.0))
    ns = [fit_tail_window(τs, Ψ[1, q, :], (1200.0, 3000.0)) for q in 1:3]
    js = sample_idx(length(τs), [0.05, 0.15, 0.3, 0.5, 0.75, 1.0])
    open(joinpath(OUT, "td_kerrtail.jl"), "w") do io
        header(io, "(c) KerrModes m=0 coupled Kerr tail a=0.9 ℓs=[0,2,4] (2,0)-seed: per-mode ψ(r*=100) samples + fit_tail_window exponents.", cfg)
        println(io, "const TD_KERRTAIL = (")
        println(io, "    nsamp = $(length(τs)),")
        println(io, "    idx = $(ivec_str(js)),")
        println(io, "    tau = $(vec_str(τs[js])),")
        println(io, "    psi_l0 = $(vec_str(Ψ[1, 1, js])),")
        println(io, "    psi_l2 = $(vec_str(Ψ[1, 2, js])),")
        println(io, "    psi_l4 = $(vec_str(Ψ[1, 3, js])),")
        println(io, "    n_l0 = $(repr(ns[1])),")
        println(io, "    n_l2 = $(repr(ns[2])),")
        println(io, "    n_l4 = $(repr(ns[3])),")
        println(io, "    # construction probes (grounds coupled-A/Sinv/V build via Spacetime tortoise):")
        println(io, "    probe_rho = $(vec_str(p.probe_ρ)),")
        println(io, "    Aflux_5_40 = $(vec_str(vec(p.Aflux[:, :, 5, 40]))),")
        println(io, "    Aflux_scri = $(vec_str(vec(p.Aflux[:, :, end, end]))),")
        println(io, "    Sinv_5_40 = $(vec_str(vec(p.Sinv[:, :, 5, 40]))),")
        println(io, "    Vmode_5_40 = $(vec_str(p.Vmode[:, 5, 40])),")
        println(io, "    # mode-coupling closed forms (App A):")
        println(io, "    cminus_20 = $(repr(cminus(2, 0))),")
        println(io, "    coup_l2 = $(vec_str(collect(mode_coupling(2, 0)))),")
        println(io, "    coup_l4 = $(vec_str(collect(mode_coupling(4, 0)))),")
        println(io, ")")
    end
    @printf("[kerrtail] n=[%.3f, %.3f, %.3f]  done in %.1f s\n", ns..., time() - t0); flush(stdout)
end

# ---------------------------------------------------------------------------------------------
# (d) KerrFlux a=0.9 circular-orbit Ė^{11}, ℓs=[1,3] (T1.3/F1.7 0.14%-class config shrunk)
# ---------------------------------------------------------------------------------------------
let t0 = time()
    cfg = "KerrFluxProblem([1,3], 1, 16.0914363989845; N=8, K=120, ρL=-300.0, R=100.0, s=300.0, P=4, M=1, a=0.9, q=1, ramp_τ=200.0, ramp_δ=0.001); solve_kerr_flux(tfinal=1400.0, dt=0.02, record_every=200, target=1)"
    rp = 16.0914363989845
    p = KerrFluxProblem([1, 3], 1, rp; N=8, K=120, ρL=-300.0, R=100.0, s=300.0, P=4,
                        a=0.9, ramp_τ=200.0, ramp_δ=0.001)
    τs, Es = solve_kerr_flux(p; tfinal=1400.0, dt=0.02, record_every=200, target=1)
    js = sample_idx(length(τs), [0.3, 0.5, 0.7, 0.9, 1.0])
    op = orbital_params(rp; M=1.0, a=0.9)
    open(joinpath(OUT, "td_kerrflux.jl"), "w") do io
        header(io, "(d) KerrFlux a=0.9 circular-orbit Ė^{11} (ℓs=[1,3], m=1, rp=16.0914363989845): scri flux series + final.", cfg)
        println(io, "const TD_KERRFLUX = (")
        println(io, "    nsamp = $(length(τs)),")
        println(io, "    idx = $(ivec_str(js)),")
        println(io, "    tau = $(vec_str(τs[js])),")
        println(io, "    Escri = $(vec_str(Es[js])),")
        println(io, "    Escri_end = $(repr(Es[end])),")
        println(io, "    # construction probes (grounds Kerr orbital_params + eigen-injection + Ylm):")
        println(io, "    kp = $(p.kp),")
        println(io, "    Omega = $(repr(p.Ω)),")
        println(io, "    E = $(repr(op.E)),")
        println(io, "    Lz = $(repr(op.Lz)),")
        println(io, "    ut = $(repr(op.ut)),")
        println(io, "    Gpre_re = $(vec_str(real.(p.Gpre))),")
        println(io, "    Gpre_im = $(vec_str(imag.(p.Gpre))),")
        println(io, "    Minjp = $(vec_str(vec(p.Minjp))),")
        println(io, "    Minjm = $(vec_str(vec(p.Minjm))),")
        println(io, "    mu_5_40_im = $(repr(imag(p.μmode[1, 5, 40]))),")
        println(io, ")")
    end
    @printf("[kerrflux] Escri_end=%.8e  done in %.1f s\n", Es[end], time() - t0); flush(stdout)
end

println("ALL CAPTURES DONE")
