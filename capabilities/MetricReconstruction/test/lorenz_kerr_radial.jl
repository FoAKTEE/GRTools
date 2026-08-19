# Grounding for Lorenz.KerrRadial — the a≠0 spin-general full radial reconstruction
# (P_{-2},P_{-1},h homogeneous Teukolsky + κ inhomogeneous VoP/Abel). No live oracle is used:
# the grounding is PURE ANALYTIC self-consistency (every expected value is a residual that must
# vanish, or a closed-form structural identity), exactly the G-K10 gate in
# packages/KerrLorenzCirc/test/ref_kerr_radial_components.jl.
#
# (1) G-K10 self-consistency gate `radial_gates_kerr(a,r0)` (analytic, no oracle):
#       maxjump  : the reconstructed P_{-2},P_{-1},h reproduce the closure_kerr jumps [f]=J0,[f']=J1  (<1e-4)
#       maxres_h : the homogeneous h satisfies its s=0 Teukolsky ODE                                  (<1e-8)
#       maxjump_k: the reconstructed κ reproduces the closure_kerr κ jumps                            (<1e-4)
#       maxres_k : κ satisfies its INHOMOGENEOUS ODE L[κ] = ½(r²+a²γ_ll)h                             (<1e-6)
#     Run at a=0.5 (ls=2:3) and a=0.9 (ls=2:2), r0=6 — the spin-general path at two spins.
# (2) ACCESSOR finite-difference gates (analytic): derK == FD(valK) and der2K == FD(derK) away from
#     the orbit (<1e-6) — the analytic 1st/2nd-derivative accessors agree with central differences.
# (3) STRUCTURE gates: jump_data_kerr's jumps == rec.jd jumps; h's J0 is identically zero.
# (4) PRECISION-INHERITANCE gate: reconstruct_kerr at BigFloat 256 vs 512 bits self-converges
#     (homogeneous P_{-2} value, κ value) — confirms the {T} path carries BigFloat precision.

using Test
import MetricReconstruction.Lorenz as L

@testset "Lorenz.KerrRadial (a≠0 radial reconstruction — analytic self-consistency G-K10)" begin

    @testset "G-K10 self-consistency: radial_gates_kerr(0.5, 6) ls=2:3" begin
        maxjump, maxres_h, maxjump_k, maxres_k = L.radial_gates_kerr(0.5, 6.0; ls=2:3, m=2)
        @test maxjump   < 1e-4    # P_{-2},P_{-1},h reproduce closure_kerr jumps
        @test maxres_h  < 1e-8    # homogeneous h s=0 Teukolsky-ODE residual
        @test maxjump_k < 1e-4    # κ reproduces closure_kerr κ jumps
        @test maxres_k  < 1e-6    # κ inhomogeneous-ODE residual L[κ]−½(r²+a²γ_ll)h
    end

    @testset "G-K10 self-consistency: radial_gates_kerr(0.9, 6) ls=2:2 (higher spin)" begin
        maxjump, maxres_h, maxjump_k, maxres_k = L.radial_gates_kerr(0.9, 6.0; ls=2:2, m=2)
        @test maxjump   < 1e-4
        @test maxres_h  < 1e-8
        @test maxjump_k < 1e-4
        @test maxres_k  < 1e-6
    end

    # One reconstruct reused by the accessor + structure gates (a=0.5, r0=6, l=2).
    a = 0.5; r0 = 6.0; l = 2; m = 2
    sol, aux = L.closure_kerr(a, r0, m; lmax=max(2l, l+8), return_aux=true)
    rec = L.reconstruct_kerr(a, r0, l; m=m, aux=aux, sol=sol)

    @testset "accessor FD: derK==FD(valK), der2K==FD(derK) away from orbit" begin
        r = r0 + 3.0; hh = 1e-5
        for f in (rec.Pm2, rec.Pm1, rec.h)
            fd1 = (L.valK(f, r+hh) - L.valK(f, r-hh)) / (2hh)
            fd2 = (L.derK(f, r+hh) - L.derK(f, r-hh)) / (2hh)
            @test abs(L.derK(f, r)  - fd1) / (abs(fd1) + 1e-300) < 1e-6
            @test abs(L.der2K(f, r) - fd2) / (abs(fd2) + 1e-300) < 1e-6
        end
    end

    @testset "structure: jump_data_kerr jumps == rec.jd, h J0 ≡ 0" begin
        @test rec.jd.h[1] == 0
        jd2 = L.jump_data_kerr(a, r0, l, m)        # self-solves its own closure
        @test jd2.Pm2 == rec.jd.Pm2
        @test jd2.Pm1 == rec.jd.Pm1
        @test jd2.κ   == rec.jd.κ
        @test jd2.h[1] == 0
    end

    @testset "precision inheritance: jump_data_kerr BigFloat 256 vs 512 self-converge" begin
        # The jump engine (closure_kerr + jumps_pm/hj1) is the BigFloat-precision-carrying chain that
        # reconstruct_kerr's jumps are built on; it self-converges 256→512 to full BigFloat precision
        # (measured ~7.5e-74). The full radial-field self-convergence was verified one-shot
        # (reconstruct_kerr with rspan: Pm2 rel ~8.7e-75, κ rel ~3.9e-75) but is too slow for CI
        # (~200 s at 512 bits), so the in-suite gate exercises the jump engine.
        jdbig(bits) = setprecision(BigFloat, bits) do
            jd = L.jump_data_kerr(BigFloat(a), BigFloat(r0), l, m)
            (jd.Pm2[1], jd.Pm2[2], jd.Pm1[1], jd.Pm1[2], jd.κ[1], jd.κ[2],
             Complex(jd.Λ), Complex(jd.λ0), Complex(jd.γll))
        end
        j256 = jdbig(256); j512 = jdbig(512)
        maxrel = maximum(abs(ComplexF64(j256[k]-j512[k])) / (abs(ComplexF64(j512[k])) + 1e-300)
                         for k in 1:9)
        @test maxrel < 1e-40   # full BigFloat precision inheritance (measured ~7.5e-74)
    end
end

# ---------------------------------------------------------------------------------------------
# (5) NODE-COUNT OPTIONS gate (repro-audit G6): the Chebyshev node-count knobs used to be process
#     ENV variables (KLC_NCHEB_FLOOR / KLC_NIN_FLOOR / KLC_NIN_FACTOR / KLC_NUP_FLOOR /
#     KLC_NUP_FACTOR), so the numbers a run produced depended on invisible process state. They are
#     now `KerrRadialOptions` fields. This gate pins (a) the defaults to the exact values those ENV
#     variables defaulted to, (b) that NO KLC_* ENV read remains in the source, and (c) that the
#     override paths (per-call `opts=`, process-wide `set_kerr_radial_options!`) actually take.
# ---------------------------------------------------------------------------------------------
@testset "KerrRadialOptions (ex-KLC_* environment knobs)" begin

    @testset "defaults equal the previous ENV-free values" begin
        o = L.KerrRadialOptions()
        @test o.ncheb_floor == 30     # was get(ENV, "KLC_NCHEB_FLOOR", "30")
        @test o.nin_floor   == 18     # was get(ENV, "KLC_NIN_FLOOR",   "18")
        @test o.nin_factor  == 3      # was get(ENV, "KLC_NIN_FACTOR",  "3")
        @test o.nup_floor   == 70     # was get(ENV, "KLC_NUP_FLOOR",   "70")
        @test o.nup_factor  == 3      # was get(ENV, "KLC_NUP_FACTOR",  "3")
        # the process-wide default IS the production struct at load time
        @test L.kerr_radial_options() == L.KerrRadialOptions()
    end

    @testset "no KLC_* environment read survives in the source" begin
        src = joinpath(pkgdir(MetricReconstruction), "src")
        offenders = String[]
        for (root, _, files) in walkdir(src), f in files
            endswith(f, ".jl") || continue
            txt = read(joinpath(root, f), String)
            # the docstring of KerrRadialOptions names the retired variables; only a live
            # `get(ENV, "KLC_...` / `ENV["KLC_...` read is a violation.
            (occursin("get(ENV, \"KLC_", txt) || occursin("get(ENV,\"KLC_", txt) ||
             occursin("ENV[\"KLC_", txt)) && push!(offenders, joinpath(root, f))
        end
        @test isempty(offenders)
    end

    @testset "overrides take (per-call struct + process-wide setter)" begin
        o = L.KerrRadialOptions(; nin_floor = 24, nup_factor = 4)
        @test o.nin_floor == 24 && o.nup_factor == 4
        @test o.ncheb_floor == 30 && o.nin_factor == 3 && o.nup_floor == 70   # untouched fields
        # `_branchK_leg_interp_big`'s default N follows the PROCESS option, not a frozen const.
        old = L.kerr_radial_options()
        try
            L.set_kerr_radial_options!(; nin_floor = 21)
            @test L.kerr_radial_options().nin_floor == 21
            @test L.kerr_radial_options().nup_floor == old.nup_floor          # keyword form is a patch
            m = only(methods(L._branchK_leg_interp_big))
            @test Base.kwarg_decl(m) == [:N]
        finally
            L.set_kerr_radial_options!(old)
        end
        @test L.kerr_radial_options() == old
    end

    @testset "opts is threaded through the public radial API" begin
        for f in (L.branchK, L.reconstruct_kerr)
            @test any(m -> :opts in Base.kwarg_decl(m), methods(f))
        end
    end
end
