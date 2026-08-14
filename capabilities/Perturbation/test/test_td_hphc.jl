# DG-D3 (Move D): Perturbation.TimeDomain HPHC spectral-collocation Teukolsky EVOLUTION backend —
# characterization vs the ORIGINAL TeukEvolution.jl `linear_field` run (Ripley,
# JLRipley314/TeukEvolution.jl @ 3c8e47f4). Grounding: test/fixtures/td_hphc_upstream.jl —
#   (A) HPHC_COEFF : Evo_lin_f coefficient arrays baked from the ORIGINAL constructor (pointwise
#       oracle for the {T}-generic closed forms);
#   (B) HPHC_EVO   : subsampled (f=ψ, p=∂_Tψ) snapshots at steps 0/200/400/600 of the a=0.7, s=−2,
#       m=2, nx=48, ny=24 run (the DECISIVE end-to-end oracle).
#
# GATES (each with pasted decisive lines in the D3 report):
#   G-D3a  DECISIVE fixture: solve the exact config, compare ALL 160 subsampled values across the
#          4 snapshots; pass = rel ≤ 1e-10 OR abs ≤ 1e-13·snapshot_scale per value.
#   G-D3b  spatial+temporal convergence: nx 48→96 vs the nx=192 reference at t≈1.24 (formula dt,
#          time-interpolated to a common target); 4th-order-consistent error drop (factor ≥ 10).
#   G-D3c  BigFloat-256 smoke: nx=24, ny=22, nt=20 runs end-to-end in Complex{BigFloat}, eltype
#          preserved, matches the same-config Float64 run to rel ≤ 1e-12. (ny raised 16→22: the
#          ported num_l keeps upstream max_s=2,max_m=6 ⇒ ny>20; 16 is below the collocation floor.)
#   G-D3d  a=0 sanity: short a=0 run stays bounded and the scri value is nonzero after propagation.
#
# Standalone-runnable (Move C G-C1 pattern): `julia --project include("test/test_td_hphc.jl")`.
using Test
using Printf
using Perturbation
using Perturbation.TimeDomain
const TDH = Perturbation.TimeDomain

include(joinpath(@__DIR__, "fixtures", "td_hphc_upstream.jl"))

_cabs(re, im) = hypot(re, im)

@testset "TimeDomain HPHC (Move D: TeukEvolution linear_field port)" begin

    # ---- G-D3 step 2: coefficient closed forms vs the baked ORIGINAL Evo_lin_f -------------------
    @testset "G-D3 coefficients (pointwise vs baked constructor)" begin
        C = HPHC_COEFF
        bhm, bhs, cl, spin, mv = C.bhm, C.bhs, C.cl, C.spin, C.mv
        worst = 0.0
        for r in C.rows
            R, cy, sy = r.R, r.cy, r.sy
            App = TDH.hphc_A_pp(R, bhm, bhs, cl)
            Apq = TDH.hphc_A_pq(R, bhm, bhs, cl)
            Bpp = TDH.hphc_B_pp(R, cy, bhm, bhs, cl, mv, spin)
            Bpq = TDH.hphc_B_pq(R, bhm, bhs, cl, mv, spin)
            Bpf = TDH.hphc_B_pf(R, bhm, bhs, cl, mv, spin)
            pre = TDH.hphc_pre(R, sy, bhm, bhs, cl)
            rel(a, b) = (m = abs(b); m == 0 ? abs(a - b) : abs(a - b) / m)
            for (a, b) in ((App, r.App), (Apq, r.Apq), (pre, r.pre),
                           (Bpp, complex(r.Bpp...)), (Bpq, complex(r.Bpq...)), (Bpf, complex(r.Bpf...)))
                worst = max(worst, rel(a, b))
            end
        end
        println("  [G-D3 coeff] worst rel deviation vs baked constructor = ", worst)
        @test worst <= 1e-12
    end

    # ---- G-D3a: DECISIVE end-to-end fixture ------------------------------------------------------
    @testset "G-D3a DECISIVE end-to-end vs upstream run" begin
        E = HPHC_EVO
        prob = HPHCTeukolskyProblem(; nx=E.nx, ny=E.ny, mv=E.mv, spin=E.spin,
            bhm=E.bhm, bhs=E.bhs, cl=E.cl, cfl=E.cfl,
            id_l_ang=E.id_l_ang, id_ru=E.id_ru, id_rl=E.id_rl, id_width=E.id_width,
            id_amp=E.id_amp, ingoing=E.ingoing, T=Float64)
        res = solve_hphc(prob; nt=E.nt, ts=E.ts)

        # per-(field,ti) scale = max modulus of the baked snapshot values
        scale = Dict{Tuple{Symbol,Int},Float64}()
        for r in E.rows
            k = (r.fld, r.ti); s = _cabs(r.re, r.im)
            scale[k] = max(get(scale, k, 0.0), s)
        end
        # worst deviations per snapshot
        worst = Dict{Tuple{Symbol,Int},NTuple{5,Float64}}()  # (relmax, absmax, i, j, scale)
        allpass = true
        for r in E.rows
            arr = r.fld === :lin_f ? res.f[r.ti] : res.p[r.ti]
            got = arr[r.i, r.j]
            baked = complex(r.re, r.im)
            ad = abs(got - baked)
            bm = abs(baked)
            rd = bm == 0 ? 0.0 : ad / bm
            k = (r.fld, r.ti); sc = scale[k]
            ok = (rd <= 1e-10) || (ad <= 1e-13 * sc)
            allpass &= ok
            w = get(worst, k, (0.0, 0.0, 0.0, 0.0, sc))
            if rd > w[1]
                worst[k] = (rd, ad, Float64(r.i), Float64(r.j), sc)
            end
        end
        println("  [G-D3a] worst per-snapshot deviations (pass = rel≤1e-10 OR abs≤1e-13·scale):")
        for fld in (:lin_f, :lin_p), ti in 1:4
            w = worst[(fld, ti)]
            @printf("    %-6s ti=%d t=%.6f  worst rel=%.3e  abs=%.3e at (i=%d,j=%d)  scale=%.3e\n",
                    string(fld), ti, E.times[ti], w[1], w[2], Int(w[3]), Int(w[4]), w[5])
        end
        @test allpass
    end

    # ---- G-D3b: convergence (4th-order-consistent error drop) ------------------------------------
    @testset "G-D3b convergence nx 48→96 vs nx=192 reference" begin
        ttar = 1.2412381986353238   # fixture ti=2 time
        function scri_at(nx)
            p = HPHCTeukolskyProblem(; nx=nx, ny=24, mv=2, spin=-2, bhm=1.0, bhs=0.7, cl=1.0,
                cfl=0.5, id_l_ang=2, id_ru=4.0, id_rl=2.5, id_width=1.5,
                id_amp=complex(0.1, 0.0), ingoing=true, T=Float64)
            k = floor(Int, ttar / p.dt)
            r = solve_hphc(p; nt=k + 1, ts=1)          # dense record; postfilter runs every step
            θ = (ttar - k * p.dt) / p.dt                # linear time-interp to the common target
            f0 = r.f[k + 1]; f1 = r.f[k + 2]           # snapshot idx = step+1
            return [(1 - θ) * f0[1, j] + θ * f1[1, j] for j in 1:24]   # scri (i=1), all j
        end
        s48 = scri_at(48); s96 = scri_at(96); s192 = scri_at(192)
        err48 = maximum(abs.(s48 .- s192))
        err96 = maximum(abs.(s96 .- s192))
        ratio = err48 / err96
        @printf("  [G-D3b] err(48 vs 192)=%.3e  err(96 vs 192)=%.3e  ratio=%.2f (4th-order≈16)\n",
                err48, err96, ratio)
        @test ratio >= 10
    end

    # ---- G-D3c: BigFloat-256 smoke ---------------------------------------------------------------
    @testset "G-D3c BigFloat-256 smoke (precision inherited)" begin
        cfg = (; nx=24, ny=22, mv=2, spin=-2, bhm=1.0, bhs=0.7, cl=1.0, cfl=0.5,
               id_l_ang=2, id_ru=4.0, id_rl=2.5, id_width=1.5, id_amp=complex(0.1, 0.0),
               ingoing=true)
        rf = solve_hphc(HPHCTeukolskyProblem(; cfg..., T=Float64); nt=20, ts=20)
        rb = setprecision(BigFloat, 256) do
            solve_hphc(HPHCTeukolskyProblem(; cfg..., T=BigFloat); nt=20, ts=20)
        end
        @test eltype(rb.f[end]) === Complex{BigFloat}
        ffb = rb.f[end]; ffd = rf.f[end]; ppb = rb.p[end]; ppd = rf.p[end]
        sc = max(maximum(abs, ffb), maximum(abs, ppb))
        devf = maximum(abs.(ComplexF64.(ffb) .- ffd))
        devp = maximum(abs.(ComplexF64.(ppb) .- ppd))
        dev = max(devf, devp) / Float64(sc)
        @printf("  [G-D3c] max |BigFloat−Float64| / scale = %.3e  (scale=%.3e)\n", dev, Float64(sc))
        @test dev <= 1e-12
    end

    # ---- G-D3d: a=0 sanity -----------------------------------------------------------------------
    @testset "G-D3d a=0 sanity (bounded + scri nonzero)" begin
        p = HPHCTeukolskyProblem(; nx=48, ny=24, mv=2, spin=-2, bhm=1.0, bhs=0.0, cl=1.0,
            cfl=0.5, id_l_ang=2, id_ru=4.0, id_rl=2.5, id_width=1.5,
            id_amp=complex(0.1, 0.0), ingoing=true, T=Float64)
        r = solve_hphc(p; nt=200, ts=200)
        ff = r.f[end]
        @test all(isfinite, ff)
        @test all(isfinite, r.p[end])
        @test maximum(abs, ff) < 10.0                       # bounded (initial peak amp = 0.1)
        scri = maximum(abs, @view ff[1, :])                 # scri = R=0 (i=1)
        @printf("  [G-D3d] a=0: max|f|=%.3e  scri|f(i=1)|=%.3e  (nonzero ⇒ propagated)\n",
                maximum(abs, ff), scri)
        @test scri > 1e-10
    end
end
