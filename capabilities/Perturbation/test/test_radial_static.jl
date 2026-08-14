# Grounding for Perturbation.Radial static (ω=0) TeukolskyRadial closed-form 2F1 path
# (teukolsky_radial_static.jl + the minimal TeukolskyRadialFunction carrier).
# Oracle: ref_teuk_radial_static.jl = BHPToolkit static TeukolskyRadial table (90 ok rows, 6 r each)
# frozen as data. The legacy Teukolsky test asserts TeukolskyRadial(s,l,m,a,0.0) radial values +
# static amplitudes against exactly this table at rtol=1e-9/atol=1e-12; reproduced for the rewritten
# form (Pochhammer/Γ delegated to BigSpecial/SpecialFunctions). Plus the Schwarzschild a=0 physical
# identity and a live-Teukolsky cross-check.

include("ref_teuk_radial_static.jl")  # const REF_RADSTATIC
include("fixtures/leg_radial_static.jl")  # const REF_RADSTATIC_LEG (baked Teukolsky In/Up, Phase 7-B)

@testset "Radial: static (ω=0) TeukolskyRadial vs BHPToolkit oracle (90 rows, closed-form 2F1)" begin
    nchecked = 0; maxRerr = 0.0
    for r in REF_RADSTATIC
        r.ok == 1 || continue
        R = Perturbation.TeukolskyRadial(r.s, r.l, r.m, r.a, 0.0)
        @test R.In.bc === :In && R.Up.bc === :Up
        @test R.In.omega == 0.0 && R.Up.omega == 0.0
        @test R.In.lambda == r.l * (r.l + 1) - r.s * (r.s + 1)   # λ = l(l+1)-s(s+1)
        for j in eachindex(r.r)
            rin = R.In(r.r[j]); rup = R.Up(r.r[j])
            rin_o = complex(r.inre[j], r.inim[j]); rup_o = complex(r.upre[j], r.upim[j])
            @test isfinite(rin) && isfinite(rup)                # no Γ-pole NaNs
            @test isapprox(rin, rin_o; rtol=1e-9, atol=1e-12)
            @test isapprox(rup, rup_o; rtol=1e-9, atol=1e-12)
            maxRerr = max(maxRerr, abs(rin - rin_o) / max(abs(rin_o), 1e-300),
                                   abs(rup - rup_o) / max(abs(rup_o), 1e-300))
            nchecked += 2
        end
        # static amplitudes: assert each MMA key flagged present (present-flag p==1)
        inA = R.In.amplitudes; upA = R.Up.amplitudes
        for (amps, pfx) in ((inA, "in"), (upA, "up"))
            for (fld, key) in ((:H, "H"), (:I, "I"), (:Hm, "Hm"), (:Hp, "Hp"), (:Im, "Im"), (:Ip, "Ip"))
                getfield(r, Symbol("$(pfx)_$(key)_p")) == 1 || continue
                @test haskey(amps, fld)
                o = complex(getfield(r, Symbol("$(pfx)_$(key)_re")),
                            getfield(r, Symbol("$(pfx)_$(key)_im")))
                @test isapprox(getfield(amps, fld), o; rtol=1e-9, atol=1e-12)
            end
        end
    end
    nok = count(r -> r.ok == 1, REF_RADSTATIC)
    @test nok == 90
    @test nchecked == 12 * nok                                  # 6 r × 2 BC per ok row
    @info "static radial: rows=$nok radial_checks=$nchecked maxRerr=$maxRerr"
end

@testset "Radial: Schwarzschild a=0 static physical identities" begin
    # τ=0: R^In regular & growing outward, R^Up decaying at large r; both finite (no Γ-pole NaN).
    Sch = Perturbation.TeukolskyRadial(-2, 2, 0, 0.0, 0.0)
    rr = (3.0, 10.0, 100.0)
    inv = [real(Sch.In(r)) for r in rr]; upv = [real(Sch.Up(r)) for r in rr]
    @test all(isfinite, inv) && all(isfinite, upv)
    @test inv[1] < inv[2] < inv[3]          # In grows outward (regular at horizon side)
    @test upv[1] > upv[2] > upv[3] > 0      # Up decays at large r
end

@testset "Radial: new static radial == legacy Teukolsky (cross-check, baked exact from Teukolsky.TeukolskyRadial)" begin
    for r in REF_RADSTATIC
        r.ok == 1 || continue
        Rnew = Perturbation.TeukolskyRadial(r.s, r.l, r.m, r.a, 0.0)
        leg = REF_RADSTATIC_LEG[(r.s, r.l, r.m, r.a)]
        for j in eachindex(r.r)
            @test Rnew.In(r.r[j]) == leg.legin[j]
            @test Rnew.Up(r.r[j]) == leg.legup[j]
        end
    end
end
