# Grounding for Perturbation.Radial raw (norm=1) MST homogeneous solutions R^in/R^up + derivs.
# Oracles: ref_teuk_radial.jl (REF_RADIAL, 5 modes × 4 r) and ref_teuk_radial_deriv.jl (REF_RADDERIV)
# = BHPToolkit MST R^in/R^up and dR^in/dr, dR^up/dr frozen as data. The legacy Teukolsky test
# asserts its _mst_radial_in/_up(...; norm=1) and the deriv variants against exactly these tables at
# rtol=1e-8/atol=1e-9; reproduced here for the rewritten form (which delegates the f_Up Pochhammer
# rescale to BigSpecial.pochhammer). Plus an independent finite-difference cross-check of the derivs
# vs the values and a live-Teukolsky cross-check.

include("ref_teuk_radial.jl")        # const REF_RADIAL
include("ref_teuk_radial_deriv.jl")  # const REF_RADDERIV
include("fixtures/leg_radial_radial.jl")  # const REF_RADRAW_LEG (baked Teukolsky._mst_radial_*, Phase 7-B)

@testset "Radial: R^in(r) vs BHPToolkit oracle (5 modes × 4 r, norm=1)" begin
    for r in REF_RADIAL
        for j in eachindex(r.r)
            Rin = Perturbation._mst_radial_in(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, r.r[j]; norm=1.0)
            @test isapprox(real(Rin), r.inre[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(Rin), r.inim[j]; rtol=1e-8, atol=1e-9)
        end
    end
end

@testset "Radial: R^up(r) vs BHPToolkit oracle (5 modes × 4 r, norm=1)" begin
    for r in REF_RADIAL
        for j in eachindex(r.r)
            Rup = Perturbation._mst_radial_up(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, r.r[j]; norm=1.0)
            @test isapprox(real(Rup), r.upre[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(Rup), r.upim[j]; rtol=1e-8, atol=1e-9)
        end
    end
end

@testset "Radial: dR^in/dr vs BHPToolkit oracle (5 modes × 4 r, norm=1)" begin
    for r in REF_RADDERIV
        for j in eachindex(r.r)
            dRin = Perturbation._mst_radial_in_deriv(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, r.r[j]; norm=1.0)
            @test isapprox(real(dRin), r.dinre[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(dRin), r.dinim[j]; rtol=1e-8, atol=1e-9)
        end
    end
end

@testset "Radial: dR^up/dr vs BHPToolkit oracle (5 modes × 4 r, norm=1)" begin
    for r in REF_RADDERIV
        for j in eachindex(r.r)
            dRup = Perturbation._mst_radial_up_deriv(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, r.r[j]; norm=1.0)
            @test isapprox(real(dRup), r.dupre[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(dRup), r.dupim[j]; rtol=1e-8, atol=1e-9)
        end
    end
end

@testset "Radial: analytic finite-difference cross-check of derivs vs values" begin
    # dR/dr ≈ [R(r+h)−R(r−h)]/(2h): an independent (oracle-free) consistency lens.
    h = 1e-5
    for r in REF_RADDERIV
        for rr in r.r
            Rin  = Perturbation._mst_radial_in
            Rup  = Perturbation._mst_radial_up
            dInFD = (Rin(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr + h) -
                     Rin(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr - h)) / (2h)
            dUpFD = (Rup(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr + h) -
                     Rup(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr - h)) / (2h)
            dIn = Perturbation._mst_radial_in_deriv(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr)
            dUp = Perturbation._mst_radial_up_deriv(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr)
            @test isapprox(dIn, dInFD; rtol=1e-6, atol=1e-7)
            @test isapprox(dUp, dUpFD; rtol=1e-6, atol=1e-7)
        end
    end
end

@testset "Radial: new R^in/R^up/derivs ≈ legacy Teukolsky (cross-check, baked, captured once)" begin
    # rtol=1e-12: the f_Up Pochhammer rescale associates its product differently from the legacy
    # inline form (BigSpecial.pochhammer vs the local copy), a ULP-level reorder — agreement to
    # 1e-12 confirms identical math; the BHPToolkit oracle at 1e-8 above is the load-bearing ground.
    for r in REF_RADDERIV
        leg = REF_RADRAW_LEG[(r.s, r.l, r.m, r.a, r.w)]
        for (j, rr) in enumerate(r.r)
            @test isapprox(Perturbation._mst_radial_in(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr; norm=1.0),
                  leg.rin[j]; rtol=1e-12)
            @test isapprox(Perturbation._mst_radial_up(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr; norm=1.0),
                  leg.rup[j]; rtol=1e-12)
            @test isapprox(Perturbation._mst_radial_in_deriv(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr; norm=1.0),
                  leg.rind[j]; rtol=1e-12)
            @test isapprox(Perturbation._mst_radial_up_deriv(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, rr; norm=1.0),
                  leg.rupd[j]; rtol=1e-12)
        end
    end
end
