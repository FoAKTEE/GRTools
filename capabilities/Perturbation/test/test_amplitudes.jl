# Module 1 (AmplitudesFluxes leaf): MST asymptotic amplitudes _mst_amplitudes.
# GROUNDING: LIVE Teukolsky oracle. For each (s,l,m,a,ω) we build the legacy radial solution
# R = Teukolsky.TeukolskyRadial(...), read its renormalized ν + eigenvalue λ, feed BOTH into the
# new Perturbation._mst_amplitudes, and assert the normalized ratios (Inc=InInc/InTrans,
# Ref=InRef/InTrans, …) match the legacy radial amplitude accessors R.In/Up.amplitudes at the
# legacy tolerance (Teukolsky/test/runtests.jl lines 243-300: In rtol 1e-7, Up rtol 1e-7).
# Legacy oracle values (ν,λ inputs + amplitude RHS) BAKED once into fixtures/leg_amplitudes.jl
# (Phase 7-B); the GrTool side + tolerances are unchanged.
include("fixtures/leg_amplitudes.jl")  # REF_AMP1 / REF_AMP2 / REF_AMP3

@testset "AmplitudesFluxes :: _mst_amplitudes vs Teukolsky radial amplitudes (baked, captured once)" begin
    grid = [(-2, 2, 2, 0.0, 0.1), (-2, 2, 2, 0.5, 0.2), (-2, 2, 1, 0.9, 0.15),
            (-2, 3, 3, 0.7, 0.25), (0, 2, 2, 0.5, 0.2), (-1, 2, 2, 0.6, 0.2),
            (-1, 2, 1, 0.9, 0.2)]
    for (s, l, m, a, ω) in grid
        leg = REF_AMP1[(s, l, m, a, ω)]            # legacy MST oracle (baked)
        ν = leg.nu; λ = leg.lam
        amp = Perturbation._mst_amplitudes(s, m, a, ω, ν, λ; Nmax=80)
        # In side: Incidence = Binc/Btrans, Reflection = Bref/Btrans
        @test isapprox(amp.InInc / amp.InTrans, leg.inInc; rtol=1e-7)
        @test isapprox(amp.InRef / amp.InTrans, leg.inRef; rtol=1e-7)
        # Up side: Incidence = Cinc/Ctrans, Reflection = Cref/Ctrans
        @test isapprox(amp.UpInc / amp.UpTrans, leg.upInc; rtol=1e-7)
        @test isapprox(amp.UpRef / amp.UpTrans, leg.upRef; rtol=1e-7)
    end
end

@testset "AmplitudesFluxes :: _mst_amplitudes matches legacy _mst_amplitudes (full set, baked)" begin
    # Direct function-level characterization at shared (ν,λ): every amplitude key, tight tol.
    for (s, l, m, a, ω) in [(-2, 2, 2, 0.5, 0.2), (0, 2, 2, 0.5, 0.2), (-1, 2, 1, 0.9, 0.2)]
        leg = REF_AMP2[(s, l, m, a, ω)]
        λ = leg.lam; ν = leg.nu
        new = Perturbation._mst_amplitudes(s, m, a, ω, ν, λ; Nmax=80)
        for k in (:InTrans, :InInc, :InRef, :UpTrans, :UpInc, :UpRef, :Aplus, :Aminus)
            @test isapprox(getfield(new, k), getfield(leg, k); rtol=1e-10)
        end
    end
end

@testset "AmplitudesFluxes :: BigFloat amplitude path matches Float64 branch" begin
    # Preserve BigFloat inheritance: Complex{BigFloat}→ComplexF64 matches the Float64 complex branch.
    s, l, m, a, ω = -2, 2, 2, 0.5, 0.2
    λ64 = REF_AMP3.lam64
    ν64 = REF_AMP3.nu64
    amp64 = Perturbation._mst_amplitudes(s, m, a, ω, ν64, λ64; Nmax=80)
    ampbig = Perturbation._mst_amplitudes(s, m, big(a), big(ω),
                                          Complex{BigFloat}(ν64), Complex{BigFloat}(λ64); Nmax=80)
    for k in (:InTrans, :InInc, :InRef, :UpTrans, :UpInc, :UpRef)
        @test isapprox(ComplexF64(getfield(ampbig, k)), getfield(amp64, k); rtol=1e-10)
    end
end
