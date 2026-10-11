# Fresh-process dispatch-isolation gate for the end-to-end eccentric solvers.
# Started only by `isolation_run(kind)` in test_eccentric_modes.jl: it redefines
# kernels to throw, so it must never run inside a shared test process.
#   circular  e=0 adiabatic solve with every circular source, master-entry and
#             scalar-Doppler selector method poisoned (ECC::gr-mode-solver)
#   li        OPTION A solves with the adiabatic and circular methods poisoned
#             (ECC::gr-li-mode-solver)
using Test, LinearAlgebra, Printf
const ECC_MODES_HELPERS_ONLY = true
include(joinpath(@__DIR__, "test_eccentric_modes.jl"))
using AccretionDisc, Spacetime
using AccretionDisc.ForcingModes
const SW = AccretionDisc.SpiralWaves
const KIND = only(ARGS)
KIND in ("circular", "li") || error("expected circular or li")

function poison_circular()
    @eval SW begin
        _master_entries(x::MasterInputs) = error("POISON circular generated entries")
        master_inputs(bg::DiscBackground, m::Integer, r; kw...) = error("POISON circular master_inputs")
        singular_field(bg::DiscBackground, r, Ψ; kw...) = error("POISON circular singular_field")
        singular_field_mode(bg::DiscBackground, m::Integer, r; kw...) = error("POISON circular singular_field_mode")
        source_modes(bg::DiscBackground, m::Integer, rgrid::AbstractVector; kw...) = error("POISON circular source_modes")
        _kernel_integrals(K::SingularFieldCoeffs, m::Integer, r::S, N::Integer) where {S<:Real} = error("POISON circular kernel")
        _singular_field_mode(K::SingularFieldCoeffs, m::Integer, r::S, N::Union{Nothing,Integer}, mE::S) where {S<:Real} =
            error("POISON circular mode kernel")
        wkb_split(C::AbstractMatrix{TC}, ωd::Real; kw...) where {TC<:Number} = error("POISON scalar-Doppler selector")
        _away_is_first(k1, k2, ωd::Real; kw...) = error("POISON scalar-Doppler branch test")
    end
end

if KIND == "circular"
    cases = map(Iterators.product((-0.4, 0.4), (2, 3))) do (a, m)
        bg = SW.DiscBackground(χ=a, aEMRI=20.0, h0=0.08, 𝔭=1.5, 𝔮=0.75, b=0.3, εp=1e-3)
        d = SW.DiscModel(bg; closure=Adiabatic, evolution=:fixed)
        p = SW.EccentricPerturber(EquatorialOrbit(a, 20.0, 0.0); mass=1e-5,
                                  softening=FixedSoftening(SW.softening_length(bg)))
        cg = SW.build_grid(bg, m; nL=1, nR=1, nC=65, n_anc=64, extent=:fiducial, r_in=8.0, r_out=45.0)
        circ = SW.solve_mode(bg, m; grid=cg, mEMRI=p.mass, linsolve=:banded)
        (; a, m, bg, d, p, cg, circ, du=SW.delta_u_phi(bg, circ))
    end
    poison_circular()
    @testset "ECC e=0 eccentric path with circular kernels poisoned" begin
        x = first(cases)
        @test_throws ErrorException SW.solve_mode(x.bg, x.m; grid=x.cg, mEMRI=x.p.mass, linsolve=:banded)
        @test_throws ErrorException SW.source_modes(x.bg, x.m, [20.0])
        @test_throws ErrorException SW.singular_field_mode(x.bg, x.m, 20.5)
        @test_throws ErrorException SW.outgoing_branch(ComplexF64[0 1; -1 0], 1.0; outward=true)
        for x in cases
            c = SW.ModeContext(x.d, x.p, x.m, 0; closure=Adiabatic, domain=(first(x.cg.r), last(x.cg.r)))
            eg = SW.EccentricGrid(copy(x.cg.r), c.r_c, 1:0, SW.mode_cache_key(c), (;circular_gate=true))
            ecc = SW.solve_mode(c, x.p; grid=eg, bank_options=(;Nq=2))
            err = norm(ecc.y - x.circ.y)/norm(x.circ.y)
            # Independent reconstruction route: circular route-R2 delta u_phi.
            duerr = norm([f.du_phi for f in SW.reconstruct_fluid(ecc)] - x.du)/norm(x.du)
            @test err < 2e-8 && duerr < 2e-8
            @printf("ISOLATED_CIRCULAR,%.2f,%d,%.12e,%.12e\n", x.a, x.m, err, duerr)
        end
    end
else
    poison_circular()
    @eval SW begin
        _context_operator(::Val{Adiabatic}, x) = error("POISON adiabatic operator")
        _context_source(::Val{Adiabatic}, x) = error("POISON adiabatic source")
        _context_reconstruction(::Val{Adiabatic}, x) = error("POISON adiabatic reconstruction")
    end
    c, _ = EccentricModeTestSupport.fixture()
    @test_throws ErrorException SW.master_operator(c, 30.0; landau=false)
    ENV["ECC_LI_DISPATCH_ONLY"] = "true"
    include(joinpath(@__DIR__, "test_eccentric_li_modes.jl"))
end
