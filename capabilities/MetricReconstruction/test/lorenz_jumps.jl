# Grounding for Lorenz.Jumps — App-H spin-general source jumps (the ORIGINAL gold-standard oracle).
#
# (1) PUBLISHED App-H upper P_{-2} jump table (gate G-K1): jumps_pm(bg,l,−1) J0/J1 at a=0 & a=0.99,
#     r0=6, m=2, l=2..10, rel < 1e-4. Literals baked verbatim from packages/KerrLorenzCirc/test/
#     ref_apph_jumps.jl (the published App-H table).
# (2) PUBLISHED App-H trace J1(h) table (gate G-K5): hj1(bg,l), baked from ref_apph_trace.jl.
# (3) SWSH S0 sub-call gated against Harmonics (the already-built oracle).

using Test
import MetricReconstruction.Lorenz as L
import Harmonics

# App-H upper P_{-2} jumps (published; ref_apph_jumps.jl appH_upper)
const appH_upper = Dict(
 0.0 => Dict(
   2=>(-0.9341652-32.035193im,-14.012478-10.106341im), 3=>(2.2106383+43.319488im,-19.895745+13.086095im),
   4=>(-1.2533141+42.979722im,77.705477+15.094069im),   5=>(-1.3855909-88.243759im,37.410955-28.566089im),
   6=>(1.6004396-44.73709im,-181.03797-15.892428im),    7=>(1.0112635+128.80814im,-53.596965+42.316776im),
   8=>(-1.7225155+45.357304im,322.32572+16.173922im),   9=>(-0.79669766-167.82924im,69.312696-55.455203im),
   10=>(1.7798579-45.64908im,-501.39643-16.306296im)),
 0.99 => Dict(
   2=>(-0.21782443-16.92872im,-15.55513-7.3009438im),   5=>(-0.32786372-47.812584im,10.89643-20.957857im),
   8=>(-0.43162509+19.345375im,248.99391+9.1424079im),  10=>(0.44368914-19.35722im,-383.60261-9.1618631im)))

# App-H trace jump J1(h) (published; ref_apph_trace.jl appHh)
const appHh = Dict(0.0 =>Dict(2=>-0.57205702,4=>0.49541591,6=>-0.48234195,8=>0.47768153,10=>-0.47548123),
                   0.99=>Dict(2=>-0.58007105,4=>0.50246972,6=>-0.48919802,8=>0.48446434,10=>-0.48222893))

@testset "Lorenz.Jumps (App-H gold-standard oracle)" begin

    @testset "G-K1: App-H upper P_{-2} jumps J0/J1 (rel < 1e-4)" begin
        mx = 0.0
        for a in (0.0, 0.99)
            bg = L.orbit_background(a, 6.0, 2)
            ls = a == 0.0 ? (2:10) : (2,5,8,10)
            for l in ls
                J0, J1 = L.jumps_pm(bg, l, -1); (r0, r1) = appH_upper[a][l]
                e0 = abs(J0 - r0)/abs(r0); e1 = abs(J1 - r1)/abs(r1)
                @test e0 < 1e-4
                @test e1 < 1e-4
                mx = max(mx, e0, e1)
            end
        end
        @test mx < 1e-4
    end

    @testset "G-K5: App-H trace jump J1(h) (rel < 1e-4)" begin
        for a in (0.0, 0.99)
            bg = L.orbit_background(a, 6.0, 2)
            for l in (2,4,6,8,10)
                J = L.hj1(bg, l); ref = appHh[a][l]
                @test abs(real(J) - ref)/abs(ref) < 1e-4
            end
            # l+m odd ⇒ J1(h) ≈ 0
            @test abs(L.hj1(bg, 3)) < 1e-8
        end
    end

    @testset "jumps_eigenvalue == Harmonics SWSH s=−2 eigenvalue" begin
        bg = L.orbit_background(0.99, 6.0, 2)
        for l in (2,3,5)
            ref = real(Harmonics.SpinWeightedSpheroidalEigenvalue(-2, l, 2, Float64(bg.aω)))
            @test L.jumps_eigenvalue(bg, l) == ref
        end
    end
end
