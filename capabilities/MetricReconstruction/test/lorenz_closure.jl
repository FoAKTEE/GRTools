# Grounding for Lorenz.Closure — the App-H LOWER jump-closure system (gold-standard oracle).
#
# (1) PUBLISHED App-H LOWER (P_{-1},κ) jump tables (gates G-K2/G-K3): closure_kerr(0.0,6,2;lmax=16) vs
#     ref_lower_a0 rel<1e-4; closure_kerr(0.99,6,2;lmax=16) vs ref_lower_a099 rel<1e-3. Literals baked
#     verbatim from packages/KerrLorenzCirc/test/ref_apph_jumps.jl (the published App-H LOWER table).
# (2) a=0 REGRESSION (G-K4, internal identity): closure_kerr(0.0,…) == closure_a0(…) max-rel<1e-4.
# (3) PARALLEL gate: probe_parallel=true == serial bit-identical (BLAS.set_num_threads(1) guarded).

using Test
import MetricReconstruction.Lorenz as L

const ref_lower_a0 = Dict(
 2=>[-269.03958+109.83495im, 67.259895-36.611649im, -30.891079, 32.178207],
 3=>[62.91587-102.74119im, 31.457935-158.39266im, 0.0, 0.0],
 4=>[350.94586-52.09929im, -143.56876+62.953309im, 26.752459, -82.486749],
 5=>[-65.986722+107.75587im, -32.993361+408.57432im, 0.0, 0.0],
 6=>[-471.56933+35.003159im, 214.3497-90.424827im, -26.046465, 163.87568],
 7=>[66.881148-109.21646im, 33.440574-769.06589im, 0.0, 0.0],
 8=>[599.8828-26.475798im, -283.72835+118.03793im, 25.794803, -275.14456],
 9=>[-67.263179+109.84031im, -33.631589+1240.2802im, 0.0, 0.0],
 10=>[-731.16863+21.32131im, 352.52773-145.69562im, -25.675987, 416.16495])
const ref_lower_a099 = Dict(
 2=>[-337.06643+85.15083im, 66.553306-18.416882im, -26.1142, 36.026256],
 3=>[19.119854-82.676796im, 15.022065-101.32406im, 0.0, 0.0],
 4=>[524.94212-38.915946im, -126.27358+48.912959im, 22.620668, -107.28865],
 5=>[-20.504489+85.504877im, -15.347091+242.82248im, 0.0, 0.0],
 6=>[-738.33137+25.959912im, 184.5521-75.107992im, -22.02319, 220.87866],
 7=>[20.918971-86.31449im, 15.435462-446.38189im, 0.0, 0.0],
 8=>[956.254-19.584755im, -242.41322+100.36684im, 21.810085, -375.96945],
 9=>[-21.097679+86.65854im, -15.472351+712.39424im, 0.0, 0.0],
 10=>[-1175.8634+15.752413im, 300.09629-125.2532im, -21.709449, 572.47272])

maxrel(sol, ref) = maximum(
    abs(ref[l][k]) > 1e-6 ? abs(sol[l][k]-ref[l][k])/abs(ref[l][k]) : abs(sol[l][k])
    for l in 2:10 for k in 1:4)

@testset "Lorenz.Closure (App-H LOWER gold-standard oracle)" begin

    @testset "G-K2: closure_kerr(0.0,6,2) vs App-H LOWER a=0 (rel<1e-4)" begin
        sol = L.closure_kerr(0.0, 6.0, 2; lmax=16)
        @test maxrel(sol, ref_lower_a0) < 1e-4
    end

    @testset "G-K3: closure_kerr(0.99,6,2) vs App-H LOWER a=0.99 (rel<1e-3)" begin
        sol = L.closure_kerr(0.99, 6.0, 2; lmax=16)
        @test maxrel(sol, ref_lower_a099) < 1e-3
    end

    @testset "G-K4: a=0 regression closure_kerr(0.0)==closure_a0 (max-rel<1e-4)" begin
        solK = L.closure_kerr(0.0, 6.0, 2; lmax=16)
        solA = L.closure_a0(6.0, 2; ls=2:10)
        mx = maximum(
            abs(solA[l][k]) > 1e-6 ? abs(solK[l][k]-solA[l][k])/abs(solA[l][k]) : abs(solK[l][k])
            for l in 2:10 for k in 1:4)
        @test mx < 1e-4
    end

    @testset "PARALLEL: probe_parallel == serial (bit-identical)" begin
        ss = L.closure_kerr(0.99, 6.0, 2; lmax=12, probe_parallel=false)
        sp = L.closure_kerr(0.99, 6.0, 2; lmax=12, probe_parallel=true)
        for l in keys(ss), k in 1:4
            @test sp[l][k] == ss[l][k]   # bit-identical
        end
    end
end
