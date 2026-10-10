using Test
using FRDerivation
# Reuse the NumericsCore Gauss rule already imported by the runtime disc grid.
using AccretionDisc.SpiralWaves: gl_nodes_weights

"Test-only orbit integral: runtime Kepler orbit and runtime softened Laplace kernel."
function factorised_phi(r::T, orbit, m, l, epsilon, Gmu; N=512) where {T<:AbstractFloat}
    step = 2T(pi)/N
    value = zero(T)
    for i in 0:N-1
        M = (i+one(T)/2)*step
        p = FRDerivation.Modes.planet_orbit(M,orbit)
        b = FRDerivation.Modes.laplace_b_soft(one(T)/2,m,r/p.r_p,epsilon/p.r_p;nphi=N)
        value += cos(m*p.phi_p-l*M)*b/p.r_p
    end
    -Gmu*value/N
end

"""
Independent double integral: composite Gauss-Legendre in eccentric anomaly E
and absolute azimuth phi, Cartesian separation, and dt=(1-e*cos(E))/n dE.
No Kepler solve, true-anomaly function, Laplace kernel, or azimuth factorisation
is used. The physical-time normalization n/(2pi^2) is retained explicitly.
"""
function double_phi(r::T, orbit, m, l, epsilon, Gmu; panels=16, order=48) where {T<:AbstractFloat}
    nodes, weights = gl_nodes_weights(T,order)
    width = 2T(pi)/panels
    angles = [(j+(x+1)/2)*width for j in 0:panels-1 for x in nodes]
    ws = [w*width/2 for _ in 1:panels for w in weights]
    cs, ss = cos.(angles), sin.(angles)
    cm, sm = cos.(m.*angles), sin.(m.*angles)
    a,e,n = orbit.a,orbit.e,orbit.n
    total = zero(T)
    for i in eachindex(angles)
        E = angles[i]
        M = E-e*ss[i]
        xp,yp = a*(cs[i]-e),a*sqrt(1-e^2)*ss[i]
        time_weight = ws[i]*(1-e*cs[i])/n
        cl,sl = cos(l*M),sin(l*M)
        inner = zero(T)
        for j in eachindex(angles)
            distance = sqrt((r*cs[j]-xp)^2+(r*ss[j]-yp)^2+epsilon^2)
            inner += ws[j]*(cm[j]*cl+sm[j]*sl)/distance
        end
        total += time_weight*inner
    end
    -Gmu*n/(2T(pi)^2)*total
end

"Fixed precision, cases, resolutions and absolute tolerances for the numerical evidence."
function check_phi_factorisation()
    cases = ((2,2,0//1,4//5), (2,3,0//1,6//5),
             (3,2,1//5,7//10), (2,3,1//2,6//5))
    println("NUMERICAL CONFIG a=7/5, n=7/3, GMp=3/7, epsilon=a/5; ",
        "factorised midpoint (NM,Nphi)=(512,512),(1024,1024); independent double ",
        "Gauss-Legendre order=48, panels=(16,16),(32,32) in (E,phi); ",
        "absolute tolerances Float64=2e-12, BigFloat-192=1e-28 (including refinement).")
    println("QUADRATURE\ttype\tbits\tm\tl\te\tr/a\tfactorised\tdouble\tabs_difference\tfactor_refinement\tdouble_refinement\ttolerance")
    @testset "FR22 independent potential quadrature Float64/BigFloat" begin
        setprecision(192) do
            for T in (Float64,BigFloat), (m,l,e0,r0) in cases
                a,n,Gmu = T(7)/5,T(7)/3,T(3)/7
                epsilon = a/5
                r,e = a*T(r0),T(e0)
                orbit = FRDerivation.Modes.EccentricOrbit(a,e,n)
                coarse_f = factorised_phi(r,orbit,m,l,epsilon,Gmu;N=512)
                fine_f = factorised_phi(r,orbit,m,l,epsilon,Gmu;N=1024)
                coarse_d = double_phi(r,orbit,m,l,epsilon,Gmu;panels=16)
                fine_d = double_phi(r,orbit,m,l,epsilon,Gmu;panels=32)
                tolerance = T === Float64 ? T(2)/10^12 : T(10)^(-28)
                difference = abs(fine_f-fine_d)
                refine_f,refine_d = abs(fine_f-coarse_f),abs(fine_d-coarse_d)
                @test fine_f isa T && fine_d isa T
                @test precision(fine_f) == precision(T)
                @test precision(fine_d) == precision(T)
                @test difference <= tolerance
                @test refine_f <= tolerance
                @test refine_d <= tolerance
                if iszero(e)
                    z,o = zero(T),one(T)
                    d = FRDerivation.Modes.NewtonianDisc{T}(o/10,z,z,o,z,epsilon/a,0,z,1024)
                    circular = m == l ? Gmu/a*FRDerivation.Modes.Phi_m(d,r/a,m) : z
                    @test abs(fine_f-circular) <= tolerance
                end
                if abs(fine_f) > 100tolerance
                    # Omitting the corrected frequency gives amplitude/n.
                    @test abs(fine_d/n-fine_f) > 100tolerance
                end
                println(join(("QUADRATURE",T,precision(T),m,l,e,r/a,fine_f,fine_d,
                    difference,refine_f,refine_d,tolerance),'\t'))
            end
        end
    end
end
