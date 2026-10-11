# Attribute adaptive-step sensitivity without replacing the shooting solver.
using AccretionDisc, TOML
const P=AccretionDisc.Modes
const F=AccretionDisc.ForcingModes
mutable struct TraceRHS{G,D,M}
    rhs::G
    d::D
    mode::M
    calls::Int
    min_D::Float64
    near_calls::Int
    max_rhs::Float64
    min_D_radius::Float64
    max_near_rhs::Float64
end
function (g::TraceRHS)(out,r,y)
    g.rhs(out,r,y)
    D=abs(P.Dfreq(g.d,r,g.mode))
    g.calls+=1
    if D<g.min_D
        g.min_D=D;g.min_D_radius=r
    end
    g.near_calls+=D<1e-8
    D<1e-8 && (g.max_near_rhs=max(g.max_near_rhs,maximum(abs,out)))
    g.max_rhs=max(g.max_rhs,maximum(abs,out))
    out
end
function main(out)
    d=NewtonianDisc(h0=.06,q=0.,p=0.,gamma=1.,beta=0.,soft=.018,lap_nphi=800)
    mode=F.NewtonianMode(40,40,1.);x=loggrid(d,.05,5.,4000)
    c=PotentialCache(d,40,.05,5.;n=1000)
    rows=Dict{String,Any}[]
    for rt in (1e-12,nextfloat(1e-12)), ih in (0,1)
        rhs=TraceRHS(P.CachedMasterRHS(d,mode,ih,c),d,mode,0,Inf,0,0.,NaN,0.)
        bc=ih==0 ? (1.0+0im,0.0+0im) : (0.0+0im,0.0+0im)
        f=P._integrate_master(d,x,bc,mode,rhs;rtol=rt,atol=rt/100,dense=false)
        push!(rows,Dict("rtol"=>rt,"ih"=>ih,"calls"=>rhs.calls,"minimum_abs_D"=>rhs.min_D,
            "calls_abs_D_below_1em8"=>rhs.near_calls,"maximum_rhs"=>rhs.max_rhs,"minimum_D_radius"=>rhs.min_D_radius,"maximum_rhs_near_Dzero"=>rhs.max_near_rhs,
            "inner_h_real"=>real(f[1,1]),"inner_h_imag"=>imag(f[1,1]),
            "outer_h_real"=>real(f[end,1]),"outer_h_imag"=>imag(f[end,1])))
        open(out,"w") do io;TOML.print(io,Dict("rows"=>rows));end
        println(last(rows));flush(stdout)
    end
end
main(ARGS[1])
