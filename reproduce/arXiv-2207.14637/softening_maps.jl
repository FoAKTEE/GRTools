# Render numerical arrays only; softening_figures.py provides the paper layout/read-back.
include("softening_bank.jl")
using DelimitedFiles
function export_maps(paper,rule,e,mmax,bankdir,out)
    s=setup(paper,rule,e,mmax)
    rows=[deserialize(joinpath(bankdir,"response-$(m.m)-$(m.m+m.k).jls")) for m in s.modes]
    all(y -> length(y.r)==4000,rows) || error("incomplete published radial grid")
    r=rows[1].r
    all(y -> y.r==r,rows) || error("common published radial grid required")
    w=M.FR22Wake(r,s.modes,hcat((y.sigma for y in rows)...))
    mkpath(out)
    radii=10.0 .^ range(-.5,.5,length=315)
    polar_r=collect(range(.1,2.,length=384))
    phi=collect(range(0,2pi,length=769))[1:end-1]
    line_phi=collect(range(-pi/2,pi/2,length=8193))
    fractions=(0.,.25,.5,.75)
    planets=[M.planet_orbit(2pi*t,s.o) for t in fractions]
    for (j,fraction) in enumerate(fractions)
        t=2pi*fraction
        rp=planets[j].r_p
        for (label,rr,frame) in (("unwrap",radii,:planet),("polar",polar_r,:guiding))
            data=M.fr22_wake(w,s.d,s.o,rr,phi,t;frame)
            open(joinpath(out,"$label-$j.bin"),"w") do io;write(io,data);end
        end
        line=vec(M.fr22_wake(w,s.d,s.o,[rp+.1],line_phi,t;frame=:planet))
        writedlm(joinpath(out,"line-$j.csv"),hcat(line_phi,line),',')
    end
    writedlm(joinpath(out,"modes.csv"),hcat([m.m for m in s.modes],[m.m+m.k for m in s.modes],
        [y.torque for y in rows],[y.cpu_s for y in rows]),',')
    open(joinpath(out,"maps.toml"),"w") do io
        TOML.print(io,Dict("e"=>e,"rule"=>rule,"mmax"=>mmax,"offset"=>s.K,
            "nphi"=>length(phi),"nr"=>length(radii),"npolar_r"=>length(polar_r),
            "planet_radii"=>[p.r_p for p in planets],
            "guiding_azimuths"=>[p.phi_p-2pi*t for (p,t) in zip(planets,fractions)],
            "normalization"=>"h_p^3 deltaSigma(Gmu=1) / Sigma","layout"=>"Float64 column-major little-endian"))
    end
    println("MAP ARRAYS $rule e=$e; $(length(rows)) modes")
end
export_maps(ARGS[1],ARGS[2],parse(Float64,ARGS[3]),parse(Int,ARGS[4]),ARGS[5],ARGS[6])
