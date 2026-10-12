# Continuous local balance probe, independent of a numerical BVP grid.
using AccretionDisc, Spacetime, LinearAlgebra, TOML, Printf
using AccretionDisc.ForcingModes
const SW=AccretionDisc.SpiralWaves
d=SW.DiscModel(;χ=0.0,r_ref=30.0,h_ref=0.06,pSigma=0.5,qTemp=0.0,closure=Adiabatic,evolution=:fixed)
p=SW.EccentricPerturber(EquatorialOrbit(0.0,30.0,0.3);mass=1e-5,softening=FixedSoftening(0.9))
c=SW.ModeContext(d,p,2,-4;closure=Adiabatic,domain=(9.0,150.0))
rows=Dict{String,Any}[]
for r in (15.0,30.0,60.0), δ in (0.02,0.01,0.005)
    radii=[r-δ,r,r+δ]
    bank=SW.source_azimuth_bank(d,p,radii;Nq=512,mmax=2)
    source=SW.mode_source(bank,c.mode)
    for kind in ("full","no_radial_metric","htt_only","zero")
        function src(x)
            jet=source(x)
            kind=="zero" && return SW._zero_source(c)
            kind=="htt_only" && return merge(SW._zero_source(c),(;htt=jet.htt,htt_r=jet.htt_r))
            kind=="no_radial_metric" && return merge(jet,(htr=zero(jet.htr),hrphi=zero(jet.hrphi)))
            jet
        end
        h=src(r); y=ComplexF64[1e-6+3e-7im,-2e-5+1e-5im]
        C,S=SW.master_matrices(c,r;hmodes=h,landau=false)
        yp=C*y+S
        fp=SW.amf_mode(c,r+δ,y+δ*yp;hmodes=src(r+δ),landau=false)
        fm=SW.amf_mode(c,r-δ,y-δ*yp;hmodes=src(r-δ),landau=false)
        lhs=(fp-fm)/(2δ)
        rhs=SW.torque_density_mode(c,r,y;hmodes=h,landau=false)
        row=Dict{String,Any}("r"=>r,"delta"=>δ,"source"=>kind,"lhs"=>lhs,"rhs"=>rhs,
            "relative"=>abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-30))
        push!(rows,row);println("LOCAL_BALANCE ",row);flush(stdout)
    end
end
open(io->TOML.print(io,Dict("rows"=>rows)),ARGS[1],"w")
