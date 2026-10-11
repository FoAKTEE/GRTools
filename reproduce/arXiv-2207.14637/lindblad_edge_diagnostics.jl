# Postprocess the saved edge sweep, with the same arguments as the audit:
# regular-edge-diagnostics <prior-evidence-directory> <output-directory>.
# Tests far-zone AMF plateaus and attributes their remaining variation by the
# directly integrated torque source, without rerunning the mode solver.
include("lindblad_crossing_audit.jl")

d=disc((;h=.1,q=0.,p=0.)); b=deserialize(joinpath(OUT,"edge-bank.jls"))
rows=[]
for (ei,e) in enumerate((.1,.25,.5)), m in (6,7,8,9), rout in (5.,6.,8.,10.,12.)
    mode=FM.NewtonianMode(m,1,1.); c=MO.harmonic_cache(b,ei,mode)
    z=deserialize(joinpath(OUT,"edge-$ei-$m-$rout.jls"))
    x=z.x; f=z.s.mode; flux=z.v.flux
    common=findall(r->r<=5.,x)
    row=Dict{String,Any}("m"=>m,"l"=>1,"e"=>e,"rout"=>rout,
        "torque_common_01_5"=>regular_observable(d,x[common],mode,c,f[common,:]).torque)
    windows=[]
    for (lo,hi) in ((8.,10.),(10.,12.))
        hi>rout && continue
        inds=findall(r->lo<=r<=hi,x)
        fs=flux[inds]; v=regular_observable(d,x[inds],mode,c,f[inds,:])
        delta=last(fs)-first(fs); scale=maximum(abs,fs)
        push!(windows,Dict("window"=>[lo,hi],"actual_window"=>[x[first(inds)],x[last(inds)]],
            "relative_flux_spread"=>(maximum(fs)-minimum(fs))/scale,
            "flux_change"=>delta,"torque_in_window"=>v.torque,
            "relative_balance_residual"=>abs(delta-v.torque)/scale))
    end
    row["windows"]=windows; push!(rows,row)
end
save("edge-diagnostics",Dict("rows"=>rows))
@assert length(rows)==60
for r in rows
    r["rout"]==12. || continue
    println("EDGE FLUX e=$(r["e"]) m=$(r["m"]) ",r["windows"])
end
println("EDGE DIAGNOSTICS COMPLETE: 60 profiles")
