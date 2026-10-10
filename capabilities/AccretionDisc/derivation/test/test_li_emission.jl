@testset "OPTION A LI generated master and stress jets" begin
    here = normpath(joinpath(@__DIR__,".."))
    target = joinpath(here,"..","src","SpiralWaves","generated")
    paths = [joinpath(target,"eccentric_$(tag).jl") for tag in ("li","li_stress")]
    mktempdir() do out
        for (tag,ex) in (("li",R3.li_response_exports()),("li_stress",R3.li_stress_exports()))
            certs = R3.emission_certificates(joinpath(target,"eccentric_$(tag).jl"),
                joinpath(here,"out","ECC_$(tag).jls.txt"),ex;tag,npts=NPTS,accounting_dir=out)
            @test all(c.pass for c in certs)
        end
        certs = R3.li_emitted_physics(paths...;npts=NPTS,accounting_dir=out)
        @test all(c.pass for c in certs)
        # A Theta'-blind generated preamble must fail at nonzero q, while
        # retaining the q=0 overlap. This is the closure review's R3 control.
        mutant = joinpath(out,"blind-gradient.jl")
        write(mutant,replace(read(paths[1],String),"p_r = inp.p_r"=>"p_r = inp.p*inp.eps_r/inp.eps"))
        bad = R3.li_emitted_physics(mutant,paths[2];npts=1,accounting_dir=out,verbose=false)
        @test only(c.pass for c in bad if occursin("fluid solution q=0//1",c.name))
        @test all(!c.pass for c in bad if occursin("fluid solution",c.name) && !occursin("q=0//1",c.name))
    end
end
