@testset "General complex adiabatic emission" begin
    here = normpath(joinpath(@__DIR__,".."))
    art = joinpath(here,"out","ECC_adiabatic.jls.txt")
    target = joinpath(here,"..","src","SpiralWaves","generated","eccentric_adiabatic.jl")
    mktempdir() do out
        certs = R3.emission_certificates(target,art,R3.general_response_exports();
            tag="adiabatic",npts=NPTS,accounting_dir=out)
        @test all(c.pass for c in certs)
        # Preamble corruption must be caught by executing the compiled target.
        mutant = joinpath(out,"mutant.jl")
        write(mutant,replace(read(target,String),"htr = inp.htr"=>"htr = zero(inp.htr)"))
        bad = R3.emission_certificates(mutant,art,R3.general_response_exports();
            tag="adiabatic",npts=1,accounting_dir=out,verbose=false)
        @test any(c -> !c.pass && occursin("compiled",c.name),bad)
    end
end
