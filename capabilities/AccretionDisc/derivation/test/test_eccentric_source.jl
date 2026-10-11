@testset "C-SOURCE symbolic local model" begin
    certs = R3.source_certificates()
    @test !isempty(certs)
    @test all(c.pass for c in certs)
end

@testset "C-SOURCE code-linked symbolic mutation controls" begin
    path = joinpath(@__DIR__,"..","src","eccentric_source.jl")
    original = read(path,String)
    for (tag,old,new,target) in (
        ("sy5","B = ulow[2]*ulow[3]","B = zero(ulow[2])","runtime distance equals chart"),
        ("sy6","K = 2mass*(g+2ulow*transpose(ulow))",
            "K = 2mass*(g+ulow*transpose(ulow))","factor four"),
        ("sy7","3j.s2_r^2/(4s^5)","j.s2_r^2/(4s^5)","analytic field derivatives"),
        ("zero-htr","h = K ./ s","h = K ./ s; h[1,2] = h[2,1] = zero(s)","h_tr contraction"),
        ("zero-hrphi","h = K ./ s","h = K ./ s; h[2,3] = h[3,2] = zero(s)","h_rphi contraction"))
        @test length(findall(old,original)) == 1
        # A private module contains each mutant. Production bytes and the
        # control module are never overwritten; exceptions are not kills.
        scope = Module(gensym(:SourceMutation))
        Base.include_string(scope,replace(original,old=>new),path)
        certs = Base.invokelatest(R3.source_certificates;
            source=Base.invokelatest(getproperty,scope,:EccentricSource),gate=target,verbose=false)
        @test !isempty(certs)
        @test all(c -> !c.pass && c.modality == :symbolic,certs)
        println("[PASS] MUTANT $tag killed symbolically by ",join((c.name for c in certs),"; "))
    end
    # The untouched control must still pass after every isolated mutant.
    @test all(c.pass for c in R3.source_certificates(;verbose=false))
end
