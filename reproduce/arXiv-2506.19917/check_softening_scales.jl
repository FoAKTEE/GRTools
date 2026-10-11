using AccretionDisc, Test, TOML
const F=AccretionDisc.ForcingModes
const M=AccretionDisc.Modes
# Transcribed from 2_linear_theory.tex:189 and 5_discussion.tex:40-51.
s=F.FixedSoftening(.3*.06)
xs=sqrt(2e-6/(3s.eps))
periods=4/(3xs)
@testset "FD25 constant softening and printed scales" begin
    @test isapprox(xs,6e-3;rtol=.02)
    @test isapprox(periods,220;rtol=.01)
    @test round(xs;sigdigits=1)==.006
    @test M._bank_softlength(s,.06,1.5,1.,.5,1.)==s.eps
    @test M._bank_softlength(s,.06,1.5,1.,1.5,1.)==s.eps
end
out=Dict("softening_length"=>s.eps,"xs_over_a"=>xs,"libration_orbital_periods"=>periods,
    "fr24_peri_to_fixed_e012"=>(1-.12)^1.5,"fr24_apo_to_fixed_e012"=>(1+.12)^1.5,
    "policy_owner"=>"AccretionDisc.ForcingModes.FixedSoftening")
open(ARGS[1],"w") do io;TOML.print(io,out);end
println(out)
