# Data JLD2 weakdep-ext activation (design §5.4): Hybrid declares JLD2 in [deps] as the ext
# trigger and loads it from src/Hybrid.jl, so `using Hybrid` (which also loads Data) must leave
# Data's DataJLD2Ext live — the declared .jld2 routing contract is real, not nominal.

import Data

@testset "Data JLD2 weakdep-ext trigger (§5.4)" begin
    @test Base.get_extension(Data, :DataJLD2Ext) !== nothing
end
