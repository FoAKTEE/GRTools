# RT-A2 velocity-reconstruction operator Q^{-1}/Q (Modes.Velocity).
# GROUNDING (self-contained, NON-CIRCULAR — no oracle value needed):
#   q_operator·qinv_operator = I (LinearAlgebra.I) to <1e-10 away from corotation. This is verbatim
#   the legacy DiscModes cov_rel_harmonics `q_resid < 1e-10` assertion (build_Q·build_Qinv = I).
#   Plus structural checks: finiteness of Q^{-1}, BigFloat element-type inheritance (Complex{BigFloat}).
using Test
using AccretionDisc
using AccretionDisc.Modes: RelativisticDisc, qinv_operator, q_operator
using LinearAlgebra

@testset "Modes.Velocity (RT-A2 Q-operator)" begin
    d = RelativisticDisc(a = 0.6, h0 = 0.1)

    # (1) Q·Q^{-1} = I away from corotation — the legacy q_resid<1e-10 identity, on the same points
    #     as verify_rel_Qoperator.jl: (r, r_p, m).
    mx = 0.0
    for (r, rp, m) in [(12.0, 15.0, 2), (20.0, 15.0, 2), (8.0, 10.0, 5), (25.0, 15.0, 3)]
        Qi = qinv_operator(d, r, rp, m)
        Q  = q_operator(d, r, rp, m)
        @test all(isfinite, Qi)                     # Q^{-1} entries finite
        @test eltype(Qi) == Complex{Float64}        # Complex{T} assembly
        mx = max(mx, maximum(abs.(Q * Qi - I)))
    end
    @test mx < 1e-10                                # build_Q·build_Qinv = I (legacy cov tolerance)

    # (2) BigFloat element-type inheritance: Complex{BigFloat} flows through the whole assembly.
    db = RelativisticDisc(a = big"0.6", h0 = big"0.1", r_I = big"1.0", eden_I = big"1.0")
    Qib = qinv_operator(db, big"12.0", big"15.0", 2)
    Qb  = q_operator(db, big"12.0", big"15.0", 2)
    @test eltype(Qib) == Complex{BigFloat}
    @test maximum(abs.(Qb * Qib - I)) < 1e-30       # BigFloat identity far tighter than Float64
end
