"Construct the existing runtime disc with a single inherited scalar type."
function certificate_disc(h, p, q)
    h, p, q = promote(h, p, q)
    T = typeof(h)
    Modes.NewtonianDisc{T}(h, q, p, one(h), zero(h), zero(h), 0, zero(h), 8)
end

"Small runtime smoke gate; physics certificates are separate entry points."
function certify_environment()
    @variables s h p q
    d = certificate_disc(h, p, q)
    c = Certificate[]
    push!(c, certify_symbolic("environment: Num through Sigma", "D-FR24-1",
                              Modes.Sigma(d, s) - 1/s^p))
    push!(c, certify_symbolic("environment: Num through LTinv", "D-FR24-1",
                              Modes.LTinv(d, s) + q/s))
    # Exercise the actual generated master kernel using Rational{BigInt}, not a
    # symbolic transcription. The LI algebra comparison belongs to the HOLE node.
    z, o = RatB(0), RatB(1)
    args = (2o, 7o/5, o/10, o, -3o/4, o, -o/2, z, z, z, z,
            3o, o/5, z, -8o, o, z, z, o/100, o, o/3, o/4, 2o)
    value = Modes.master_equation_cool(args...)
    push!(c, certify_bool("environment: exact generated-code arithmetic", "D-FR24-1",
        all(x -> x isa RatB, value), "runtime master_equation_cool returns only Rational{BigInt}"))
    c
end
