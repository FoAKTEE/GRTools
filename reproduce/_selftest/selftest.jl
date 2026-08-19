# Harness self-test: a toy "result" exercising every golden value type and the
# check/regen flow.  THIN driver: zero physics here — the numbers come from
# GrTools.Spacetime (Bardeen ISCO, Kerr fundamental frequencies).
#
#     julia --project=reproduce reproduce/_selftest/selftest.jl --regen
#     julia --project=reproduce reproduce/_selftest/selftest.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()

using GrTools
using GrTools.Spacetime: kerr_isco, kerr_frequencies

# scalars
@repro_result "isco_a0.9" kerr_isco(0.9)
@repro_result "isco_a0" kerr_isco(0.0)                       # exact 6
@repro_result "isco_a0_is_six" kerr_isco(0.0) == 6.0        # Bool
@repro_result "n_spins" 6                                    # Int
@repro_result "label" "bardeen-isco"                         # String

# BigFloat inheritance (256 bits): the closed form evaluated in extended precision
isco_big = setprecision(BigFloat, 256) do
    kerr_isco(big(9) / 10)
end
@repro_result "isco_a0.9_big256" isco_big
@repro_result "isco_big_minus_f64" Float64(isco_big - kerr_isco(0.9))   # tiny, atol-gated below

# small arrays
spins = [-0.99, -0.5, 0.0, 0.5, 0.9, 0.99]
@repro_result "spins" spins
@repro_result "isco_vs_spin" kerr_isco.(spins)
freq = kerr_frequencies(0.9, 8.0, 0.3, 0.8)
@repro_result "freq_matrix" [freq.Omegar freq.Omegath freq.Omegaph; 2freq.Omegar 2freq.Omegath 2freq.Omegaph]
@repro_result "complex_mode" complex(freq.Omegaph, -freq.Omegar)
@repro_result "nonfinite_vec" [1.0, Inf, NaN]

# serializer round trip is part of the self-test: encode → text → parse → decode must be identical
let sample = Dict{String,Any}(
        "f" => 1.5e-7, "i" => -3, "b" => true, "s" => "a\"b\\c\né", "n" => nothing,
        "big" => big"1.234567890123456789012345678901234567890", "z" => 2.0 - 3.5im,
        "v" => [1.0, 2.0, 3.0], "m" => [1 2; 3 4], "nf" => [NaN, -Inf, 0.0],
        "nested" => Any[Any[1, 2], Any[3.0]])
    enc = Dict{String,Any}(k => ReproHarness.encode_value(v) for (k, v) in sample)
    back = Dict{String,Any}(k => ReproHarness.decode_value(v) for (k, v) in ReproHarness.json_read(ReproHarness.json_write(enc)))
    ok = all(k -> isequal(back[k], sample[k]) || (k == "nested"), keys(sample)) &&
         back["nested"] == [[1, 2], [3.0]] && back["big"] isa BigFloat && back["m"] isa Matrix{Int}
    @repro_result "json_roundtrip_ok" ok
end

exit(ReproHarness.conclude!(opts; paper = "_selftest", id = "harness-selftest", script = @__FILE__,
                            rtol = 1e-12, tolerances = Dict("isco_big_minus_f64" => (rtol = 0.0, atol = 1e-14))))
