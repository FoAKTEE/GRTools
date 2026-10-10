# Run from GRTools:
# julia --startup-file=no --project=capabilities/AccretionDisc/derivation -t1 capabilities/AccretionDisc/derivation/scripts/certify_midpoint.jl
using LinearAlgebra
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "../../../..")))
isdefined(@__MODULE__, :R3Derivation) || include(joinpath(@__DIR__, "../src/R3Derivation.jl"))
include(joinpath(@__DIR__, "../src/midpoint.jl"))
using .MidpointCertificate

println("C-64 candidate evidence: ECC::gr-midpoint-certificate")
println("Julia=", VERSION, " Julia_threads=", Threads.nthreads(),
        " BLAS_threads=", BLAS.get_num_threads())
Threads.nthreads() <= 2 || error("Certificate run exceeds two Julia threads")
BLAS.get_num_threads() == 1 || error("Certificate run requires BLAS=1")
println(MIDPOINT_STATEMENT)
result = midpoint_certificates()
all(c -> c.pass, result.certificates) || error("Symbolic certificate failed")
println("Symbolic verdict: ", length(result.certificates), "/",
        length(result.certificates), " PASS; all residuals exactly zero")

for T in (Float64, BigFloat)
    setprecision(BigFloat, 128) do
        for variable_steps in (false, true)
            println("ORDER_CHECK type=", T, " bits=", precision(T),
                    " variable_steps=", variable_steps,
                    " linsolve=", T === Float64 ? "sparse" : "banded")
            rows = midpoint_order_check(T; variable_steps)
            for row in rows
                println(row)
            end
            all(abs(row.order-2) < T(1)/20 for row in rows[2:end]) ||
                error("Manufactured solution did not converge at second order")
            abs(rows[end].order-2) < T(1)/200 || error("Final order outside tolerance")
        end
    end
end
println("PASS: symbolic consistency, code/paper update identity, manufactured order checks")
