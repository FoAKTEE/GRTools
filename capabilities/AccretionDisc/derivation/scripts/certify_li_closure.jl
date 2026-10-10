# julia --project=capabilities/AccretionDisc/derivation -t2 \
#   capabilities/AccretionDisc/derivation/scripts/certify_li_closure.jl OUTDIR [NPTS]
include(joinpath(@__DIR__, "..", "src", "R3Derivation.jl"))
const R3 = R3Derivation
length(ARGS) in (1, 2) || error("usage: certify_li_closure.jl OUTDIR [NPTS]")
npts = length(ARGS) == 2 ? parse(Int, ARGS[2]) : 50
println("Option A: p=Theta(r)*eps; delta Theta=0; G26 singular-field/Boyer-Lindquist gauge")
println("Unchanged energy/Euler equations; conserved baryons; zero bath; retained relation (19)")
println("Delta p=Theta*Delta eps + eps*Theta_prime*xi^r; Theta=p/eps")
println("Julia ", VERSION, "; threads=", Threads.nthreads(), "; npts=", npts, "; seed=191101428")
println("Measure: existing nonuniform rational_point; no probability-of-identity claim")
println("Degree bound: <=1 in response/metric amplitudes; no global background numerator bound")
println("Exclusions: nonzero p, eps, eps+p, enu, epsi, emu, W, omega_d and all solve pivots")
certs = R3.LocallyIsothermal.closure_certificates(ARGS[1]; npts)
all(c.pass for c in certs) || error("LI closure certificate FAILED")
println("[PASS] LI closure: ", length(certs), "/", length(certs), " certificates")
