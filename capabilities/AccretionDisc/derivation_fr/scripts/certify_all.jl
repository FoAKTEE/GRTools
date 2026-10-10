using FRDerivation
println("Runtime: Julia ", VERSION, "; Julia threads=", Threads.nthreads())

groups = Dict("environment" => FRDerivation.certify_environment,
              "kepler" => FRDerivation.certify_fr22_kepler,
              "background" => FRDerivation.certify_fr22_background,
              "fr24-background" => FRDerivation.certify_fr24_background,
              "normalisation" => FRDerivation.certify_fr24_normalisation,
              "fourier" => FRDerivation.certify_fr22_fourier,
              "second-order-averages" => FRDerivation.certify_fr24_second_order_averages,
              "linearised" => FRDerivation.certify_fr22_linearised,
              "lindblad" => FRDerivation.certify_fr22_lindblad,
              "phi-factorisation" => FRDerivation.certify_fr22_phi_factorisation)
all(name -> haskey(groups, name), ARGS) || error("unknown certificate group")
order = ["environment", "kepler", "background", "normalisation", "fourier", "second-order-averages",
         "linearised", "fr24-background", "lindblad", "phi-factorisation"]
# Fourier reality uses the landed Kepler parity certificate even in a focused run.
requested = copy(ARGS)
any(in(requested), ["linearised","lindblad"]) && append!(requested,["background","fourier"])
"phi-factorisation" in requested && push!(requested,"fourier")
"fourier" in requested && push!(requested,"kepler")
selected = isempty(ARGS) ? order : filter(in(requested), order)
for name in selected
    certs = groups[name]()
    for c in certs
        println("CERTIFICATE\t",name,'\t',c.name,'\t',c.eqref,'\t',c.modality,'\t',c.pass)
    end
    println("Certificates: ", count(c -> c.pass, certs), "/", length(certs), " PASS (", name, ")")
    flush(stdout)
    all(c -> c.pass, certs) || exit(1)
end
