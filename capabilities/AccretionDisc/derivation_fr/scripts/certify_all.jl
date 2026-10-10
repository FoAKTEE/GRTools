using FRDerivation

groups = Dict("environment" => FRDerivation.certify_environment,
              "kepler" => FRDerivation.certify_fr22_kepler,
              "background" => FRDerivation.certify_fr22_background,
              "normalisation" => FRDerivation.certify_fr24_normalisation,
              "fourier" => FRDerivation.certify_fr22_fourier,
              "second-order-averages" => FRDerivation.certify_fr24_second_order_averages)
all(name -> haskey(groups, name), ARGS) || error("unknown certificate group")
order = ["environment", "kepler", "background", "normalisation", "fourier", "second-order-averages"]
# Fourier reality uses the landed Kepler parity certificate even in a focused run.
requested = "fourier" in ARGS ? vcat(ARGS,["kepler"]) : ARGS
selected = isempty(ARGS) ? order : filter(in(requested), order)
for name in selected
    certs = groups[name]()
    for c in certs
        println("CERTIFICATE\t",name,'\t',c.name,'\t',c.eqref,'\t',c.modality,'\t',c.pass)
    end
    println("Certificates: ", count(c -> c.pass, certs), "/", length(certs), " PASS (", name, ")")
    all(c -> c.pass, certs) || exit(1)
end
