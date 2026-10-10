using FRDerivation

groups = Dict("environment" => FRDerivation.certify_environment,
              "kepler" => FRDerivation.certify_fr22_kepler,
              "background" => FRDerivation.certify_fr22_background,
              "normalisation" => FRDerivation.certify_fr24_normalisation)
all(name -> haskey(groups, name), ARGS) || error("unknown certificate group")
order = ["environment", "kepler", "background", "normalisation"]
selected = isempty(ARGS) ? order : filter(in(ARGS), order)
for name in selected
    certs = groups[name]()
    println("Certificates: ", count(c -> c.pass, certs), "/", length(certs), " PASS (", name, ")")
    all(c -> c.pass, certs) || exit(1)
end
