"FR24 v1 (40),(51), and the thermal mass in section 2.2 line 79."
function certify_fr24_normalisation()
    @variables Sigma a Omega h Mp Ms G
    flux = Sigma*a^4*Omega^2/h^3*(Mp/Ms)^2
    rate = flux/(Mp*Omega*a^2)
    c = Certificate[]
    push!(c, certify_symbolic("FR24 (40) F_J0 with Kepler law", "40",
        Symbolics.substitute(flux,Dict(Omega^2=>G*Ms/a^3))-Sigma*G*a*Mp^2/(h^3*Ms)))
    push!(c, certify_symbolic("FR24 (51) inverse rate unit", "51",
                              rate-Sigma*a^2*Omega*Mp/(h^3*Ms^2)))
    push!(c, certify_thermal_mass())
    units = Dict(a=>1,Omega=>1,Ms=>1,Sigma=>1)
    push!(c, certify_symbolic("FR24 GRTools flux/rate/mass units", "40,51; section 2.2",
        [Symbolics.substitute(flux,units)-Mp^2/h^3,
         Symbolics.substitute(rate,units)-Mp/h^3,
         Symbolics.substitute(h^3*Ms,units)-h^3]))
    # These identities fix conversion factors for runtime quantities as well.
    @variables s p q
    d = certificate_disc(h,p,q)
    push!(c, certify_symbolic("FR24 runtime profile unit map", "2–4; a=Omega_p=Ms=Sigmap=1",
        [Modes.Sigma(d,s)-1/s^p, Modes.csiso(d,s)-h/s^(q/2)]))
    append!(c, certify_fr24_dimensions())
    c
end
