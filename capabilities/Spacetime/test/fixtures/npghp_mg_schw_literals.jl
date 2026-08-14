# npghp_mg_schw_literals.jl — baked ModifiedGravity.schw_spincoeff literals (G-B3).
#
# PROVENANCE: the equatorial-Schwarzschild GHP spin-coefficient reduction rules quoted verbatim in
# R_surface.md §2 Site A (from ModifiedGravity/src/ghp_operators.jl:32-37, `schw_spincoeff`):
#   ρ  = −1/r
#   ρ′ = (r − 2M)/(2 r²)
#   ε′ = −M/(2 r²)          (ε′ = ε′†)
#   ψ₂ = −M/r³
#   κ = σ = τ = τ′ = ε = ε† = 0
# with M = 1 (Spacetime M=1 canon). BAKED as literals rather than a live ModifiedGravity dep:
# adding MG to Spacetime's test env would INVERT the Tier-2 dependency ordering and break the §7.7
# hygiene gate. These are the a=0 / equator reduction target for Spacetime.NPGHP.np_spincoeff and
# its primed view (ρ′ = −μ, ε′ = −γ; the [AXIOM-4] GHP prime relations).
module NPGHPMGSchwLiterals

# key r -> (ρ, ρ′, ε′, ψ₂) evaluated from the closed forms above (M = 1); all real.
const SCHW = Dict{Float64,NamedTuple}(
    3.0  => (ρ = -1 / 3,  ρp = (3 - 2) / (2 * 3^2),  εp = -1 / (2 * 3^2),  ψ2 = -1 / 3^3),
    5.0  => (ρ = -1 / 5,  ρp = (5 - 2) / (2 * 5^2),  εp = -1 / (2 * 5^2),  ψ2 = -1 / 5^3),
    10.0 => (ρ = -1 / 10, ρp = (10 - 2) / (2 * 10^2), εp = -1 / (2 * 10^2), ψ2 = -1 / 10^3),
)

end # module NPGHPMGSchwLiterals
