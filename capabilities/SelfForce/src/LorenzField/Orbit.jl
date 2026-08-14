# Circular-orbit parameters + frequency-domain BLS source coefficients R1..R10 that the h1-corpus
# generator OWNS. Schwarzschild geometry (lapse/tortoise) and the circular-orbit constants are
# CONSUMED from Spacetime; the normalized Legendre factors from NumericsCore. M = 1 throughout.

# Schwarzschild geometry helpers — consumed from Spacetime, not re-owned.
const _f = Spacetime.schwarzschild_lapse          # f(r) = 1 - 2M/r
const _df = Spacetime.schwarzschild_dlapse         # f'(r) = 2M/r²
const _Vl = Spacetime.bls_potential                # Vₗ(r) = ¼ f (f'/r + l(l+1)/r²)
const rs_of_r = Spacetime.tortoise_schw            # r*(r)
const r_of_rs = Spacetime.tortoise_schw_inverse    # r(r*)

# Circular-orbit constants of motion — re-exported from Spacetime (libSch.c circular sector).
const orbit_energy = Spacetime.schw_orbit_energy
const orbit_ang_mom = Spacetime.schw_orbit_ang_mom
const orbit_phi_freq = Spacetime.schw_orbit_phi_freq

"""
    OrbitalParams{T}

Circular Schwarzschild orbit at radius `p` (= `r0`). Holds the constants of motion `E`, `L`, `Ω_φ`
and the radial grid plus the (0-based, as in the `.h5` attrs) important indices. Mirrors the
populated fields of the H1Lorenz `orbital_params` struct used by the circular path. All scalars are
type `T`; `e = Ω_r = 0` for the circular case.
"""
struct OrbitalParams{T}
    p::T
    e::T
    E::T
    L::T
    Omega_phi::T
    Omega_r::T
    grid::Vector{T}
    ra_index::Int
    r0_index::Int
    rb_index::Int
end

"""
    OrbitalParams(p, grid, ra, r0, rb) -> OrbitalParams{T}

Construct circular orbital params at radius `p` with the supplied radial grid and the (0-based)
ra/r0/rb important indices. `T = typeof(float(p))` (BigFloat inherited).
"""
function OrbitalParams(p, grid::AbstractVector, ra::Integer, r0::Integer, rb::Integer)
    T = typeof(float(p))
    pp = T(p)
    OrbitalParams{T}(pp, zero(T), orbit_energy(pp), orbit_ang_mom(pp), orbit_phi_freq(pp),
                     zero(T), T.(grid), Int(ra), Int(r0), Int(rb))
end

# Frequency-domain BLS source coefficients (port of ref-code/h1Lorenz/src/sources.c). Circular orbit
# (e=0): the scaling-coefficient integral collapses to a single evaluation at χ=π/2, where r_p=p,
# t_p=φ_p=0 ⇒ source_cos=1, source_sin=0, contra_ur=√(E²-f(p)(L²/p²+1))=0. Hence R5_source and
# R9_source (∝ source_sin / contra_ur) vanish identically. Sources return SOURCE_TOTAL at χ=π/2.
# `_Plm`/`_Plm_θ` are the normalized associated Legendre value / θ-derivative at the equator.
_Plm(l, m, ::Type{T}) where {T} = assoc_legendre_P(l, m, T(π) / 2; normalized = true)
_Plm_θ(l, m, ::Type{T}) where {T} = legendre_P_thetaderiv(l, m, T(π) / 2; normalized = true)

"R1 source (SOURCE_TOTAL at χ=π/2). Uses source_cos=1."
function R1_source(l, m, orbit::OrbitalParams{T}) where {T}
    p = orbit.p; fp = _f(p); E = orbit.E; L = orbit.L
    coeff = -4 * 8 * T(π) * _Plm(l, m, T)
    Complex{T}(coeff * fp^2 * (2 * E^2 * p^2 - fp * p^2 - L^2 * fp) / (E * p^3))
end

"R3 source (SOURCE_TOTAL at χ=π/2). Uses source_cos=1."
function R3_source(l, m, orbit::OrbitalParams{T}) where {T}
    p = orbit.p; fp = _f(p); E = orbit.E; L = orbit.L
    coeff = -4 * 8 * T(π) * _Plm(l, m, T)
    Complex{T}(coeff * fp^2 * (p^2 + L^2) / (E * p^3))
end

"R5 source (SOURCE_TOTAL at χ=π/2). ∝ source_sin = 0 ⇒ vanishes for circular orbits."
R5_source(l, m, orbit::OrbitalParams{T}) where {T} = zero(Complex{T})

"R6 source (SOURCE_TOTAL at χ=π/2). Uses source_cos=1."
function R6_source(l, m, orbit::OrbitalParams{T}) where {T}
    p = orbit.p; fp = _f(p); E = orbit.E; L = orbit.L
    coeff = -4 * 8 * T(π) * L^2 / E * _Plm(l, m, T)
    Complex{T}(coeff * fp^2 / p^3)
end

"R7 source (SOURCE_TOTAL at χ=π/2) = (l(l+1) - 2m²)·R6_source. Uses source_cos=1."
R7_source(l, m, orbit::OrbitalParams{T}) where {T} =
    (l * (l + 1) - 2 * m^2) * R6_source(l, m, orbit)

"R9 source (SOURCE_TOTAL at χ=π/2). ∝ contra_ur · source_sin = 0 ⇒ vanishes for circular orbits."
R9_source(l, m, orbit::OrbitalParams{T}) where {T} = zero(Complex{T})

"R10 source (SOURCE_TOTAL at χ=π/2). Carries explicit I; uses source_cos=1."
function R10_source(l, m, orbit::OrbitalParams{T}) where {T}
    p = orbit.p; fp = _f(p); E = orbit.E; L = orbit.L
    coeff = -4 * 16 * T(π) * m * L^2 / E * _Plm_θ(l, m, T)
    im * Complex{T}(coeff * fp^2 / p^3)
end
