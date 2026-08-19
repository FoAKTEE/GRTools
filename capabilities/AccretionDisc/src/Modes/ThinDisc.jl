# Relativistic Novikov–Thorne / Shakura–Sunyaev thin-disc STRUCTURE and the total Lindblad-torque
# assembly of arXiv:2510.02433 (RT4 fig:discProfile, RT5 fig:TorqueStrengths, tab:FitParameters,
# RT6 fig:TorqueGWComp). PROMOTED verbatim-in-physics from the legacy figure-local library
# `results/emri-disk/figures/2510.02433/lindblad_torque_lib.jl` + `repro_discProfile.jl`, which had
# no package home. Everything kinematic is DELEGATED — the Kerr circular frequencies to
# `Background` (disc_OmegaK/kappar2), the resonance radii to `Resonances.lindblad_locations`, the
# per-mode strength to `Forcing.lindblad_strength`, the ISCO to `Spacetime.kerr_isco` — so no
# geodesic formula is written twice.
#
# Structure equations: Novikov & Thorne (1973) / Page & Thorne (1974) relativistic correction
# factors 𝒜,ℬ,𝒞,𝒟,ℰ,𝒬 in the Abramowicz & Fragile (2013, Living Rev.) notation, radiation-pressure
# / electron-scattering INNER region (their Eqs. 99), y = √(r/M):
#   𝒜 = 1 + a²y⁻⁴ + 2a²y⁻⁶     ℬ = 1 + a y⁻³            𝒞 = 1 − 3y⁻² + 2a y⁻³
#   𝒟 = 1 − 2y⁻² + a²y⁻⁴       ℰ = 1 + 4a²y⁻⁴ − 4a²y⁻⁶ + 3a⁴y⁻⁸
#   𝒬 = 𝒬₀[ y − y₀ − (3/2)a ln(y/y₀) − Σᵢ 3(yᵢ−a)²/(yᵢ(yᵢ−yⱼ)(yᵢ−yₖ)) ln((y−yᵢ)/(y₀−yᵢ)) ],
#   𝒬₀ = ℬ/(y√𝒞), y₀ = √(r_ISCO), yᵢ = roots of y³ − 3y + 2a = 0 (Page–Thorne closed form).
#   Σ = Σ₀ (α ṁ)⁻¹ r^{3/2} 𝒜⁻²ℬ³𝒞^{1/2}ℰ𝒬⁻¹  [g/cm²]      H = H₀ ṁ 𝒜²ℬ⁻³𝒞^{1/2}𝒟⁻¹ℰ⁻¹𝒬  [cm]
# Newtonian limit 𝒜,ℬ,𝒞,𝒟,ℰ→1, 𝒬→y: Σ ∝ r^{3/2} and h = H/r ∝ 1/r — analytically the
# Shakura–Sunyaev scalings, which is why the SS overlay is normalized to coincide with the a=0 NT
# profile at `r_norm` ("normalize all profiles to the same value in the Newtonian limit", l.573).
# The aspect ratio carries ONE free overall constant (the paper: the normalization "alters the
# overall normalisation … but not their radial dependence"), fixed by `h_anchor` at `r_norm`; the
# gravitational radius cancels identically out of that ratio, so no physical constants are needed.
# Type-generic {T}; BigFloat-inherited.

using ..Modes: RelativisticDisc, lindblad_locations, lindblad_strength, kerr_isco

"""
    NovikovThorneDisc(; a, mdot, alpha, h_anchor, r_norm, sigma0, height0)

Relativistic Novikov–Thorne thin disc around a Kerr primary (M = 1, geometrized radii). `a` = spin,
`mdot` = Ṁ/Ṁ_Edd, `alpha` = Shakura–Sunyaev viscosity; `h_anchor` = the published aspect ratio at
`r_norm` that fixes the single free normalization of h = H/r (arXiv:2510.02433 fig:discProfile:
h(a=0, 1000M) ≈ 0.0033); `sigma0` [g/cm²] and `height0` [cm] are the Abramowicz–Fragile inner-region
prefactors. Only `a` changes the SHAPE; `mdot`/`alpha`/`sigma0` set the overall Σ scale.
"""
Base.@kwdef struct NovikovThorneDisc{T<:Real}
    a::T        = 0.0
    mdot::T     = 0.1
    alpha::T    = 0.1
    h_anchor::T = 0.0033
    r_norm::T   = 1000.0
    sigma0::T   = 5.0
    height0::T  = 1.0e5
end

"Kerr equatorial prograde ISCO of the disc — delegates to `Spacetime.kerr_isco`."
nt_isco(d::NovikovThorneDisc) = kerr_isco(d.a)

"The same disc at zero spin — the reference curve every normalization (h anchor, SS overlay) uses."
_zero_spin(d::NovikovThorneDisc) = NovikovThorneDisc(a = zero(d.a), mdot = d.mdot, alpha = d.alpha,
                                                    h_anchor = d.h_anchor, r_norm = d.r_norm,
                                                    sigma0 = d.sigma0, height0 = d.height0)

"""
    nt_factors(a, r) -> (; A, B, C, D, E, Q)

Novikov–Thorne / Page–Thorne relativistic correction factors at radius `r` (M = 1) for spin `a`;
see the file header for the closed forms. `Q` is the Page–Thorne flux factor (the same 𝒬 that
carries the disc-flux profile F ∝ ℬ⁻¹𝒞^{-1/2}𝒬), evaluated from the ISCO `y₀ = √r_ISCO`.
"""
function nt_factors(a::Real, r::Real)
    T = float(promote_type(typeof(a), typeof(r)))
    aa = T(a); y = sqrt(T(r))
    A = 1 + aa^2*y^-4 + 2*aa^2*y^-6
    B = 1 + aa*y^-3
    C = 1 - 3*y^-2 + 2*aa*y^-3
    D = 1 - 2*y^-2 + aa^2*y^-4
    E = 1 + 4*aa^2*y^-4 - 4*aa^2*y^-6 + 3*aa^4*y^-8
    y0 = sqrt(kerr_isco(aa))
    Q0 = B / (y * sqrt(C))
    ac = acos(aa)
    y1 =  2*cos((ac - T(π))/3)
    y2 =  2*cos((ac + T(π))/3)
    y3 = -2*cos(ac/3)
    term(yi, yj, yk) = 3*(yi - aa)^2 / (yi*(yi - yj)*(yi - yk)) * log((y - yi)/(y0 - yi))
    Q = Q0 * (y - y0 - (3//2)*aa*log(y/y0) - term(y1, y2, y3) - term(y2, y1, y3) - term(y3, y1, y2))
    return (A = A, B = B, C = C, D = D, E = E, Q = Q)
end

"""
    nt_surface_density(d::NovikovThorneDisc, r) -> Σ(r)  [g/cm²]

Novikov–Thorne inner-region surface density Σ = Σ₀(α ṁ)⁻¹ r^{3/2} 𝒜⁻²ℬ³𝒞^{1/2}ℰ𝒬⁻¹
(Abramowicz & Fragile Eqs. 99). Absolute cgs value — the paper's Σ normalization needs no rescale.
"""
function nt_surface_density(d::NovikovThorneDisc, r)
    f = nt_factors(d.a, r)
    d.sigma0 * (d.alpha^-1 * d.mdot^-1) * r^(3//2) * f.A^-2 * f.B^3 * sqrt(f.C) * f.E / f.Q
end

"""
    nt_scale_height(d::NovikovThorneDisc, r) -> H(r)  [cm]

Novikov–Thorne inner-region half-thickness H = H₀ ṁ 𝒜²ℬ⁻³𝒞^{1/2}𝒟⁻¹ℰ⁻¹𝒬 (radiation-pressure
supported, hence M- and r-independent in cgs up to the correction factors).
"""
function nt_scale_height(d::NovikovThorneDisc, r)
    f = nt_factors(d.a, r)
    d.height0 * d.mdot * f.A^2 * f.B^-3 * sqrt(f.C) / f.D / f.E * f.Q
end

"""
    nt_aspect_ratio(d::NovikovThorneDisc, r) -> h(r) = H/r

Dimensionless aspect ratio with the EXACT relativistic radial+spin shape and the single free overall
constant fixed by `d.h_anchor` at `d.r_norm` on the a=0 curve:

    h(r) = h_anchor · [H_NT(a, r)/r] / [H_NT(0, r_norm)/r_norm].

The gravitational radius cancels identically, so no physical constants enter.
"""
function nt_aspect_ratio(d::NovikovThorneDisc, r)
    ref = _zero_spin(d)
    d.h_anchor * (nt_scale_height(d, r)/r) / (nt_scale_height(ref, d.r_norm)/d.r_norm)
end

"""
    ss_surface_density(d::NovikovThorneDisc, r) -> Σ_SS(r)  [g/cm²]

Shakura–Sunyaev (Newtonian) overlay Σ_SS = Σ_NT(a=0, r_norm)(r/r_norm)^{3/2}: the same power law the
NT profile approaches as r → ∞, normalized to coincide with the a=0 NT curve at `d.r_norm`.
"""
function ss_surface_density(d::NovikovThorneDisc, r)
    ref = _zero_spin(d)
    nt_surface_density(ref, d.r_norm) * (r/d.r_norm)^(3//2)
end

"""
    ss_aspect_ratio(d::NovikovThorneDisc, r) -> h_SS(r)

Shakura–Sunyaev aspect ratio h_SS = h_NT(a=0, r_norm)·r_norm/r (∝ 1/r), the Newtonian limit of
[`nt_aspect_ratio`](@ref) with the same anchor.
"""
function ss_aspect_ratio(d::NovikovThorneDisc, r)
    ref = _zero_spin(d)
    nt_aspect_ratio(ref, d.r_norm) * d.r_norm / r
end

"""
    pressure_cutoff(ξ) -> f(ξ)

Pressure regularization of the resonant torque (arXiv:2510.02433 eq:TPressureCutoff),

    f(ξ) = [2/(1+e^ξ)] · [1/(1+4ξ²)],   ξ = m h(r_LR),

which suppresses modes whose azimuthal wavelength falls below the disc scale height (m ≳ 1/h) and
so makes the mode sum converge. f(0)=1, f decays faster than any power.
"""
pressure_cutoff(ξ::T) where {T<:Real} = (2/(1 + exp(ξ))) * (1/(1 + 4*ξ^2))

"""
    total_lindblad_torque(d::NovikovThorneDisc, r_s; profile=:nt, mmax=nothing, mcap=600,
                          enhancement=nothing) -> T_LR

SIGNED total Lindblad torque on a circular equatorial secondary at `r_s` (per q², in the Σ cgs
scale), summed over azimuthal modes and over the inner/outer resonances (arXiv:2510.02433 RT5):

    T_LR(r_s) = Σ_m Σ_{kind∈{ILR,OLR}} 2π r_LR Σ(r_LR) N^{(m)}(r_s) f(m h(r_LR)),

with the signed per-mode strength `N` from [`lindblad_strength`](@ref) (ILR < 0 ⇒ outward migration,
OLR > 0 ⇒ inward) and the pressure cutoff [`pressure_cutoff`](@ref). `profile` selects the
Novikov–Thorne (`:nt`) or the Shakura–Sunyaev overlay (`:ss`, Newtonian Σ,h with the a=0 relativistic
kinematics). The mode ceiling defaults to `min(mcap, max(60, ⌈8/h(r_s)⌉))` — the cutoff kills the
tail well below it. `enhancement(m, r_s, kind)` optionally multiplies each mode (the relativistic
amplitude enhancement |N_Teukolsky/N_Newtonian|, cf. [`lindblad_torque`](@ref)).
"""
function total_lindblad_torque(d::NovikovThorneDisc{T}, r_s; profile::Symbol = :nt,
                               mmax::Union{Nothing,Integer} = nothing, mcap::Integer = 600,
                               enhancement = nothing) where {T}
    profile in (:nt, :ss) || throw(ArgumentError("profile must be :nt or :ss, got $(profile)"))
    Σof = profile === :nt ? (r -> nt_surface_density(d, r)) : (r -> ss_surface_density(d, r))
    hof = profile === :nt ? (r -> nt_aspect_ratio(d, r))    : (r -> ss_aspect_ratio(d, r))
    kin = RelativisticDisc(a = profile === :nt ? d.a : zero(d.a))
    if mmax === nothing
        hmin = max(hof(r_s), T(1e-3))
        mmax = min(mcap, max(60, ceil(Int, 8/hmin)))
    end
    Tnet = zero(float(T))
    for m in 1:mmax
        rI, rO = lindblad_locations(kin, r_s, m)
        for (r1, kind) in ((rI, :ILR), (rO, :OLR))
            isnan(r1) && continue
            N = lindblad_strength(kin, r_s, m; kind = kind)
            isfinite(N) || continue
            e = enhancement === nothing ? one(N) : enhancement(m, r_s, kind)
            Tnet += 2*T(π)*r1*Σof(r1)*N*e*pressure_cutoff(m*hof(r1))
        end
    end
    return Tnet
end

"""
    powerlaw_fit(x, y; x0=10) -> (; A, n)

Least-squares single power law `y = A (x/x0)^n` on POSITIVE data, solved linearly in log–log
(arXiv:2510.02433 tab:FitParameters, the "All r" rows with x0 = 10M).
"""
function powerlaw_fit(x::AbstractVector, y::AbstractVector; x0 = 10)
    length(x) == length(y) || throw(DimensionMismatch("powerlaw_fit: x and y must have equal length"))
    length(x) >= 2 || throw(ArgumentError("powerlaw_fit needs at least two points"))
    lx = log.(x ./ x0); ly = log.(abs.(y))
    n = length(lx); sx = sum(lx); sy = sum(ly); sxx = sum(lx .^ 2); sxy = sum(lx .* ly)
    den = n*sxx - sx^2
    abs(den) > 0 || throw(ArgumentError("powerlaw_fit: degenerate abscissae"))
    slope = (n*sxy - sx*sy)/den
    (A = exp((sy - slope*sx)/n), n = slope)
end
