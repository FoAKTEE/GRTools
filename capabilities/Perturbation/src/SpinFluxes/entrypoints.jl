# Public TeukolskySpinMode* entrypoints (BHPToolkit names) — faithful 1:1 of the wrapper layer in
# TeukolskySpinFluxes/src/TeukolskySpinFluxes.jl (the 10 entrypoints missing from the original
# circular/eccentric SpinFluxes slice; TeukolskySpinModeCircularCorrectionAnalytical and
# TeukolskySpinModeEquatorialCorrectionAnalyticalPhases live in circular.jl/eccentric.jl).
# Together the 12 exported TeukolskySpinMode* symbols complete the design §4.3 SpinFluxes surface.

"""
    TeukolskySpinMode(l, m, n, k, orbit) -> NamedTuple
    TeukolskySpinMode(l, m, n, k, a, p, e, x, nmax, kmax, s) -> NamedTuple

Full finite-s linear-in-spin s=−2 Teukolsky mode from a materialized numeric spinning orbit — the
BHPToolkit `TeukolskySpinMode` (L1210-1404). `orbit` is the `KerrSpinOrbitNum` object (spin baked
into the trajectory + frequencies); the wr×wθ Mino double-sum evaluates the source at the
PERTURBED point rp=r0+s·r1 with the radial function Taylor-expanded about r0 and the full
A0/A1/A2/B1/B2/B3 source assembly (cf. `TeukolskySpinModeFromCorrection`, which strictly
linearises in s). The convenience method builds `orbit = KerrSpinOrbitNum(KerrSpinOrbitCorrection(
a,p,e,x,nmax,kmax), s)`. Returns `(; l,m,n,k,omega, α, Amplitudes=(I,H), S,
Fluxes=Energy/AngularMomentum→(I,H), stepsr, stepsθ)`.
"""
function TeukolskySpinMode(l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real, nmax::Integer, kmax::Integer, s::Real)
    orbit = KerrSpinOrbitNum(KerrSpinOrbitCorrection(a, p, e, x, nmax, kmax), s)
    return TeukolskySpinMode(l, m, n, k, orbit)
end

function TeukolskySpinMode(l::Integer, m::Integer, n::Integer, k::Integer, orbit::NamedTuple)
    md = _spin_mode_num(orbit, l, m, n, k)
    return (l = md.l, m = md.m, n = md.n, k = md.k, omega = md.omega, α = md.α,
            Amplitudes = (I = md.CPlus, H = md.CMinus), S = md.Spi2,
            Fluxes = (Energy = (I = md.EnergyI, H = md.EnergyH),
                      AngularMomentum = (I = md.AngMomI, H = md.AngMomH)),
            stepsr = md.stepsr, stepsθ = md.stepsθ)
end

"""
    TeukolskySpinModeFromCorrection(l, m, n, k, a, p, e, x, nmax, kmax, s) -> NamedTuple

Linear-in-spin s=−2 Teukolsky mode from a NUMERIC orbit correction — the BHPToolkit
`TeukolskySpinModeFromCorrection` (L1656-1872), the "DH-trajectory" flux path. The numeric
`KerrSpinOrbitCorrection(a,p,e,x,nmax,kmax)` is materialized internally, then the wr×wθ Mino
double-sum runs over the spinning orbit: the spin enters via the numeric trajectory correction
(rS,zS,UrS,UzS,ΔtS,ΔφS) and the spin-dipole source projection. Returns `(; l,m,n,k,omega, α,
Amplitudes=(I,H), S, Fluxes=Energy/AngularMomentum/CarterConstantK→(I,H), stepsr, stepsθ)`.
`TeukolskySpinModeFromCorrection(l, m, n, k, orbitCorrection, s)` is also accepted.
"""
function TeukolskySpinModeFromCorrection(l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real, nmax::Integer, kmax::Integer, s::Real)
    oc = KerrSpinOrbitCorrection(a, p, e, x, nmax, kmax)
    return TeukolskySpinModeFromCorrection(l, m, n, k, oc, s)
end

function TeukolskySpinModeFromCorrection(l::Integer, m::Integer, n::Integer, k::Integer,
        orbitCorrection::NamedTuple, s::Real)
    md = _from_correction_mode(orbitCorrection, l, m, n, k, s)
    return (l = md.l, m = md.m, n = md.n, k = md.k, omega = md.omega, α = md.α,
            Amplitudes = (I = md.CPlus, H = md.CMinus), S = md.Spi2,
            Fluxes = (Energy = (I = md.EnergyI, H = md.EnergyH),
                      AngularMomentum = (I = md.AngMomI, H = md.AngMomH),
                      CarterConstantK = (I = md.CarterKI, H = md.CarterKH)),
            stepsr = md.stepsr, stepsθ = md.stepsθ)
end

"""
    TeukolskySpinModeAnalytical(l, m, n, k, a, p, e, x, spar; gauge="CTilde", Nmax=80) -> NamedTuple

Generic (inclined) linear-in-spin s=−2 Teukolsky mode — the BHPToolkit `TeukolskySpinModeAnalytical`
(L1879-2013). Materializes the inclined `KerrSpinOrbit` internally, then performs the wr×wz Mino
double-sum → `(; l,m,n,k,omega, α, Amplitudes=(I,H),
Fluxes=Energy/AngularMomentum/CarterConstantK→(I,H), stepsr, stepsz)`. The Drasco–Hughes P5 path.
"""
function TeukolskySpinModeAnalytical(l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real, spar::Real; gauge::AbstractString="CTilde", Nmax::Integer=80)
    orbit = KerrSpinOrbit(a, p, e, x, spar; gauge=gauge)
    md = _inclined_spin_mode(orbit, l, m, n, k; Nmax=Nmax)
    return (l=md.l, m=md.m, n=md.n, k=md.k, omega=md.omega, α=md.α,
            Amplitudes = (I = md.CPlus, H = md.CMinus),
            Fluxes = (Energy = (I = md.EnergyI, H = md.EnergyH),
                      AngularMomentum = (I = md.AngMomI, H = md.AngMomH),
                      CarterConstantK = (I = md.CarterKI, H = md.CarterKH)),
            stepsr = md.stepsr, stepsz = md.stepsz)
end

"""
    TeukolskySpinModeCorrection(l, m, n, k, a, p, e, x, nmax, kmax) -> NamedTuple
    TeukolskySpinModeCorrection(l, m, n, k, orbitCorrection::NamedTuple) -> NamedTuple

Numeric-orbit-correction linear-in-spin Teukolsky mode + flux CORRECTION — the BHPToolkit
`TeukolskySpinModeCorrection` (L2609-2860). Returns BOTH the geodesic monopole (Amplitudes/Fluxes,
s⁰) AND its linear-in-spin correction (AmplitudesCorrection/FluxesCorrection = dE¹/dLz¹) for the
NUMERIC inclined orbit `KerrSpinOrbitCorrection(a,p,e,x,nmax,kmax)` (built internally by the first
method). The .wl `{angparNew,RCorrection}` ω-derivative callbacks are replaced by internal central
FD over ω. Returns `(; l,m,n,k,omega, Amplitudes=(I,H), AmplitudesCorrection=(I,H), α, S,
Fluxes=Energy/AngularMomentum→(I,H), FluxesCorrection=Energy/AngularMomentum→(I,H), stepsr,stepsθ)`.
"""
function TeukolskySpinModeCorrection(l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real, nmax::Integer, kmax::Integer)
    oc = KerrSpinOrbitCorrection(a, p, e, x, nmax, kmax)
    return TeukolskySpinModeCorrection(l, m, n, k, oc)
end

function TeukolskySpinModeCorrection(l::Integer, m::Integer, n::Integer, k::Integer,
        orbitCorrection::NamedTuple)
    md = _spin_mode_correction(orbitCorrection, l, m, n, k)
    return (l = md.l, m = md.m, n = md.n, k = md.k, omega = md.omega,
            Amplitudes = (I = md.CPlus0, H = md.CMinus0),
            AmplitudesCorrection = (I = md.CPlus1, H = md.CMinus1),
            α = md.α, S = md.Spi2,
            Fluxes = (Energy = (I = md.EnergyI, H = md.EnergyH),
                      AngularMomentum = (I = md.AngMomI, H = md.AngMomH)),
            FluxesCorrection = (Energy = (I = md.dE1I, H = md.dE1H),
                                AngularMomentum = (I = md.dL1I, H = md.dL1H)),
            stepsr = md.stepsr, stepsθ = md.stepsθ)
end

"""
    TeukolskySpinModeCorrectionAnalytical(l, m, n, k, a, p, e, x; Nmax=80, h=3e-4) -> NamedTuple

Generic (inclined) linear-in-spin flux CORRECTION — the BHPToolkit
`TeukolskySpinModeCorrectionAnalytical` (L2864-3090). The geodesic orbit + its O(s) correction via
the wr×wz Mino double-sum; the `angparNew`/`RCorrection` ω-derivative callbacks are materialized
internally (central FD). Returns `(; l,m,n,k,omega, Amplitudes=(I,H), AmplitudesCorrection=(I,H),
FluxesCorrection=Energy/AngularMomentum→(I,H), stepsr, stepsz)`. The Drasco–Hughes P5 path:
`dE_I(2,2,0,−10, 0.9,12,0.2,0.5)≈3.818e-15`.
"""
function TeukolskySpinModeCorrectionAnalytical(l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real; Nmax::Integer=80, h::Real=3e-4)
    sf = inclined_spin_flux(l, m, n, k, a, p, e, x; Nmax=Nmax, h=h)
    return (l=l, m=m, n=n, k=k, omega=sf.omega,
            Amplitudes = (I = sf.CPlus0, H = sf.CMinus0),
            AmplitudesCorrection = (I = sf.CPlus1, H = sf.CMinus1),   # linear-in-spin C^{+,1}/C^{-,1}
            FluxesCorrection = (Energy = (I = sf.dE1I, H = sf.dE1H),
                                AngularMomentum = (I = sf.dL1I, H = sf.dL1H)),
            stepsr = sf.stepsr, stepsz = sf.stepsz)
end

"""
    TeukolskySpinModeSpherical(l, m, k, a, p, e, x, kmax, s) -> NamedTuple
    TeukolskySpinModeSpherical(l, m, k, orbitCorrection::NamedTuple, s) -> NamedTuple

Spherical-orbit linear-in-spin s=−2 Teukolsky mode — the BHPToolkit `TeukolskySpinModeSpherical`
(L3331-3513). The θ-only analog of `TeukolskySpinModeFromCorrection`: a SINGLE wz (polar Mino) sum
at fixed r≡p (Ur=0), ω = m·Ωφ + k·Ωθ. The convenience method builds
`oc=KerrSpinOrbitCorrectionSpherical(a,p,e,x,kmax)` internally (e must be 0). Returns
`(; l,m,k,n=0,omega, Amplitudes=(I,H), α, S, Fluxes=Energy/AngularMomentum→(I,H), stepsθ)`.
"""
function TeukolskySpinModeSpherical(l::Integer, m::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real, kmax::Integer, s::Real)
    oc = KerrSpinOrbitCorrectionSpherical(a, p, e, x, kmax)
    return TeukolskySpinModeSpherical(l, m, k, oc, s)
end

function TeukolskySpinModeSpherical(l::Integer, m::Integer, k::Integer,
        orbitCorrection::NamedTuple, s::Real)
    md = _spin_mode_spherical(orbitCorrection, l, m, k, s)
    return (l = md.l, m = md.m, k = md.k, n = md.n, omega = md.omega,
            Amplitudes = (I = md.CPlus, H = md.CMinus), α = md.α, S = md.Spi2,
            Fluxes = (Energy = (I = md.EnergyI, H = md.EnergyH),
                      AngularMomentum = (I = md.AngMomI, H = md.AngMomH)),
            stepsθ = md.stepsθ)
end

"""
    TeukolskySpinModeSphericalCorrectionNum(l, m, k, a, p, e, x, kmax, δω) -> NamedTuple
    TeukolskySpinModeSphericalCorrectionNum(l, m, k, orbitCorrection::NamedTuple, δω) -> NamedTuple

Spherical-orbit linear-in-spin Teukolsky mode + flux CORRECTION (numerical ω-derivative) — the
BHPToolkit `TeukolskySpinModeSphericalCorrectionNum` (L3516-3754). Returns BOTH the geodesic
monopole AND its linear-in-spin correction (dE¹/dLz¹) for the SPHERICAL orbit
`KerrSpinOrbitCorrectionSpherical(a,p,e,x,kmax)` (built internally by the first method; e must be
0). The ω-derivatives of R/SWSH use EXPLICIT central FD with the user step `δω` plus the analytic
fab-basis ω-derivative. Returns `(; l,m,k,n=0,omega, Amplitudes=(I,H), AmplitudesCorrection=(I,H),
α, S, Fluxes=Energy/AngularMomentum→(I,H), FluxesCorrection=Energy/AngularMomentum→(I,H), stepsθ)`.
"""
function TeukolskySpinModeSphericalCorrectionNum(l::Integer, m::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real, kmax::Integer, δω::Real)
    oc = KerrSpinOrbitCorrectionSpherical(a, p, e, x, kmax)
    return TeukolskySpinModeSphericalCorrectionNum(l, m, k, oc, δω)
end

function TeukolskySpinModeSphericalCorrectionNum(l::Integer, m::Integer, k::Integer,
        orbitCorrection::NamedTuple, δω::Real)
    md = _spin_mode_spherical_correction_num(orbitCorrection, l, m, k, δω)
    return (l = md.l, m = md.m, k = md.k, n = md.n, omega = md.omega,
            Amplitudes = (I = md.CPlus0, H = md.CMinus0),
            AmplitudesCorrection = (I = md.CPlus1, H = md.CMinus1),
            α = md.α, S = md.Spi2,
            Fluxes = (Energy = (I = md.EnergyI, H = md.EnergyH),
                      AngularMomentum = (I = md.AngMomI, H = md.AngMomH)),
            FluxesCorrection = (Energy = (I = md.dE1I, H = md.dE1H),
                                AngularMomentum = (I = md.dL1I, H = md.dL1H)),
            stepsθ = md.stepsθ)
end

"""
    TeukolskySpinModeSphericalCorrection(l, m, k, a, p, e, x, kmax; Nmax=80, hh=1e-6, h=1e-5) -> NamedTuple
    TeukolskySpinModeSphericalCorrection(l, m, k, orbitCorrection, orbitDerivatives; Nmax=80, hh=1e-6) -> NamedTuple

Spherical-orbit linear-in-spin Teukolsky mode + flux CORRECTION + r/x DERIVATIVES — the BHPToolkit
`TeukolskySpinModeSphericalCorrection` (L3759-4110), the largest flux function. Returns the
geodesic monopole, its linear-in-spin correction (same physics as `…CorrectionNum`), AND the
∂r/∂x derivatives of the flux (the inspiral-evolution content). The .wl
`{angparNew,TeukolskySolverHS1spin}` callbacks are materialized internally (ω-derivatives via
central FD of SWSH/TeukolskyRadial at ω±hh + the analytic `_dfabidω`); `orbitDerivatives` is the
FD `_orbit_derivatives_spherical` object (built internally by the first method; e must be 0).
Returns `(; l,m,k,n=0,omega, Amplitudes=(I,H), AmplitudesCorrection=(I,H),
AmplitudesDerivatives=(r=(I,H), x=(I,H)), α, S, Fluxes=…→(I,H), FluxesCorrection=…→(I,H),
FluxesDerivatives=(r=…,x=…), stepsθ)`.
"""
function TeukolskySpinModeSphericalCorrection(l::Integer, m::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real, kmax::Integer;
        Nmax::Integer=80, hh::Real=1e-6, h::Real=1e-5)
    oc = KerrSpinOrbitCorrectionSpherical(a, p, e, x, kmax)
    odv = _orbit_derivatives_spherical(a, p, e, x, kmax; h=h)
    return TeukolskySpinModeSphericalCorrection(l, m, k, oc, odv; Nmax=Nmax, hh=hh)
end

function TeukolskySpinModeSphericalCorrection(l::Integer, m::Integer, k::Integer,
        orbitCorrection::NamedTuple, orbitDerivatives::NamedTuple; Nmax::Integer=80, hh::Real=1e-6)
    md = _spin_mode_spherical_correction(orbitCorrection, orbitDerivatives, l, m, k; Nmax=Nmax, hh=hh)
    return (l = md.l, m = md.m, k = md.k, n = md.n, omega = md.omega,
            Amplitudes = (I = md.CPlus0, H = md.CMinus0),
            AmplitudesCorrection = (I = md.CPlus1, H = md.CMinus1),
            AmplitudesDerivatives = (r = (I = md.dCPlusdr, H = md.dCMinusdr),
                                     x = (I = md.dCPlusdx, H = md.dCMinusdx)),
            α = md.α, S = md.Spi2,
            Fluxes = (Energy = (I = md.EnergyI, H = md.EnergyH),
                      AngularMomentum = (I = md.AngMomI, H = md.AngMomH)),
            FluxesCorrection = (Energy = (I = md.dE1I, H = md.dE1H),
                                AngularMomentum = (I = md.dL1I, H = md.dL1H)),
            FluxesDerivatives = (r = (Energy = (I = md.dEdrI, H = md.dEdrH),
                                      AngularMomentum = (I = md.dLdrI, H = md.dLdrH)),
                                 x = (Energy = (I = md.dEdxI, H = md.dEdxH),
                                      AngularMomentum = (I = md.dLdxI, H = md.dLdxH))),
            stepsθ = md.stepsθ)
end

"""
    TeukolskySpinModeSphericalCorrectionAnalytical(l, m, k, a, p, x; Nmax=80, h=3e-4) -> NamedTuple

Spherical-orbit (e=0) linear-in-spin Teukolsky mode + flux CORRECTION — the BHPToolkit
`TeukolskySpinModeSphericalCorrectionAnalytical` (L4472-4672). The θ-only restriction of
`TeukolskySpinModeCorrectionAnalytical`: a SINGLE wz (polar Mino) sum at fixed r≡p (Ur=0). Returns
the geodesic monopole (Amplitudes/Fluxes incl. CarterConstantK) AND its linear-in-spin correction
(dE¹/dLz¹/dK¹). The .wl `{angparNew,TeukolskySolverHS1spin}` analytic ω-derivative callbacks are
replaced by internal central FD over ω (`h`); the geodesic orbit is built internally. Returns
`(; l,m,k,n=0,omega, Amplitudes=(I,H), AmplitudesCorrection=(I,H), α, S,
Fluxes=Energy/AngularMomentum/CarterConstantK→(I,H), FluxesCorrection=…→(I,H), stepsθ)`.
"""
function TeukolskySpinModeSphericalCorrectionAnalytical(l::Integer, m::Integer, k::Integer,
        a::Real, p::Real, x::Real; Nmax::Integer=80, h::Real=3e-4)
    md = _spin_mode_spherical_correction_analytical(l, m, k, a, p, x; Nmax=Nmax, h=h)
    return (l = md.l, m = md.m, k = md.k, n = md.n, omega = md.omega,
            Amplitudes = (I = md.CPlus0, H = md.CMinus0),
            AmplitudesCorrection = (I = md.CPlus1, H = md.CMinus1),
            α = md.α, S = md.Spi2,
            Fluxes = (Energy = (I = md.EnergyI, H = md.EnergyH),
                      AngularMomentum = (I = md.AngMomI, H = md.AngMomH),
                      CarterConstantK = (I = md.CarterKI, H = md.CarterKH)),
            FluxesCorrection = (Energy = (I = md.dE1I, H = md.dE1H),
                                AngularMomentum = (I = md.dL1I, H = md.dL1H),
                                CarterConstantK = (I = md.dCarterKI, H = md.dCarterKH)),
            stepsθ = md.stepsz)
end

"""
    TeukolskySpinModeEquatorialCorrectionAnalyticalDarwin(l, m, n, a, p, e; Nmax=80, h=3e-4) -> NamedTuple

Eccentric-EQUATORIAL s=-2 spinning-body mode in the DARWIN χr parametrization (BHPToolkit name,
.wl L1004-1199). The Darwin-reparametrized sibling of
`TeukolskySpinModeEquatorialCorrectionAnalyticalPhases`: the SAME equatorial flux via a χr-grid sum
(Jacobian dλ/dχr, prefactor Ωr) instead of the Mino-phase grid. Returns `(; l,m,k=0,n,omega,
Amplitudes=(I,H), AmplitudesCorrection=(I,H), Fluxes=Energy/AngularMomentum→(I,H),
FluxesCorrection=Energy/AngularMomentum→(I,H), stepsr)`.
"""
function TeukolskySpinModeEquatorialCorrectionAnalyticalDarwin(l::Integer, m::Integer, n::Integer,
        a::Real, p::Real, e::Real; Nmax::Integer=80, h::Real=3e-4)
    sf = _spin_mode_eq_correction_darwin(l, m, n, a, p, e; Nmax=Nmax, h=h)
    λ = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * sf.omega)
    eg = energy_flux(-2, l, m, a, sf.omega, λ, sf.CPlus0, sf.CMinus0)
    jg = angular_momentum_flux(-2, l, m, a, sf.omega, λ, sf.CPlus0, sf.CMinus0)
    return (l = l, m = m, k = 0, n = n, omega = sf.omega,
            Amplitudes = (I = sf.CPlus0, H = sf.CMinus0),
            AmplitudesCorrection = (I = sf.CPlus1, H = sf.CMinus1),
            Fluxes = (Energy = (I = eg.I, H = eg.H), AngularMomentum = (I = jg.I, H = jg.H)),
            FluxesCorrection = (Energy = (I = sf.dE1I, H = sf.dE1H),
                                AngularMomentum = (I = sf.dL1I, H = sf.dL1H)),
            stepsr = sf.stepsr)
end
