"""
    MetricReconstruction.GHZ

Sourced Green--Hollands--Zimmerman reconstruction.  This first bounded capability solves the
ingoing-radiation-gauge corrector head on outgoing Schwarzschild rays,

    -(d²/dr² + 2/r d/dr) x_mmb = S_ll,

for an already-normalized Einstein-equation source.  It is deliberately separate from vacuum
Hertz/CCK reconstruction (`Radiation`) and the Li--Yunes static transport specializations
(`Transport`).
"""
module GHZ

import QuadGK: quadgk

abstract type SourceNormalization end
"`S_ab` in `delta G_ab[h] = S_ab`; no additional factor of `8pi` is applied."
struct EinsteinRHS <: SourceNormalization end
"Physical stress energy `T_ab`; rejected until explicitly mapped to an Einstein RHS."
struct PhysicalStressEnergy <: SourceNormalization end

"Kinnersley tetrad in the project's `(-,+,+,+)` convention."
struct KinnersleyMinusPlus end
"Outgoing Schwarzschild coordinates, with the principal-null ray directed toward larger `r`."
struct SchwarzschildOutgoing end
"Separated mode phase `exp(-i*omega*u + i*m*phi)`."
struct ExpMinusIomegaUPlusImPhi end
"The sourced GHZ ingoing-radiation-gauge construction."
struct GHZIRG end
"Vacuum Hertz ingoing-radiation-gauge reconstruction (a distinct construction marker)."
struct HertzIRG end

"""
    TetradSourceMode{N}(Sll; Slm=0, Sln=0, M, l, m, omega, source_fingerprint,
                        normalization=EinsteinRHS(), tetrad=KinnersleyMinusPlus(),
                        slicing=SchwarzschildOutgoing(),
                        phase=ExpMinusIomegaUPlusImPhi())

A separated mode of the three l-contracted source components `S_ll(r)`, `S_lm(r)`, and
`S_ln(r)`.  The positive perturbative order `N`, tetrad, slicing, normalization,
harmonic/frequency labels, background mass, phase convention, and content fingerprint are
carried with the callable sources.  `source_fingerprint` is a stable `UInt128` supplied by the
producer and covers all three components.
"""
struct ZeroRadialSource end
(::ZeroRadialSource)(r) = zero(r)

struct TetradSourceMode{N,T<:AbstractFloat,Norm<:SourceNormalization,Fll,Flm,Fln}
    Sll::Fll
    Slm::Flm
    Sln::Fln
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    source_fingerprint::UInt128
    normalization::Norm
    tetrad::KinnersleyMinusPlus
    slicing::SchwarzschildOutgoing
    phase::ExpMinusIomegaUPlusImPhi
end

_finite(z::Real) = isfinite(z)
_finite(z::Complex) = isfinite(real(z)) && isfinite(imag(z))

function _require_nonbool_real(value, label)
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("$label must be a non-Boolean real scalar"))
    value
end

function _require_nonbool_number(value, label)
    value isa Number || throw(ArgumentError("$label must be a numeric scalar"))
    components = value isa Real ? (value,) : (real(value), imag(value))
    any(component -> component isa Bool, components) &&
        throw(ArgumentError("$label must have non-Boolean scalar components"))
    value
end

function _checked_fingerprint(value, label)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$label must be a non-Boolean integer"))
    fp = try
        UInt128(value)
    catch
        throw(ArgumentError("$label must be representable as UInt128"))
    end
    fp != 0 || throw(ArgumentError("$label must be nonzero"))
    fp
end

function _checked_mode_index(value::Integer, label)
    value isa Bool &&
        throw(ArgumentError("$label must be a non-Boolean integer representable as Int"))
    try
        Int(value)
    catch
        throw(ArgumentError("$label must be a non-Boolean integer representable as Int"))
    end
end

function _checked_real(::Type{T}, value, label) where {T<:AbstractFloat}
    _require_nonbool_real(value, label)
    T(value)
end

function _checked_complex(::Type{T}, value, label) where {T<:AbstractFloat}
    _require_nonbool_number(value, label)
    Complex{T}(value)
end

function TetradSourceMode{N}(Sll; Slm=ZeroRadialSource(), Sln=ZeroRadialSource(), M,
        l::Integer, m::Integer, omega=zero(M),
        source_fingerprint, normalization::SourceNormalization=EinsteinRHS(),
        tetrad=KinnersleyMinusPlus(), slicing=SchwarzschildOutgoing(),
        phase=ExpMinusIomegaUPlusImPhi()) where {N}
    N isa Integer && !(N isa Bool) && N > 0 ||
        throw(ArgumentError("perturbative order N must be a positive integer"))
    tetrad isa KinnersleyMinusPlus ||
        throw(ArgumentError("GHZ head requires the Kinnersley (-,+,+,+) tetrad"))
    slicing isa SchwarzschildOutgoing ||
        throw(ArgumentError("GHZ head requires outgoing Schwarzschild slicing"))
    phase isa ExpMinusIomegaUPlusImPhi || throw(ArgumentError(
        "GHZ source phase must be exp(-i*omega*u + i*m*phi)"))
    lI = _checked_mode_index(l, "l")
    mI = _checked_mode_index(m, "m")
    lI >= 0 || throw(ArgumentError("l must be nonnegative"))
    -lI <= mI <= lI || throw(ArgumentError("mode requires |m| <= l"))
    _require_nonbool_number(M, "M")
    _require_nonbool_number(omega, "omega")
    T = promote_type(typeof(float(M)), typeof(float(real(omega))), typeof(float(imag(omega))))
    T <: AbstractFloat || throw(ArgumentError("M and omega must promote to a floating type"))
    MT, omegaT = T(M), Complex{T}(omega)
    _finite(MT) && MT > zero(T) || throw(ArgumentError("M must be finite and positive"))
    _finite(omegaT) || throw(ArgumentError("omega must be finite"))
    fp = _checked_fingerprint(source_fingerprint, "source_fingerprint")
    TetradSourceMode{N,T,typeof(normalization),typeof(Sll),typeof(Slm),typeof(Sln)}(
        Sll, Slm, Sln, MT, lI, mI, omegaT, fp, normalization, tetrad, slicing, phase)
end

"Return the compile-time perturbative order carried by a GHZ source, boundary, or solution."
perturbative_order(::TetradSourceMode{N}) where {N} = N

"""
    InnerRayIVP(source; rb, x_mmb=0, dr_x_mmb=0)

Inner-ray data at `rb >= 2M`.  All mode labels and the source fingerprint are copied into the
boundary carrier so that accidental cross-mode solves are rejected.
"""
struct InnerRayIVP{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    source_fingerprint::UInt128
    rb::T
    x0::Complex{T}
    dx0::Complex{T}
end

function InnerRayIVP(src::TetradSourceMode{N,T}; rb, x_mmb=zero(T),
        dr_x_mmb=zero(T)) where {N,T}
    rbT = _checked_real(T, rb, "inner-ray boundary rb")
    x0 = _checked_complex(T, x_mmb, "x_mmb boundary value")
    dx0 = _checked_complex(T, dr_x_mmb, "dr_x_mmb boundary value")
    _finite(rbT) && rbT >= 2src.M ||
        throw(ArgumentError("inner-ray boundary rb must be finite and satisfy rb >= 2M"))
    _finite(x0) && _finite(dx0) || throw(ArgumentError("boundary data must be finite"))
    InnerRayIVP{N,T}(src.M, src.l, src.m, src.omega, src.source_fingerprint, rbT, x0, dx0)
end
perturbative_order(::InnerRayIVP{N}) where {N} = N

"Keyed inner-ray data for the `x_nm` middle corrector equation."
struct MiddleRayIVP{N,T<:AbstractFloat}
    M::T
    l::Int
    m::Int
    omega::Complex{T}
    source_fingerprint::UInt128
    rb::T
    x0::Complex{T}
    dx0::Complex{T}
end

function _require_middle_spin_domain(l::Integer)
    lI = _checked_mode_index(l, "middle mode l")
    lI >= 1 || throw(ArgumentError("GHZ middle corrector requires l >= 1"))
    lI
end

function MiddleRayIVP(src::TetradSourceMode{N,T}; rb, x_nm=zero(T),
        dr_x_nm=zero(T)) where {N,T}
    _require_middle_spin_domain(src.l)
    rbT = _checked_real(T, rb, "middle boundary rb")
    x0 = _checked_complex(T, x_nm, "x_nm boundary value")
    dx0 = _checked_complex(T, dr_x_nm, "dr_x_nm boundary value")
    _finite(rbT) && rbT >= 2src.M ||
        throw(ArgumentError("middle boundary rb must be finite and satisfy rb >= 2M"))
    _finite(x0) && _finite(dx0) || throw(ArgumentError("middle boundary data must be finite"))
    MiddleRayIVP{N,T}(src.M, src.l, src.m, src.omega, src.source_fingerprint,
                       rbT, x0, dx0)
end
perturbative_order(::MiddleRayIVP{N}) where {N} = N

"A solved GHZ IRG corrector-head mode, evaluated by its one-integral Green formula."
struct IRGCorrectorHeadMode{N,T<:AbstractFloat,Fll,Flm,Fln}
    source::TetradSourceMode{N,T,EinsteinRHS,Fll,Flm,Fln}
    boundary::InnerRayIVP{N,T}
    gauge::GHZIRG
end
perturbative_order(::IRGCorrectorHeadMode{N}) where {N} = N

function _same_mode(src::TetradSourceMode, bc::Union{InnerRayIVP,MiddleRayIVP})
    src.M == bc.M && src.l == bc.l && src.m == bc.m && src.omega == bc.omega &&
        src.source_fingerprint == bc.source_fingerprint
end

"""
    solve_irg_corrector_head(source, boundary) -> IRGCorrectorHeadMode

Validate the complete mode key and return the analytic inner-ray solution.  A
`PhysicalStressEnergy` source is rejected: callers must explicitly construct the corresponding
`EinsteinRHS`, fixing trace reversal and the `8pi` normalization before reconstruction.
"""
function solve_irg_corrector_head(src::TetradSourceMode{N,T,EinsteinRHS,Fll,Flm,Fln},
        bc::InnerRayIVP{N,T}) where {N,T<:AbstractFloat,Fll,Flm,Fln}
    _same_mode(src, bc) || throw(ArgumentError("source and inner-ray boundary mode keys differ"))
    value = _checked_complex(T, src.Sll(bc.rb), "S_ll")
    _finite(value) || throw(ArgumentError("S_ll must be finite at the inner boundary"))
    IRGCorrectorHeadMode{N,T,Fll,Flm,Fln}(src, bc, GHZIRG())
end
function solve_irg_corrector_head(src::TetradSourceMode{N,T,PhysicalStressEnergy},
        bc::InnerRayIVP) where {N,T<:AbstractFloat}
    throw(ArgumentError("PhysicalStressEnergy is not an Einstein RHS; explicitly map T_ab to S_ab"))
end
function solve_irg_corrector_head(src::TetradSourceMode, bc::InnerRayIVP)
    src.normalization isa EinsteinRHS ||
        throw(ArgumentError("unsupported source normalization; GHZ requires EinsteinRHS"))
    throw(ArgumentError("source and inner-ray boundary order, precision, or mode keys differ"))
end

function _radius(sol::IRGCorrectorHeadMode{N,T}, r) where {N,T}
    rT = _checked_real(T, r, "evaluation radius")
    _finite(rT) && rT >= sol.boundary.rb && rT >= 2sol.source.M ||
        throw(ArgumentError("evaluation radius must be finite, outside the horizon, and r >= rb"))
    rT
end

function _Sll(sol::IRGCorrectorHeadMode{N,T}, r::T) where {N,T}
    value = _checked_complex(T, sol.source.Sll(r), "S_ll")
    _finite(value) || throw(ArgumentError("S_ll must be finite along the integration ray"))
    value
end

function _tolerances(::Type{T}, atol, rtol) where {T<:AbstractFloat}
    at = _checked_real(T, atol, "atol")
    rt = _checked_real(T, rtol, "rtol")
    _finite(at) && _finite(rt) && at >= zero(T) && rt >= zero(T) &&
        (at > zero(T) || rt > zero(T)) || throw(ArgumentError(
            "atol and rtol must be finite, nonnegative, and not both zero"))
    at, rt
end

function _quad(f, a::T, b::T; atol=zero(T), rtol=sqrt(eps(T))) where {T}
    at, rt = _tolerances(T, atol, rtol)
    a == b && return zero(Complex{T})
    first(quadgk(f, a, b; atol=at, rtol=rt))
end

"Evaluate `x_mmb(r)` using the single-integral Green solution of the corrector-head equation."
function x_mmb(sol::IRGCorrectorHeadMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    rT = _radius(sol, r)
    bc = sol.boundary
    at, rt = _tolerances(T, atol, rtol)
    rT == bc.rb && return bc.x0
    ratio = bc.rb / rT
    homogeneous = bc.x0 + bc.rb * bc.dx0 * (one(T) - ratio)
    sourced = _quad(q -> q * (one(T) - q / rT) * _Sll(sol, q), bc.rb, rT;
                     atol=at, rtol=rt)
    homogeneous - sourced
end

"Evaluate `d x_mmb/dr` from the differentiated single-integral Green solution."
function dr_x_mmb(sol::IRGCorrectorHeadMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    rT = _radius(sol, r)
    bc = sol.boundary
    at, rt = _tolerances(T, atol, rtol)
    rT == bc.rb && return bc.dx0
    ratio = bc.rb / rT
    moment = _quad(q -> (q / rT)^2 * _Sll(sol, q), bc.rb, rT;
                   atol=at, rtol=rt)
    ratio^2 * bc.dx0 - moment
end

"Pure residual `-(x'' + 2x'/r) - S_ll` for independently supplied local values."
function irg_head_residual(r, dx, d2x, Sll)
    _require_nonbool_real(r, "head residual radius")
    _require_nonbool_number(dx, "dr_x_mmb")
    _require_nonbool_number(d2x, "d2r_x_mmb")
    _require_nonbool_number(Sll, "S_ll")
    _finite(r) && r > zero(r) ||
        throw(ArgumentError("head residual radius must be finite and positive"))
    all(_finite, (dx, d2x, Sll)) ||
        throw(ArgumentError("head residual inputs must be finite"))
    -(d2x + 2dx / r) - Sll
end

"A solved outgoing-Schwarzschild IRG `x_nm` corrector mode."
struct IRGCorrectorMiddleMode{N,T<:AbstractFloat,Fll,Flm,Fln}
    head::IRGCorrectorHeadMode{N,T,Fll,Flm,Fln}
    boundary::MiddleRayIVP{N,T}
    gauge::GHZIRG
end
perturbative_order(::IRGCorrectorMiddleMode{N}) where {N} = N

function MiddleRayIVP(head::IRGCorrectorHeadMode; rb=head.boundary.rb, kwargs...)
    MiddleRayIVP(head.source; rb=rb, kwargs...)
end

function solve_irg_corrector_middle(head::IRGCorrectorHeadMode{N,T,Fll,Flm,Fln},
        bc::MiddleRayIVP{N,T}) where {N,T<:AbstractFloat,Fll,Flm,Fln}
    _require_middle_spin_domain(head.source.l)
    _require_middle_spin_domain(bc.l)
    _same_mode(head.source, bc) ||
        throw(ArgumentError("head and middle boundary mode keys differ"))
    bc.rb == head.boundary.rb ||
        throw(ArgumentError("head and middle correctors must share the same inner-ray boundary"))
    value = _checked_complex(T, head.source.Slm(bc.rb), "S_lm")
    _finite(value) || throw(ArgumentError("S_lm must be finite at the inner boundary"))
    IRGCorrectorMiddleMode{N,T,Fll,Flm,Fln}(head, bc, GHZIRG())
end
function solve_irg_corrector_middle(head::IRGCorrectorHeadMode, bc::MiddleRayIVP)
    throw(ArgumentError("head and middle boundary order, precision, or mode keys differ"))
end

function _radius(sol::IRGCorrectorMiddleMode{N,T}, r) where {N,T}
    _require_middle_spin_domain(sol.head.source.l)
    _require_middle_spin_domain(sol.boundary.l)
    rT = _checked_real(T, r, "middle evaluation radius")
    _finite(rT) && rT >= sol.boundary.rb && rT >= 2sol.head.source.M ||
        throw(ArgumentError("evaluation radius must be finite, outside the horizon, and r >= rb"))
    rT
end

function _Slm(sol::IRGCorrectorMiddleMode{N,T}, r::T) where {N,T}
    value = _checked_complex(T, sol.head.source.Slm(r), "S_lm")
    _finite(value) || throw(ArgumentError("S_lm must be finite along the integration ray"))
    value
end

function _middle_coupling(::Type{T}, l::Integer) where {T}
    lI = _require_middle_spin_domain(l)
    c = sqrt(T(lI) / T(2)) * sqrt(T(lI) + one(T))
    _finite(c) || throw(ArgumentError("middle angular coupling is not finite in $T for l=$lI"))
    c
end

"Evaluate `x_nm` with the one-integral Green formula for `z=r^2 x_nm`."
function x_nm(sol::IRGCorrectorMiddleMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    rT = _radius(sol, r)
    bc = sol.boundary
    at, rt = _tolerances(T, atol, rtol)
    rT == bc.rb && return bc.x0
    three = T(3)
    A = bc.dx0 + T(2) * bc.x0 / bc.rb
    ratio = bc.rb / rT
    homogeneous = ratio^2 * bc.x0 + A * rT * (one(T) - ratio^3) / three
    c = _middle_coupling(T, sol.head.source.l)
    k_slm(q) = one(T) - (q / rT)^3
    k_sll(q) = inv(T(2)) + (q / rT)^3 - T(3) * (q / rT)^2 / T(2)
    slm = _quad(q -> k_slm(q) * _Slm(sol, q), bc.rb, rT; atol=at, rtol=rt)
    sll = _quad(q -> k_sll(q) * _Sll(sol.head, q), bc.rb, rT; atol=at, rtol=rt)
    sourced = rT * (-T(2) * slm +
                    c * (sol.head.boundary.dx0 * k_sll(bc.rb) - sll)) / three
    homogeneous + sourced
end

"Evaluate `d x_nm/dr` by differentiating the one-integral Green formula."
function dr_x_nm(sol::IRGCorrectorMiddleMode{N,T}, r; atol=zero(T),
        rtol=sqrt(eps(T))) where {N,T}
    rT = _radius(sol, r)
    bc = sol.boundary
    at, rt = _tolerances(T, atol, rtol)
    rT == bc.rb && return bc.dx0
    three = T(3)
    A = bc.dx0 + T(2) * bc.x0 / bc.rb
    ratio = bc.rb / rT
    xhom = ratio^2 * bc.x0 + A * rT * (one(T) - ratio^3) / three
    c = _middle_coupling(T, sol.head.source.l)
    k_slm(q) = (one(T) + T(2) * (q / rT)^3) / three
    k_sll(q) = inv(T(6)) + (q / rT)^2 / T(2) - T(2) * (q / rT)^3 / three
    slm = _quad(q -> k_slm(q) * _Slm(sol, q), bc.rb, rT; atol=at, rtol=rt)
    sll = _quad(q -> k_sll(q) * _Sll(sol.head, q), bc.rb, rT; atol=at, rtol=rt)
    sourced = -T(2) * slm + c * (sol.head.boundary.dx0 * k_sll(bc.rb) - sll)
    A - T(2) * xhom / rT + sourced
end

"Pure residual of the outgoing-Schwarzschild IRG middle equation."
function irg_middle_residual(r, x, dx, d2x, Slm, dx_mmb, l)
    _require_nonbool_real(r, "middle residual radius")
    _require_nonbool_number(x, "x_nm")
    _require_nonbool_number(dx, "dr_x_nm")
    _require_nonbool_number(d2x, "d2r_x_nm")
    _require_nonbool_number(Slm, "S_lm")
    _require_nonbool_number(dx_mmb, "dr_x_mmb")
    l isa Integer || throw(ArgumentError("middle mode l must be an integer"))
    _finite(r) && r > zero(r) ||
        throw(ArgumentError("middle residual radius must be finite and positive"))
    all(_finite, (x, dx, d2x, Slm, dx_mmb)) ||
        throw(ArgumentError("middle residual inputs must be finite"))
    lI = _require_middle_spin_domain(l)
    T = typeof(float(r))
    c = _middle_coupling(T, lI)
    d2x + T(2) * dx / r - T(2) * x / r^2 + T(2) * Slm - c * dx_mmb / r
end

include("GHZ/CorrectorPair.jl")
include("GHZ/CorrectorTail.jl")
include("GHZ/CorrectorMonopole.jl")
include("GHZ/RayMomentCache.jl")
include("GHZ/AdjointSource.jl")
include("GHZ/HertzInversion.jl")
include("GHZ/HertzMetric.jl")

export SourceNormalization, EinsteinRHS, PhysicalStressEnergy
export KinnersleyMinusPlus, SchwarzschildOutgoing, ExpMinusIomegaUPlusImPhi, GHZIRG, HertzIRG
export TetradSourceMode, InnerRayIVP, MiddleRayIVP
export IRGCorrectorHeadMode, IRGCorrectorMiddleMode, perturbative_order
export solve_irg_corrector_head, x_mmb, dr_x_mmb, irg_head_residual
export solve_irg_corrector_middle, x_nm, dr_x_nm, irg_middle_residual
export PhysicalCorrectorMiddlePair, x_nmb, dr_x_nmb
export FiniteRayFreeIVP, TailRayIVP, IRGCorrectorTailMode
export solve_irg_corrector_tail, x_nn, dr_x_nn, irg_tail_residual
export IRGCorrectorTetradComponents, corrector_tetrad_components
export MonopoleSpinOneComponentsAbsent, MonopoleTailRayIVP
export IRGCorrectorMonopoleTailMode, IRGCorrectorMonopoleCoefficients
export solve_irg_corrector_monopole_tail, corrector_monopole_coefficients
export RayMomentCache, RayMomentProvenance, RayMomentComponents
export RayMomentQuadratureDiagnostics, RayMomentDiagnosticPoint
export materialize_grid, ray_moment_components, ray_quadrature_diagnostics
export TwoRealT0Dagger, ReducedAdjointSourceMode, EtaRayIVP, IRGEtaMode
export solve_irg_eta, eta, dr_eta, d2r_eta, irg_eta_residual
export FixedOutgoingURadialDerivative, FiniteHertzRadialIVP
export ConjugateHertzRadialJet, PhysicalHertzPair, IRGHertzMetricMode
export conjugate_hertz_radial_jet, physical_hertz_jet, hertz_metric, static_hertz_metric

end # module GHZ
