# Eccentric disc contexts. The circular background remains the sole owner of
# the equilibrium arithmetic; its reference radius is independent of any orbit.
using ..ForcingModes: HarmonicMode, ClosureTag, Adiabatic, LocallyIsothermalOptionA,
                      SofteningPolicy, FixedSoftening

_closure_tag(c::ClosureTag) = c in (Adiabatic, LocallyIsothermalOptionA) ? c :
    throw(ArgumentError("SpiralWaves supports adiabatic or locally isothermal OPTION A"))
function _closure_tag(c::Symbol)
    c === :adiabatic && return Adiabatic
    c === :locally_isothermal && return LocallyIsothermalOptionA
    throw(ArgumentError("closure must be :adiabatic or :locally_isothermal"))
end

"""
    DiscModel(; χ, r_ref, h_ref, pSigma=0, qTemp=0, epsilon_ref=1,
              closure, evolution, εp=1e-5, σε=nothing)
    DiscModel(bg::DiscBackground; closure, evolution)

Power-law disc in M=1 units, normalized at the independently fixed `r_ref`:
`c_T(r_ref)=h_ref*r_ref*Omega_ref(model)`. The equilibrium is shared by the
adiabatic and locally isothermal OPTION A closures. Both `closure` and
`evolution=:fixed|:moving` are mandatory. The latter records the requested
evolution policy; changing the normalization requires constructing a new model
in the evolution layer. It never silently follows a perturber.

`reference` stores the frozen circular profile evaluator: its `aEMRI` slot is
used ONLY as r_ref and its softening is unused. No orbit or perturber is stored.
The read-only aliases χ, r_ref, h_ref, pSigma, qTemp, epsilon_ref, εp and σε
expose the physical parameters. BigFloat inputs inherit precision.
"""
struct DiscModel{T<:AbstractFloat}
    reference::DiscBackground{T}
    closure::ClosureTag
    evolution::Symbol
    function DiscModel(bg::DiscBackground{T}; closure, evolution::Symbol) where {T<:AbstractFloat}
        evolution in (:fixed, :moving) || throw(ArgumentError("evolution must be :fixed or :moving"))
        all(isfinite, (bg.χ,bg.aEMRI,bg.h0,bg.𝔭,bg.𝔮,bg.ε0,bg.εp,bg.σε)) ||
            throw(ArgumentError("disc parameters must be finite"))
        new{T}(bg, _closure_tag(closure), evolution)
    end
end

function DiscModel(; χ, r_ref, h_ref, pSigma=0, qTemp=0, epsilon_ref=1,
                   closure, evolution, εp=1e-5, σε=nothing)
    bg = DiscBackground(; χ, aEMRI=r_ref, h0=h_ref, 𝔭=pSigma, 𝔮=qTemp,
                        ε0=epsilon_ref, b=0, εp, σε)
    DiscModel(bg; closure, evolution)
end

function Base.getproperty(d::DiscModel, name::Symbol)
    bg = getfield(d, :reference)
    name === :χ && return bg.χ
    name === :r_ref && return bg.aEMRI
    name === :h_ref && return bg.h0
    name === :pSigma && return bg.𝔭
    name === :qTemp && return bg.𝔮
    name === :epsilon_ref && return bg.ε0
    name === :εp && return bg.εp
    name === :σε && return bg.σε
    getfield(d, name)
end
Base.propertynames(::DiscModel) =
    (:reference,:closure,:evolution,:χ,:r_ref,:h_ref,:pSigma,:qTemp,:epsilon_ref,:εp,:σε)

"Circular Kerr reference frequency, in inverse M; independent of perturber frequencies."
Omega_ref(d::DiscModel) = disc_OmegaK(d, d.r_ref)

# These methods deliberately call the frozen evaluator, preserving its operation
# order for circular regression, instead of repeating the equilibrium formulas.
for f in (:disc_OmegaK, :kappar2, :kappa2, :eps_profile, :eden0, :csound,
          :pressure, :gamma1, :dlogp_dr_over_enthalpy, :rotation_profile,
          :Omega, :dOmega_dr, :dv_dr, :rotation_derivs)
    @eval $f(d::DiscModel, r) = $f(d.reference, r)
end
r_isco(d::DiscModel) = r_isco(d.reference)
F_J0(d::DiscModel; q=1) = F_J0(d.reference; q)
F_J1(d::DiscModel; q=1) = F_J1(d.reference; q)

"""
    EccentricPerturber(orbit::EquatorialOrbit; mass=1, softening)

Spacetime-owned equatorial worldline and a shared softening policy. The mass
and policy must fit the orbit's working type without lowering precision.
Neither the orbit nor its evolution changes the disc reference normalization.
"""
struct EccentricPerturber{T,O<:EquatorialOrbit{T},S<:SofteningPolicy}
    orbit::O
    mass::T
    softening::S
    function EccentricPerturber(orbit::EquatorialOrbit{T}; mass=1,
                               softening::SofteningPolicy{U}) where {T,U}
        promote_type(T,typeof(float(mass)),U) === T ||
            throw(ArgumentError("perturber mass/softening must not exceed orbit precision type"))
        isfinite(mass) && mass > 0 || throw(ArgumentError("perturber mass must be finite and positive"))
        new{T,typeof(orbit),typeof(softening)}(orbit,T(mass),softening)
    end
end

"""
    ModeContext(disc, mode::HarmonicMode; closure, domain=(r_isco(disc),5disc.r_ref))
    ModeContext(disc, perturber, m, k; closure, domain=...)
    ModeContext(disc, perturber, mode::HarmonicMode; closure, domain=...)

Frequency-aware operator context. `mode.omega` is the total coordinate-time
frequency in inverse M, never a pattern speed. Its type must match the disc.
The mandatory closure must agree with the disc. `r_c` is the root in the
physical `domain`, including endpoints, or `nothing`; no reference radius is
substituted for an absent root. The static axisymmetric (0,0-frequency) operator
is singular everywhere and is rejected. Other signed/zero frequencies are valid.

The perturber overload checks spin/type and ω=mΩφ+kΩr through Spacetime's sole
frequency implementation. The direct mode overload also supports arbitrary
prescribed frequencies. Source jets are supplied at evaluation time.
"""
struct ModeContext{T<:AbstractFloat}
    disc::DiscModel{T}
    mode::HarmonicMode{T}
    closure::ClosureTag
    domain::NTuple{2,T}
    r_c::Union{Nothing,T}
    function ModeContext(d::DiscModel{T}, mode::HarmonicMode; closure,
                         domain=(r_isco(d),5d.r_ref)) where {T}
        typeof(mode.omega) === T || throw(ArgumentError("disc and mode frequency types must match"))
        tag = _closure_tag(closure)
        tag == d.closure || throw(ArgumentError("context closure must match DiscModel"))
        mode.m == 0 && iszero(mode.omega) &&
            throw(ArgumentError("static axisymmetric operator is singular everywhere"))
        lo,hi = T.(domain)
        all(isfinite,(lo,hi)) && r_isco(d) <= lo < hi ||
            throw(ArgumentError("domain must be finite, ordered and outside the ISCO"))
        rc = corotation_radius(d,mode;domain=(lo,hi))
        new{T}(d,mode,tag,(lo,hi),rc)
    end
end

function ModeContext(d::DiscModel{T}, p::EccentricPerturber, mode::HarmonicMode;
                     closure, kwargs...) where {T}
    typeof(p.orbit.a) === T && p.orbit.a == d.χ ||
        throw(ArgumentError("disc and orbit spin/type must match"))
    ω = orbit_frequency(p.orbit,mode.m,mode.k)
    scale = max(abs(ω),abs(mode.m*p.orbit.Ωφ),abs(mode.k*p.orbit.Ωr))
    abs(mode.omega-ω) <= 64eps(T)*scale ||
        throw(ArgumentError("mode frequency does not match mΩφ+kΩr"))
    ModeContext(d,mode;closure,kwargs...)
end
function ModeContext(d::DiscModel, p::EccentricPerturber, m::Integer, k::Integer; kwargs...)
    ModeContext(d,p,HarmonicMode(m,k,orbit_frequency(p.orbit,m,k));kwargs...)
end

omega_d(c::ModeContext, r) = c.mode.omega-c.mode.m*Omega(c.disc,r)
corotation_radius(c::ModeContext) = c.r_c

"""
Complex source inputs for the emitted eccentric kernels (R3 alphabet).
The geometry, equilibrium, and thermodynamics are real; omega_d, all metric
amplitudes (including htr/hrphi), and the response y0/y1 are Complex{T}.
No circular zeroing or real projection is applied.
"""
struct EccentricMasterInputs{T<:AbstractFloat}
    m::Int
    enu::T; epsi::T; emu::T; psi_r::T; mu_r::T; omt::T; omt_r::T
    W::T; v::T; v_r::T; p::T; p_r::T; eps::T; eps_r::T; Gamma1::T
    omega_d::Complex{T}
    htt::Complex{T}; htt_r::Complex{T}; htphi::Complex{T}; htphi_r::Complex{T}
    hphiphi::Complex{T}; hphiphi_r::Complex{T}; hrr::Complex{T}; hzz::Complex{T}
    htr::Complex{T}; hrphi::Complex{T}
    y0::Complex{T}; y1::Complex{T}
end

"Zero covariant source jets in the context's working type (operator evaluation)."
function _zero_source(c::ModeContext{T}) where {T}
    z = zero(Complex{T})
    (;htt=z,htt_r=z,htphi=z,htphi_r=z,hphiphi=z,hphiphi_r=z,hrr=z,hzz=z,htr=z,hrphi=z)
end

"""
    mode_cache_key(context)

Value key for operator/grid work only, including closure, frequency, domain,
disc normalization and stored/current precision. A response/source cache must
add its own source identity. No hash of the source is inferred from (m,k).
"""
function mode_cache_key(c::ModeContext{T}) where {T}
    d = c.disc
    (; closure=c.closure, type=T, precision=precision(T), stored_precision=precision(d.r_ref),
       disc=(d.χ,d.r_ref,d.h_ref,d.pSigma,d.qTemp,d.epsilon_ref,d.εp,d.σε,d.evolution),
       mode=c.mode, domain=c.domain, r_c=c.r_c)
end

"Context solution; closure is taken from the validated context, never defaulted."
struct EccentricModeSolution{T,F,P}
    context::ModeContext{T}
    closure::ClosureTag
    r::Vector{T}
    y::Matrix{Complex{T}}
    source::F
    residual_abs::T
    residual::T
    bc_residual::T
    params::P
    seconds::Float64
    function EccentricModeSolution(c::ModeContext{T},r,y,source,ra,rr,bcres,params,seconds) where {T}
        new{T,typeof(source),typeof(params)}(c,c.closure,r,y,source,ra,rr,bcres,params,seconds)
    end
end
Base.length(sol::EccentricModeSolution) = length(sol.r)

"Caller-owned operator cache, bound to closure/context, grid values and Landau setting."
struct EccentricOperatorCache{T,K}
    key::K
    r::Vector{T}
    Cmid::Vector{Matrix{Complex{T}}}
    Cedges::NTuple{2,Matrix{Complex{T}}}
end
