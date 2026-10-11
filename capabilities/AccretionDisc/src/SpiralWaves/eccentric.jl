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
