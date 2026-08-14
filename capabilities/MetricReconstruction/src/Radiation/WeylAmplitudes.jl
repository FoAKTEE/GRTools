# =============================================================================================================
# WeylAmplitudes.jl — Teukolsky–Starobinsky constants (calD/calC/Dhat/DhatPrime), the radial-asymptotic
# constants Γ/Γ̃ and in/out normalizations, and the Weyl-scalar amplitude presets that turn the Teukolsky
# radial/angular data into the CCK Hertz-potential reconstruction amplitudes {A, Ā, B, B̄}.
# Transcription of ExampleUsage.nb cell 7/17/19. λ₊₂ = teuk_lambda(2,…) (Modes.jl); horizons/Ω₊ delegated
# to Spacetime. Type-generic; the |m|≥2 Dhat branches form their integer factor exactly, then convert
# it to the caller's arithmetic type.  In particular, the m≤−2 formulas in the source notebook are
# ratios of factorials at negative integers: they must be evaluated through the finite Gamma-ratio
# continuation of arXiv:2403.20311 Eq. `ConstantsD`, not by calling `factorial` on singular endpoints.
#
# The {A,Ā,B,B̄} preset picks WHICH reconstruction (ψ0/ψ4, IRG/ORG) the metric represents; for the vacuum/
# IRG↔ORG-consistency gates the perturbation is vacuum for ANY choice, so those gates are amplitude-blind.
# The ψ0-IRG / ψ4-ORG presets additionally use the αin/αout incidence/transmission split through Γ,Γ̃; their
# numeric values match the .nb cell-21 oracle. The PRESET CHOICE is now CLOSED by the H10 physics gate in
# Radiation.WeylScalars: ψ4[h_ORG(:psi4_ORG)]/(ρ⁴·R^up_{−2}·S_{−2}) is (r,θ)-constant to <1e-9 (the regenerated
# gauge-invariant Weyl scalar has EXACTLY the Teukolsky s=−2 shape), so :psi4_ORG/:psi0_IRG are physics-validated.
# =============================================================================================================

using Spacetime: horizons, omega_plus

"Teukolsky–Starobinsky angular constant 𝒟 (cell 7). λ₊₂ = λ[2,ω,l,m]."
function calD(ω, l::Integer, m::Integer, a)
    λ2 = teuk_lambda(2, l, m, a, real(ω))
    t1 = (λ2 + 4)^2 * (λ2 + 6)^2
    t2 = 8*a*ω*(m - a*ω)*(λ2 + 4)*(5*λ2 + 26)
    t3 = 48*(a*ω)^2*(2*λ2 + 8 + 3*(m - a*ω)^2)
    t1 + t2 + t3
end

"Teukolsky–Starobinsky constant 𝒞 = 𝒟 + (12 ω M)² (M=1)."
calC(ω, l::Integer, m::Integer, a) = calD(ω, l, m, a) + (12*ω)^2

# angular polynomial helpers (cell 7)
_p1(ω, x, a) = x + 6*a*ω + 4
_p2(ω, x, a) = x^2 + 2*(4*a*ω + 5)*x + 4*(3*(a*ω)^2 + 8*a*ω + 6)
_p3(ω, x, a) = x^3 - 2*(3*a*ω - 8)*x^2 - 4*((a*ω)^2 + 8*a*ω - 21)*x +
               8*(3*(a*ω)^3 - 8*(a*ω)^2 - a*ω + 18)

# ExampleUsage.nb cell 7 and arXiv:2403.20311 define the |m|≥2 branches through
# one four-factor magnitude F(q)=(q-1)q(q+1)(q+2), q=|m|.  Widen before `abs`
# so even typemin(Int) is safe.  Keep the internal product exact, but convert it
# before arithmetic: mixing BigInt with Float32/Float64 otherwise widens public
# amplitudes to ambient-precision BigFloat.
function _ts_factor(m::Integer)
    q = abs(big(m))
    return (q - 1) * q * (q + 1) * (q + 2)
end
_factor_like(x, m::Integer) = oftype(x, _ts_factor(m))
_factor_seed(ω, a) = float(ω + zero(a))
function _validate_ts_mode(l::Integer, m::Integer)
    (l < 2 || abs(big(m)) > l) &&
        throw(ArgumentError("invalid (s,l,m)=(2,$l,$m): need l>=|s| and |m|<=l"))
    return nothing
end

"𝒟̂(ω,l,m): m-dependent Teukolsky–Starobinsky factor (cell 7 `Dhat`)."
function Dhat(ω, l::Integer, m::Integer, a::Real)
    if m >= 2
        d = calD(ω, l, m, a)
        return d / _factor_like(d, m)
    elseif m <= -2
        _validate_ts_mode(l, m)
        return _factor_like(_factor_seed(ω, a), m)
    end

    λ2 = teuk_lambda(2, l, m, a, real(ω))
    if m == 1
        return (-one(ω)/6) * _p3(ω, λ2, a)
    elseif m == 0
        return _p2(-ω, λ2, a)
    else  # m == -1
        return -6 * _p1(-ω, λ2, a)
    end
end

"𝒟̂′(ω,l,m): the ψ4 partner factor (cell 7 `DhatPrime`)."
function DhatPrime(ω, l::Integer, m::Integer, a::Real)
    if m >= 2
        _validate_ts_mode(l, m)
        return _factor_like(_factor_seed(ω, a), m)
    elseif m <= -2
        d = calD(ω, l, m, a)
        return d / _factor_like(d, m)
    end

    λ2 = teuk_lambda(2, l, m, a, real(ω))
    if m == 1
        return -6 * _p1(ω, λ2, a)
    elseif m == 0
        return _p2(ω, λ2, a)
    else  # m == -1
        return (-one(ω)/6) * _p3(-ω, λ2, a)
    end
end

# ── Radial-asymptotic constants Γ, Γ̃ and in/out normalizations (cell 17) ─────────────────────────────────
# Self-contained algebraic functions of (ω,m) through the outer-horizon data (r₊, r₋, Ω₊):
#   σ2 = r₊ − r₋,  w = 4(ω − m Ω₊) r₊,
#   Γ  = (w+2iσ2)(w+iσ2) w (w−iσ2),   Γ̃ = (w−2iσ2)(w−iσ2) w (w+iσ2),
#   in = Γ,  inPrime = 𝒞/Γ,  out = 𝒞/Γ̃,  outPrime = Γ̃.

"Outer/inner horizon radii r₊, r₋ (M=1) — delegated to Spacetime."
_horizons(a) = horizons(a)
"Horizon angular velocity Ω₊ = a/(2 M r₊) (M=1) — delegated to Spacetime."
_Omega_plus(a) = omega_plus(a)

"Γ(ω,m) radial-asymptotic constant (cell 17). M=1."
function _Gamma(ω, m::Integer, a)
    rp, rm = _horizons(a)
    σ2 = rp - rm
    w  = 4*(ω - m*_Omega_plus(a))*rp
    (w + 2*im*σ2)*(w + im*σ2)*w*(w - im*σ2)
end

"Γ̃(ω,m) radial-asymptotic constant (cell 17). M=1."
function _GammaTilde(ω, m::Integer, a)
    rp, rm = _horizons(a)
    σ2 = rp - rm
    w  = 4*(ω - m*_Omega_plus(a))*rp
    (w - 2*im*σ2)*(w - im*σ2)*w*(w + im*σ2)
end

"in[ω,l,m] = Γ(ω,m) (cell 17)."
radial_in(ω, l::Integer, m::Integer, a) = _Gamma(ω, m, a)
"out[ω,l,m] = 𝒞/Γ̃ (cell 17)."
radial_out(ω, l::Integer, m::Integer, a) = calC(ω, l, m, a) / _GammaTilde(ω, m, a)
"inPrime[ω,l,m] = 𝒞/Γ (cell 17)."
radial_inPrime(ω, l::Integer, m::Integer, a) = calC(ω, l, m, a) / _Gamma(ω, m, a)
"outPrime[ω,l,m] = Γ̃(ω,m) (cell 17)."
radial_outPrime(ω, l::Integer, m::Integer, a) = _GammaTilde(ω, m, a)

"""
    weyl_amplitudes(which, ω, l, m, a; M=1, αin=1, αout=0) -> NamedTuple{(:A,:Abar,:B,:Bbar)}

Weyl-scalar amplitude preset (cell 19). `which`:
  `:psi0_ORG` → A = 192iωM/𝒞,  B = 16 conj(𝒟̂)/conj(𝒞)
  `:psi4_IRG` → A = −192iωM/𝒞, B = 16 conj(𝒟̂′)/conj(𝒞)
  `:psi0_IRG` → A = 0, B = 4/(αin·in + αout·out)                       (in=Γ, out=𝒞/Γ̃)   [H10 CLOSED]
  `:psi4_ORG` → A = 0, B = 64/conj(αin·inPrime + αout·outPrime)         (inPrime=𝒞/Γ, outPrime=Γ̃) [H10 CLOSED]
  `:unit`     → A = B = 1 (generic non-trivial amplitude for the vacuum/consistency gates)
Ā,B̄ are the conjugates. `αin,αout` select the incidence/transmission mix (.nb cell-21 uses 1,0).
"""
function weyl_amplitudes(which::Symbol, ω, l::Integer, m::Integer, a; M=1, αin=1, αout=0)
    cC = calC(ω, l, m, a)
    if which === :psi0_ORG
        A = 192*im*ω*M / cC
        B = 16 * conj(Dhat(ω, l, m, a)) / conj(cC)
        return (A=A, Abar=conj(A), B=B, Bbar=conj(B))
    elseif which === :psi4_IRG
        A = -192*im*ω*M / cC
        B = 16 * conj(DhatPrime(ω, l, m, a)) / conj(cC)
        return (A=A, Abar=conj(A), B=B, Bbar=conj(B))
    elseif which === :unit
        A = complex(one(float(real(ω))))
        B = complex(one(float(real(ω))))
        return (A=A, Abar=conj(A), B=B, Bbar=conj(B))
    elseif which === :psi0_IRG
        # [H10 CLOSED] preset choice physics-validated by Radiation.WeylScalars (ψ0[h_IRG]∝R^up_{+2}·S_{+2}).
        B = 4 / (αin*radial_in(ω, l, m, a) + αout*radial_out(ω, l, m, a))
        return (A=zero(B), Abar=zero(B), B=B, Bbar=conj(B))
    elseif which === :psi4_ORG
        # [H10 CLOSED] preset choice physics-validated by Radiation.WeylScalars (ψ4[h_ORG]==Teukolsky s=−2 shape).
        # The .nb conjugates the (αin·inPrime+αout·outPrime) sum, not αin/αout separately.
        B = 64 / conj(αin*radial_inPrime(ω, l, m, a) + αout*radial_outPrime(ω, l, m, a))
        return (A=zero(B), Abar=zero(B), B=B, Bbar=conj(B))
    else
        error("unknown Weyl preset $which")
    end
end
