# MST scattering-phase pipeline (numerical branch). A strict dep-SUBSET of the radial machinery:
# REUSES the already-built ν/MST/amplitude tier (RenormalizedAngularMomentum, _mst_fn[_array],
# _mst_amplitudes) + the Harmonics spheroidal eigenvalue, rather than re-porting them. Notebook
# convention: ε = 2Mω (so ω = ε/2), q = a (M=1). The analytic ε-expansion branch is the closed-form
# low-frequency series (firewall-safe: a plain Julia polynomial in ε with γE/π/ζ3/log, NO Symbolics).

"""
    starobinsky_constant(s, l, m, q, ε; P=+1) -> Complex

Teukolsky–Starobinsky spin-normalization constant C* for spin weight `s`, multipole `(l, m)`,
Kerr spin `q = a`, `ε = 2Mω`. `P = ±1` selects the sign of the imaginary Teukolsky–Starobinsky
term for `s = −2`. Closed-form in the spheroidal eigenvalues `λ = SWSH(s,…)`, `negλ = SWSH(−s,…)`
(no continued fraction). Defined for `s ∈ {0, −1, −2}`.
"""
function starobinsky_constant(s::Integer, l::Integer, m::Integer, q::Real, ε::Real; P::Integer=+1)
    ω = ε / 2
    aω = q * ω
    if s == 0
        return complex(one(float(ω)))
    elseif s == -1
        negλ = SpinWeightedSpheroidalEigenvalue(-s, l, m, aω)
        radicand = (negλ + 2 + q^2 * ω^2 - 2 * m * q * ω)^2 + 4 * m * q * ω - 4 * q^2 * ω^2
        return sqrt(complex(radicand))
    elseif s == -2
        λ = SpinWeightedSpheroidalEigenvalue(s, l, m, aω)
        radicand = ((λ + 2)^2 + 4 * q * m * ω - 4 * q^2 * ω^2) *
                   (λ^2 + 36 * q * m * ω - 36 * q^2 * ω^2) +
                   (2 * λ + 3) * (96 * q^2 * ω^2 - 48 * q * ω * m) -
                   144 * ω^2 * q^2
        return sqrt(complex(radicand)) + 12 * P * im * ω
    else
        throw(ArgumentError("starobinsky_constant: spin weight s=$s unsupported (defined for s ∈ {0, −1, −2})"))
    end
end

"""
    _nu_lambda_omega(s, l, m, q, ε) -> (ν, λ, ω)

Convention shim: map the notebook `ε = 2Mω`, `q = a` to `(a=q, ω=ε/2)`, and return the spheroidal
eigenvalue λ and renormalized angular momentum ν. REUSES the Radial ν-solve + Harmonics eigenvalue.
"""
function _nu_lambda_omega(s::Integer, l::Integer, m::Integer, q::Real, ε::Real)
    ω = ε / 2
    λ = SpinWeightedSpheroidalEigenvalue(s, l, m, q * ω)
    ν = RenormalizedAngularMomentum(s, l, m, q, ω, λ)
    return ν, λ, ω
end

"""
    _brefl_over_binc(s, l, m, q, ε; n=10, nUp=5) -> Complex

MST reflection/incidence ratio `Bref/Binc` of the Teukolsky "In" solution = `_mst_amplitudes`
`InRef / InInc`. The notebook prefactors `ω^{−2s} e^{−iε(1−κ)} e^{2iε log ε}` are exactly the
folded ratio of the package's UpTrans/InInc prefactors, so no extra prefactor is applied.
"""
function _brefl_over_binc(s::Integer, l::Integer, m::Integer, q::Real, ε::Real;
                          n::Integer=10, nUp::Integer=5)
    ν, λ, ω = _nu_lambda_omega(s, l, m, q, ε)
    amp = _mst_amplitudes(s, m, q, ω, ν, λ; Nmax=max(n, nUp))
    return amp.InRef / amp.InInc
end

"""
    mst_coefficient(nn, s, l, m, q, ε; nUp=5) -> Complex

The `nn`-th renormalized MST coefficient `f_nn` of the "In" solution (`f_0 = 1`, MST gauge). Thin
accessor over `_mst_fn` fed `(ν, λ, ω)` from `_nu_lambda_omega`. `nUp` is accepted for parity with
the notebook; the self-converging continued fraction makes the result independent of it.
"""
function mst_coefficient(nn::Integer, s::Integer, l::Integer, m::Integer, q::Real, ε::Real;
                         nUp::Integer=5)
    ν, λ, ω = _nu_lambda_omega(s, l, m, q, ε)
    return _mst_fn(s, m, q, ω, ν, λ, nn)
end

"""
    mst_convergence(s, l, m, q, ε; N=20) -> NamedTuple

MST-series convergence diagnostic. Returns `(a0, tail, ratio, monotone)`: `a0` the n=0 coefficient,
`tail = max(|f_N|,|f_{−N}|)/|f_0|`, `ratio = |f_N|/|f_{N−1}|`, `monotone` iff `|f_n|` is
non-increasing for n=1…N.
"""
function mst_convergence(s::Integer, l::Integer, m::Integer, q::Real, ε::Real; N::Integer=20)
    ν, λ, ω = _nu_lambda_omega(s, l, m, q, ε)
    farr, _, _, _ = _mst_fn_array(s, m, q, ω, ν, λ, N)
    a0 = farr[0]
    tail = max(abs(farr[N]), abs(farr[-N])) / abs(a0)
    ratio = abs(farr[N]) / abs(farr[N - 1])
    monotone = all(abs(farr[n + 1]) <= abs(farr[n]) for n in 1:(N - 1))
    return (a0=a0, tail=tail, ratio=ratio, monotone=monotone)
end

"""
    scattering_phase_numerical(s, l, m, q, ε; n=10, nUp=5, P=+1) -> Complex

Numerical MST scattering phase factor `(−1)^{l+1} · C* · (2ω)^{2s} · Bref/Binc` (ω = ε/2). `C*` is
the Teukolsky–Starobinsky constant (`P=±1` its branch); `Bref/Binc` the MST reflection/incidence
ratio (truncation `Nmax = max(n, nUp)`). The scattering phase is `δ = log(value) / (2i)`.
"""
function scattering_phase_numerical(s::Integer, l::Integer, m::Integer, q::Real, ε::Real;
                                    n::Integer=10, nUp::Integer=5, P::Integer=+1)
    ω = ε / 2
    Cstar = starobinsky_constant(s, l, m, q, ε; P=P)
    BoB = _brefl_over_binc(s, l, m, q, ε; n=n, nUp=nUp)
    return (-1)^(l + 1) * Cstar * (2ω)^(2s) * BoB
end

# ζ(3) (Apéry's constant) to ~100 digits — parsed at the working precision T (BigFloat-safe). Not a
# Julia Irrational, so it is carried as a literal rather than `T(...)` of a built-in constant.
const _ZETA3_STR = "1.2020569031595942853997381615114499907649862923404988817922715553418382057863130901864558736093352581"
_zeta3(::Type{T}) where {T<:AbstractFloat} = parse(T, _ZETA3_STR)

# Published low-frequency ε-series of the scattering phase δ (= log(e^{2iδ})/2i) for the modes
# tabulated in NumericalMST.nb (q = a, ε = 2Mω; exact through ε³, O(ε⁴) remainder). Transcribed
# VERBATIM from the notebook cells. Each is type-generic so BigFloat precision is inherited.
function _phase000(q, ε)               # NumericalMST.nb In[10]: s=0, l=0, m=0
    T = float(real(promote_type(typeof(q), typeof(ε))))
    γ = T(MathConstants.eulergamma); π_ = T(pi); ζ3 = _zeta3(T); qT = T(q); εT = T(ε)
    κ = sqrt(1 - qT^2)
    return (-one(T) + 2γ + 2log(T(2)) + 2log(εT)) / 2 * εT +
        (6im - 6im*qT^2 + 6im*κ + 11π_*κ) * εT^2 / (12κ) +
        (1 / (36κ)) * (-33 + 18im*π_ + 33qT^2 - 18im*π_*qT^2 + 51κ - 36γ*κ + 18im*π_*κ +
            11π_^2*κ + 3qT^2*κ - 36κ*log(T(2)) - 36κ*log(κ*εT) - 12κ*ζ3) * εT^3
end

function _phase010(q, ε)               # NumericalMST.nb In[8]: s=0, l=1, m=0
    T = float(real(promote_type(typeof(q), typeof(ε))))
    γ = T(MathConstants.eulergamma); π_ = T(pi); ζ3 = _zeta3(T); qT = T(q); εT = T(ε)
    return complex((-3 + 2γ + 2log(T(2)) + 2log(εT)) / 2 * εT + 19π_ * εT^2 / 60 +
        (19π_^2 - 9qT^2 - 60ζ3) * εT^3 / 180)
end

function _phase110(q, ε)               # NumericalMST.nb In[9]: s_phys=1 (input s=-1), l=1, m=0
    T = float(real(promote_type(typeof(q), typeof(ε))))
    γ = T(MathConstants.eulergamma); π_ = T(pi); ζ3 = _zeta3(T); qT = T(q); εT = T(ε)
    return (-5 + 4γ + 4log(T(2)) + 4log(εT)) / 4 * εT + (94π_ - 15im*qT^2) * εT^2 / 240 +
        (105 + 188π_^2 - 180qT^2 - 480ζ3) * εT^3 / 1440
end

function _phase220(q, ε)               # NumericalMST.nb In[5]: s_phys=2 (input s=-2, P=-1), l=2, m=0
    T = float(real(promote_type(typeof(q), typeof(ε))))
    γ = T(MathConstants.eulergamma); π_ = T(pi); ζ3 = _zeta3(T); qT = T(q); εT = T(ε)
    return complex((-5 + 3γ + 3log(T(2)) + 3log(εT)) / 3 * εT + 107π_ * εT^2 / 420 +
        (-1015 + 1926π_^2 - 945qT^2 - 7560ζ3) * εT^3 / 22680)
end

"""
    scattering_phase_analytic(s, l, m, q, ε) -> Complex

Closed-form low-frequency ε-series of the scattering phase `δ` (= log(e^{2iδ})/2i) for the spin
weight `s`, multipole `(l, m)`, Kerr spin `q = a`, `ε = 2Mω`. Exact through `ε³` (O(ε⁴) remainder).
The published `phase{slm}Expansion` coefficients (NumericalMST.nb) are transcribed as explicit
polynomials in `ε` with `γE, π, ζ(3), log` — NO symbolic-Series dependency (the light-dep firewall).
Defined for the tabulated modes `(0,0,0), (0,1,0), (-1,1,0), (-2,2,0)`; the `s=-2` series matches
the numerical branch at `P=-1`. Type-generic, so BigFloat precision is inherited.
"""
function scattering_phase_analytic(s::Integer, l::Integer, m::Integer, q::Real, ε::Real)
    key = (Int(s), Int(l), Int(m))
    key == (0, 0, 0)  && return _phase000(q, ε)
    key == (0, 1, 0)  && return _phase010(q, ε)
    key == (-1, 1, 0) && return _phase110(q, ε)
    key == (-2, 2, 0) && return _phase220(q, ε)
    throw(ArgumentError("scattering_phase_analytic: no published ε-series for (s,l,m)=$key " *
        "(available: (0,0,0), (0,1,0), (-1,1,0), (-2,2,0))"))
end

"""
    scattering_phase(s, l, m, q, ε; mode=:numerical, kwargs...) -> Complex

`mode=:numerical` returns the scattering phase FACTOR `e^{2iδ}` (MST series, kwargs `(n, nUp, P)`
forwarded); `mode=:analytic` returns the scattering PHASE `δ` itself via the closed-form
low-frequency ε-series. Spin weight `s`, multipole `(l, m)`, Kerr spin `q = a`, `ε = 2Mω`.
"""
function scattering_phase(s, l, m, q, ε; mode::Symbol=:numerical, kwargs...)
    mode === :analytic  && return scattering_phase_analytic(s, l, m, q, ε)
    mode === :numerical && return scattering_phase_numerical(s, l, m, q, ε; kwargs...)
    throw(ArgumentError("mode must be :analytic or :numerical, got $mode"))
end

"""
    scattering_phase_table(s, l, m, q, εvals; n=10, nUp=5, P=+1) -> Vector

Scattering phase factor over a grid of `ε` values. Batch wrapper over
`scattering_phase_numerical` — each grid point is the identical per-point call, so the result
equals the elementwise serial broadcast; only the loop is parallelized (`Threads.@threads`). The
untyped buffer + final comprehension reproduce the serial element-type widening exactly (the table
is byte-identical to serial).
"""
function scattering_phase_table(s::Integer, l::Integer, m::Integer, q::Real, εvals;
                                n::Integer=10, nUp::Integer=5, P::Integer=+1)
    εv = collect(εvals)
    buf = Vector{Any}(undef, length(εv))
    Threads.@threads for i in eachindex(εv)
        buf[i] = scattering_phase_numerical(s, l, m, q, εv[i]; n=n, nUp=nUp, P=P)
    end
    return [x for x in buf]
end
