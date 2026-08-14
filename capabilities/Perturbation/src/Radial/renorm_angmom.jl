# MST renormalized angular momentum ν (Sasaki–Tagoshi). The first ingredient every MST radial /
# amplitude needs. Conventions: q=a, ε=2ω, κ=√(1−q²), τ=(ε−mq)/κ; λ in SWSH convention. Default
# method = :Auto: BigFloat inputs → precision-inherited :Monodromy, Float64 → machine-precision
# :FindRoot (real-ν branch PLUS the ½+iν_i and iν_i COMPLEX branches, |cos2πν|>1). :Series =
# closed-form O(ω⁴) seed. :Monodromy = confluent-Heun monodromy in Complex{BigFloat}.
#
# Consolidation vs the legacy Teukolsky form: the CF tails of the characteristic function reuse the
# promoted NumericsCore.cf_finite primitive (no local _cf_finite copy); the 5-arg auto-λ form pulls
# Harmonics.SpinWeightedSpheroidalEigenvalue. Type-generic throughout (BigFloat 256..2048-bit
# inherited through the Complex{BigFloat} monodromy path).

"Sasaki–Tagoshi (LRR 2003) Eq.173 small-ω Cos[2πν] series seed (depends on ω,s,l)."
function _cos2pinu_series(omega::Real, s::Integer, l::Integer)
    if l == 0
        return 1 - (8 * ((-11 + 15l + 15l^2)^2 * pi^2) * omega^4) /
                   ((1 - 2l)^2 * (1 + 2l)^2 * (3 + 2l)^2)
    else
        return 1 - (8 * (pi^2 * (30l^3 + 15l^4 + 3s^2 * (-2 + s^2) +
                                 l * (-11 + 6s^2) + l^2 * (4 + 6s^2))^2) * omega^4) /
                   ((1 - 2l)^2 * l^2 * (1 + l)^2 * (1 + 2l)^2 * (3 + 2l)^2)
    end
end

"MST recurrence coefficient closures (αγ, β) for the ν characteristic function (ST §2.1)."
function _mst_alphagamma_beta(s::Integer, m::Integer, a::Real, omega::Real, lambda::Real)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa
    alphagamma = (n, nu) -> eps^2 * kappa^2 * (n + nu) * (2 + n + nu) *
        ((1 + n + nu - s)^2 + eps^2) * ((1 + n + nu + s)^2 + eps^2) *
        (-1 + 2n + 2nu) * (5 + 2n + 2nu) * ((1 + n + nu)^2 + tau^2)
    beta = (n, nu) -> (2n + 2nu + 3) * (2n + 2nu - 1) *
        ((-lambda - s * (s + 1) + (n + nu) * (n + nu + 1) + eps^2 + eps * (eps - m * q)) *
             ((n + nu) * (n + nu + 1)) +
         eps * (eps - m * q) * (s^2 + eps^2))
    return alphagamma, beta
end

"MST characteristic function β[0,ν] + R[1,ν] + L[-1,ν] (finite nmax tails via NumericsCore.cf_finite)."
function _nu_characteristic(nu, alphagamma, beta, nmax::Integer)
    R1  = cf_finite(i -> -alphagamma(i - 1, nu), i -> beta(i, nu), 1, nmax)        # upward tail
    Lm1 = cf_finite(i -> -alphagamma(-2 - i, nu), i -> beta(-2 - i, nu), -1, nmax) # downward tail
    return beta(0, nu) + R1 + Lm1
end

"""
    _nu_rootfind_real(nu0, alphagamma, beta, nmax) -> Real

Solve `real(char(ν)) = 0` near `nu0` (real branch). The characteristic function has POLES (where
some β(n,ν)=0), so a bare secant can diverge: scan a window, pick the sign-change with the
SMALLEST magnitude (a root has |f|→0, a pole |f|→∞), then bracketed secant + bisection.
"""
function _nu_rootfind_real(nu0::Real, alphagamma, beta, nmax::Integer; tol::Real=1e-13)
    f = nu -> real(_nu_characteristic(nu, alphagamma, beta, nmax))
    lo, hi, N = nu0 - 0.45, nu0 + 0.45, 256        # window < 1 contains exactly one root (ν mod 1)
    h = (hi - lo) / (N - 1)
    bracket = nothing; bestmag = Inf
    xprev = lo; fprev = f(lo)
    for k in 1:(N - 1)
        xcur = lo + k * h; fcur = f(xcur)
        if isfinite(fprev) && isfinite(fcur) && signbit(fprev) != signbit(fcur)
            mag = max(abs(fprev), abs(fcur))
            mag < bestmag && (bestmag = mag; bracket = (xprev, fprev, xcur, fcur))
        end
        xprev, fprev = xcur, fcur
    end
    bracket === nothing && return NaN
    a, fa, b, fb = bracket
    for _ in 1:200                                  # bracketed secant, bisection fallback
        m = abs(fb - fa) > 1e-300 ? b - fb * (b - a) / (fb - fa) : (a + b) / 2
        (m <= min(a, b) || m >= max(a, b)) && (m = (a + b) / 2)
        fm = f(m)
        isfinite(fm) || (m = (a + b) / 2; fm = f(m))
        abs(b - a) <= tol * (1 + abs(m)) && return m
        if signbit(fm) == signbit(fa); a, fa = m, fm else b, fb = m, fm end
    end
    return (a + b) / 2
end

"""
    _nu_rootfind_imag(nui0, repart, alphagamma, beta, nmax) -> Real

Solve `real(char(repart + i·νi)) = 0` for `νi>0` — the COMPLEX-ν branches: `repart=1/2` when
`Cos[2πν]<−1` (ν=½+i νi), `repart=0` when `Cos[2πν]>1` (ν=i νi). Same pole-robust strategy.
"""
function _nu_rootfind_imag(nui0::Real, repart::Real, alphagamma, beta, nmax::Integer; tol::Real=1e-13)
    f = nui -> real(_nu_characteristic(complex(repart, nui), alphagamma, beta, nmax))
    lo, hi, N = max(1e-6, nui0 - 0.6), nui0 + 0.6, 256
    h = (hi - lo) / (N - 1)
    bracket = nothing; bestmag = Inf
    xprev = lo; fprev = f(lo)
    for k in 1:(N - 1)
        xcur = lo + k * h; fcur = f(xcur)
        if isfinite(fprev) && isfinite(fcur) && signbit(fprev) != signbit(fcur)
            mag = max(abs(fprev), abs(fcur))
            mag < bestmag && (bestmag = mag; bracket = (xprev, fprev, xcur, fcur))
        end
        xprev, fprev = xcur, fcur
    end
    bracket === nothing && return NaN
    a, fa, b, fb = bracket
    for _ in 1:200
        m = abs(fb - fa) > 1e-300 ? b - fb * (b - a) / (fb - fa) : (a + b) / 2
        (m <= min(a, b) || m >= max(a, b)) && (m = (a + b) / 2)
        fm = f(m); isfinite(fm) || (m = (a + b) / 2; fm = f(m))
        abs(b - a) <= tol * (1 + abs(m)) && return m
        if signbit(fm) == signbit(fa); a, fa = m, fm else b, fb = m, fm end
    end
    return (a + b) / 2
end

"""
    _nu_monodromy(s, l, m, a, omega, lambda; nmax=:auto, bits=256) -> Real|Complex

ν via monodromy of the confluent Heun equation (Casals–Ottewill). Builds the two-term a1/a2
recurrences + Pochhammer products in `Complex{BigFloat}` (the a1/a2 cancellation needs arbitrary
precision), assembles `Cos[2πν]`, then branch-selects ν. Reflection `Γ[w]Γ[-w]=-π/(w sin πw)`
(w=μ1C−μ2C) avoids complex BigFloat Γ. Returns the same ν as `:FindRoot` except at degenerate
boundary modes where the two BHPToolkit methods themselves disagree (Re[Cos2πν]≈−1).
"""
function _nu_monodromy(s::Integer, l::Integer, m::Integer, a::Real, omega::Real,
        lambda::Real; nmax::Union{Symbol,Integer}=:auto, bits::Integer=max(256, precision(BigFloat)))
    setprecision(BigFloat, bits) do
        CT = Complex{BigFloat}
        q = BigFloat(a); eps = 2 * BigFloat(omega); kappa = sqrt(1 - q^2)
        tau = (eps - m * q) / kappa; lam = BigFloat(lambda); im1 = im
        gammaCH = 1 - s - im1*eps - im1*tau                       # confluent-Heun params
        deltaCH = 1 + s + im1*eps - im1*tau
        epsCH   = 2*im1*eps*kappa
        aeCH    = 1 - s + im1*(eps - tau)                          # αCHεCH (single symbol)
        qCH     = -(-s*(1+s) + eps^2 + im1*(-1+2s)*eps*kappa - lam - tau*(im1+tau))
        w = (aeCH - (gammaCH+deltaCH)) - (-aeCH)                   # μ1C - μ2C

        build = function (nm::Integer)                            # a1[n],a2[n], index 1→n=0
            a1 = Vector{CT}(undef, nm+1); a2 = Vector{CT}(undef, nm+1)
            a1[1] = one(CT); a2[1] = one(CT)
            for n in 1:nm
                a1m2 = n == 1 ? zero(CT) : a1[n-1]                 # a1[-1]=0
                a1[n+1] = ((aeCH-(-1+n+deltaCH))*(aeCH-(-2+n+gammaCH+deltaCH))*epsCH * a1m2)/n -
                          ((aeCH^2 + aeCH*(1-2n-gammaCH-deltaCH+epsCH) +
                            (n^2 - qCH + n*(-1+gammaCH+deltaCH-epsCH) + epsCH - deltaCH*epsCH)) * a1[n])/n
                a2m2 = n == 1 ? zero(CT) : a2[n-1]
                a2[n+1] = -(((aeCH+(-2+n))*(aeCH+(-1+n-gammaCH))*epsCH * a2m2)/n) +
                          ((aeCH^2 + (n^2 - qCH + gammaCH + deltaCH - n*(1+gammaCH+deltaCH-epsCH) - epsCH) +
                            aeCH*(-1+2n-gammaCH-deltaCH+epsCH)) * a2[n])/n
            end
            a1, a2
        end
        cos2pinu = function (nm::Integer, a1, a2)                  # Pochhammer sums + assembly
            S1 = zero(CT); S2 = zero(CT)
            poch_p = Vector{CT}(undef, nm+1); poch_m = Vector{CT}(undef, nm+1)
            poch_p[1] = one(CT); poch_m[1] = one(CT)
            for i in 1:nm
                poch_p[i+1] = (w + i - 1) * poch_p[i]; poch_m[i+1] = (-w + i - 1) * poch_m[i]
            end
            for j in 0:cld(nm, 2)
                S1 += a1[j+1] * poch_p[nm-j+1]; S2 += (-1)^j * a2[j+1] * poch_m[nm-j+1]
            end
            pref = -2*BigFloat(pi)*w*sin(BigFloat(pi)*w) / (S1 * S2)
            cos(BigFloat(pi)*w) + pref * CT((-1)^(nm-1)) * a1[nm+1] * a2[nm+1]
        end

        local C
        if nmax isa Integer
            a1, a2 = build(nmax); C = cos2pinu(nmax, a1, a2)
        else
            # precision-INHERITED tol: the C-series locks by n≈128-256 and |ΔC| floors at the a1/a2
            # Pochhammer-cancellation roundoff floor ≈2^(22−bits); set tol to that floor + 2-bit margin.
            tolC = big(2.0)^(24 - bits)
            n = 8; a1, a2 = build(n); Cprev = cos2pinu(n, a1, a2); C = Cprev
            converged = false
            while n < max(400, 2*bits)             # cap scales with precision
                n = max(n + 4, n * 11 ÷ 10); a1, a2 = build(n); C = cos2pinu(n, a1, a2)
                if abs(C - Cprev) <= tolC * (1 + abs(C)); converged = true; break end
                Cprev = C
            end
            converged || @warn "_nu_monodromy: nmax cap hit without convergence" s l m a omega C
        end

        if imag(BigFloat(omega)) != 0                             # branch select
            return acos(C) / (2*BigFloat(pi))                     # Complex{BigFloat}, full precision
        end
        ReC = real(C)
        if ReC < -1
            return complex(BigFloat(1//2), abs(imag(acos(complex(ReC)) / (2*BigFloat(pi)))))
        elseif ReC <= 1
            return l - acos(ReC) / (2*BigFloat(pi))               # BigFloat, full precision
        else
            return complex(zero(BigFloat), abs(imag(acos(complex(ReC)) / (2*BigFloat(pi)))))
        end
    end
end

"""
    RenormalizedAngularMomentum(s, l, m, a, omega, lambda; method=:Auto, nmax=25)
    RenormalizedAngularMomentum(s, l, m, a, omega; ...)  # auto λ via Harmonics SWSH eigenvalue

MST renormalized angular momentum `ν` (Sasaki–Tagoshi). The 5-arg form auto-computes
`λ = SpinWeightedSpheroidalEigenvalue(s,l,m,a*ω)`. Methods: `:FindRoot` (real-ν via root-find +
the ½+iν_i / iν_i COMPLEX branches), `:Series` (closed-form small-ω seed, no root-find),
`:Monodromy` (confluent-Heun monodromy in `Complex{BigFloat}`; BHPToolkit's default).
"""
function RenormalizedAngularMomentum(s::Integer, l::Integer, m::Integer,
        a::Real, omega::Real, lambda::Real; method::Symbol=:Auto, nmax::Integer=25)
    if method === :Auto                       # BigFloat inputs → precision-inherited monodromy; Float64 → FindRoot
        method = (a isa BigFloat || omega isa BigFloat) ? :Monodromy : :FindRoot
    end
    l < abs(s) && return 0.0                  # unphysical (case A)
    iszero(omega) && return float(l)          # static spheroidal limit (case B)
    c2pn = _cos2pinu_series(omega, s, l)
    if method === :Series                     # closed-form seed, no root-find
        return -1 <= c2pn <= 1 ? l - acos(c2pn) / (2pi) :
               complex(c2pn < -1 ? 0.5 : 0.0, abs(imag(acos(complex(c2pn)))) / (2pi))
    elseif method === :Monodromy              # confluent-Heun monodromy (Complex{BigFloat})
        return _nu_monodromy(s, l, m, a, omega, lambda)
    elseif method !== :FindRoot
        throw(ArgumentError("RenormalizedAngularMomentum: method must be :FindRoot, :Series, or :Monodromy; got :$method"))
    end
    alphagamma, beta = _mst_alphagamma_beta(s, m, a, omega, lambda)
    if -1 <= c2pn <= 1                        # real-ν branch (the usual flux case)
        nu0 = l - acos(c2pn) / (2 * pi)
        return _nu_rootfind_real(nu0, alphagamma, beta, nmax)
    end
    # |cos2πν_series| > 1: the O(ω⁴) seed is unreliable at moderate ω / high ℓ and can FALSELY flag
    # the complex branch. Defer to the EXACT confluent-Heun monodromy, which branch-selects robustly.
    return _nu_monodromy(s, l, m, a, omega, lambda)
end

RenormalizedAngularMomentum(s::Integer, l::Integer, m::Integer, a::Real, omega::Real; kwargs...) =
    RenormalizedAngularMomentum(s, l, m, a, omega,
        SpinWeightedSpheroidalEigenvalue(s, l, m, a * omega); kwargs...)
