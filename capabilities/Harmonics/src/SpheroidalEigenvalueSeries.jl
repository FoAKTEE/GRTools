# SpheroidalEigenvalueSeries — asymptotic-series eigenvalue backend (design §4.2): small-γ Taylor
# (Stewart/Press-Teukolsky spheroidal expansion) + large-γ Laurent (Casals-Ottewill-Warburton,
# arXiv:1810.00432). Pure recurrences, independent of the spectral/Leaver dispatcher (a sibling
# backend). The shared small-γ d-tables it builds (`_swsh_series_tables`) feed the harmonic Series
# in SpheroidalHarmonics.
#
# FIDELITY NOTE: the recurrences call `_swsh_alpha` at SHIFTED angular indices `l+j+k` that can be
# zero/negative even for physical `l`. Mathematica evaluates `Sqrt[-n]` as the finite complex
# `I·Sqrt[n]`, so `0·alpha = 0` there. We mirror this with COMPLEX sqrt — the (possibly complex)
# value is always multiplied by a structurally-zero `d`-coefficient and contributes 0, so the final
# real eigenvalue is reproduced exactly. A real `sqrt(-n)` would give NaN and poison the sum.

"""
    _swsh_alpha(s, l, m) -> ComplexF64

Spheroidal off-diagonal coupling `alpha[s,l,m]`. Dispatch precedence mirrors the Mathematica rule
ordering: the `s==0` and `m==0` forms beat the `l==0` literal, which beats the general form.
"""
function _swsh_alpha(s::Integer, l::Integer, m::Integer)
    (2l == 1 || 2l == -1) && return 0.0 + 0.0im            # alpha[s,±1/2,m]=0
    den = sqrt(complex(2l - 1)) * sqrt(complex(2l + 1))
    if s == 0
        return sqrt(complex(l^2 - m^2)) / den
    elseif m == 0
        return sqrt(complex(l^2 - s^2)) / den
    elseif l == 0
        return 0.0 + 0.0im
    else
        return (1 / l) * sqrt(complex(l^2 - m^2)) * sqrt(complex(l^2 - s^2)) / den
    end
end

"Spheroidal diagonal coupling `beta[s,l,m] = -(m s)/(l(l+1))`; 0 at `l in {0,-1}`."
_swsh_beta(s::Integer, l::Integer, m::Integer) = (l == 0 || l == -1) ? 0.0 : -(m * s) / (l * (l + 1))

"""
    _swsh_series_tables(s, l, m, order) -> (Lam, Dd)

Build the small-γ series tables up to `order`: `Lam[i+1]` are the eigenvalue Taylor coefficients
`slmLambda[s,l,m][i]` and `Dd[(i,j)]` the harmonic mixing coefficients `d[s,l,m][i,j]`. Mutually
recursive but strictly lower-order, so one bottom-up pass `i=0..order` suffices.
"""
function _swsh_series_tables(s::Integer, l::Integer, m::Integer, order::Integer)
    Dd  = Dict{Tuple{Int,Int},ComplexF64}()
    Lam = Vector{ComplexF64}(undef, order + 1)
    dget(i, j) = get(Dd, (i, j), 0.0 + 0.0im)             # 0 for i<0, |j|>i, or unfilled
    α(ll) = _swsh_alpha(s, ll, m)
    β(ll) = _swsh_beta(s, ll, m)

    Dd[(0, 0)] = 1.0 + 0.0im
    Lam[1]     = l * (l + 1) - s * (s + 1)                # slmLambda[0]
    for i in 1:order
        if i == 1                                         # explicit small cases
            Dd[(1, 1)]  = -(s * α(l + 1)) / (l + 1)
            Dd[(1, 0)]  = 0.0 + 0.0im
            Dd[(1, -1)] = l == 0 ? 0.0 + 0.0im : (s * α(l)) / l
        else
            for j in -i:i
                if 1 + j + 2l == 0                        # overrides general
                    Dd[(i, j)] = 0.0 + 0.0im
                elseif j == 0
                    acc = 0.0 + 0.0im
                    for k in 1:(i - 1), jj in -i:i
                        acc += dget(i - k, jj) * dget(k, jj)
                    end
                    Dd[(i, 0)] = -acc / (2 * Dd[(0, 0)])
                else                                      # general j≠0
                    acc = 0.0 + 0.0im
                    for k in 1:(i - 1)
                        acc += dget(i - k, j) * Lam[k + 1]
                    end
                    acc += α(j + l - 1) * α(j + l) * dget(i - 2, j - 2)
                    acc += α(j + l) * (β(j + l - 1) + β(j + l)) * dget(i - 2, j - 1)
                    acc += (-1 + α(j + l)^2 + α(1 + j + l)^2 + β(j + l)^2) * dget(i - 2, j)
                    acc += α(1 + j + l) * (β(j + l) + β(1 + j + l)) * dget(i - 2, j + 1)
                    acc += α(1 + j + l) * α(2 + j + l) * dget(i - 2, j + 2)
                    acc += -2 * s * α(j + l) * dget(i - 1, j - 1)
                    acc += 2 * (m - s * β(j + l)) * dget(i - 1, j)
                    acc += -2 * s * α(1 + j + l) * dget(i - 1, j + 1)
                    Dd[(i, j)] = acc / (j * (1 + j + 2l))
                end
            end
        end
        # slmLambda[i] — uses only orders < i (the d[i,0]·slmLambda[0] term is commented out upstream)
        acc = 0.0 + 0.0im
        for k in 1:(i - 1)
            acc += dget(i - k, 0) * Lam[k + 1]
        end
        acc += α(l - 1) * α(l) * dget(i - 2, -2)
        acc += α(l) * (β(l - 1) + β(l)) * dget(i - 2, -1)
        acc += (-1 + α(l)^2 + α(l + 1)^2 + β(l)^2) * dget(i - 2, 0)
        acc += α(l + 1) * (β(l) + β(l + 1)) * dget(i - 2, 1)
        acc += α(l + 1) * α(l + 2) * dget(i - 2, 2)
        acc += -2 * s * α(l) * dget(i - 1, -1)
        acc += 2 * (m - s * β(l)) * dget(i - 1, 0)
        acc += -2 * s * α(l + 1) * dget(i - 1, 1)
        Lam[i + 1] = -acc / Dd[(0, 0)]
    end
    return Lam, Dd
end

"""
    SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=6) -> Vector{Float64}

Small-γ Taylor coefficients `[c_0,…,c_order]` of the spheroidal eigenvalue about `γ=0`, i.e.
`λ(γ) = Σ_i c_i γ^i + O(γ^{order+1})`. `c_0 = l(l+1)-s(s+1)`. The coefficients are real; the
complex series machinery is verified to carry negligible imaginary part (guard `<1e-9`).
"""
function SpinWeightedSpheroidalEigenvalueSeries(s::Integer, l::Integer, m::Integer; order::Integer=6)
    order >= 0 || throw(ArgumentError("order must be ≥ 0; got $order"))
    Lam, _ = _swsh_series_tables(s, l, m, order)
    maximum(abs ∘ imag, Lam) < 1e-9 ||
        throw(ErrorException("SpinWeightedSpheroidalEigenvalueSeries: unexpected complex coefficient (s=$s,l=$l,m=$m)"))
    return real.(Lam)
end

# ============================== Large-γ asymptotic Series ==============================
# Casals-Ottewill-Warburton, arXiv:1810.00432. Mathematica builds the recurrence RecRelgen
# symbolically, expands it in c=γ about ∞, and Solves order-by-order for δgen[p]/aFgen[k,p]. Each
# Solve is LINEAR in its single unknown, so we reproduce it exactly with truncated Laurent-in-c
# arithmetic + a finite-difference solve (X=-B/A where B=eqn|_{X=0}, A=eqn|_{X=1}-B) — no symbolics.

"`SWSHq[s,l,m]` integer helper for the large-γ asymptotic series."
function _swsh_q(s::Integer, l::Integer, m::Integer)
    slm = abs(m + abs(s)) + abs(s)
    z0  = iseven(l + m) ? 0 : 1
    return l >= slm ? l + 1 - z0 : 2 * l + 1 - slm
end

"""
    SpinWeightedSpheroidalEigenvalueSeriesLargeGamma(s, l, m; order=4) -> Vector{Float64}

Large-γ asymptotic Laurent coefficients of the spheroidal eigenvalue about `γ=∞`
(Casals-Ottewill-Warburton). Returns the coefficient list `cl` of `λ(γ)·γ^order` as an
ascending-power polynomial in γ: `cl[i] =` coefficient of `γ^(i-1)` in `λ(γ)·γ^order`. Length
`order+2`; the `±γ²` terms cancel so the leading entry `cl[end] = 2(q-m)` (may be 0).
"""
function SpinWeightedSpheroidalEigenvalueSeriesLargeGamma(s::Integer, l::Integer, m::Integer; order::Integer=4)
    order >= 0 || throw(ArgumentError("order must be ≥ 0; got $order"))
    q  = _swsh_q(s, l, m)
    ν  = 0.5 * (-1 + q - s - abs(m + s))
    K  = order + 1                       # highest c⁻ power kept
    kmin = -(K)                          # truncate Laurent products below this
    ams, ad = abs(m + s), abs(m - s)

    aF = Dict{Tuple{Int,Int},Float64}()  # aFgen[(r,k)]; seed aFgen[0,0]=1
    aF[(0, 0)] = 1.0
    δ = Dict{Int,Float64}()              # δgen[k], k = 1..K

    function AF(r::Integer)               # AFgen[r] as a Laurent series in c (Dict: power => coeff)
        d = Dict{Int,Float64}()
        if r == 0
            d[0] = 1.0
        else
            for k in abs(r):K
                v = get(aF, (r, k), 0.0)
                v != 0.0 && (d[-k] = get(d, -k, 0.0) + v)
            end
        end
        return d
    end
    addto!(d, p, v) = (v != 0.0 && (d[p] = get(d, p, 0.0) + v))
    function smuls(d, a)                  # scalar * series
        out = Dict{Int,Float64}()
        a != 0.0 && for (p, v) in d; addto!(out, p, a * v); end
        return out
    end
    function smul(d1, d2)                 # series * series, truncated to powers ≥ kmin
        out = Dict{Int,Float64}()
        for (p1, v1) in d1, (p2, v2) in d2
            p = p1 + p2
            p >= kmin && addto!(out, p, v1 * v2)
        end
        return out
    end

    function serngen(nn::Integer, p::Integer)   # coeff of c^(-(p-1)) in RecRelgen(n=nn)
        pref1 = 0.5 * (1 + nn + ν + ams) * (2 * (1 + nn + ν + s) - ad + ams)
        pref3 = 0.5 * (nn + ν) * (2 * (nn + ν + s) + ad + ams)
        Bm = Dict{Int,Float64}()
        addto!(Bm, 1, 4.0 * nn)
        addto!(Bm, 0, -1 - m - m^2 - 2 * (nn + ν) - 2 * (nn + ν)^2 - s - m * s
                       - 2 * (nn + ν) * s - (1 + 2 * nn + 2 * ν + s) * ams)
        for j in 1:K                     # -4 c * Asgen = Σ_j -4 δ[j] c^(1-j)
            addto!(Bm, 1 - j, -4.0 * get(δ, j, 0.0))
        end
        rec = smuls(AF(nn + 1), pref1)
        for (pp, vv) in smul(Bm, AF(nn)); addto!(rec, pp, vv); end
        for (pp, vv) in smuls(AF(nn - 1), pref3); addto!(rec, pp, vv); end
        return get(rec, -(p - 1), 0.0)
    end

    for p in 1:K                         # order-by-order elimination, FD linear solve
        B = serngen(0, p)                # δ[p]
        δ[p] = 1.0; A = serngen(0, p) - B; δ[p] = -B / A
        for k in -p:p
            k == 0 && continue           # aFgen[k,p]
            B = serngen(k, p)
            aF[(k, p)] = 1.0; A = serngen(k, p) - B; aF[(k, p)] = -B / A
        end
    end

    # slmLambda = 2(q-m) γ + C0 + Σ_{k=2}^{K} 4 δ[k] γ^(1-k)  (the ±γ² cancel)
    C0 = -0.5 * (q^2 - m^2 - 2 * s^2 + 1) - s * (s + 1)
    cl = zeros(Float64, order + 2)       # cl[i] = coeff of γ^(i-1) in λ·γ^order
    cl[order + 2] = 2.0 * (q - m)        # γ^1 coeff of λ
    cl[order + 1] = C0                   # γ^0 coeff of λ
    for k in 2:K
        cl[order + 2 - k] = 4.0 * δ[k]
    end
    return cl
end
