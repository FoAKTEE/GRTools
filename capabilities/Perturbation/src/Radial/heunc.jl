# Confluent-Heun radial backend (:HeunC) for the homogeneous Teukolsky solution. Leaf special-
# function port — NO ODE/MST dependency: a z=0 Frobenius series + an analytic-continuation marcher
# (Taylor re-expansion about regular interior points), wrapped by the two BHPToolkit
# TeukolskyRadialHeunC radial prefactor expressions. HeunC(q,α,γ,δ,ε,z) is the local Frobenius
# solution at z=0 (=1 there) of  w'' + (γ/z + δ/(z-1) + ε) w' + (α z - q)/(z(z-1)) w = 0.
#
# The :HeunC branch is wired into the carrier (R(r) callable), the TeukolskyRadial dispatcher, and
# teukolsky_radial_deriv. The In/Up VALUES reuse the same unit-transmission amplitudes as :MST (same
# physical solution); the DERIVATIVE is the genuine confluent-Heun derivative (analytic-continued
# (w,w') + product rule on the raw prefactor), not an MST fall-through. Type-generic {T}: every
# literal base is an Integer so BigFloat precision is inherited (2^Complex{BigFloat} → Complex{BigFloat}).

"""
    _heunc_series(q, α, γ, δ, ε, z; Nmax=2000, tol=1e-16) -> Complex

z=0 Frobenius series Σ_{n≥0} c_n z^n of HeunC, normalized c_0=1. Three-term recurrence
c_1=-q/γ and (k+1)(k+γ) c_{k+1} = [k(k-1+γ+δ-ε)-q] c_k + [α+ε(k-1)] c_{k-1}. Converges for
|z|<1 (distance to the z=1 singularity).
"""
function _heunc_series(q, α, γ, δ, ε, z; Nmax::Integer=2000, tol::Real=1e-16)
    T = complex(float(promote_type(typeof(q), typeof(α), typeof(γ), typeof(δ), typeof(ε), typeof(z))))
    zc = T(z)
    cprev = one(T)                       # c_0
    w = cprev                            # running w = Σ c_n z^n
    zpow = one(T)
    c = -T(q) / T(γ)                     # c_1
    zpow *= zc
    w += c * zpow
    ckm1 = cprev; ck = c
    for k in 1:Nmax
        ckp1 = ((k * (k - 1 + γ + δ - ε) - q) * ck + (α + ε * (k - 1)) * ckm1) /
               ((k + 1) * (k + γ))
        zpow *= zc
        term = ckp1 * zpow
        w += term
        (abs(term) <= tol * (abs(w) + tol)) && (k > 2) && return w
        ckm1 = ck; ck = ckp1
    end
    return w
end

"""
    _heunc_series_deriv(q, α, γ, δ, ε, z; Nmax=2000, tol=1e-16) -> Complex

w'(z) from the z=0 series (term-by-term derivative), same recurrence as `_heunc_series`. Seeds
the analytic-continuation marcher at an interior point.
"""
function _heunc_series_deriv(q, α, γ, δ, ε, z; Nmax::Integer=2000, tol::Real=1e-16)
    T = complex(float(promote_type(typeof(q), typeof(α), typeof(γ), typeof(δ), typeof(ε), typeof(z))))
    zc = T(z)
    dw = zero(T)
    c1 = -T(q) / T(γ)
    dw += c1
    ckm1 = one(T); ck = c1
    zpow = one(T)
    for k in 1:Nmax
        ckp1 = ((k * (k - 1 + γ + δ - ε) - q) * ck + (α + ε * (k - 1)) * ckm1) /
               ((k + 1) * (k + γ))
        zpow *= zc
        term = (k + 1) * ckp1 * zpow
        dw += term
        (abs(term) <= tol * (abs(dw) + tol)) && k > 2 && return dw
        ckm1 = ck; ck = ckp1
    end
    return dw
end

"""
    _heunc_local(q, α, γ, δ, ε, z0, w0, dw0, z; Nmax=400, tol=1e-16) -> (w, dw)

Taylor re-expansion of the HeunC ODE solution about a REGULAR point z0 (z0∉{0,1}), given
(w0,dw0)=(w(z0),w'(z0)). Returns (w(z),w'(z)) at nearby z. The polynomial form
z(z-1)w'' + [γ(z-1)+δz+ε z(z-1)]w' + (αz-q)w = 0 gives an exact recurrence for the local Taylor
coefficients. Valid for |z-z0| < dist(z0,{0,1}).
"""
function _heunc_local(q, α, γ, δ, ε, z0::T, w0::T, dw0::T, z::T;
                      Nmax::Integer=400, tol::Real=1e-16) where {T<:Complex}
    h = z - z0
    p2_0 = z0 * (z0 - 1);            p2_1 = 2 * z0 - 1;                    p2_2 = one(T)
    p1_0 = γ * (z0 - 1) + δ * z0 + ε * z0 * (z0 - 1)
    p1_1 = γ + δ + ε * (2 * z0 - 1); p1_2 = ε
    p0_0 = α * z0 - q;               p0_1 = α
    a = Vector{T}(undef, Nmax + 3)
    a[1] = w0                                   # a_0
    a[2] = dw0                                  # a_1
    w = w0 + dw0 * h
    dw = dw0
    hpow = h
    for n in 0:Nmax
        an   = a[n + 1]
        anm1 = n >= 1 ? a[n] : zero(T)
        anp1 = a[n + 2]
        rhs = p2_1 * (n + 1) * n * anp1 + p2_2 * n * (n - 1) * an +
              p1_0 * (n + 1) * anp1 + p1_1 * n * an + p1_2 * (n - 1) * anm1 +
              p0_0 * an + p0_1 * anm1
        anp2 = -rhs / (p2_0 * (n + 2) * (n + 1))
        a[n + 3] = anp2
        hpow *= h                                # h^{n+2}
        wterm = anp2 * hpow
        w += wterm
        dw += (n + 2) * anp2 * (hpow / h)        # (n+2) a_{n+2} h^{n+1}
        (abs(wterm) <= tol * (abs(w) + tol)) && n > 2 && break
    end
    return (w, dw)
end

"""
    _heunc(q, α, γ, δ, ε, z; rfrac=0.5) -> Complex

Confluent Heun function HeunC(q,α,γ,δ,ε,z) = local solution at z=0 (=1 there). For |z|≤0.5 use the
z=0 series directly; otherwise analytically continue along the straight segment [0,z] by Taylor
re-expansion, stepping ≤ rfrac × (distance to the nearest singularity z∈{0,1}).
"""
function _heunc(q, α, γ, δ, ε, z; rfrac::Real=0.5)
    return _heunc_wderiv(q, α, γ, δ, ε, z; rfrac=rfrac)[1]
end

"""
    _heunc_wderiv(q, α, γ, δ, ε, z; rfrac=0.5) -> (w, dw)

HeunC value and its z-derivative together. Same series/continuation path as `_heunc`; the marcher
already carries (w,w'), so returning both is free. Used by the HeunC radial derivative.
"""
function _heunc_wderiv(q, α, γ, δ, ε, z; rfrac::Real=0.5)
    T = complex(float(promote_type(typeof(q), typeof(α), typeof(γ), typeof(δ), typeof(ε), typeof(z))))
    qc, αc, γc, δc, εc, zc = T(q), T(α), T(γ), T(δ), T(ε), T(z)
    if abs(zc) <= 0.5
        return (_heunc_series(qc, αc, γc, δc, εc, zc), _heunc_series_deriv(qc, αc, γc, δc, εc, zc))
    end
    dir = zc / abs(zc)                           # unit vector toward target along the ray
    z0 = T(35 // 100) * dir                       # seed point inside |·|<1
    w = _heunc_series(qc, αc, γc, δc, εc, z0)
    dw = _heunc_series_deriv(qc, αc, γc, δc, εc, z0)
    target = abs(zc)
    pos = abs(z0)
    zj = z0
    maxsteps = 100000; step = 0
    while pos < target - eps(real(T)) && step < maxsteps
        step += 1
        rad = min(abs(zj), abs(zj - 1))          # distance to nearest singularity
        h = rfrac * rad
        newpos = min(pos + h, target)
        znext = newpos * dir
        (w, dw) = _heunc_local(qc, αc, γc, δc, εc, zj, w, dw, znext)
        zj = znext
        pos = newpos
    end
    return (w, dw)
end

# --- Teukolsky HeunC radial solutions (BHPToolkit Teukolsky TeukolskyRadialHeunC) ---
# The In solution = prefactor(powers·exp) × HeunC[6 params from s,m,a,ω,λ at z=(1+κ-r)/(2κ)]. The Up
# solution = Reflection·In_raw + Incidence·(a conjugate HeunC with ωc=conj(ω), λc=conj(λ), -s, and a
# (r²-2r+a²)^{-s} factor), amplitudes already rescaled to unit transmission. ω real here (ωc=ω) but
# kept general per the source.

"""
    _heunc_in_params(s, m, a, omega, lambda, r) -> NamedTuple

Shared building blocks of the raw HeunC "In" radial solution: prefactor `pref`, HeunC parameters
(q,α,γ,δ,ε), argument `z`, its r-derivative `zp`, and the logarithmic r-derivative `ld` of the
prefactor. `pref * HeunC(q,α,γ,δ,ε,z)` is `solFuncs["In"]`; the raw value and derivative both build
from this pack so their constants can never drift.
"""
function _heunc_in_params(s::Integer, m::Integer, a::Real, omega, lambda, r::Real)
    am = a * m; w = omega
    kappa = sqrt(complex(1 - a^2))
    z = (1 + kappa - r) / (2 * kappa)
    zp = -1 / (2 * kappa)                        # dz/dr
    C = (im * am) / (2 * kappa) - s - im * (1 + 1 / kappa) * w   # exponent of (r-1-κ)/κ
    D = (im * (am + 2 * (-1 + kappa) * w)) / (2 * kappa)         # exponent of (r-1+κ)/κ
    pref = 2^((-im * am - kappa * s + 2im * w) / kappa) *
           complex(1 - a^2)^((im * am) / (2 * (1 + kappa)) - s - im * w) *
           exp(0.5im * (am - 2 * r * w)) *
           ((-1 - kappa + r) / kappa)^C *
           ((-1 + kappa + r) / kappa)^D
    q = s + s^2 + lambda - (am - 2 * w)^2 / (-1 + a^2) - 4 * w^2 +
        (im * (-am - 4 * (-1 + s) * w + 2 * a^2 * (-1 + 2 * s) * w)) / kappa
    α = -4 * w * (-im * kappa + am + im * kappa * s - 2 * w + 2 * kappa * w)
    γ = 1 - s + (im * (am - 2 * w)) / kappa - 2im * w
    δ = 1 + s + (im * (am - 2 * w)) / kappa + 2im * w
    ε = 4im * kappa * w
    ld = -im * w + C / (r - 1 - kappa) + D / (r - 1 + kappa)     # d/dr log(pref)
    return (pref=pref, q=q, α=α, γ=γ, δ=δ, ε=ε, z=z, zp=zp, ld=ld)
end

"""
    _heunc_up_conj_params(s, m, a, omega, lambda, r) -> NamedTuple

Shared building blocks of the second (conjugate) HeunC term of the "Up" solution
(`solFuncs["Up"]`): the (r²-2r+a²)^{-s} prefactor × a HeunC built with ωc=conj(ω), λc=conj(λ),
spin -s, and q using (λc+2s). Same fields as [`_heunc_in_params`](@ref).
"""
function _heunc_up_conj_params(s::Integer, m::Integer, a::Real, omega, lambda, r::Real)
    am = a * m; wc = conj(omega); lamc = conj(lambda)
    kappa = sqrt(complex(1 - a^2))
    z = (1 + kappa - r) / (2 * kappa)
    zp = -1 / (2 * kappa)
    Δ = r^2 - 2 * r + a^2
    C = (-im * am) / (2 * kappa) - (-s) + im * (1 + 1 / kappa) * wc
    D = (-im * (am + 2 * (-1 + kappa) * wc)) / (2 * kappa)
    pref = complex(Δ)^(-s) *
           2^((im * am - kappa * (-s) - 2im * wc) / kappa) *
           complex(1 - a^2)^((-im * am) / (2 * (1 + kappa)) - (-s) + im * wc) *
           exp(-0.5im * (am - 2 * r * wc)) *
           ((-1 - kappa + r) / kappa)^C *
           ((-1 + kappa + r) / kappa)^D
    q = (-s) + (-s)^2 + (lamc + 2 * s) - (am - 2 * wc)^2 / (-1 + a^2) - 4 * wc^2 +
        (-im * (-am - 4 * (-1 + (-s)) * wc + 2 * a^2 * (-1 + 2 * (-s)) * wc)) / kappa
    α = -4 * wc * (im * kappa + am - im * kappa * (-s) - 2 * wc + 2 * kappa * wc)
    γ = 1 - (-s) + (-im * (am - 2 * wc)) / kappa + 2im * wc
    δ = 1 + (-s) + (-im * (am - 2 * wc)) / kappa - 2im * wc
    ε = -4im * kappa * wc
    ld = -s * (2 * r - 2) / Δ + im * wc + C / (r - 1 - kappa) + D / (r - 1 + kappa)
    return (pref=pref, q=q, α=α, γ=γ, δ=δ, ε=ε, z=z, zp=zp, ld=ld)
end

"Raw HeunC \"In\" radial value (prefactor × HeunC), BHPToolkit `solFuncs[\"In\"]`."
function _heunc_radial_in_raw(s::Integer, m::Integer, a::Real, omega, lambda, r::Real)
    p = _heunc_in_params(s, m, a, omega, lambda, r)
    return p.pref * _heunc(p.q, p.α, p.γ, p.δ, p.ε, p.z)
end

"Second (conjugate) HeunC term of the \"Up\" radial value, BHPToolkit `solFuncs[\"Up\"]`."
function _heunc_radial_up_conj_raw(s::Integer, m::Integer, a::Real, omega, lambda, r::Real)
    p = _heunc_up_conj_params(s, m, a, omega, lambda, r)
    return p.pref * _heunc(p.q, p.α, p.γ, p.δ, p.ε, p.z)
end

"Raw HeunC \"In\" radial derivative dR/dr = pref·(w·ld + w'·z'), analytic (no FD, no MST)."
function _heunc_radial_in_raw_deriv(s::Integer, m::Integer, a::Real, omega, lambda, r::Real)
    p = _heunc_in_params(s, m, a, omega, lambda, r)
    w, dw = _heunc_wderiv(p.q, p.α, p.γ, p.δ, p.ε, p.z)
    return p.pref * (w * p.ld + dw * p.zp)
end

"Raw \"Up\"-conjugate HeunC term derivative dR/dr, analytic (product rule on the raw prefactor)."
function _heunc_radial_up_conj_raw_deriv(s::Integer, m::Integer, a::Real, omega, lambda, r::Real)
    p = _heunc_up_conj_params(s, m, a, omega, lambda, r)
    w, dw = _heunc_wderiv(p.q, p.α, p.γ, p.δ, p.ε, p.z)
    return p.pref * (w * p.ld + dw * p.zp)
end

"""
    _heunc_radial_value(R, r) -> Complex

Normalized HeunC radial value. In = raw HeunC In (already unit-transmission). Up =
`Reflection·In_raw + Incidence·conj_raw`, with Reflection/Incidence the Up amplitudes rescaled to
unit transmission (= `R.amplitudes`, shared with the :MST backend).
"""
function _heunc_radial_value(R::TeukolskyRadialFunction, r::Real)
    inraw = _heunc_radial_in_raw(R.s, R.m, R.a, R.omega, R.lambda, r)
    R.bc === :In && return inraw
    conjraw = _heunc_radial_up_conj_raw(R.s, R.m, R.a, R.omega, R.lambda, r)
    return R.amplitudes.Reflection * inraw + R.amplitudes.Incidence * conjraw
end

"""
    _heunc_radial_deriv(R, r) -> Complex

Normalized HeunC radial first derivative dR/dr, from the genuine confluent-Heun derivative (In: raw
In deriv; Up: `Reflection·In_deriv + Incidence·conj_deriv`). Same unit-transmission amplitudes as
the value — the two are the same physical solution as :MST, so this agrees with the MST derivative
to the shared conditioning tolerance.
"""
function _heunc_radial_deriv(R::TeukolskyRadialFunction, r::Real)
    inderiv = _heunc_radial_in_raw_deriv(R.s, R.m, R.a, R.omega, R.lambda, r)
    R.bc === :In && return inderiv
    conjderiv = _heunc_radial_up_conj_raw_deriv(R.s, R.m, R.a, R.omega, R.lambda, r)
    return R.amplitudes.Reflection * inderiv + R.amplitudes.Incidence * conjderiv
end
