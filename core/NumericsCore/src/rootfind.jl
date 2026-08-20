# Vendored 1D root-finders: secant, bisection, Brent, and Alefeld-Potra-Shi A42 bracketing.
#
# Consolidates the scattered derivative-free solvers (Teukolsky/H1 `_secant`, FEW `_brentq`,
# KerrBackwards Brent, the FEW separatrix A42) into one type-generic `{T}` family so BigFloat
# precision is inherited (a hard Float64 floor would cap the resolved root at ~14 digits no matter
# the mantissa). Vendoring A42 lets `Spacetime.Separatrix` drop the external `Roots` dependency.
#
# Every solver promotes its inputs to `T = float(promote_type(...))`, so the default `atol`/`rtol`
# (zero / `eps(T)`) automatically tighten at higher BigFloat precision.

@inline _solver_T(xs...) = float(promote_type(map(typeof, xs)...))

"""
    secant(f, x0, x1; atol=0, rtol=eps, maxiter=100, throw_on_fail=true) -> root

Derivative-free secant iteration from the two seeds `x0`, `x1` (no bracketing required).
Converges superlinearly (order φ≈1.618) for a simple root. Stops when the step
`|xₙ₊₁−xₙ| ≤ atol + rtol·|xₙ₊₁|`. (Teukolsky/H1 `_secant` pattern, generic `{T}`.)

COMPLEX-CAPABLE: `x0`/`x1` may be complex (Leaver-CF eigenvalues, QNM frequencies). The iterate
type is `T = float(promote_type(...))` while the tolerances live in `RT = real(T)` — the step test
compares the REAL magnitudes `|c − b|` and `|c|`, which is what a complex `T(rtol)` could not
express and what forced the two private complex copies (`Harmonics._swsh_secant`,
`Perturbation.QNM._qnm_secant`) that this method replaces.

`throw_on_fail=false` makes the two non-convergence exits RETURN the current iterate `b` instead of
erroring (a zero secant slope, or `maxiter` exhausted) — the behaviour the spectral-seeded
eigenvalue/QNM solvers rely on, where the seed itself is already an acceptable answer and the
caller judges the residual.
"""
function secant(f, x0, x1; atol=nothing, rtol=nothing, maxiter::Integer=100,
                throw_on_fail::Bool=true)
    T = _solver_T(x0, x1)
    RT = real(T)
    a = T(x0); b = T(x1)
    at = atol === nothing ? zero(RT) : RT(atol)
    rt = rtol === nothing ? eps(RT) : RT(rtol)
    (at >= zero(RT) && rt >= zero(RT)) || throw(ArgumentError("secant: tolerances must be nonnegative"))
    fa = f(a); fb = f(b)
    iszero(fa) && return a
    iszero(fb) && return b
    for _ in 1:maxiter
        denom = fb - fa
        if iszero(denom)
            throw_on_fail || return b
            error("secant failed: zero secant slope before convergence")
        end
        c = b - fb * (b - a) / denom
        if abs(c - b) <= at + rt * abs(c)
            return c
        end
        a, fa = b, fb
        b, fb = c, f(c)
    end
    throw_on_fail || return b
    error("secant did not converge in $maxiter iterations (last residual $(abs(fb)))")
end

"""
    bisect(f, a, b; atol=0, maxiter=10000) -> root

Bracketing bisection: requires `f(a)·f(b) ≤ 0`. Each step halves the bracket, so the
half-width `(b−a)/2` is the exact error bound; stops when it falls below `atol` (default a
precision-inherited `eps`-scaled bound). Linear convergence, but unconditionally robust.
"""
function bisect(f, a, b; atol=nothing, maxiter::Integer=10_000)
    T = _solver_T(a, b)
    lo = T(a); hi = T(b)
    lo > hi && ((lo, hi) = (hi, lo))
    at = atol === nothing ? eps(T) * max(one(T), abs(lo), abs(hi)) : T(atol)
    at >= zero(T) || throw(ArgumentError("bisect: atol must be nonnegative"))
    flo = f(lo); fhi = f(hi)
    iszero(flo) && return lo
    iszero(fhi) && return hi
    sign(flo) * sign(fhi) < 0 || throw(ArgumentError("bisect: f(a)=$flo and f(b)=$fhi must bracket a root"))
    c = (lo + hi) / 2
    for _ in 1:maxiter
        c = (lo + hi) / 2
        (hi - lo) / 2 <= at && return c
        fc = f(c)
        iszero(fc) && return c
        if sign(flo) * sign(fc) < 0
            hi, fhi = c, fc
        else
            lo, flo = c, fc
        end
    end
    error("bisect did not converge in $maxiter iterations (last bracket [$lo,$hi])")
end

"""
    brent(f, a, b; atol=0, rtol=eps, maxiter=200) -> root

Brent's method: inverse-quadratic / secant interpolation with a bisection safeguard and the
bracket invariant `f(a)·f(b) < 0` maintained every step. Superlinear in the smooth regime,
bisection-robust in the worst case. Replaces the `Roots` Brent used in FEW `_brentq` /
KerrBackwards. Requires `f(a)·f(b) ≤ 0`.
"""
function brent(f, a, b; atol=nothing, rtol=nothing, maxiter::Integer=200)
    T = _solver_T(a, b)
    aa = T(a); bb = T(b)
    aa > bb && ((aa, bb) = (bb, aa))
    at = atol === nothing ? zero(T) : T(atol)
    rt = rtol === nothing ? eps(T) : T(rtol)
    (at >= zero(T) && rt >= zero(T)) || throw(ArgumentError("brent: tolerances must be nonnegative"))
    fa = f(aa); fb = f(bb)
    iszero(fa) && return aa
    iszero(fb) && return bb
    sign(fa) * sign(fb) < 0 || throw(ArgumentError("brent: f(a)=$fa and f(b)=$fb must bracket a root"))
    if abs(fa) < abs(fb)
        aa, bb = bb, aa; fa, fb = fb, fa
    end
    c = aa; fc = fa
    d = aa
    mflag = true
    for _ in 1:maxiter
        tol = 2 * rt * abs(bb) + at
        if abs(bb - aa) <= tol || iszero(fb)
            return bb
        end
        if fa != fc && fb != fc
            # inverse quadratic interpolation
            s = aa * fb * fc / ((fa - fb) * (fa - fc)) +
                bb * fa * fc / ((fb - fa) * (fb - fc)) +
                c * fa * fb / ((fc - fa) * (fc - fb))
        else
            # secant
            s = bb - fb * (bb - aa) / (fb - fa)
        end
        # safeguard conditions → fall back to bisection
        lo = (3 * aa + bb) / 4
        cond1 = !((lo < s < bb) || (bb < s < lo))
        cond2 = mflag && abs(s - bb) >= abs(bb - c) / 2
        cond3 = !mflag && abs(s - bb) >= abs(c - d) / 2
        cond4 = mflag && abs(bb - c) < tol
        cond5 = !mflag && abs(c - d) < tol
        if cond1 || cond2 || cond3 || cond4 || cond5
            s = (aa + bb) / 2
            mflag = true
        else
            mflag = false
        end
        fs = f(s)
        d = c
        c, fc = bb, fb
        if sign(fa) * sign(fs) < 0
            bb, fb = s, fs
        else
            aa, fa = s, fs
        end
        if abs(fa) < abs(fb)
            aa, bb = bb, aa; fa, fb = fb, fa
        end
    end
    error("brent did not converge in $maxiter iterations (last residual $(abs(fb)))")
end

# --- Alefeld-Potra-Shi A42 helpers ------------------------------------------------------------------
# Newton's method on the quadratic interpolating (a,fa),(b,fb),(d,fd): k inner steps.
@inline function _newton_quadratic(a::T, b::T, d::T, fa::T, fb::T, fd::T, k::Int) where {T}
    fab = (fb - fa) / (b - a)
    fbd = (fd - fb) / (d - b)
    A = (fbd - fab) / (d - a)            # second divided difference f[a,b,d]
    B = fab
    iszero(A) && return a - fa / B
    r = (sign(A) == sign(fa)) ? a : b    # start where A·f(r)>0 so Newton stays in bracket
    @inbounds for _ in 1:k
        P = fa + B * (r - a) + A * (r - a) * (r - b)
        Pp = B + A * (2r - a - b)
        iszero(Pp) && break
        r -= P / Pp
    end
    return r
end

# Inverse cubic (Neville) interpolation through four points — the A42 fast step.
@inline function _ipzero(a::T, b::T, c::T, d::T, fa::T, fb::T, fc::T, fd::T) where {T}
    Q11 = (c - d) * fc / (fd - fc)
    Q21 = (b - c) * fb / (fc - fb)
    Q31 = (a - b) * fa / (fb - fa)
    D21 = (b - c) * fc / (fc - fb)
    D31 = (a - b) * fb / (fb - fa)
    Q22 = (D21 - Q11) * fb / (fd - fb)
    Q32 = (D31 - Q21) * fa / (fc - fa)
    D32 = (D31 - Q21) * fc / (fc - fa)
    Q33 = (D32 - Q22) * fa / (fd - fa)
    return a + (Q31 + Q32 + Q33)
end

# Given a bracket [a,b] (f(a)·f(b)<0) and an interior c, return the new bracket plus the discarded
# endpoint d (and their f-values). The bracket invariant f(a')·f(b')<0 always holds on return.
@inline function _bracket(a::T, b::T, c::T, fa::T, fb::T, fc::T) where {T}
    if sign(fa) * sign(fc) < 0
        return (a, c, b, fa, fc, fb)
    else
        return (c, b, a, fc, fb, fa)
    end
end

"""
    a42(f, a, b; atol=0, rtol=eps, maxiter=100) -> root

Alefeld-Potra-Shi Algorithm 4.2: enhanced bracketing root-finder combining inverse-cubic /
Newton-quadratic interpolation, a double-length secant step, and a bisection safeguard that
forces the bracket width to at least halve every iteration. Maintains `f(a)·f(b) < 0` at every
step and converges with high order (asymptotic ≈ 1.7 per function evaluation). This is the
solver behind FEW's separatrix root; vendoring it characterization-gated lets `Spacetime` drop
`Roots`. Requires `f(a)·f(b) ≤ 0`.
"""
function a42(f, a, b; atol=nothing, rtol=nothing, maxiter::Integer=100)
    T = _solver_T(a, b)
    aa = T(a); bb = T(b)
    at = atol === nothing ? zero(T) : T(atol)
    rt = rtol === nothing ? eps(T) : T(rtol)
    (at >= zero(T) && rt >= zero(T)) || throw(ArgumentError("a42: tolerances must be nonnegative"))
    if aa > bb
        aa, bb = bb, aa
    end
    fa = f(aa); fb = f(bb)
    iszero(fa) && return aa
    iszero(fb) && return bb
    sign(fa) * sign(fb) < 0 || throw(ArgumentError("a42: f(a)=$fa and f(b)=$fb must bracket a root"))

    convtol(lo, hi) = at + 2 * rt * max(abs(lo), abs(hi))
    interior(c, lo, hi) = (lo < c < hi) ? c : (lo + hi) / 2

    # initial secant step → first interior point + bracket
    c = aa - fa * (bb - aa) / (fb - fa)
    c = interior(c, aa, bb)
    fc = f(c)
    iszero(fc) && return c
    aa, bb, d, fa, fb, fd = _bracket(aa, bb, c, fa, fb, fc)
    e, fe = d, fd

    for _ in 1:maxiter
        (bb - aa) <= convtol(aa, bb) && return (aa + bb) / 2
        width0 = bb - aa

        # --- interpolation step: inverse cubic if 4 distinct points, else Newton-quadratic ---
        distinct = (fa != fb != fc) && fa != fd != fe &&
                   (e != aa) && (e != bb) && (d != aa) && (d != bb)
        c = distinct ? _ipzero(aa, bb, d, e, fa, fb, fd, fe) :
            _newton_quadratic(aa, bb, d, fa, fb, fd, 2)
        c = interior(c, aa, bb)
        fc = f(c)
        iszero(fc) && return c
        e, fe = d, fd
        aa, bb, d, fa, fb, fd = _bracket(aa, bb, c, fa, fb, fc)

        # --- double-length secant step from the smaller-|f| endpoint ---
        if abs(fa) < abs(fb)
            u, fu = aa, fa
        else
            u, fu = bb, fb
        end
        c = u - 2 * fu * (bb - aa) / (fb - fa)
        if abs(c - u) > (bb - aa) / 2
            c = (aa + bb) / 2
        end
        c = interior(c, aa, bb)
        fc = f(c)
        iszero(fc) && return c
        e, fe = d, fd
        aa, bb, d, fa, fb, fd = _bracket(aa, bb, c, fa, fb, fc)

        # --- bisection safeguard: force the width to at least halve ---
        if (bb - aa) > width0 / 2
            m = (aa + bb) / 2
            fm = f(m)
            iszero(fm) && return m
            e, fe = d, fd
            aa, bb, d, fa, fb, fd = _bracket(aa, bb, m, fa, fb, fm)
        end
    end
    error("a42 did not converge in $maxiter iterations (last bracket [$aa,$bb])")
end
