# =============================================================================================================
# RADIAL_KERR — a≠0 (spin-general) FULL RADIAL RECONSTRUCTION: P_{-2}(r),P_{-1}(r),h(r),κ(r) at spin a.
#
# Generalizes radial.jl from a=0 to Kerr.  The ONLY structural changes vs the a=0 path are:
#   1. spin-carrying geometry: K=ω(r²+a²)−am, Δ=r²−2r+a², c=iK/Δ (the s=−2/−1/0 ODEs and their R'' inherit this);
#   2. spheroidal eigenvalues Λ,λ,λ0 = SpinWeightedSpheroidalEigenvalue(s,l,m,aω) (NOT the integer a=0 forms);
#   3. the a≠0 jumps {J0,J1}(P_{-1}),{J0,J1}(κ) come from the BANDED closure_kerr (App-H a=0.99 lower 6.85e-8);
#   4. the κ source is the spin-general −½(r²+a²γ_ll)/Δ·h ⇒ inhomogeneous ODE  L[κ]=½(r²+a²γ_ll)·h  (a=0 ⇒ ½r²h);
#      L = 𝒟Δ𝒟†−2iωr−λ0, the same s=0 operator.  γ_ll is the DRESSED diagonal (b₀(c²)⁰⁰b₀ᵀ)_{ll} from closure_kerr.
#
# The homogeneous solutions are TeukolskyRadial(s,l,m,a,ω) (a≠0 R^up/R^in, MST/Monodromy — BigFloat-inheritable);
# derivatives via teukolsky_radial_deriv; 2nd derivatives from the homogeneous ODE itself (FD-free).  The κ VoP +
# Abel-summed tails are identical in form to a=0, with the horizon at rh=1+√(1−a²) and the spin-general source.
#
# GATE (a): the reconstructed P_{-2},P_{-1},h,κ satisfy their spin-general ODEs (residual <1e-8 for the homogeneous
# ones; <1e-6 for the inhomogeneous κ) AND reproduce the closure_kerr jumps ([f]=J0,[f']=J1) to <1e-4.
# =============================================================================================================

# Teukolsky homogeneous radial surface + the low-level MST internals, REUSED from the GrTools
# Perturbation capability (the migrated, verified Teukolsky port). _d2R_ode/_mst_*/_hu_*/_prefac_up
# are flat unexported Perturbation internals; pochhammer is owned by BigSpecial (used by the
# high-magnitude node sampler to build the r-INDEPENDENT MST f_n array ONCE per leg).
import Perturbation: TeukolskyRadial, teukolsky_radial_deriv,
                     _d2R_ode, _mst_fn_array, _mst_sum, _prefac_up, _hu_array, _dhu_array
import BigSpecial: pochhammer
import NumericsCore: bary_weights, bary_eval
import QuadGK: quadgk
import Printf: @printf
import Base.Threads

export reconstruct_kerr, HomFieldK, KappaFieldK, valK, derK, der2K, κvalK, κderK, radial_gates_kerr, jump_data_kerr

# spin-general c=iK/Δ and c'(r)
cfunK(r,a,ω,m)  = (Δ=r^2-2r+a^2; K=ω*(r^2+a^2)-a*m; im*K/Δ)
cpfunK(r,a,ω,m) = (Δ=r^2-2r+a^2; K=ω*(r^2+a^2)-a*m; Kp=2ω*r; Δp=2r-2; im*(Kp*Δ-K*Δp)/Δ^2)

# NB: the spin-general homogeneous 2nd derivative R'' is taken from the package's EXACT Teukolsky radial ODE
# `Teukolsky._d2R_ode` (valid for ALL spins incl. s=+2, package-tested) inside `branchK` below — NOT a hand-coded
# GHP form. cfunK/cpfunK above are still used by the residual checks (resK_*) which verify the GHP-form ODEs.

# R^up Float64 precision ceiling: the MST sum's terms are each ~|R^up|, so a Float64 eval keeps only
# 16−log10|R^up| relative digits — corrupt above |R^up|~_BIGMAG. Fix: build R^up/R^up'/R^up'' in BigFloat
# at ~26 Chebyshev nodes and serve Float64 barycentric interpolants (rel ~8e-15); triggered DYNAMICALLY
# on the span-max magnitude, not a fixed l. Full diagnosis: the F5 research note (progress/ dossiers).
const _BIGMAG = 1e8          # |R^up| above which Float64 loses ≳8 relative digits ⇒ use the BigFloat-interp branch
const _NCHEB  = parse(Int, get(ENV,"KLC_NCHEB_FLOOR","30"))  # Cheb nodes/leg floor (env-override; default 30 = F8-byte-identical)
# bits sized to the |R^up| magnitude + working headroom; accuracy is N-limited, not bits-limited (measured).
_bigprec_for_mag(magpow) = clamp(128 + 6*ceil(Int, magpow), 160, 320)   # magpow 13→206→bits 224; 33→bits 320

# Chebyshev–Gauss nodes on [lo,hi] (open; avoids endpoints), SORTED ascending. NOT NumericsCore.cheb_nodes
# (which is Gauss–LOBATTO, endpoints included) — the interpolant panels need the open Gauss abscissae here.
# The barycentric weights/eval feeding these nodes ARE the canonical NumericsCore.bary_weights/bary_eval.
_cheb_nodes(lo, hi, N) = sort!([0.5*(lo+hi) + 0.5*(hi-lo)*cospi((2k-1)/(2N)) for k in 1:N])

# Log-domain Chebyshev-barycentric interpolant for a smooth NONZERO R(r): |R^up| spans 5–8 decades, so
# interpolate L=log R (phase-unwrapped), Lp=R'/R, Lpp=R''/R−(R'/R)² and reconstruct R=e^L, R'=R·Lp,
# R''=R(Lpp+Lp²) — relative precision holds across the dynamic range (a linear-domain panel fails at ~1e-1).
# `outer` = exact evaluator for r outside [lo,hi] (the κ Abel tail; exponentially damped there).
struct LogChebInterp{F}
    lo::Float64; hi::Float64
    x::Vector{Float64}; w::Vector{Float64}
    L::Vector{ComplexF64}; Lp::Vector{ComplexF64}; Lpp::Vector{ComplexF64}
    which::Symbol     # :val | :der | :d2  (selects which reconstruction this closure returns)
    outer::F
end
@inline function (ci::LogChebInterp)(r::Real)
    rf = Float64(r)
    if !(ci.lo - 1e-9 <= rf <= ci.hi + 1e-9)
        return ci.outer(rf)
    end
    Lr = bary_eval(ci.x, ci.w, ci.L, rf)
    ci.which === :val && return exp(Lr)
    Lpr = bary_eval(ci.x, ci.w, ci.Lp, rf)
    Rv = exp(Lr)
    ci.which === :der && return Rv * Lpr
    Lppr = bary_eval(ci.x, ci.w, ci.Lpp, rf)
    return Rv * (Lppr + Lpr*Lpr)     # :d2
end

# Build a spin-general homogeneous solution branch (value/deriv/2nd-deriv) for spin s. Pass eig=nothing to let
# TeukolskyRadial resolve its own SWSH eigenvalue (needed for s=+2, where the s=−2 Λ is the wrong radial eigenvalue
# and MST may not converge if forced); the R'' ODE then uses the resolved λ (read back from R.Up.lambda).
#
# `rspan=(rlo,rhi)` is the radial domain the branch will be sampled over (the mode grid range). When the Float64
# |R^up(rhi)| would exceed _BIGMAG, the up/inn legs are replaced by BigFloat-built Chebyshev interpolants (FAST +
# full relative precision); otherwise the native Float64 MST closures are used unchanged.
function branchK(s, l, m, a, ω, eig; rspan=nothing, r0=nothing)
    T = typeof(float(real(a)))
    R = eig === nothing ? TeukolskyRadial(s, l, m, a, ω) : TeukolskyRadial(s, l, m, a, ω; lambda=eig)
    # 2nd derivative from the package's EXACT Teukolsky radial ODE (_d2R_ode), valid for ALL spins, using the
    # resolved SWSH spheroidal eigenvalue R.Up.lambda (NOT the GHP Λ — the package convention).
    λteuk = real(R.Up.lambda)
    d2(Rf,r) = _d2R_ode(Rf(r), teukolsky_radial_deriv(Rf,r), s, m, a, r, ω, λteuk)
    up(r)  = R.Up(r);  upd(r)=teukolsky_radial_deriv(R.Up,r); updd(r)=d2(R.Up,r)
    inn(r) = R.In(r);  ind(r)=teukolsky_radial_deriv(R.In,r); indd(r)=d2(R.In,r)
    nat = (R=R, up=up, upd=upd, updd=updd, inn=inn, ind=ind, indd=indd, eig=(eig===nothing ? λteuk : eig), big_at=nothing)
    # BigFloat-type fast path: native-BigFloat MST evals are seconds each — replace the Up leg with log-r
    # BigFloat barycentric panels (sampled once); the In leg is O(1)-conditioned so Float64-promote suffices
    # below l=10. Panel spans: s=0 must reach the κ Abel-tail cap; s≠0 covers the grid+pad only (a too-wide
    # panel loses accuracy). `big_at` serves all four legs in BigFloat for the cancellation-prone jump solve.
    if (T !== Float64) && rspan !== nothing && (s == 0 || s == -2 || s == -1 || s == 2)
        rh0 = one(T) + sqrt(one(T) - real(a)^2)
        rlo_b = max(Float64(rh0) + 0.05, Float64(real(rspan[1])) - 0.5)
        hi_b = s == 0 ? _KVOP_SPAN_HI : max(Float64(real(rspan[2])) + 0.5, rlo_b + 1.0)
        upI_b, updI_b = _branchK_up_interp_big(s, l, m, BigFloat(real(a)), BigFloat(real(ω)), R.Up.lambda,
                                               R.Up.norm, R.Up.nu, rlo_b, hi_b, up, upd)
        d2b(r) = _d2R_ode(upI_b(r), updI_b(r), s, m, a, BigFloat(real(r)), ω, real(R.Up.lambda))
        # In leg: Float64 fast path with a lazy BigLogRInterp fallback where Float64 |R^in| overflows
        # (high ℓ, r<r0); byte-identical where Float64 is finite.
        Rf = TeukolskyRadial(s, l, m, Float64(real(a)), Float64(real(ω)); lambda=Float64(real(R.Up.lambda)))
        λf = real(Rf.Up.lambda)
        _isfinC(z) = isfinite(real(z)) && isfinite(imag(z))
        # In-panel span constraints: UPPER edge anchored to the PHYSICAL reach Rout_hi=r0+11 (+pad), not the
        # grid — grid-independent panel (else interp error is amplified ×1e3-1e5 by the LHS source cancellation);
        # LOWER edge at rh+0.05 so the whole κ inner Abel integral stays on the cheap panel (out-of-span native
        # BigFloat MST in a quadgk loop is the 30+min/mode blocker). Up-panel edges are separate and tighter.
        in_hi_b = max(Float64(real(r0 === nothing ? rspan[2] : r0)) + 11.5, rlo_b + 1.0)
        in_lo_b = min(rlo_b, Float64(rh0) + 0.05)
        _inI = Ref{Any}(nothing); _indI = Ref{Any}(nothing)
        # N_in = 3ℓ (floor 18): measured spectral law — N=18 leaves 1e-4..1e-3 inter-node error at ℓ≈50,
        # amplified ×1e4-1e5 by the tetrad-sector cancellations; 3ℓ holds ≤1e-12 through ℓ≈56.
        _N_in = max(_NCHEB_IN, parse(Int, get(ENV, "KLC_NIN_FACTOR", "3"))*l)   # env knob (default 3 = production)
        # √-span node bump (capped ×1.6): the rh-extension region is smooth so needs fewer nodes/unit;
        # reference span anchored to r0 so _N_in is grid-invariant.
        _span_ref = log(in_hi_b) - log(r0 === nothing ? max(rlo_b, Float64(rh0)+0.05) : max(Float64(real(r0))-3.0, Float64(rh0)+0.05))
        _span_act = log(in_hi_b) - log(in_lo_b)
        _N_in = max(_N_in, ceil(Int, _N_in * min(1.6, sqrt(_span_ref > 0 ? _span_act/_span_ref : 1.0))))
        _ensure_in_interp() = begin
            if _inI[] === nothing
                vI, dI = _branchK_leg_interp_big(R.In, in_lo_b, in_hi_b; N=_N_in)  # log-r Cheb panel of R.In (val,der)
                _inI[] = vI; _indI[] = dI
            end
        end
        # l≥10 ⇒ always the BigFloat panel: the Float64 In MST is relative-precision-corrupt from l≳10,
        # long before it overflows (finite-but-wrong values would pass a finiteness gate).
        _in_force_interp = l >= 10
        innF(r) = begin
            if !_in_force_interp
                v = ComplexF64(Rf.In(Float64(real(r))))
                _isfinC(v) && return Complex{BigFloat}(v)
            end
            _ensure_in_interp(); _inI[](BigFloat(real(r)))
        end
        indF(r) = begin
            if !_in_force_interp
                v = ComplexF64(teukolsky_radial_deriv(Rf.In, Float64(real(r))))
                _isfinC(v) && return Complex{BigFloat}(v)
            end
            _ensure_in_interp(); _indI[](BigFloat(real(r)))
        end
        inddF(r) = begin
            rf = Float64(real(r))
            if !_in_force_interp
                vF = ComplexF64(Rf.In(rf)); dF = ComplexF64(teukolsky_radial_deriv(Rf.In, rf))
                if _isfinC(vF) && _isfinC(dF)
                    return Complex{BigFloat}(_d2R_ode(vF, dF, s, m, Float64(real(a)), rf, Float64(real(ω)), λf))
                end
            end
            rB = BigFloat(real(r)); _ensure_in_interp()
            _d2R_ode(_inI[](rB), _indI[](rB), s, m, a, rB, ω, real(R.Up.lambda))
        end
        # s=0: big_at=nothing (the one-point _big_legs_at fallback covers it); s≠0: four legs from the SAME R.
        big_atI = s == 0 ? nothing : function (r)
            rB = BigFloat(real(r))
            (Complex{BigFloat}(R.Up(rB)), Complex{BigFloat}(teukolsky_radial_deriv(R.Up, rB)),
             Complex{BigFloat}(R.In(rB)), Complex{BigFloat}(teukolsky_radial_deriv(R.In, rB)))
        end
        return (R=R, up=upI_b, upd=updI_b, updd=d2b, inn=innF, ind=indF, indd=inddF,
                eig=(eig===nothing ? λteuk : eig), big_at=big_atI)
    end
    # Float64 dynamic magnitude gate.
    (T === Float64 && rspan !== nothing) || return nat
    rlo = Float64(rspan[1]); rhi = Float64(rspan[2])
    # gate on the MAX over BOTH endpoints (|R^up| can peak at either end depending on s, l).
    magpow = log10(max(abs(R.Up(rlo)), abs(R.Up(rhi)), 1.0))
    magpow < log10(_BIGMAG) && return nat   # Float64 keeps enough relative digits ⇒ native path
    # Only the Up leg is precision-poisoned; the In leg is O(1)-conditioned and stays native Float64.
    # s=0: native-Float64 out-of-span fallback (the κ Abel tail must stay off per-point BigFloat) and the
    # resolved λ passed for VoP consistency; s≠0: BigFloat fallback + self-resolved λ (no Abel tail).
    outer_native = s == 0 ? (up, upd, updd) : nothing
    eig_interp   = s == 0 ? λteuk : nothing
    upI, updI, upddI, big_atI = _branchK_up_interp(s, l, m, Float64(a), Float64(ω), rlo, rhi, magpow;
                                                   outer_native=outer_native, eig=eig_interp)
    (R=R, up=upI, upd=updI, updd=upddI, inn=inn, ind=ind, indd=indd, eig=(eig===nothing ? λteuk : eig), big_at=big_atI)
end

# Shared-f_n MST evaluator: Perturbation's _mst_radial_up formulas verbatim, but with a precomputed
# f_n array and value+derivative in ONE pass (the HU array is the dominant node cost). norm=1.
function _up_val_der(farr, s, m, eps, kappa, tau, nu, zhat, Nmax, atol, rtol)
    HU  = _hu_array(s, nu, tau, eps, zhat, Nmax)
    dHU = _dhu_array(s, nu, tau, eps, zhat, HU, Nmax)
    a1 = nu + 1 + s - im*eps; a2 = nu + 1 - s + im*eps
    prefac = _prefac_up(s, eps, kappa, tau, nu, zhat)
    dprefac = prefac * (im + (nu + im*(eps + tau)/2)/zhat + (-im*(eps + tau)/2 - s)/(zhat - eps*kappa))
    fUp = n -> (iseven(n) ? 1 : -1) * pochhammer(a1, n) / pochhammer(a2, n) * farr[n]
    v = _mst_sum(n -> prefac * fUp(n) * HU[n], Nmax, atol, rtol)
    d = (eps/2) * _mst_sum(n -> fUp(n) * (dprefac * HU[n] + prefac * dHU[n]), Nmax, atol, rtol)
    (v, d)
end

# Up-leg high-magnitude path: one BigFloat R + one f_n array, sample L/Lp/Lpp at _NCHEB Chebyshev
# nodes over [rlo,rhi], phase-unwrap, return Float64 log-domain interpolant closures + BigFloat fallback.
function _branchK_up_interp(s, l, m, a::Float64, ω::Float64, rlo::Float64, rhi::Float64, magpow::Real;
                            outer_native=nothing, eig=nothing)
    bits = _bigprec_for_mag(magpow)
    ab = BigFloat(a; precision=bits); ωb = BigFloat(ω; precision=bits)
    R = setprecision(BigFloat, bits) do
        eig === nothing ? TeukolskyRadial(s, l, m, ab, ωb) :
                          TeukolskyRadial(s, l, m, ab, ωb; lambda=BigFloat(Float64(eig); precision=bits))
    end
    λteuk = R.Up.lambda                                   # BigFloat λ (keep full precision for the ODE 2nd-deriv)
    nu = R.Up.nu; norm_up = R.Up.norm
    rh = 1.0 + sqrt(1.0 - a^2)
    lo = max(rlo, rh + 0.05); hi = rhi
    Nmax = 80; atol = BigFloat(1e-30; precision=bits); rtol = BigFloat(1e-30; precision=bits)
    eps = 2ωb; kappa = sqrt(1 - ab^2); tau = (eps - m*ab)/kappa; rm = 1 - kappa
    farr_up = setprecision(BigFloat, bits) do; first(_mst_fn_array(s, m, ab, ωb, nu, λteuk, Nmax)); end
    λr = real(λteuk)
    vd_at(rB) = (zhat = eps*rB/2 - eps*rm/2;
                 vd = _up_val_der(farr_up, s, m, eps, kappa, tau, nu, zhat, Nmax, atol, rtol);
                 (vd[1]/norm_up, vd[2]/norm_up))
    # exact BigFloat point-evaluator (out-of-span fallback): value/deriv/2nd-deriv at one r.
    function bigfn(which::Symbol)
        r -> setprecision(BigFloat, bits) do
            rB = BigFloat(r; precision=bits); v, d = vd_at(rB)
            which === :val && return ComplexF64(v)
            which === :der && return ComplexF64(d)
            ComplexF64(_d2R_ode(v, d, s, m, ab, rB, ωb, λr))
        end
    end
    x = _cheb_nodes(lo, hi, _NCHEB)           # ascending (⇒ phase-unwrap is valid)
    L = Vector{ComplexF64}(undef, _NCHEB); Lp = similar(L); Lpp = similar(L)
    setprecision(BigFloat, bits) do
        @inbounds for k in 1:_NCHEB
            rB = BigFloat(x[k]; precision=bits)
            v, d = vd_at(rB)
            dd = _d2R_ode(v, d, s, m, ab, rB, ωb, λr)
            LR = d / v
            L[k]   = ComplexF64(log(v))
            Lp[k]  = ComplexF64(LR)
            Lpp[k] = ComplexF64(dd / v - LR*LR)
        end
    end
    # phase-unwrap the imaginary part of log along ascending r so the interpolant sees a smooth branch.
    @inbounds for k in 2:_NCHEB
        d = imag(L[k]) - imag(L[k-1])
        while d >  π; L[k] -= 2π*im; d = imag(L[k]) - imag(L[k-1]); end
        while d < -π; L[k] += 2π*im; d = imag(L[k]) - imag(L[k-1]); end
    end
    w = bary_weights(x)
    # det = Rup·Rind − Rin·Rupd is a difference of huge near-equal terms ⇒ ALL four legs must be
    # evaluated in BigFloat at the cancellation point (the In leg is Float64-safe only for SAMPLING).
    big_at = function (r)
        setprecision(BigFloat, bits) do
            rB = BigFloat(r; precision=bits)
            (Complex{BigFloat}(R.Up(rB)), Complex{BigFloat}(teukolsky_radial_deriv(R.Up, rB)),
             Complex{BigFloat}(R.In(rB)), Complex{BigFloat}(teukolsky_radial_deriv(R.In, rB)))
        end
    end
    # Out-of-span fallback: s≠0 exact BigFloat (never hot); s=0 native Float64 (the κ Abel tail is hot).
    outer_val, outer_der, outer_d2 = outer_native === nothing ?
        (bigfn(:val), bigfn(:der), bigfn(:d2)) : outer_native
    (LogChebInterp(lo, hi, x, w, L, Lp, Lpp, :val, outer_val),
     LogChebInterp(lo, hi, x, w, L, Lp, Lpp, :der, outer_der),
     LogChebInterp(lo, hi, x, w, L, Lp, Lpp, :d2,  outer_d2),
     big_at)
end

# =============================================================================================================
# BigFloat-type s=0 fast interp leg: per-node native BigFloat MST in the κ VoP quadgk is hours/mode —
# sample the smooth Up leg ONCE on a log-r Chebyshev panel instead. Up span must NOT reach rh (a too-wide
# log-r panel loses the accuracy the W̃ Wronskian cancellation needs; the In leg covers the rh region).
const _NCHEB_BIG = parse(Int, get(ENV,"KLC_NUP_FLOOR","70"))  # BigFloat s=0 Up-leg Cheb-node floor (env-override; default 70 = F8-safe)
const _KVOP_SPAN_HI = 120.0  # κ integration upper edge for the Up interp (Abel integrand <1e-24 past r∼100)

# LOG-r barycentric interpolant returning Complex{BigFloat} value/derivative for a SMOOTH nonzero R(r).
struct BigLogRInterp{F}
    lo::BigFloat; hi::BigFloat
    u::Vector{BigFloat}; w::Vector{BigFloat}
    L::Vector{Complex{BigFloat}}; Lp::Vector{Complex{BigFloat}}   # L=log R(r), Lp=dlogR/dr=R'/R sampled at nodes
    which::Symbol; outer::F
end
@inline function (ci::BigLogRInterp)(r::Real)
    rB = BigFloat(r)
    (ci.lo - BigFloat(1e-9) <= rB <= ci.hi + BigFloat(1e-9)) || return ci.outer(rB)
    uq = log(rB)
    Lr = bary_eval(ci.u, ci.w, ci.L, uq)
    ci.which === :val && return exp(Lr)
    Lpr = bary_eval(ci.u, ci.w, ci.Lp, uq)
    exp(Lr) * Lpr        # :der  (R' = R·(dlogR/dr); Lp samples already store R'/R)
end

# Build the BigFloat s=0 Up-leg fast interp closures (value/deriv) over [lo,hi] in LOG-r. Shared f_n; value+deriv
# in one pass (reuses _up_val_der). `nat_up/nat_upd` are the native BigFloat Up closures used out-of-span.
function _branchK_up_interp_big(s, l, m, ab::BigFloat, ωb::BigFloat, λteuk, norm_up, nu,
                                lo::Float64, hi::Float64, nat_up, nat_upd)
    bits = precision(ab)
    kappa = sqrt(1 - ab^2); eps = 2ωb; tau = (eps - m*ab)/kappa; rm = 1 - kappa
    Nmax = 80; atol = BigFloat(1e-30; precision=bits); rtol = BigFloat(1e-30; precision=bits)
    farr = first(_mst_fn_array(s, m, ab, ωb, nu, λteuk, Nmax))
    vd_at(rB) = (zhat = eps*rB/2 - eps*rm/2;
                 vd = _up_val_der(farr, s, m, eps, kappa, tau, nu, zhat, Nmax, atol, rtol);
                 (vd[1]/norm_up, vd[2]/norm_up))
    # N = 3l (floor 70): at high l the potential turning point r_t≈l/ω sits inside the span and a fixed N
    # leaves oscillatory inter-node error, amplified ×1e2-5e3 by the κ-sector cancellations (measured law).
    N = max(_NCHEB_BIG, parse(Int, get(ENV, "KLC_NUP_FACTOR", "3"))*l)          # env knob (default 3 = production)
    uF = _cheb_nodes(log(lo), log(hi), N)        # log-r Cheb nodes (ascending ⇒ phase-unwrap valid)
    uB = BigFloat.(uF; precision=bits); xB = exp.(uB)
    L = Vector{Complex{BigFloat}}(undef, N); Lp = similar(L)
    # Threaded node sampling (independent nodes; byte-identical to serial); phase-unwrap AFTER the barrier.
    Threads.@threads for k in 1:N
        v, d = vd_at(xB[k]); L[k] = log(v); Lp[k] = d/v
    end
    @inbounds for k in 2:N
        dd = imag(L[k]) - imag(L[k-1])
        while dd >  BigFloat(π); L[k] -= 2π*im; dd = imag(L[k]) - imag(L[k-1]); end
        while dd < -BigFloat(π); L[k] += 2π*im; dd = imag(L[k]) - imag(L[k-1]); end
    end
    w = bary_weights(uB)
    loB = BigFloat(lo; precision=bits); hiB = BigFloat(hi; precision=bits)
    (BigLogRInterp(loB, hiB, uB, w, L, Lp, :val, nat_up),
     BigLogRInterp(loB, hiB, uB, w, L, Lp, :der, nat_upd))
end

# In-leg node count: the In leg is smooth/monotone so 18 nodes reach full precision (measured); each
# node costs ~2 native BigFloat 2F1 evals, so N directly sets the panel build time.
const _NCHEB_IN = parse(Int, get(ENV,"KLC_NIN_FLOOR","18"))  # In-leg Cheb-node floor (env-override; default 18 = F8-byte-identical; raise for low-ℓ/low-m anti-scallop)

# Generic log-r Chebyshev BigFloat interpolant of a smooth nonzero radial function (the In leg at high ℓ
# where Float64 overflows): sample once (threaded), phase-unwrap after the barrier, serve barycentric.
function _branchK_leg_interp_big(Rfun, lo::Float64, hi::Float64; N::Int=_NCHEB_IN)
    bits = precision(BigFloat)
    uF = _cheb_nodes(log(lo), log(hi), N)
    uB = BigFloat.(uF; precision=bits); xB = exp.(uB)
    L = Vector{Complex{BigFloat}}(undef, N); Lp = similar(L)
    Threads.@threads for k in 1:N
        v = Complex{BigFloat}(Rfun(xB[k])); d = Complex{BigFloat}(teukolsky_radial_deriv(Rfun, xB[k]))
        L[k] = log(v); Lp[k] = d/v
    end
    @inbounds for k in 2:N
        dd = imag(L[k]) - imag(L[k-1])
        while dd >  BigFloat(π); L[k] -= 2π*im; dd = imag(L[k]) - imag(L[k-1]); end
        while dd < -BigFloat(π); L[k] += 2π*im; dd = imag(L[k]) - imag(L[k-1]); end
    end
    w = bary_weights(uB)
    loB = BigFloat(lo; precision=bits); hiB = BigFloat(hi; precision=bits)
    nat_v(r)  = Complex{BigFloat}(Rfun(BigFloat(real(r))))
    nat_d(r)  = Complex{BigFloat}(teukolsky_radial_deriv(Rfun, BigFloat(real(r))))
    (BigLogRInterp(loB, hiB, uB, w, L, Lp, :val, nat_v),
     BigLogRInterp(loB, hiB, uB, w, L, Lp, :der, nat_d))
end

"""
    HomFieldK{T}

a≠0 (spin-general) jumped homogeneous radial field: ``A^{\\rm in} R^{\\rm in}(r)`` for ``r < r_0``
and ``A^{\\rm up} R^{\\rm up}(r)`` for ``r ≥ r_0`` on the `branchK` Teukolsky legs (native Float64,
or the log-Chebyshev-barycentric BigFloat interpolants where the Float64 ``R^{\\rm up}`` magnitude
overflows). Evaluate with [`valK`](@ref) / [`derK`](@ref) / [`der2K`](@ref). The a=0 counterpart is
[`HomField`](@ref). Legacy: `KerrLorenzCirc/radial_kerr.jl`.
"""
struct HomFieldK{T}
    r0::T; Aup::Complex{T}; Ain::Complex{T}; br
end
"""
    valK(f, r)

Value of an a≠0 radial field ([`HomFieldK`](@ref) or [`ScaledHomFieldK`](@ref)): the In leg
``A^{\\rm in} R^{\\rm in}(r)`` for ``r < r_0``, else the Up leg ``A^{\\rm up} R^{\\rm up}(r)``
(times ``Δ^n`` for the scaled field).
"""
valK(f::HomFieldK,r)  = r<f.r0 ? f.Ain*f.br.inn(r)  : f.Aup*f.br.up(r)
"""
    derK(f, r)

First radial derivative of [`valK`](@ref) (analytic leg derivative; product rule for
[`ScaledHomFieldK`](@ref)).
"""
derK(f::HomFieldK,r)  = r<f.r0 ? f.Ain*f.br.ind(r)  : f.Aup*f.br.upd(r)
"""
    der2K(f, r)

Second radial derivative of [`valK`](@ref), closed by the field's homogeneous Teukolsky ODE
(product rule for [`ScaledHomFieldK`](@ref)).
"""
der2K(f::HomFieldK,r) = r<f.r0 ? f.Ain*f.br.indd(r) : f.Aup*f.br.updd(r)

# A homogeneous field times Δ(r)ⁿ (n usually 2 for P_{+2}=Δ²R_{+2}).  value/1st/2nd via product rule (analytic).
"""
    ScaledHomFieldK{T}

An a≠0 homogeneous field times ``Δ(r)^n`` — typically ``n=2`` for ``P_{+2} = Δ^2 R_{+2}``.
`valK`/`derK`/`der2K` apply the analytic product rule on top of the base [`HomFieldK`](@ref).
"""
struct ScaledHomFieldK{T}
    base::HomFieldK{T}; a::T; n::Int
end
function _Δpow(f::ScaledHomFieldK, r)
    Δ=r^2-2r+f.a^2; Δp=2r-2; Δpp=2.0
    A=Δ^f.n; Ap=f.n*Δ^(f.n-1)*Δp; App=f.n*(f.n-1)*Δ^(f.n-2)*Δp^2 + f.n*Δ^(f.n-1)*Δpp
    (A,Ap,App)
end
function valK(f::ScaledHomFieldK,r); A,_,_=_Δpow(f,r); A*valK(f.base,r); end
function derK(f::ScaledHomFieldK,r); A,Ap,_=_Δpow(f,r); Ap*valK(f.base,r)+A*derK(f.base,r); end
function der2K(f::ScaledHomFieldK,r); A,Ap,App=_Δpow(f,r); App*valK(f.base,r)+2Ap*derK(f.base,r)+A*der2K(f.base,r); end

# spin-general κ source coefficient: 𝒮(r)=½(r²+a²γ_ll)·h(r)
"""
    KappaFieldK{T}

a≠0 (spin-general) inhomogeneous κ field: the jumped homogeneous piece plus the Green's-function
variation-of-parameters particular integral on the s=0 legs, sourced by
``𝒮(r) = \\tfrac12 (r^2 + a^2 γ_{ll})\\, h(r)``; the ``r→∞``/horizon tails are Abel-regularized and
cached (`Iup_inf`/`Iin_hor`). Evaluate with [`κvalK`](@ref) / [`κderK`](@ref). The a=0 counterpart
is [`KappaField`](@ref). Legacy: `KerrLorenzCirc/radial_kerr.jl`.
"""
struct KappaFieldK{T}
    r0::T; a::T; γll::T; hom::HomFieldK{T}; brs0; hfield::HomFieldK{T}
    W̃::Complex{T}; Rin_lo::T; Rout_hi::T
    Tup_tail::Complex{T}; Tin_head::Complex{T}; Iup_inf::Complex{T}; Iin_hor::Complex{T}; rtol::Float64
end
srcK(kf::KappaFieldK, s) = 0.5*(s^2 + kf.a^2*kf.γll)*valK(kf.hfield, s)
function kappa_IK(kf::KappaFieldK, r)
    brs0=kf.brs0; W̃=kf.W̃; r0=kf.r0; rt=kf.rtol
    integ_up(s) = brs0.up(s) * srcK(kf,s) / W̃
    integ_in(s) = brs0.inn(s) * srcK(kf,s) / W̃
    if r ≤ r0
        I1 = first(quadgk(integ_up, r, r0; rtol=rt, atol=1e-12)) + kf.Tup_tail
        if r ≥ kf.Rin_lo
            I2 = kf.Iin_hor + first(quadgk(integ_in, kf.Rin_lo, r; rtol=rt, atol=1e-12))
        else
            # r < Rin_lo: DIRECT horizon→r Abel-regularized forward integral — the difference form
            # (Iin_hor + ∫_{Rin_lo}^r) cancels catastrophically at high l (the tiny true value is a
            # difference of pieces ~(Rin_lo/r)^{2l} larger). See F5_SMALLR_DIVERGENCE_DIAGNOSIS.md.
            ε1 = 4e-3; ε2 = 2e-3
            rh = 1 + sqrt(1 - kf.a^2)
            lo = rh + oftype(float(real(kf.r0)), 1e-6)
            fI(s, ε) = integ_in(s) * exp(-ε*(s - rh))
            t1 = first(quadgk(s -> fI(s, ε1), lo, r; rtol=1e-8, atol=0.0))
            t2 = first(quadgk(s -> fI(s, ε2), lo, r; rtol=1e-8, atol=0.0))
            I2 = 2t2 - t1
        end
    else
        # outer branch: no cancellation here — the difference form is fine (a direct r→∞ Abel mirror
        # was trialed and reverted: ~100× slower, no correctness gain).
        I1 = first(quadgk(integ_up, r, kf.Rout_hi; rtol=rt, atol=1e-12)) + kf.Iup_inf
        I2 = kf.Tin_head + first(quadgk(integ_in, r0, r; rtol=rt, atol=1e-12))
    end
    (I1, I2)
end
kappa_pK(kf::KappaFieldK, r)       = let (I1,I2)=kappa_IK(kf,r); kf.brs0.inn(r)*I1 + kf.brs0.up(r)*I2 end
kappa_p_derivK(kf::KappaFieldK, r) = let (I1,I2)=kappa_IK(kf,r); kf.brs0.ind(r)*I1 + kf.brs0.upd(r)*I2 end
"""
    κvalK(kf::KappaFieldK, r)

Spin-general ``κ(r)``: the jumped homogeneous piece plus the Green's-function particular integral
(VoP on the s=0 legs, Abel-regularized tails).
"""
κvalK(kf::KappaFieldK, r) = valK(kf.hom, r) + kappa_pK(kf, r)
"""
    κderK(kf::KappaFieldK, r)

First radial derivative ``κ'(r)`` (analytic: the VoP integrals differentiate onto the leg
derivatives only).
"""
κderK(kf::KappaFieldK, r) = derK(kf.hom, r) + kappa_p_derivK(kf, r)

# Four legs in BigFloat for the cancellation-prone jump/Wronskian solves: prefer br.big_at, else a fresh
# one-point BigFloat TeukolskyRadial when |R^up(r)| is a cancellation risk; `nothing` = Float64 is safe.
function _big_legs_at(br, r)
    if hasproperty(br, :big_at) && br.big_at !== nothing
        return br.big_at(r)
    end
    upr = br.up(r)
    (upr isa Complex{Float64} && abs(upr) ≥ _BIGMAG) || return nothing
    # re-wrapped branches without a stored R degrade to Float64 (carry big_at if the BigFloat det is needed).
    hasproperty(br, :R) || return nothing
    up0 = br.R.Up
    s = up0.s; l = up0.l; m = up0.m; a = Float64(up0.a); ω = Float64(up0.omega); λ = up0.lambda
    magpow = log10(max(abs(br.up(r)), 1.0)); bits = _bigprec_for_mag(magpow)
    setprecision(BigFloat, bits) do
        ab = BigFloat(a; precision=bits); ωb = BigFloat(ω; precision=bits)
        R = TeukolskyRadial(s, l, m, ab, ωb; lambda=BigFloat(Float64(λ); precision=bits))
        rB = BigFloat(Float64(r); precision=bits)
        (Complex{BigFloat}(R.Up(rB)), Complex{BigFloat}(teukolsky_radial_deriv(R.Up, rB)),
         Complex{BigFloat}(R.In(rB)), Complex{BigFloat}(teukolsky_radial_deriv(R.In, rB)))
    end
end

# spin-general 2×2 jump solve ⇒ (A_up, A_in). The det is O(1) while the legs reach ~1e18 at high l —
# a catastrophic Float64 cancellation ⇒ form det/Aup/Ain in BigFloat when at risk, down-convert after.
function jump_ampsK(br, r0, J0, J1)
    legs = _big_legs_at(br, r0)
    if legs !== nothing
        Rup, Rupd, Rin, Rind = legs
        J0b = Complex{BigFloat}(J0); J1b = Complex{BigFloat}(J1)
        det = Rup*(-Rind) - (-Rin)*Rupd
        Aup = (J0b*(-Rind) - (-Rin)*J1b)/det
        Ain = (Rup*J1b - Rupd*J0b)/det
        return (ComplexF64(Aup), ComplexF64(Ain))
    end
    Rup=br.up(r0); Rin=br.inn(r0); Rupd=br.upd(r0); Rind=br.ind(r0)
    det = Rup*(-Rind) - (-Rin)*Rupd
    Aup = (J0*(-Rind) - (-Rin)*J1)/det
    Ain = (Rup*J1 - Rupd*J0)/det
    (Aup, Ain)
end

# κ normalization W̃ = ΔW·(R^in·R^up' − R^up·R^in')|_{rW} (Wronskian × Δ ⇒ r-independent; rW=r0).
# Same cancellation class as jump_ampsK ⇒ BigFloat legs when at risk.
function _W̃_kerr(br0, ΔW, rW)
    legs = _big_legs_at(br0, rW)
    legs === nothing && return ΔW*(br0.inn(rW)*br0.upd(rW) - br0.up(rW)*br0.ind(rW))
    Rup, Rupd, Rin, Rind = legs
    ΔB = BigFloat(Float64(ΔW); precision=precision(real(Rup)))
    ComplexF64(ΔB*(Rin*Rupd - Rup*Rind))
end

"""
    jump_data_kerr(a, r0, l, m; lmax=auto, aux=nothing) -> NamedTuple

a≠0 jumps {J0,J1} of P_{-2},P_{-1},h,κ at (a,r0,l) and the spin-general spectrum (Λ,λ,λ0,γll,ω). P_{-2} from
`jumps_pm`, h from `hj1`; P_{-1},κ from the banded `closure_kerr` (App-H a=0.99 lower-validated). Pass a precomputed
`aux` (from `closure_kerr(...;return_aux=true)`) to avoid re-solving the closure across l.
"""
function jump_data_kerr(a, r0, l, m; lmax=nothing, aux=nothing)
    if aux === nothing
        lmx = lmax === nothing ? max(2*l, l+8) : lmax
        sol, aux = closure_kerr(a, r0, m; lmax=lmx, return_aux=true)
        clo = sol[l]
    else
        # caller supplied aux AND the per-l closure separately is unusual; recompute solution sans aux is wasteful,
        # so require aux to carry it — but closure_kerr returns sol separately. Handle by re-solving once if needed.
        error("jump_data_kerr: pass aux=nothing or use reconstruct_kerr which threads sol+aux together")
    end
    bg = aux.bg
    J0m2,J1m2 = jumps_pm(bg, l, -1)
    J1h = hj1(bg, l)
    (ω=bg.ω, Δ0=bg.Δ0, a=bg.a,
     Pm2=(Complex{typeof(bg.r0)}(J0m2), Complex{typeof(bg.r0)}(J1m2)),
     Pm1=(clo[1], clo[2]), h=(zero(Complex{typeof(bg.r0)}), Complex{typeof(bg.r0)}(J1h)),
     κ=(clo[3], clo[4]),
     Λ=aux.Λ[l], λ=aux.λ[l], λ0=aux.λ0[l], γll=aux.γll[l], m=m)
end

# Module-level Abel-tail ε scale (default 1 = F8-validated). The FIGURE (log plot, ~1e-3 κ suffices) can set this
# >1 to shorten the κ Abel tail; callers that don't thread `eps_scale` (kerr_modes→reconstruct_kerr) pick it up here.
const _RECON_EPS_SCALE = Ref{Int}(1)
# Module-level default radial span for the high-magnitude BigFloat-interp branch gate. `kerr_modes`/`kerr_mode`
# set this to the actual mode grid range so branchK knows where |R^up| is sampled; standalone callers (gates,
# tests, a single reconstruct_kerr) leave it `nothing` ⇒ the native Float64 path (the magnitude gate is inert
# without a span, so those byte-unchanged). `reconstruct_kerr(...; rspan=(lo,hi))` overrides per-call.
const _RECON_RSPAN = Ref{Union{Nothing,Tuple{Float64,Float64}}}(nothing)
"""
    reconstruct_kerr(a, r0, l; m=2, Rin_lo=nothing, Rout_hi=nothing, aux=nothing, sol=nothing) -> NamedTuple

a≠0 reconstructed radial fields at (a,r0,l): `rec.Pm2,rec.Pm1,rec.h` are `HomFieldK` (use `valK/derK/der2K`),
`rec.κ` a `KappaFieldK` (use `κvalK/κderK`). `rec.jd` carries ω,Δ0,Λ,λ,λ0,γll and the jumps. Pass `(sol,aux)`
from `closure_kerr(...;return_aux=true)` to reuse one banded solve across l (and threads).
"""
function reconstruct_kerr(a, r0, l; m=2, Rin_lo=nothing, Rout_hi=nothing, aux=nothing, sol=nothing,
                          eps_scale=_RECON_EPS_SCALE[], rspan=_RECON_RSPAN[])
    if aux === nothing || sol === nothing
        lmx = max(2*l, l+8)
        sol, aux = closure_kerr(a, r0, m; lmax=lmx, return_aux=true)
    end
    clo = sol[l]
    bg = aux.bg
    T = typeof(bg.r0); CT = Complex{T}
    ω = bg.ω; aT = bg.a
    # ω-guard: assert the circular-equatorial frequency against the RAW args — a corrupted-run-state
    # incident once produced a wrong bg.ω (F5_SMALLR_DIVERGENCE_DIAGNOSIS.md); abort loudly.
    let ωref = m/(sqrt(float(real(r0))^3) + float(real(a)))
        abs(Float64(real(ω)) - ωref) ≤ 1e-10*abs(ωref) ||
            error("reconstruct_kerr ω-guard: bg.ω=$(Float64(real(ω))) ≠ m/(√r0³+a)=$ωref (corrupted run state?)")
    end
    Λ = aux.Λ[l]; λ = aux.λ[l]; λ0 = aux.λ0[l]; γll = aux.γll[l]
    J0m2,J1m2 = jumps_pm(bg, l, -1); J1h = hj1(bg, l)
    jd = (ω=ω, Δ0=bg.Δ0, a=aT, Pm2=(CT(J0m2),CT(J1m2)), Pm1=(CT(clo[1]),CT(clo[2])),
          h=(zero(CT),CT(J1h)), κ=(CT(clo[3]),CT(clo[4])), Λ=Λ, λ=λ, λ0=λ0, γll=γll, m=m)

    brm2 = branchK(-2, l, m, aT, ω, Λ; rspan=rspan, r0=r0); brm1 = branchK(-1, l, m, aT, ω, λ; rspan=rspan, r0=r0)
    br0  = branchK(0, l, m, aT, ω, λ0; rspan=rspan, r0=r0)
    Pm2 = HomFieldK{T}(r0, jump_ampsK(brm2,r0,jd.Pm2...)..., brm2)
    Pm1 = HomFieldK{T}(r0, jump_ampsK(brm1,r0,jd.Pm1...)..., brm1)
    h   = HomFieldK{T}(r0, jump_ampsK(br0, r0,jd.h...)...,  br0)
    κhom = HomFieldK{T}(r0, jump_ampsK(br0, r0, jd.κ...)..., br0)

    rW=T(r0); ΔW=rW^2-2rW+aT^2
    W̃ = _W̃_kerr(br0, ΔW, rW)
    rh = one(T) + sqrt(one(T) - aT^2)
    Rin_lo = Rin_lo===nothing ? T(r0)-min(T(2.6), (T(r0)-rh)/2) : T(Rin_lo)
    Rout_hi = Rout_hi===nothing ? T(r0)+T(11) : T(Rout_hi)
    rtol=1e-9
    srcf(s) = 0.5*(s^2 + aT^2*γll)*valK(h,s)
    iup(s)=br0.up(s) *srcf(s)/W̃
    iin(s)=br0.inn(s)*srcf(s)/W̃
    # Abel-tail upper limit: BigFloat caps at the Up-interp edge (integrand <1e-24 past r~100; keeps the
    # quadgk off the slow out-of-span fallback); Float64 keeps the full lo+30/ε.
    _abel_hi(lo, ε) = (T === Float64) ? lo + 30.0/ε : min(lo + 30/ε, T(_KVOP_SPAN_HI))
    abel_tail(integ, lo, ε) = first(quadgk(s->integ(s)*exp(-ε*(s-lo)), lo, _abel_hi(lo, ε); rtol=1e-8, atol=1e-14))
    function abel0(integ, lo)
        ε1=4e-3*eps_scale; ε2=2e-3*eps_scale     # eps_scale>1 shortens the tail (figure-only; gates use 1)
        t1=abel_tail(integ,lo,ε1); t2=abel_tail(integ,lo,ε2)
        2t2 - t1
    end
    Iup_inf = abel0(iup, Rout_hi)
    function abel0_in(integ, b)
        ε1=4e-3*eps_scale; ε2=2e-3*eps_scale
        f(s,ε)=integ(s)*exp(-ε*(s-rh))
        lo = rh + T(1e-6)
        t1=first(quadgk(s->f(s,ε1), lo, b; rtol=1e-8, atol=1e-14))
        t2=first(quadgk(s->f(s,ε2), lo, b; rtol=1e-8, atol=1e-14))
        2t2 - t1
    end
    Iin_hor = abel0_in(iin, Rin_lo)
    Tup_tail=first(quadgk(iup, T(r0), Rout_hi; rtol=rtol, atol=1e-12)) + Iup_inf
    Tin_head=first(quadgk(iin, Rin_lo, T(r0); rtol=rtol, atol=1e-12)) + Iin_hor
    κf = KappaFieldK{T}(r0, aT, γll, κhom, br0, h, W̃, Rin_lo, Rout_hi, Tup_tail, Tin_head, Iup_inf, Iin_hor, rtol)
    (jd=jd, Pm2=Pm2, Pm1=Pm1, h=h, κ=κf, aux=aux, sol=sol, rspan=rspan)
end

# ---- spin-general residuals ----
resK_sm2(r,a,ω,m,Λ,R,Rp,Rpp)=(Δ=r^2-2r+a^2;Δp=2r-2;c=cfunK(r,a,ω,m);cp=cpfunK(r,a,ω,m); Δ*(Rpp-cp*R-c^2*R-(Δp/Δ)*(Rp-c*R))-6im*ω*r*R-Λ*R)
resK_sm1(r,a,ω,m,λ,R,Rp,Rpp)=(Δ=r^2-2r+a^2;c=cfunK(r,a,ω,m);cp=cpfunK(r,a,ω,m); Δ*(Rpp-cp*R-c^2*R)-2im*ω*r*R-λ*R)
resK_s0(r,a,ω,m,λ0,R,Rp,Rpp)=(Δ=r^2-2r+a^2;Δp=2r-2;K=ω*(r^2+a^2)-a*m;Kp=2ω*r; Δ*Rpp+Δp*Rp+(im*Kp+K^2/Δ)*R-2im*ω*r*R-λ0*R)
function resK_kappa(kf::KappaFieldK, r, ω, m, λ0)
    brs0=kf.brs0; aT=kf.a
    local I1,I2; I1,I2 = kappa_IK(kf, r)
    S = srcK(kf, r)
    κ  = valK(kf.hom,r)  + brs0.inn(r)*I1  + brs0.up(r)*I2
    κp = derK(kf.hom,r)  + brs0.ind(r)*I1  + brs0.upd(r)*I2
    κpp = der2K(kf.hom,r) + brs0.indd(r)*I1 + brs0.updd(r)*I2 + (-brs0.ind(r)*brs0.up(r) + brs0.upd(r)*brs0.inn(r))*S/kf.W̃
    Δ=r^2-2r+aT^2; Δp=2r-2; K=ω*(r^2+aT^2)-aT*m; Kp=2ω*r
    Δ*κpp + Δp*κp + (im*Kp + K^2/Δ)*κ - 2im*ω*r*κ - λ0*κ - S
end

"""
    radial_gates_kerr(a, r0; ls=2:6, m=2) -> (maxjump, maxres_h, maxjump_k, maxres_k)

Spin-general self-consistency gates (GATE a): P_{-2},P_{-1},h jumps+continuity, the h s=0-ODE residual, and the
κ jumps + inhomogeneous-ODE residual L[κ]−½(r²+a²γ_ll)h. One banded `closure_kerr` solve is shared across l.
"""
function radial_gates_kerr(a, r0; ls=2:6, m=2, lmax=nothing)
    lsv = collect(ls)
    lmx = lmax === nothing ? max(2*maximum(lsv), maximum(lsv)+8) : lmax
    sol, aux = closure_kerr(a, r0, m; lmax=lmx, return_aux=true)
    println("\n############################################################################")
    @printf("#######  RADIAL_KERR GATES — a=%.3f, m=%d, r0=%.1f  #######\n", Float64(a), m, Float64(r0))
    println("############################################################################")
    ε=1e-6
    rh = 1.0 + sqrt(1.0 - Float64(a)^2)
    grid = Float64[r0+1.5, r0+3.0, r0+6.0, r0-1.0, r0-2.0]
    grid = filter(r->r>rh+0.2, grid)
    # rspan covers the jump points + residual grid so high-l gates sample the accurate interp leg.
    span = (max(rh + 0.1, minimum(grid) - 0.5), maximum(grid) + 0.5)
    maxjump=0.0; maxres_h=0.0; maxres_k=0.0; maxjump_k=0.0
    for l in lsv
        rec = reconstruct_kerr(a, r0, l; m=m, aux=aux, sol=sol, rspan=span); jd=rec.jd; ω=Float64(jd.ω)
        flush(stdout)
        @printf("\n=== l=%d  (ω=%.5f, Λ=%.4f, λ=%.4f, λ0=%.4f, γll=%.4f) ===\n",
                l, ω, Float64(jd.Λ), Float64(jd.λ), Float64(jd.λ0), Float64(jd.γll))
        for (nm,f,J) in (("P_{-2}",rec.Pm2,jd.Pm2), ("P_{-1}",rec.Pm1,jd.Pm1), ("h     ",rec.h,jd.h))
            jv = valK(f,r0+ε)-valK(f,r0-ε); jd1=derK(f,r0+ε)-derK(f,r0-ε)
            e0 = abs(J[1])>1e-8 ? abs(jv-J[1])/abs(J[1]) : abs(jv)
            e1 = abs(J[2])>1e-8 ? abs(jd1-J[2])/abs(J[2]) : abs(jd1)
            maxjump=max(maxjump,e0,e1)
            @printf("  GATE1 %s: [f]=%+.4e%+.4ei (J0 rel %.1e)  [f']=%+.4e%+.4ei (J1 rel %.1e)\n",
                    nm, real(jv),imag(jv),e0, real(jd1),imag(jd1),e1)
        end
        relh=maximum(let R=valK(rec.h,r),Rp=derK(rec.h,r),Rpp=der2K(rec.h,r);
                         abs(resK_s0(r,Float64(a),ω,m,Float64(jd.λ0),R,Rp,Rpp))/(abs(Float64(jd.λ0)*R)+1e-300) end for r in grid)
        maxres_h=max(maxres_h,relh)
        @printf("  GATE2 h ODE residual (max over grid) = %.2e\n", relh)
        kf=rec.κ
        jv = κvalK(kf,r0+ε)-κvalK(kf,r0-ε); jd1=κderK(kf,r0+ε)-κderK(kf,r0-ε)
        e0 = abs(jd.κ[1])>1e-8 ? abs(jv-jd.κ[1])/abs(jd.κ[1]) : abs(jv)
        e1 = abs(jd.κ[2])>1e-8 ? abs(jd1-jd.κ[2])/abs(jd.κ[2]) : abs(jd1)
        maxjump_k=max(maxjump_k,e0,e1)
        @printf("  GATE3 κ jumps: [κ]=%+.4e%+.4ei (rel %.1e)  [κ']=%+.4e%+.4ei (rel %.1e)\n",
                real(jv),imag(jv),e0, real(jd1),imag(jd1),e1)
        relk=0.0; rworst=0.0
        for r in grid
            res=resK_kappa(kf,r,ω,m,Float64(jd.λ0)); S=Float64(abs(srcK(kf,r))); scale=S+abs(Float64(jd.λ0)*κvalK(kf,r))+1e-300
            rel=abs(res)/scale; if rel>relk; relk=rel; rworst=r; end
            maxres_k=max(maxres_k,rel)
        end
        @printf("  GATE3 κ inhomog ODE residual |L[κ]−½(r²+a²γ_ll)h|/scale (max) = %.2e (at r=%.1f)\n", relk, rworst)
    end
    @printf("\n---- a=%.3f r0=%.1f SUMMARY ----\n", Float64(a), Float64(r0))
    @printf("  GATE1 max jump rel.err (P_{-2},P_{-1},h)      = %.2e  ⟹ %s\n", maxjump, maxjump<1e-4 ? "PASS" : "CHECK")
    @printf("  GATE2 max h s=0-ODE residual                  = %.2e  ⟹ %s\n", maxres_h, maxres_h<1e-8 ? "PASS" : "CHECK")
    @printf("  GATE3 max κ jump rel.err                      = %.2e  ⟹ %s\n", maxjump_k, maxjump_k<1e-4 ? "PASS" : "CHECK")
    @printf("  GATE3 max κ inhomog-ODE residual              = %.2e  ⟹ %s\n", maxres_k, maxres_k<1e-6 ? "PASS" : "CHECK")
    (maxjump,maxres_h,maxjump_k,maxres_k)
end
