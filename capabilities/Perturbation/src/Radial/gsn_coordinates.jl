# :GSN backend — Kerr primitives + tortoise coordinate map (GSN-2 absorption of the vendored
# GeneralizedSasakiNakamura.jl v0.8.5 commit 982b087, arXiv:2306.16469; plan:
# progress/emri-disk/GSN_PLAN.md). M=1 throughout (upstream convention, same as this capability).
#
# Local `_gsn_*` Kerr one-liners mirror upstream src/Homogeneous/Kerr.jl VERBATIM (they are also
# definitionally identical to the Spacetime horizon/Δ forms; kept local so the GSN expression
# blobs below stay textually faithful to their Mathematica export). The tortoise map replaces
# upstream Roots.find_zero with the in-firewall NumericsCore rootfinders and — the audited
# BigFloat blocker GSN_PLAN §3c-1 — carries a {T}-typed bracket instead of the Float64 (0, 1.4).

"Kerr Δ(a,r) = r² − 2r + a² (M=1). [GSN Kerr.jl]"
_gsn_Delta(a, r) = r^2 - 2*r + a^2
"Outer horizon r₊ = 1 + √(1−a²) (M=1). [GSN Kerr.jl]"
_gsn_rp(a) = 1 + sqrt(1 - a^2)
"Inner horizon r₋ = 1 − √(1−a²) (M=1). [GSN Kerr.jl]"
_gsn_rm(a) = 1 - sqrt(1 - a^2)
"Horizon angular frequency Ω_H = a/(2r₊) (M=1). [GSN Kerr.jl]"
_gsn_Omega_h(a) = a/(2*_gsn_rp(a))
"Teukolsky K(r) = (r²+a²)ω − ma (M=1). [GSN Kerr.jl]"
_gsn_KK(m::Int, a, omega, r) = (r^2 + a^2)*omega - m*a

"""
    _gsn_rstar_from_r(a, r) -> rstar

Boyer-Lindquist r → tortoise r* (M=1): r* = r + 2r₊/(r₊−r₋)·log((r−r₊)/2) − 2r₋/(r₊−r₋)·log((r−r₋)/2)
for |a|<1, with the distinct extremal |a|=1 form. Type-generic. [GSN Coordinates.jl rstar_from_rp]
"""
function _gsn_rstar_from_r(a, r)
    rp = _gsn_rp(a)
    return _gsn_rstar_from_rp(a, r - rp)
end

function _gsn_rstar_from_rp(a, r_from_rp)
    rp = _gsn_rp(a)
    rm = _gsn_rm(a)
    if abs(a) < 1
        # [GSN Coordinates.jl l.14, verbatim algebra; 2*1 factors are the upstream M=1 export]
        return rp + r_from_rp + (2*rp)/(rp - rm) * log(r_from_rp/2) - (2*rm)/(rp - rm) * log((r_from_rp + rp - rm)/2)
    elseif abs(a) == 1
        return rp + r_from_rp - log(4*one(r_from_rp)) + 2*(log(r_from_rp) - 1/r_from_rp)
    else
        throw(ArgumentError("a must be in the range [-1, 1]"))
    end
end

"""
    _gsn_r_from_rstar(a, rstar) -> r

Tortoise r* → Boyer-Lindquist r (M=1), assuming r ≥ r₊. Bracketed a42 for rstar ≤ 0 (h′ = r−r₊
is bounded: rstar_from_rp(h′)=0 has h′ ≲ 1.328 as |a|→1, upstream comment), Newton with the
analytic derivative dr*/dr = (r²+a²)/Δ otherwise. The bracket is {T}-typed — this fixes the
audited upstream Float64 cap `find_zero(f, (0, 1.4))` (GSN_PLAN §3c-1). [GSN Coordinates.jl]
"""
function _gsn_r_from_rstar(a, rstar)
    T = typeof(float(a + rstar))
    rp = _gsn_rp(T(a))
    f(x) = _gsn_rstar_from_rp(T(a), abs(x)) - T(rstar)
    # Newton with the analytic derivative dr*/dr = (r²+a²)/Δ — quadratic, so a coarse seed
    # reaches full T-precision in ~log2(bits) steps. The bracketing stage below runs at a FIXED
    # coarse tolerance and the high-precision work is confined to the polish (2 log evaluations
    # per step) — an eps(T)-tight a42 would cost hundreds of BigFloat log evaluations per call,
    # which is the hot path of every ODE right-hand-side evaluation near the horizon.
    function newton_polish(x)
        for _ in 1:(8 + ceil(Int, log2(precision(one(T)) + 1)))
            fx = f(x)
            dfx = sign(x) * ((rp + abs(x))^2 + T(a)^2) / _gsn_Delta(T(a), rp + abs(x))
            dx = fx / dfx
            x -= dx
            abs(dx) <= 4*eps(T)*(1 + abs(x)) && break
        end
        return abs(x)
    end
    if rstar <= 0
        # NumericsCore.a42 on the {T}-typed bracket (0, 1.4) — the audited upstream Float64 cap
        # is gone (GSN_PLAN §3c-1) — at 1e-6 relative only, then Newton to full precision.
        seed = a42(f, zero(T), T(14)/10; rtol=T(1)/10^6)
        return rp + newton_polish(seed)
    else
        return rp + newton_polish(T(rstar))    # r* ≈ r seed for large r* [upstream Newton branch]
    end
end
