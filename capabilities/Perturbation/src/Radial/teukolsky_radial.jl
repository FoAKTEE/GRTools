# Public `TeukolskyRadial` dispatcher + `teukolsky_radial_deriv`. Resolves λ (SWSH eigenvalue) and
# ν (renormalized angular momentum), computes the MST amplitudes, and normalizes the homogeneous
# In/Up solutions to UNIT TRANSMISSION (R = raw[norm=1]/Btrans|Ctrans; amplitudes / Transmission).
# Static ω=0 modes route to the closed-form 2F1 path (the MST machinery degenerates at ε=0). Real ω.
# Needs `_mst_amplitudes` (AmplitudesFluxes), so it is included after that leaf.

"""
    TeukolskyRadial(s, l, m, a, omega; method=:MST, Nmax=80, lambda=nothing, nu=nothing,
                    gsn=NamedTuple()) -> (In=TeukolskyRadialFunction, Up=...)

Public homogeneous radial Teukolsky solutions, normalized to unit transmission. λ via
`SpinWeightedSpheroidalEigenvalue`, ν via `RenormalizedAngularMomentum` (both resolved when not
supplied). `method=:MST` is the MST series path; ω=0 uses the static closed-form 2F1 path.
`method=:GSN` is the Generalized Sasaki-Nakamura backend (arXiv:2306.16469 absorption): the full
s ∈ {0,±1,±2} × {In,Up} surface at real ω, type-generic to BigFloat; backend options pass
through the `gsn` NamedTuple (`rsin`, `rsout`, `hor_order`, `inf_order`, `tol`, `form`,
`autowiden` — see `_gsn_teukolsky_radial`). For :GSN, ν is carried for cross-backend
comparability only (never used by the GSN pipeline).
"""
function TeukolskyRadial(s::Integer, l::Integer, m::Integer, a::Real, omega::Real;
                         method::Symbol=:MST, Nmax::Integer=80, lambda=nothing, nu=nothing,
                         gsn::NamedTuple=NamedTuple())
    method in (:MST, :HeunC, :SasakiNakamura, :GSN) ||
        error("TeukolskyRadial: method must be :MST, :HeunC, :SasakiNakamura or :GSN (got :$method).")
    iszero(omega) && return _teukolsky_radial_static(s, l, m, a)
    lambda = lambda === nothing ? SpinWeightedSpheroidalEigenvalue(s, l, m, a * omega) : lambda
    nu = nu === nothing ? RenormalizedAngularMomentum(s, l, m, a, omega, lambda) : nu
    method === :GSN && return _gsn_teukolsky_radial(s, l, m, a, omega, lambda, nu; gsn...)
    amp = _mst_amplitudes(s, m, a, omega, nu, lambda; Nmax=Nmax)
    mkamp(inc, trans, ref) = (Incidence=inc / trans, Transmission=one(inc / trans), Reflection=ref / trans)
    In = TeukolskyRadialFunction(s, l, m, float(a), float(omega), nu, lambda, :In,
                                 amp.InTrans, mkamp(amp.InInc, amp.InTrans, amp.InRef), method)
    Up = TeukolskyRadialFunction(s, l, m, float(a), float(omega), nu, lambda, :Up,
                                 amp.UpTrans, mkamp(amp.UpInc, amp.UpTrans, amp.UpRef), method)
    return (In=In, Up=Up)
end

"""
    teukolsky_radial_deriv(R, r) -> Complex

Normalized first derivative dR/dr of a `TeukolskyRadialFunction`. Static (ω=0) → closed-form
`_static_radial_deriv`; `:MST` series path; or the genuine `:HeunC` confluent-Heun derivative.
"""
function teukolsky_radial_deriv(R::TeukolskyRadialFunction, r::Real)
    iszero(R.omega) &&
        return _static_radial_deriv(R.s, R.l, R.m, R.a, R.bc, r; norm=R.norm)
    R.method === :HeunC && return _heunc_radial_deriv(R, r)
    R.method === :SasakiNakamura && return _sn_value(R, r, 1)
    R.method === :GSN && return _gsn_value(R, r, 1)
    return R.bc === :In ?
        _mst_radial_in_deriv(R.s, R.m, R.a, R.omega, R.nu, R.lambda, r; norm=R.norm) :
        _mst_radial_up_deriv(R.s, R.m, R.a, R.omega, R.nu, R.lambda, r; norm=R.norm)
end
