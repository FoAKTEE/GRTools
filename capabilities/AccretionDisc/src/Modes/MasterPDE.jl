# RT-A1 core: assemble the master PDE (arXiv:2601.19123 eq:MasterEquations, l.385) at the equator
# and reduce it to the radial operator L_m by the test-function method, extracting (A_m, B_m, P).
# This is the RIGOROUS-construction route that VALIDATES the WKB B_m_full already shipped:
#
#   c_s²∇_μ(Q^μ_α prj^{αν}∇_ν h) + [R_μ Q^μ_α prj^{αν} + c_s² Q^ν_α B^α]∇_ν h
#       − [iαω̃ − (c_s²∇_μ+R_μ)Q^μ_α B^α] h = 𝔉,    prj^{μν}=g^{μν}+U^μU^ν (l.54).
#
# Apply L to radial test functions 1, (r−r₀), (r−r₀)²: P = L[(·)²] (∂² coeff), L[(·)¹], L[1] ⇒
# A_m = L[(r−r₀)]/P, B_m = L[1]/P. Equatorial leading order; the PMLA h²/Gaussian height-average is
# a [FUTURE] refinement. Reuses Velocity.q_operator + Spacetime.{kerr_metric,kerr_metric_inv}.
# Complex{T}, BigFloat-inherited. From DiscModes/relativistic_masterPDE.jl.

using ..Modes: RelativisticDisc, disc_OmegaK, csound, eden0, q_operator
import Spacetime: kerr_metric, kerr_metric_inv

# Background pieces (R_μ lower, B^α upper, α, prj, √−g) at the equator — recomputed locally for the
# assembly. Returns a NamedTuple of the per-point geometric/thermodynamic quantities.
function _bg_pieces(d::RelativisticDisc{T}, r) where {T}
    a = d.a; θ = T(π) / 2; hh = sqrt(eps(T)) * max(abs(r), one(T))
    gi = kerr_metric_inv(a, r, θ)
    Ω = disc_OmegaK(d, r)
    αof(s) = (Ωs = disc_OmegaK(d, s); gs = kerr_metric(a, s, θ);
              sqrt(-one(T) / (gs[1, 1] + 2 * Ωs * gs[1, 4] + Ωs^2 * gs[4, 4])))
    α = αof(r); U = T[α, 0, 0, α * Ω]
    prj = gi .+ U * U'                                          # prj^{μν} = g^{μν}+U^μU^ν
    pof(s) = csound(d, s)^2 * eden0(d, s); dpdr = (pof(r + hh) - pof(r - hh)) / (2hh)
    εp = eden0(d, r) + pof(r)
    Rlo = zeros(T, 4); Rlo[2] = dpdr / εp                       # R_r (lower)
    cs2of(s) = csound(d, s)^2; dlogcs2 = (log(cs2of(r + hh)) - log(cs2of(r - hh))) / (2hh)
    Blo = zeros(T, 4); Blo[2] = -dlogcs2 / (1 + cs2of(r)); Bup = gi * Blo   # B^α (upper)
    sqrtg = r^2                                                 # √−g = Σ sinθ = r² at equator
    (g = kerr_metric(a, r, θ), gi = gi, prj = prj, U = U, α = α, Ω = Ω,
     Rlo = Rlo, Bup = Bup, cs2 = cs2of(r), sqrtg = sqrtg)
end

# W^μ(r) = Q^μ_α prj^{αν} ∂_ν h for a radial test (value hval, derivative hr) and mode m
# (∂_φ = im, ∂_t = −iω with ω = mΩ_p). Returns (W, bg).
function _Wvec(d, r, r_p, m, hval, hr)
    bg = _bg_pieces(d, r); Q = q_operator(d, r, r_p, m)
    ω = m * disc_OmegaK(d, r_p)
    dh = Complex.(zeros(typeof(bg.α), 4)); dh[1] = -im * ω * hval; dh[2] = hr; dh[4] = im * m * hval
    W = Q * (bg.prj * dh)                                       # W^μ = Q^μ_α (prj^{αν}∂_ν h)
    W, bg
end

"""
    master_radial_coeffs(d::RelativisticDisc, r, r_p, m) -> (A_m, B_m, P)

Reduced radial-ODE coefficients at `r` from the equatorial leading-order reduction of the master
PDE (eq:MasterEquations). Apply the operator `L` to radial test functions `1, (r−r₀), (r−r₀)²`:
`P=L[(r−r₀)²]` (the ∂² coefficient), `A_m=L[(r−r₀)]/P`, `B_m=L[1]/P`. Cross-validated against the
independent exact `A_m=−d/dr log𝒲`. The `∂_r(√−g W^r)` divergence is taken by finite difference of
the radial flux. Complex coefficients (BigFloat-inherited).
"""
function master_radial_coeffs(d::RelativisticDisc{T}, r, r_p, m) where {T}
    δ = sqrt(eps(T)) * max(abs(r), one(T))
    ω = m * disc_OmegaK(d, r_p)
    # L acting on a radial test fn h(s)=h0+hp(s−r)+½h2(s−r)², via the ∂_r(√−g W^r) divergence
    function Lval(h0, hp, h2)
        bg = _bg_pieces(d, r)
        function WrSg(s)                                        # √−g W^r as a function of s near r
            hs = h0 + hp * (s - r) + (h2 / 2) * (s - r)^2; hps = hp + h2 * (s - r)
            s^2 * first(_Wvec(d, s, r_p, m, hs, hps))[2]        # √−g(s)=s², W^r component
        end
        divr = (WrSg(r + δ) - WrSg(r - δ)) / (2δ) / bg.sqrtg    # (1/√−g)∂_r(√−g W^r)
        W0, _ = _Wvec(d, r, r_p, m, h0, hp)
        divtφ = (-im * ω) * W0[1] + (im * m) * W0[4]            # (−iω)W^t+(im)W^φ
        princ = bg.cs2 * (divr + divtφ)
        # second term: [R_μ Q^μ_α prj^{αν} + c_s² Q^ν_α B^α] ∂_ν h
        Q = q_operator(d, r, r_p, m)
        dh = Complex.(zeros(T, 4)); dh[1] = -im * ω * h0; dh[2] = hp; dh[4] = im * m * h0
        vecν = zeros(Complex{T}, 4)
        for ν in 1:4
            s1 = sum(bg.Rlo[μ] * Q[μ, α] * bg.prj[α, ν] for μ in 1:4, α in 1:4)
            s2 = bg.cs2 * sum(Q[ν, α] * bg.Bup[α] for α in 1:4)
            vecν[ν] = s1 + s2
        end
        second = sum(vecν[ν] * dh[ν] for ν in 1:4)
        # third term: −[iαω̃ − (c_s²∇_μ+R_μ)Q^μ_α B^α] h ; ω̃=ω−mΩ ; V^μ=Q^μ_α B^α
        ω̃ = ω - m * bg.Ω
        function Vr(s)
            bg2 = _bg_pieces(d, s); Qs = q_operator(d, s, r_p, m)
            s^2 * sum(Qs[2, α] * bg2.Bup[α] for α in 1:4)       # √−g V^r
        end
        divV = (Vr(r + δ) - Vr(r - δ)) / (2δ) / bg.sqrtg
        QB = [sum(Q[μ, α] * bg.Bup[α] for α in 1:4) for μ in 1:4]
        RV = sum(bg.Rlo[μ] * QB[μ] for μ in 1:4)
        third = -(im * bg.α * ω̃ - (bg.cs2 * divV + RV)) * h0
        princ + second + third
    end
    R = Lval(one(T), zero(T), zero(T))      # L[1]
    Qc = Lval(zero(T), one(T), zero(T))     # L[(r−r₀)]
    P = Lval(zero(T), zero(T), one(T))      # L[(r−r₀)²] (h2=1 ⇒ ∂²-coefficient P)
    (A_m = Qc / P, B_m = R / P, P = P)
end
