# MModeOperator/hyperboloidal.jl — hyperboloidal slicing KINEMATICS (ESSF node
# `hyperboloidal-transformation`, W5-C; arXiv:2606.04998 Sec. II.B.2, main.tex:255-262
# tortoise chart + 267-314 height/boost functions and the operator conjugation). The
# slicing is a radially-dependent shift of the time coordinate, s = t − h(r_*)
# (main.tex:270-272), fully described by the height function h and its tortoise
# derivatives — the boost H ≡ dh/dr_* and H′ ≡ dH/dr_* (main.tex:274-279).
#
# SCOPE BOUNDARIES (predecessors own the rest; this file only adds the kinematics):
#   - ansatz.jl (node mode-ansatz-decomposition) owns the phase map Ψ_m = Ψ̃_m e^{−imΩh}
#     and its inverse (eq:toTilde, `psi_from_tilde`/`tilde_from_psi`);
#   - tsliced.jl (node elliptic-operator-tsliced) owns the flux-form coefficients with
#     caller-supplied pointwise (H, H′) slots (replacement rules main.tex:306-308).
# Composing the three layers realizes the transformed operator △_m: coefficients at
# (H(r_*(r)), H′(r_*(r))) + the e^{−imΩh} phase. The closure oracle (test file, gate (e))
# is the CONJUGATION IDENTITY △_m[Ψ̃ e^{−imΩh}] = e^{−imΩh} △̃_m[Ψ̃], exact because the
# replacement rules ARE conjugation by the phase: (∂_{r_*} + imΩH)[Ψ̃e^{−imΩh}] =
# e^{−imΩh} ∂_{r_*}Ψ̃ term by term (main.tex:291-309). The vtu piecewise instance
# (H ∈ {−1, 0, +1}, main.tex:315-330) and its junction conditions are the
# vtu-slicing-jumps node; t-slicing is `t_slicing()` here.

"""
    rstar(orbit, r) -> r_*

Kerr tortoise coordinate r_*(r) of the slicing kinematics (main.tex:255-258), a thin wrap
of `Spacetime.tortoise_kerr` (whose M = 1 log form and integration constant are exactly
the printed ones):

    r_* = r + (2M/(r₊−r₋))·[ r₊ ln((r−r₊)/2M) − r₋ ln((r−r₋)/2M) ] ,   M = 1,

defined by dr_*/dr = (r²+a²)/Δ (main.tex:260-262) — the inverse of the radial flux
coefficient: `flux_coeff_r(orbit, r)` = A^rr = Δ/(r²+a²) IS dr/dr_*. r_* → −∞
logarithmically at the horizon and ~r at infinity; t ± r_* are the null pairs that make
this the natural argument of the height function (main.tex:265). At a = 0 Spacetime
delegates to `tortoise_schw` (identical limit, r₋ = 0). Type-generic (BigFloat inherited).
"""
rstar(orbit::CircularOrbit, r) = tortoise_kerr(r, orbit.a)

"""
    HeightFunction(h, H, Hprime)

Slicing-kinematics carrier: the height function of the hyperboloidal transformation
s = t − h(r_*) (main.tex:270-272) together with its two tortoise derivatives, each a
callable of r_* (NOT of r — compose with [`rstar`](@ref) via
[`slicing_kinematics`](@ref)):

  - `h(r_*)`      — the height entering the phase map Ψ_m = Ψ̃_m e^{−imΩh} (ansatz.jl);
  - `H(r_*)`      — the boost H = dh/dr_* (main.tex:274-277), the (H, ·) slot of the
                    tsliced.jl coefficients;
  - `Hprime(r_*)` — H′ = dH/dr_*, the (·, H′) slot (β_m's −imΩH′ term).

INVARIANT (the boost/height relation): `H` must be the exact r_*-derivative of `h`, and
`Hprime` of `H` — the operator conjugation identity △_m[Ψ̃e^{−imΩh}] = e^{−imΩh}△̃_m[Ψ̃]
holds ONLY for consistent triples (FD-gated in the node tests; the constructor cannot
check closures). The carrier does not enforce the asymptotic hyperboloidal condition
H(±∞) = ±1 (main.tex:278) — single slabs like [`t_slicing`](@ref) (H ≡ 0, the source
region, main.tex:280) are valid pieces; the full vtu assembly owns the condition. The
callables must be type-generic so BigFloat precision is inherited.
"""
struct HeightFunction{Fh,FH,FHp}
    h::Fh
    H::FH
    Hprime::FHp
end

"""
    t_slicing() -> HeightFunction

The trivial slicing h ≡ 0 (so s = t, H = H′ = 0): the t-domain of vtu slicing, enforced
in the source region (main.tex:280, 315-330). Its fields are `Base.zero`, so every
evaluation returns the exact zero OF THE ARGUMENT'S TYPE — through
[`slicing_kinematics`](@ref) the tsliced.jl coefficients and the ansatz.jl phase map
reproduce their t-sliced defaults bitwise (node gate (c)).
"""
t_slicing() = HeightFunction(zero, zero, zero)

"""
    slicing_kinematics(hf, orbit, r) -> (h, H, Hprime)

Evaluate the slicing kinematics at Boyer–Lindquist radius `r`: one r_* = [`rstar`](@ref)
evaluation composed into the [`HeightFunction`](@ref) triple. Downstream operator
assembly consumes the result directly — `h` feeds `psi_from_tilde`/`tilde_from_psi`
(ansatz.jl), `(H, Hprime)` feed `beta_m`/`gamma_r_m` (tsliced.jl).
"""
function slicing_kinematics(hf::HeightFunction, orbit::CircularOrbit, r)
    rs = rstar(orbit, r)
    (h = hf.h(rs), H = hf.H(rs), Hprime = hf.Hprime(rs))
end
