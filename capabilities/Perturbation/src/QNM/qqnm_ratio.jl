# QNM/qqnm_ratio.jl — the QQNM amplitude ratio R = A_Q/(A₁A₂).
# The chain R = 2E/(ω_Q²ω₁²ω₂²Σaₙ₁Σaₙ₂) assembles Paper I:987/991/1000/1004 (scri readout,
# Kinnersley excitation, tetrad ½, Ψ₄=−ḧ strain rung); factors are fixed by those cites, not
# tunable. Mirror assembly: 2410 main.tex:310. Child basis is spheroidal at aω_Q; the fit
# tabs are SPHERICAL (2410 supp:226; identical at a=0).

"""
    qqnm_ratio(parent1, parent2, channel, l_Q; a=0, n=64, T=Complex{BigFloat}, N_leaver=400) -> NamedTuple

The single-channel amplitude ratio `R_{c₁c₂}` for the excitation `channel ∈ {:pp,:pm,:mp,:mm}` via the
frozen chain (file header): `R = 2E/(ω_Q²ω₁²ω₂²Σaₙ₁Σaₙ₂)`.  `−−` returns `R=0` structurally
(G-Q4H-5).  Returns `(; R, E, omega_Q, cond, Sigma_an1, Sigma_an2, mag, phase, resonance_ok)`.
The ++ channel is DERIVED at general a (Q5-STAGE-1); +− Kerr is [FUTURE].  {T}/BigFloat-inherited.
"""
function qqnm_ratio(parent1, parent2, channel::Symbol, l_Q::Integer;
                    a::Real=0, n::Integer=64, T::Type=Complex{BigFloat}, N_leaver::Integer=400,
                    hertz_scale::Number=1)
    C = T <: Complex ? T : Complex{T}
    sol = qqnm_solve(parent1, parent2, channel, l_Q; a=a, n=n, T=C, N_leaver=N_leaver,
                     hertz_scale=hertz_scale)
    if channel in (:mm, :minusminus)
        return (; R=zero(C), E=zero(C), omega_Q=sol.omega_Q, cond=sol.cond,
                Sigma_an1=zero(C), Sigma_an2=zero(C), mag=zero(real(C)), phase=zero(real(C)),
                resonance_ok=true, zero=qqnm_mm_is_zero())
    end
    ω1 = C(parent1.ω); ω2 = C(parent2.ω); ω_Q = sol.omega_Q
    # per-parent linear normalizations Σaₙ (s=+2 Hertz Leaver); scaled by `hertz_scale` consistently
    # with the source E (G-Q4H-2: E→c²E, Σ_j→cΣ_j ⇒ R bit-stable).
    _, Σ1 = qqnm_hertz_leaver(parent1.l, parent1.m, a, ω1; N=N_leaver, T=C)
    _, Σ2 = qqnm_hertz_leaver(parent2.l, parent2.m, a, ω2; N=N_leaver, T=C)
    Σ1 *= C(hertz_scale); Σ2 *= C(hertz_scale)
    R = 2 * sol.E / (ω_Q^2 * ω1^2 * ω2^2 * Σ1 * Σ2)
    return (; R=R, E=sol.E, omega_Q=ω_Q, cond=sol.cond, Sigma_an1=Σ1, Sigma_an2=Σ2,
            mag=abs(R), phase=angle(R), resonance_ok=sol.resonance_ok)
end

"""
    qqnm_ratio_np(parent1, parent2, l_Q; a=0, n=64, T=Complex{BigFloat}, N_leaver=400) -> NamedTuple

The nonprecessing total `R_np = R_++ + (−1)^{ℓ₂}R_+− + (−1)^{ℓ₁}R_−+` (R_−−=0; 2410 main.tex:310).
For identical parents the `−+` row equals `+−` under exchange and the certified assembly is
`R_np = R_++ + (−1)^{ℓ₂}R_+−` (owner channel-fit check).  Returns the per-channel R's + the total.
"""
function qqnm_ratio_np(parent1, parent2, l_Q::Integer;
                       a::Real=0, n::Integer=64, T::Type=Complex{BigFloat}, N_leaver::Integer=400)
    C = T <: Complex ? T : Complex{T}
    ℓ1 = parent1.l; ℓ2 = parent2.l
    rpp = qqnm_ratio(parent1, parent2, :pp, l_Q; a=a, n=n, T=C, N_leaver=N_leaver)
    rpm = qqnm_ratio(parent1, parent2, :pm, l_Q; a=a, n=n, T=C, N_leaver=N_leaver)
    identical = (parent1.l == parent2.l && parent1.m == parent2.m &&
                 parent1.n == parent2.n && parent1.ω == parent2.ω)
    Rnp = rpp.R + (-1)^ℓ2 * rpm.R
    if !identical
        rmp = qqnm_ratio(parent1, parent2, :mp, l_Q; a=a, n=n, T=C, N_leaver=N_leaver)
        Rnp += (-1)^ℓ1 * rmp.R
        return (; R_np=Rnp, R_pp=rpp.R, R_pm=rpm.R, R_mp=rmp.R,
                mag=abs(Rnp), phase=angle(Rnp), mag_pp=rpp.mag, phase_pp=rpp.phase)
    end
    return (; R_np=Rnp, R_pp=rpp.R, R_pm=rpm.R, mag=abs(Rnp), phase=angle(Rnp),
            mag_pp=rpp.mag, phase_pp=rpp.phase)
end
