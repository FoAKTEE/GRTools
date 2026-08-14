# LorenzWeyl — gauge-invariant ψ4 of the a≠0 Lorenz reconstruction:
# ψ4 = ρ⁴·R_{-2}(r)·{}_{-2}S_{lm}(θ) — the (r,θ)-independent ratio vs the CCK/Teukolsky ψ4 is the
# reconstruction-identity cross-check.

import Harmonics: SpinWeightedSpheroidalHarmonicS

export weyl_psi4_lorenz, weyl_psi4_lorenz_amp

# `jd_l` (spheroidal-l of a reconstruct_kerr result) is owned by KerrSMatrices, included earlier in this module.

"""
    weyl_psi4_lorenz(rec, r, θ; bc=:Up, l=nothing) -> Complex

The (un-normalized) gauge-invariant ψ4 radial×angular mode of the LORENZ reconstruction at (r,θ): the s=−2
Teukolsky Weyl scalar ψ4 = ρ⁴·R_{-2}(r)·{}_{-2}S_{lm}(θ), ρ=−1/(r−i a cosθ). `bc` selects the branch
(:Up r>r0, :In r<r0); `rec` is a `reconstruct_kerr` result. Pass `l` when `rec` carries no `_l` metadata.
"""
function weyl_psi4_lorenz(rec, r, θ; bc::Symbol=:Up, l=nothing)
    jd = rec.jd; a = jd.a; m = jd.m; ω = jd.ω
    l = l === nothing ? jd_l(rec) : l
    br = rec.Pm2.br
    Rval = bc===:Up ? rec.Pm2.Aup*br.up(r) : rec.Pm2.Ain*br.inn(r)
    Sm2 = real(SpinWeightedSpheroidalHarmonicS(-2, l, m, Float64(a*ω))(θ, 0.0))
    ρ = -1/(r - im*a*cos(θ))
    ρ^4 * Rval * Sm2
end

"""
    weyl_psi4_lorenz_amp(rec, r; bc=:Up) -> Complex

The branch amplitude×radial part R_{-2}(r) of ψ4 (angular factor stripped) — for forming the (r,θ)-independent
ratio against the CCK/Teukolsky ψ4 radial part.
"""
function weyl_psi4_lorenz_amp(rec, r; bc::Symbol=:Up)
    br = rec.Pm2.br
    bc===:Up ? rec.Pm2.Aup*br.up(r) : rec.Pm2.Ain*br.inn(r)
end
