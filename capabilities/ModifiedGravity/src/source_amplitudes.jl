# Per-(l,m,n,k) effective-source amplitudes (Sourcesamplitudes, Sourceincoordinates.nb SC l.1114,
# Polcar & Witzany 2025). JOIN node composing hertz_potential_mode ⊕ source_kernels ⊕ Harmonics.
#
# PRESERVE verbatim: σ=(-1)^(l+kloc) minus-mode conjugation, the (mi→-m, ω→-ω) substitution on
# the −leg, and Snorm≡0 (upstream substitutes a literal 0 for the normalisation leg).

using Harmonics: SpinWeightedSphericalHarmonicY

"""
    source_amplitudes(l, m, n, k, a, p, e, x; Nmax=80) -> (Sl, Sn, Sm, Smdag, omega, Snorm)

Tetrad-leg amplitudes of the effective ψ4 source for the `(l,m,n,k)` mode on a Kerr geodesic
`(a,p,e,x)` (M=1), porting `Sourcesamplitudes`. `Sl,Sn,Sm,Smdag,Snorm` are callables of
Boyer-Lindquist radius `r`; `omega` is the +mode frequency.

Each leg is the sum of a +ω and a −ω contribution (SC `Substplus`/`Substminus`):
`S(r) = K₊(1,r,ω,m,l, RHz(r),RHz'(r), yv) + K₋(1,r,−ω,−m,l, σ·RHz(r),σ·RHz'(r), yv)`,
with `σ = (−1)^(l+kloc)` the minus-mode conjugation sign (`kloc = 0` if `x=1` else `k`),
`RHz`/`RHz'` from [`hertz_potential_mode`](@ref), `K±` the source kernels, and
`yv(s,m) = {}_sY_{l,m}(π/2)` (zero when `|s|>l`). `Snorm ≡ 0`, matching upstream.
"""
function source_amplitudes(l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real; Nmax::Integer=80)
    kloc = x == 1 ? 0 : k
    RadHZ, dRadHZ, ω = hertz_potential_mode(l, m, n, k, a, p, e, x; Nmax=Nmax)
    σ = (-1)^(l + kloc)
    yv(s, mm) = abs(s) > l ? zero(ComplexF64) :
                complex(SpinWeightedSphericalHarmonicY(s, l, mm, π/2, 0.0))
    leg(Kplus, Kminus) = function (r::Real)
        rhz = RadHZ(r); drhz = dRadHZ(r)
        return Kplus(1, r, ω, m, l, rhz, drhz, yv) +
               Kminus(1, r, -ω, -m, l, σ * rhz, σ * drhz, yv)
    end
    return (Sl    = leg(_source_kernel_l_plus,    _source_kernel_l_minus),
            Sn    = leg(_source_kernel_n_plus,    _source_kernel_n_minus),
            Sm    = leg(_source_kernel_m_plus,    _source_kernel_m_minus),
            Smdag = leg(_source_kernel_mdag_plus, _source_kernel_mdag_minus),
            omega = ω,
            Snorm = (r::Real) -> zero(ComplexF64))
end

export source_amplitudes
