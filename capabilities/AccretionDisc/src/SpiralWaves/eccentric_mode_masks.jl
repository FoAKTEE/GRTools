"""
    eccentric_mode_mask(orbit; policy=:fr24, mmax=nothing, kmax=nothing,
                        mask=:paper, anchor_modes=nothing,
                        secular_rtol=64eps(typeof(orbit.p)))

One mode inventory for paper comparisons and physical sums. Paper membership
is the integer rule `m >= 1, m+k >= 1`, independent of the sign of the Kerr
frequency. Seeds are FR22 `(150,8)`, FR24 `(150/175,40)` with the switch at
`e > 0.06`, and FD25 `(170,40)`; circular seeds have `kmax=0`.
These are truncation masks, with no amplitude or exterior-resonance pruning.
Pass the actual Newtonian `anchor_modes` to transport additional paper-specific
cuts (especially FD25's exterior-resonance policy) by the exact integer pairs.
They must be unique m>=1,l>=1 modes inside this seed. Their Newtonian frequencies
are ignored. Omitted positive-l modes remain in `sectors.anchor_excluded`,
so an anchor cut never silently deletes a physical GR mode. Without that list
the FD25 mask is its seed mask, not a certificate of its domain-dependent cuts.
Sources (read as text under emri-disk/ref-paper):
`arXiv-2207.14637/src/Sections/3_numerical_procedure.tex:72`,
`arXiv-2407.20398/src/Sections/2_linear_framework.tex:192,204-213`,
`arXiv-2407.20398/src/Sections/4_eccentric_torques.tex:19`,
`arXiv-2506.19917/src/Sections/2_linear_theory.tex:189`.
Further domain exclusions need a domain.

`paper` preserves the anchor mask exactly. `physical` additionally contains
the excluded `m+k <= 0` modes and the independent reality representatives
`m=0,k=1:kmax`; it omits only the stationary axisymmetric background `(0,0)`.
`selected` is `paper` or `physical` according to `mask`. The five disjoint
`sectors` partition `physical`: `paper`, `excluded`, `secular`, `axisymmetric`,
`anchor_excluded` (empty unless the supplied anchor applies additional cuts).
The secular sector takes precedence over paper/excluded membership: `(1,-1)`
is always conditional per the design, as are modes whose computed frequency
is unresolved relative to its summands,
`abs(omega) <= secular_rtol*max(abs(m*Omega_phi),abs(k*Omega_r))`.
It is a routing guard, not a damping prescription or a physical smallness test.
Near-zero total frequency alone is not an operator pole for m>0; callers may
certify these modes separately. Nothing in the inventory deletes their torque.
"""
function eccentric_mode_mask(orbit::EquatorialOrbit{T}; policy::Symbol=:fr24,
                             mmax::Union{Nothing,Integer}=nothing,
                             kmax::Union{Nothing,Integer}=nothing,
                             mask::Symbol=:paper, anchor_modes=nothing,
                             secular_rtol=64eps(T)) where {T}
    policy in (:fr22,:fr24,:fd25) || throw(ArgumentError("policy must be :fr22, :fr24 or :fd25"))
    mask in (:paper,:physical) || throw(ArgumentError("mask must be :paper or :physical"))
    M = something(mmax,policy === :fd25 ? 170 :
                  (policy === :fr24 && orbit.e > T(3)/50 ? 175 : 150))
    K = something(kmax,iszero(orbit.e) ? 0 : (policy === :fr22 ? 8 : 40))
    M >= 0 && K >= 0 || throw(ArgumentError("mode truncations must be nonnegative"))
    isfinite(secular_rtol) && secular_rtol >= 0 || throw(ArgumentError("secular_rtol must be finite and nonnegative"))
    nonaxisymmetric = [HarmonicMode(m,k,orbit_frequency(orbit,m,k)) for m in 1:M for k in -K:K]
    paper = filter(s -> s.m+s.k >= 1,nonaxisymmetric)
    anchor_pairs = Set((s.m,s.k) for s in (anchor_modes === nothing ? paper : anchor_modes))
    if anchor_modes !== nothing
        length(anchor_pairs)==length(anchor_modes) &&
            all(t->1<=t[1]<=M && abs(t[2])<=K && t[1]+t[2]>=1,anchor_pairs) ||
            throw(ArgumentError("anchor modes must be unique m>=1,l>=1 pairs inside the seed"))
        filter!(s -> (s.m,s.k) in anchor_pairs,paper)
    end
    axisymmetric = [HarmonicMode(0,k,orbit_frequency(orbit,0,k)) for k in 1:K]
    physical = vcat(nonaxisymmetric,axisymmetric)
    sectors = (;paper=HarmonicMode{T}[],excluded=HarmonicMode{T}[],
                secular=HarmonicMode{T}[],axisymmetric=HarmonicMode{T}[],anchor_excluded=HarmonicMode{T}[])
    for s in physical
        scale = max(abs(s.m*orbit.Ωφ),abs(s.k*orbit.Ωr))
        sector = (s.m==1 && s.k == -1) || abs(s.omega) <= secular_rtol*scale ? :secular :
                 s.m==0 ? :axisymmetric : s.m+s.k<=0 ? :excluded :
                 (s.m,s.k) in anchor_pairs ? :paper : :anchor_excluded
        push!(getproperty(sectors,sector),s)
    end
    (;selected=mask === :paper ? paper : physical,paper,physical,sectors,
      policy,mmax=M,kmax=K,secular_rtol=T(secular_rtol))
end
