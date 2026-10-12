"""
    synthesize_eccentric(sols, rq, phiq, t; component=1)

Real field on the inertial polar grid `(r, phi)` at coordinate time `t`:
`2 sum(real(f_mk(r) * exp(im*(m*phi - omega_mk*t))))`.
Each entry carries `context.mode`, `r`, and a two-row amplitude matrix `y`;
an `EccentricModeSolution` supplies enthalpy (component 1) or radial velocity
master variable (component 2). Density amplitudes must first be reconstructed
by the closure owner. The time origin is the source's pericentre convention.

Use one representative per reality pair: m>0 with any k, or dynamic m=0,k>0.
The pair factor is the diagnostic owner's C-NORM weight divided by 2pi;
in particular a dynamic axisymmetric pair is doubled too. Static (0,0) is
excluded. Signed frequencies are supported. Radial interpolation and zero
extension reuse `mode_interp`; entries accumulate in the supplied order.
"""
function synthesize_eccentric(sols, rq::AbstractVector, phiq::AbstractVector, t::Real;
                              component::Integer=1)
    isempty(sols) && throw(ArgumentError("synthesize_eccentric: no modes given"))
    T = float(promote_type(eltype(rq),eltype(phiq),typeof(t)))
    for sol in sols
        T = promote_type(T,eltype(sol.r),typeof(sol.context.mode.omega))
    end
    out = zeros(T,length(rq),length(phiq))
    for sol in sols
        mode = sol.context.mode
        weight = reality_pair_weight(mode)/(2oftype(mode.omega,pi))
        mi = mode_interp((;m=mode.m,r=sol.r,y=sol.y);component)
        phases = [cis(T(mode.m)*T(phi)-T(mode.omega)*T(t)) for phi in phiq]
        for (i,r) in enumerate(rq)
            amplitude = mi(r)
            iszero(amplitude) && continue
            for (j,phase) in enumerate(phases)
                out[i,j] += weight*real(amplitude*phase)
            end
        end
    end
    out
end
