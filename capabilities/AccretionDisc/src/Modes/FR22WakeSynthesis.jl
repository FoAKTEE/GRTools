using NumericsCore: CubicSplineInterpolant

"""
    FR22Wake(r, modes, density)

Radial density amplitudes for FR22 (12), with `density[i,j] = δΣ_lm(r[i])`
and `modes[j] = NewtonianMode(m,l,n)`. Store a global cubic spline in physical
radius, with real/imaginary rows for each mode. Positive-m modes only; their
real-amplitude convention has weight ONE, unlike SpiralWaves.synthesize.
No mass or background normalization is applied here. Values inherit the input
floating type; no extrapolation is allowed by [`fr22_wake`](@ref).
"""
struct FR22Wake{M,S}
    modes::M
    spline::S
    function FR22Wake(r::AbstractVector, modes, density::AbstractMatrix)
        ms = collect(modes)
        size(density) == (length(r),length(ms)) ||
            throw(DimensionMismatch("density must have radius rows and mode columns"))
        !isempty(ms) && all(m -> m.m > 0,ms) && allunique(ms) ||
            throw(ArgumentError("distinct positive-m modes required"))
        all(x -> isfinite(x) && x > 0,r) && all(isfinite,density) ||
            throw(ArgumentError("finite density and finite positive radii required"))
        components = vcat(permutedims(real.(density)),permutedims(imag.(density)))
        spline = CubicSplineInterpolant(r,components)
        new{typeof(ms),typeof(spline)}(ms,spline)
    end
end

"""
    fr22_wake(wake, disc, orbit, radii, angles, t; frame=:planet, Gmu=1)

Return a `(length(radii), length(angles))` matrix of
`(Mp/Mth)^(-1) δΣ/Σ`. Input angles are radians in `:inertial`, `:guiding`
(phi' = phi - n*t), or `:planet` (phi'' = phi - phi_p(t)) coordinates.
The time phase is the TOTAL mode frequency, `m*phi - mode.omega*t`.

`Gmu` MUST match the mass used to compute the input density amplitudes.
The factor `h_p^3/Gmu` is applied exactly once here, with GMstar=1 and
`h_p = disc.h0 * orbit.a^((1-disc.q)/2)`. Thus unit-Gmu solver output gains
h_p^3; output already computed at Gmu=h_p^3 gains no further mass factor.
FR22 figure reproduction uses a=n=1. No weight-2 stationary convention enters.
Sources: ref-paper/arXiv-2207.14637/src/Sections/2_framework.tex:83-120 (thermal mass and ansatz (12)),
3_numerical_procedure.tex:79; 4_results.tex:10-33 (normalization and frames).
"""
function fr22_wake(wake::FR22Wake, d::NewtonianDisc, orbit::EccentricOrbit,
                   radii::AbstractVector, angles::AbstractVector, t::Real;
                   frame::Symbol=:planet, Gmu::Real=1)
    isfinite(Gmu) && Gmu > 0 || throw(DomainError(Gmu,"input Gmu must be positive and finite"))
    all(r -> first(wake.spline.t) <= r <= last(wake.spline.t),radii) ||
        throw(ArgumentError("requested radius outside the wake spline"))
    shift = if frame == :inertial
        zero(t)
    elseif frame == :guiding
        orbit.n*t
    elseif frame == :planet
        planet_orbit(orbit.n*t,orbit).phi_p
    else
        throw(ArgumentError("frame must be :inertial, :guiding or :planet"))
    end
    T = float(promote_type(eltype(wake.spline.t),eltype(radii),eltype(angles),
                           typeof(t),typeof(orbit.a),typeof(d.h0),typeof(Gmu),
                           mapreduce(m -> typeof(m.omega),promote_type,wake.modes)))
    azimuths = unique(m.m for m in wake.modes)
    amplitude = zeros(Complex{T},length(radii),length(azimuths))
    nm = length(wake.modes)
    for (j,mode) in enumerate(wake.modes)
        k = findfirst(==(mode.m),azimuths)
        phase = cis(T(mode.m)*T(shift)-T(mode.omega)*T(t))
        for (i,r) in enumerate(radii)
            amplitude[i,k] += complex(wake.spline(r,j),wake.spline(r,j+nm))*phase
        end
    end
    angular = [cis(T(m)*T(phi)) for m in azimuths,phi in angles]
    hp = T(d.h0)*T(orbit.a)^((1-T(d.q))/2)
    real.(amplitude*angular) .* (hp^3 ./ (T(Gmu) .* Sigma.(Ref(d),radii)))
end
