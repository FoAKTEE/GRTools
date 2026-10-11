# OPTION A radiation boundaries. Include after locally_isothermal.jl and boundary.jl.
# The operator and its thermal gradients come from generated/eccentric_li.jl;
# no adiabatic wave-flux formula or circular pattern frequency enters this file.

"""
    li_boundary(C, Cω, Cr; outward)

First-order WKB boundary data for the real, unregularized LI system `y′=C*y`,
`y=(δp/(ε+p),-iδu^r)`. `Cω` differentiates total frequency at fixed disc/radius/m;
`Cr` differentiates radius at fixed total frequency. The existing `wkb_split`
selects the outward group velocity, or exterior decay for real roots.

With `C=[a b; c d]`, normalize each eigenvector to `Y=(1,z)`,
`z=(λ-a)/b`. Differentiating its characteristic equation gives
`λ′=(tr(C)′*λ-det(C)′)/(2λ-tr(C))`. Hence
`z′=(λ′-a′-b′*z)/b`. For `gap=λ_t-λ_a`, the first transport correction is
`β=b*z′/gap²`, and the corrected state is `Y_a+β*Y_t`. The returned `row`
annihilates that state by calling `radiation_impedance_row`.
In `y≈A*Y_a*exp(∫λ_a dr)`, `dlog(A)/dr=transport=b*z′/gap` to first order;
`logderivative=λ_a+transport`. The impedance implied by `row` differs from
that truncated logarithmic derivative only at second WKB order.

This is a homogeneous, local asymptotic condition, not an exact exterior map.
Subtract a particular solution if the exterior is forced. Require separated
roots, nonzero `b` (pressure chart), and slow variation relative to the gap;
`mixing=β` exposes the expansion parameter. No uniform turning-point or
corotation certificate is implied. Grid refinement does not remove WKB error.

Newtonian limit: with `L=(log(c_T²))′`, `D=κ²-(ω-mΩ)²`, the pressure amplitude
gradient is `D′/(4D)-Σ′/(2Σ)-1/(2r)+3L/4`, including the extra `L/2` in
Modes/NewtonianPlanet.jl:114–129. It implements the closure-dependent amplitude
prescribed in MR20 src/main.tex:301 (footnote) and :404–419 (F_J/c_T² law).
"""
function li_boundary(C::AbstractMatrix, Cω::AbstractMatrix, Cr::AbstractMatrix;
                     outward::Bool)
    size(Cr) == (2,2) || throw(DimensionMismatch("radial derivative must be 2×2"))
    all(isfinite,Cr) && all(isreal,Cr) ||
        throw(ArgumentError("LI boundary needs a finite real radial derivative"))
    λa, _, λt, _ = wkb_split(C,Cω;outward)
    a,b,c,d = C[1,1],C[1,2],C[2,1],C[2,2]
    iszero(b) && throw(DomainError(b,"LI pressure chart is singular: move the edge"))
    ar,br,cr,dr = Cr[1,1],Cr[1,2],Cr[2,1],Cr[2,2]
    λr = ((ar+dr)*λa-(ar*d+a*dr-br*c-b*cr))/(λa-λt)
    za,zt = (λa-a)/b,(λt-a)/b
    zr = (λr-ar-br*za)/b
    Ya,Yt = (one(za),za),(one(zt),zt)
    row = radiation_impedance_row(λa,Ya,λt,Yt,(zero(zr),zr))
    transport = b*zr/(λt-λa)
    mixing = transport/(λt-λa)
    vector = Ya .+ mixing .* Yt
    (;row,vector,eigenvalue=λa,other_eigenvalue=λt,transport,
      logderivative=λa+transport,mixing)
end

"""
    li_boundary(context::ModeContext, r; outward, radial_step=nothing)

Evaluate the emitted OPTION A operator and its edge derivatives, then return
`li_boundary(C,Cω,Cr;outward)`. Landau regularization is always disabled here.
Frequency differentiation uses a type-inherited complex step of the emitted
arithmetic at fixed background. Radial differentiation is second order, centered
when possible and one-sided inside the context domain at its edges. The default
step is `cbrt(eps(T))*max(abs(r),1)`, bounded by one quarter of the domain width.
Override `radial_step` to audit derivative convergence. A stencil crossing
corotation is rejected. BigFloat inputs inherit the current working precision.

For `solve_mode`, pass `boundary_rows=((bin.row,0),(bout.row,0))` explicitly,
where `bin=li_boundary(context,r_in;outward=false)` and likewise for `bout`.
No orbit, corotation location, sign of m, or sign of ω sets the edge normal.
"""
function li_boundary(c::ModeContext{T}, r; outward::Bool, radial_step=nothing) where {T}
    c.closure == LocallyIsothermalOptionA ||
        throw(ArgumentError("li_boundary requires locally isothermal OPTION A"))
    x = master_inputs(c,r;hmodes=_zero_source(c),landau=false)
    S = typeof(x.enu)
    rr = S(r)
    nt = (;(n => getfield(x,n) for n in fieldnames(typeof(x)))...)
    C = _context_operator(Val(c.closure),x)
    hω = sqrt(eps(S))*abs(real(x.omega_d))
    Cω = imag.(_context_operator(Val(c.closure),merge(nt,(omega_d=x.omega_d+im*hω,))))/hω
    lo,hi = c.domain
    h = radial_step === nothing ? min(cbrt(eps(S))*max(abs(rr),one(S)),(hi-lo)/4) : S(radial_step)
    isfinite(h) && h>0 && 2h<=hi-lo ||
        throw(ArgumentError("radial_step must be finite, positive and fit inside the domain"))
    radii,weights = if rr-h < lo
        rr+2h<=hi || throw(ArgumentError("radial stencil leaves domain"))
        ((rr,rr+h,rr+2h),(-3,4,-1))
    elseif rr+h > hi
        rr-2h>=lo || throw(ArgumentError("radial stencil leaves domain"))
        ((rr,rr-h,rr-2h),(3,-4,1))
    else
        ((rr,rr-h,rr+h),(0,-1,1))
    end
    c.r_c !== nothing && minimum(radii)<=c.r_c<=maximum(radii) &&
        throw(DomainError(c.r_c,"LI boundary stencil crosses corotation"))
    Cr = (weights[1]*C+weights[2]*master_operator(c,radii[2];landau=false)+
          weights[3]*master_operator(c,radii[3];landau=false))/(2h)
    li_boundary(C,Cω,Cr;outward)
end
