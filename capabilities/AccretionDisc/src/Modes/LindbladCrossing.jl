# Regular locally-isothermal crossing of Lindblad resonances. Hand-derived from
# ref-paper/arXiv-2207.14637/src/Sections/2_framework.tex:120-134 (read as text).
# Reuses the MR shooting/refinement and NumericsCore's global C2 spline / DP5.
using NumericsCore: CubicSplineInterpolant

struct RegularPotential{S,I}
    spline::S
    integrator::I
end
RegularPotential(c::PotentialCache; integrator=integrate_grid) =
    RegularPotential(CubicSplineInterpolant(c.lg, c.P),integrator)

function potential_at(c::RegularPotential, r)
    t = log10(r); a = inv(r * log(oftype(r,10)))
    p = c.spline(t,1)
    p1 = c.spline(t,1;deriv_order=1) * a
    p2 = c.spline(t,1;deriv_order=2) * a^2 - p1/r
    p,p1,p2
end

# psi = h + Phi, F = r Sigma u_r. With w = omega-m Omega+i epsilon,
# psi' = 2m Omega psi/(rw) + L_T^-1 h - i D F/(w r Sigma),
# F' = i r Sigma (w h/c_s^2 - m^2 psi/(r^2 w)) - m kappa^2 F/(2 Omega r w).
# D is the landed real Dfreq convention; no coefficient divides by D.
struct RegularRHS{D,M,C}
    d::D
    mode::M
    ih::Int
    cache::C
end
function (g::RegularRHS)(out::AbstractVector, r, y)
    d=g.d; m=_azimuthal_number(g.mode)
    O=Omega(d,r); S=Sigma(d,r); w=omega_pat(d,g.mode)-m*O+im*d.pole_disp
    psi=complex(y[1],y[2]); F=complex(y[3],y[4])
    p=g.ih == 0 ? zero(r) : g.cache.spline(log10(r),1)
    h=psi-p
    dp=2m*O/(r*w)*psi + LTinv(d,r)*h - im*Dfreq(d,r,g.mode)/(w*r*S)*F
    df=im*r*S*(w*h/csiso(d,r)^2-m^2*psi/(r^2*w)) - m*kappa(d,r)^2/(2O*r*w)*F
    out[1]=real(dp); out[2]=imag(dp); out[3]=real(df); out[4]=imag(df)
    out
end
(g::RegularRHS)(r,y) = g(similar(y),r,y)

function integrate_master(d::NewtonianDisc{T}, x, bc, mode, ih, c::RegularPotential;
                          rtol=1e-8,atol=1e-10,dense::Bool=true) where {T}
    ico=corotation_index(d,x,mode); r=x[ico]; m=_azimuthal_number(mode)
    O=Omega(d,r); S=Sigma(d,r); w=omega_pat(d,mode)-m*O+im*d.pole_disp
    D=Dfreq(d,r,mode)
    iszero(D) && throw(ArgumentError("resolve corotation on the grid: the shooting seed lies at D=0"))
    iszero(w) && throw(ArgumentError("a corotation grid point requires nonzero pole_disp"))
    p,p1,_=potential_at(c,r); p*=ih; p1*=ih
    psi=bc[1]+p
    F=im*r*S/D*(w*(bc[2]+p1-LTinv(d,r)*bc[1])-2m*O/r*psi)
    y0=T[real(psi),imag(psi),real(F),imag(F)]
    rhs=RegularRHS(d,mode,ih,c)
    yo=c.integrator(rhs,y0,x,ico,+1;rtol,atol,dense)
    yi=c.integrator(rhs,y0,x,ico,-1;rtol,atol,dense)
    f=zeros(Complex{T},length(x),3); dy=similar(y0)
    for i in eachindex(x)
        y=@view (i>=ico ? yo : yi)[i,:]
        rhs(dy,x[i],y)
        p,p1,_=potential_at(c,x[i])
        f[i,1]=complex(y[1],y[2])-ih*p
        f[i,2]=complex(dy[1],dy[2])-ih*p1
        f[i,3]=complex(y[3],y[4])
    end
    f[ico,1]=bc[1]; f[ico,2]=bc[2]
    f
end

function _regular_vars(d::NewtonianDisc{T}, f, x, mode, cache::PotentialCache) where {T}
    iszero(d.beta) && size(f,2)==3 ||
        throw(ArgumentError("regular reconstruction requires beta=0 and (h,h′,rΣu_r) columns"))
    c=RegularPotential(cache); m=_azimuthal_number(mode)
    ur=zeros(Complex{T},length(x)); up=similar(ur); sig=similar(ur)
    for i in eachindex(x)
        r=x[i]; O=Omega(d,r); S=Sigma(d,r)
        w=omega_pat(d,mode)-m*O+im*d.pole_disp
        p=c.spline(log10(r),1)
        ur[i]=f[i,3]/(r*S)
        up[i]=(m*(f[i,1]+p)/r-im*kappa(d,r)^2*ur[i]/(2O))/w
        sig[i]=S*f[i,1]/csiso(d,r)^2
    end
    ur,up,sig
end
