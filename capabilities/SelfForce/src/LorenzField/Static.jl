# Analytic static first-order fields for the h1-corpus generator: the l=m=0 closed-form static
# monopole (no integration needed). Port of the static monopole In/Out closed forms of
# ref-code/h1Lorenz/src/static_monopole.c (legacy H1Lorenz/static.jl). Type-generic {T}, M=1; leaf
# of the LorenzField solver chain.

# Static monopole (l=0,m=0) closed forms — exact algebraic/log expressions in (r0, r). M=1.
_h1MonoIn(r0, r)  = (-4*sqrt(pi)*(-2+r)^3*(2+r)*(-4+(-3+r0)*log((-2+r0)/r0)))/(3*r^3*sqrt((-3+r0)*r0))
_h3MonoIn(r0, r)  = (8*sqrt(pi)*(-32+r^3-(-8+r^3)*(-3+r0)*log((-2+r0)/r0)))/(3*r^2*sqrt((-3+r0)*r0))
_h6MonoIn(r0, r)  = (4*sqrt(pi)*(-64+2*r^3+(16+r^3)*(-3+r0)*log((-2+r0)/r0)))/(3*r^2*sqrt((-3+r0)*r0))
_h1MonoOut(r0, r) = (4*sqrt(pi)*(-2*(-8+r0)*(-4+r0)-6*r^2*(-3+r0)+4*r^3*(-3+r0)+2*r*(20+(-8+r0)*r0)-(-4+r)*r^3*(-3+r0)*log((-2+r)/r)+16*(-1+r)*(-3+r0)*log(r0/(-2+r))))/(3*r^3*sqrt((-3+r0)*r0))
_h3MonoOut(r0, r) = (-8*sqrt(pi)*(32+3*r*(4+r)-12*r0-r*(4+r)*r0+r0^2+r^3*(-3+r0)*log((-2+r)/r)+8*(-3+r0)*log(r0/(-2+r))))/(3*r^2*sqrt((-3+r0)*r0))
_h6MonoOut(r0, r) = (4*sqrt(pi)*(-2*(-8+r0)*(-4+r0)+8*r*(-3+r0)+2*r^2*(-3+r0)+r^3*(-3+r0)*log((-2+r)/r)+16*(-3+r0)*log((-2+r)/r0)))/(3*r^2*sqrt((-3+r0)*r0))

_dh1MonoIn(r0, r) = (-4*sqrt(pi)*(-2+r)^3*(-4+(-3+r0)*log((-2+r0)/r0)))/(3*r^3*sqrt((-3+r0)*r0))+(4*sqrt(pi)*(-2+r)^3*(2+r)*(-4+(-3+r0)*log((-2+r0)/r0)))/(r^4*sqrt((-3+r0)*r0))-(4*sqrt(pi)*(-2+r)^2*(2+r)*(-4+(-3+r0)*log((-2+r0)/r0)))/(r^3*sqrt((-3+r0)*r0))
_dh3MonoIn(r0, r) = (8*sqrt(pi)*(3*r^2-3*r^2*(-3+r0)*log((-2+r0)/r0)))/(3*r^2*sqrt((-3+r0)*r0))-(16*sqrt(pi)*(-32+r^3-(-8+r^3)*(-3+r0)*log((-2+r0)/r0)))/(3*r^3*sqrt((-3+r0)*r0))
_dh6MonoIn(r0, r) = (4*sqrt(pi)*(6*r^2+3*r^2*(-3+r0)*log((-2+r0)/r0)))/(3*r^2*sqrt((-3+r0)*r0))-(8*sqrt(pi)*(-64+2*r^3+(16+r^3)*(-3+r0)*log((-2+r0)/r0)))/(3*r^3*sqrt((-3+r0)*r0))
_dh1MonoOut(r0, r) = (4*sqrt(pi)*((-16*(-1+r)*(-3+r0))/(-2+r)-12*r*(-3+r0)+12*r^2*(-3+r0)-((-((-2+r)/r^2)+1/r)*(-4+r)*r^4*(-3+r0))/(-2+r)+2*(20+(-8+r0)*r0)-3*(-4+r)*r^2*(-3+r0)*log((-2+r)/r)-r^3*(-3+r0)*log((-2+r)/r)+16*(-3+r0)*log(r0/(-2+r))))/(3*r^3*sqrt((-3+r0)*r0))-(4*sqrt(pi)*(-2*(-8+r0)*(-4+r0)-6*r^2*(-3+r0)+4*r^3*(-3+r0)+2*r*(20+(-8+r0)*r0)-(-4+r)*r^3*(-3+r0)*log((-2+r)/r)+16*(-1+r)*(-3+r0)*log(r0/(-2+r))))/(r^4*sqrt((-3+r0)*r0))
_dh3MonoOut(r0, r) = (-8*sqrt(pi)*(3*r+3*(4+r)-(8*(-3+r0))/(-2+r)+((-((-2+r)/r^2)+1/r)*r^4*(-3+r0))/(-2+r)-r*r0-(4+r)*r0+3*r^2*(-3+r0)*log((-2+r)/r)))/(3*r^2*sqrt((-3+r0)*r0))+(16*sqrt(pi)*(32+3*r*(4+r)-12*r0-r*(4+r)*r0+r0^2+r^3*(-3+r0)*log((-2+r)/r)+8*(-3+r0)*log(r0/(-2+r))))/(3*r^3*sqrt((-3+r0)*r0))
_dh6MonoOut(r0, r) = (4*sqrt(pi)*(8*(-3+r0)+(16*(-3+r0))/(-2+r)+4*r*(-3+r0)+((-((-2+r)/r^2)+1/r)*r^4*(-3+r0))/(-2+r)+3*r^2*(-3+r0)*log((-2+r)/r)))/(3*r^2*sqrt((-3+r0)*r0))-(8*sqrt(pi)*(-2*(-8+r0)*(-4+r0)+8*r*(-3+r0)+2*r^2*(-3+r0)+r^3*(-3+r0)*log((-2+r)/r)+16*(-3+r0)*log((-2+r)/r0)))/(3*r^3*sqrt((-3+r0)*r0))

"""
    construct_static_monopole(orbit::OrbitalParams{T}) -> Array{T,3}

Analytic asymptotically-flat static monopole (l=m=0) on the orbit's radial grid. Returns
`inhom[k, part, grididx]` with `k = 1:4` the three coupled monopole fields `h1,h3,h6` plus a zero
gauge slot (R2 ≡ 0); `part = 1:4 = (Re, Im, Re∂_r, Im∂_r)` (all imaginary parts 0); `grididx =
1:gridsize+1`, where the last index holds the right-limit at the worldline `r0`. The ingoing branch
(`grididx ≤ r0_index+1`) uses the In closed forms, the rest the Out forms. Port of
`ref-code/h1Lorenz/src/static_monopole.c`, type-generic `{T}` (M=1).
"""
function construct_static_monopole(orbit::OrbitalParams{T}) where {T}
    r0 = orbit.p; grid = orbit.grid; gs = length(grid); r0i = orbit.r0_index
    nfg = 3 + 1                              # 3 coupled monopole fields + 1 (zero) gauge slot
    inhom = zeros(T, nfg, 4, gs + 1)
    for i in 0:r0i
        r = grid[i+1]
        inhom[1,1,i+1] = _h1MonoIn(r0,r); inhom[1,3,i+1] = _dh1MonoIn(r0,r)
        inhom[2,1,i+1] = _h3MonoIn(r0,r); inhom[2,3,i+1] = _dh3MonoIn(r0,r)
        inhom[3,1,i+1] = _h6MonoIn(r0,r); inhom[3,3,i+1] = _dh6MonoIn(r0,r)
    end
    for i in (r0i+1):(gs-1)
        r = grid[i+1]
        inhom[1,1,i+1] = _h1MonoOut(r0,r); inhom[1,3,i+1] = _dh1MonoOut(r0,r)
        inhom[2,1,i+1] = _h3MonoOut(r0,r); inhom[2,3,i+1] = _dh3MonoOut(r0,r)
        inhom[3,1,i+1] = _h6MonoOut(r0,r); inhom[3,3,i+1] = _dh6MonoOut(r0,r)
    end
    # worldline right-limit (Out at r0) stored at the last index gs+1
    inhom[1,1,gs+1] = _h1MonoOut(r0,r0); inhom[1,3,gs+1] = _dh1MonoOut(r0,r0)
    inhom[2,1,gs+1] = _h3MonoOut(r0,r0); inhom[2,3,gs+1] = _dh3MonoOut(r0,r0)
    inhom[3,1,gs+1] = _h6MonoOut(r0,r0); inhom[3,3,gs+1] = _dh6MonoOut(r0,r0)
    return inhom
end
