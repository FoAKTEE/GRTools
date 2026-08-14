"""
    MetricReconstruction.Lorenz.LorenzRadial

Canonical home of the a=0 first-order Lorenz-gauge metric-perturbation radial field-equation RHS
systems (Barack–Lousto–Sago tensor-harmonic field equations, Schwarzschild M=1, tortoise variable
`rs = r + 2 ln(r/2 − 1)`). Pure arithmetic, HDF5-free, Kernels-free, type-generic `{T}` so BigFloat
precision is inherited.

State layout (matching the legacy C/Julia convention): for field index `k` (0-based),
`y[1+4k]=Re(R)`, `y[2+4k]=Re(R')`, `y[3+4k]=Im(R)`, `y[4+4k]=Im(R')`; the RHS `g` has the same
layout. Each system returns `g`.

These RHS are CONSUMED by `SelfForce.LorenzField` and `AccretionDisc.Force` but OWNED here. The
legacy `monopole_ode!` carried two arities — H1Lorenz's `(g,rs,y,l,ω)` with a DEAD `l`, and
DiscModes' `(g,rs,y,ω)`; their bodies are bit-identical, so this is the single canonical arity-4
primitive plus a trivial arity-5 method `monopole_ode!(g,rs,y,l,ω)=monopole_ode!(g,rs,y,ω)` kept
for the H1Lorenz consumer convention (AHA: one decision, one implementation).

Background helpers `f`, `f'`, `Vₗ`, and the tortoise inverse `r(rs)` are delegated to
`Spacetime` (`schwarzschild_lapse`/`schwarzschild_dlapse`/`bls_potential`/`tortoise_schw_inverse`),
already grounded there.
"""
module LorenzRadial

using Spacetime: schwarzschild_lapse, schwarzschild_dlapse, bls_potential, tortoise_schw_inverse

export even_ode!, even_dipole_ode!, even_static_ode!, monopole_ode!, odd_ode!, odd_dipole_ode!

# background-scalar bundle at tortoise rs: (r, f, f', Vₗ, W) with W = 4Vₗ − ω²
@inline function _bg(rs, l, ω)
    r  = tortoise_schw_inverse(rs)
    f  = schwarzschild_lapse(r)
    df = schwarzschild_dlapse(r)
    Vl = bls_potential(r, l)
    return r, f, df, Vl, 4Vl - ω^2
end

"""
    even_ode!(g, rs, y, l, ω) -> g

Even-parity (l≥2) homogeneous Lorenz-gauge system: 5 coupled fields R1,R3,R5,R6,R7 → 20 components.
"""
function even_ode!(g, rs, y, l, ω)
    r, eff, df, _, W = _bg(rs, l, ω)
    # R1 (k=0): couples R3(k=1), R5(k=2), R6(k=3)
    g[1]=y[2];  g[2]=W*y[1] + 2eff*df*y[6] + 2eff/r^2*(1-4/r)*(y[1]-y[9]-eff*y[5]) - 2eff^2/r^2*(1-6/r)*y[13]
    g[3]=y[4];  g[4]=W*y[3] + 2eff*df*y[8] + 2eff/r^2*(1-4/r)*(y[3]-y[11]-eff*y[7]) - 2eff^2/r^2*(1-6/r)*y[15]
    # R3 (k=1)
    g[5]=y[6];  g[6]=W*y[5] - 2eff/r^2*(y[1]-y[9]-(1-4/r)*(y[5]+y[13]))
    g[7]=y[8];  g[8]=W*y[7] - 2eff/r^2*(y[3]-y[11]-(1-4/r)*(y[7]+y[15]))
    # R5 (k=2)
    g[9]=y[10];   g[10]=W*y[9]  + 4eff/r^2*((1-4.5/r)*y[9]  - (one(r)/2)*l*(l+1)*(y[1]-eff*y[5]) + (one(r)/2)*(1-3/r)*(l*(l+1)*y[13]-y[17]))
    g[11]=y[12];  g[12]=W*y[11] + 4eff/r^2*((1-4.5/r)*y[11] - (one(r)/2)*l*(l+1)*(y[3]-eff*y[7]) + (one(r)/2)*(1-3/r)*(l*(l+1)*y[15]-y[19]))
    # R6 (k=3)
    g[13]=y[14];  g[14]=W*y[13] - 2eff/r^2*(y[1]-y[9]-(1-4/r)*(y[5]+y[13]))
    g[15]=y[16];  g[16]=W*y[15] - 2eff/r^2*(y[3]-y[11]-(1-4/r)*(y[7]+y[15]))
    # R7 (k=4)
    g[17]=y[18];  g[18]=W*y[17] - 2eff/r^2*(y[17]+(l+2)*(l-1)*y[9])
    g[19]=y[20];  g[20]=W*y[19] - 2eff/r^2*(y[19]+(l+2)*(l-1)*y[11])
    g
end

"""
    even_dipole_ode!(g, rs, y, l, ω) -> g

Even-parity dipole (l=1) system: 4 coupled fields R1,R3,R5,R6 → 16 components (no R7).
"""
function even_dipole_ode!(g, rs, y, l, ω)
    r, eff, df, _, W = _bg(rs, l, ω)
    g[1]=y[2];  g[2]=W*y[1] + 2eff*df*y[6] + 2eff/r^2*(1-4/r)*(y[1]-y[9]-eff*y[5]) - 2eff^2/r^2*(1-6/r)*y[13]
    g[3]=y[4];  g[4]=W*y[3] + 2eff*df*y[8] + 2eff/r^2*(1-4/r)*(y[3]-y[11]-eff*y[7]) - 2eff^2/r^2*(1-6/r)*y[15]
    g[5]=y[6];  g[6]=W*y[5] - 2eff/r^2*(y[1]-y[9]-(1-4/r)*(y[5]+y[13]))
    g[7]=y[8];  g[8]=W*y[7] - 2eff/r^2*(y[3]-y[11]-(1-4/r)*(y[7]+y[15]))
    g[9]=y[10];   g[10]=W*y[9]  + 4eff/r^2*((1-4.5/r)*y[9]  - (one(r)/2)*l*(l+1)*(y[1]-eff*y[5]) + (one(r)/2)*(1-3/r)*(l*(l+1)*y[13]))
    g[11]=y[12];  g[12]=W*y[11] + 4eff/r^2*((1-4.5/r)*y[11] - (one(r)/2)*l*(l+1)*(y[3]-eff*y[7]) + (one(r)/2)*(1-3/r)*(l*(l+1)*y[15]))
    g[13]=y[14];  g[14]=W*y[13] - 2eff/r^2*(y[1]-y[9]-(1-4/r)*(y[5]+y[13]))
    g[15]=y[16];  g[16]=W*y[15] - 2eff/r^2*(y[3]-y[11]-(1-4/r)*(y[7]+y[15]))
    g
end

"""
    even_static_ode!(g, rs, y, l, ω) -> g

Even-static special case (l+m even, m=0, n=0 ⇒ ω=0): 3 coupled fields R1,R3,R5 → 12 components.
R6 is eliminated algebraically, R6 = 1/(2f)·(r/f·R1' − r·R3' + R1 − R5 − f·R3).
"""
function even_static_ode!(g, rs, y, l, ω)
    r, eff, df, _, W = _bg(rs, l, ω)
    R6re = 1/(2eff)*(r/eff*y[2] - r*y[6] + y[1] - y[9] - eff*y[5])
    R6im = 1/(2eff)*(r/eff*y[4] - r*y[8] + y[3] - y[11] - eff*y[7])
    g[1]=y[2];  g[2]=W*y[1] + 2eff*df*y[6] + 2eff/r^2*(1-4/r)*(y[1]-y[9]-eff*y[5]) - 2eff^2/r^2*(1-6/r)*R6re
    g[3]=y[4];  g[4]=W*y[3] + 2eff*df*y[8] + 2eff/r^2*(1-4/r)*(y[3]-y[11]-eff*y[7]) - 2eff^2/r^2*(1-6/r)*R6im
    g[5]=y[6];  g[6]=W*y[5] - 2eff/r^2*(y[1]-y[9]-(1-4/r)*(y[5]+R6re))
    g[7]=y[8];  g[8]=W*y[7] - 2eff/r^2*(y[3]-y[11]-(1-4/r)*(y[7]+R6im))
    g[9]=y[10];   g[10]=W*y[9]  + 4eff/r^2*((1-4.5/r)*y[9]  - (one(r)/2)*l*(l+1)*(y[1]-eff*y[5]) - (one(r)/2)*(1-3/r)*(r/eff*y[10]+2*y[9]))
    g[11]=y[12];  g[12]=W*y[11] + 4eff/r^2*((1-4.5/r)*y[11] - (one(r)/2)*l*(l+1)*(y[3]-eff*y[7]) - (one(r)/2)*(1-3/r)*(r/eff*y[12]+2*y[11]))
    g
end

"""
    monopole_ode!(g, rs, y, ω) -> g

Monopole (l=0) special case: 3 coupled fields R1,R3,R6 → 12 components. Canonical arity-4 form;
the potential reduces to Vₗ = ¼ f f'/r (the l(l+1) term vanishes).
"""
function monopole_ode!(g, rs, y, ω)
    r  = tortoise_schw_inverse(rs)
    eff = schwarzschild_lapse(r); df = schwarzschild_dlapse(r)
    Vl = (one(r)/4)*eff*(df/r); W = 4Vl - ω^2
    g[1]=y[2];  g[2]=W*y[1] + 2eff*df*y[6] + 2eff/r^2*(1-4/r)*(y[1]-eff*y[5]) - 2eff^2/r^2*(1-6/r)*y[9]
    g[3]=y[4];  g[4]=W*y[3] + 2eff*df*y[8] + 2eff/r^2*(1-4/r)*(y[3]-eff*y[7]) - 2eff^2/r^2*(1-6/r)*y[11]
    g[5]=y[6];  g[6]=W*y[5] - 2eff/r^2*(y[1]-(1-4/r)*(y[5]+y[9]))
    g[7]=y[8];  g[8]=W*y[7] - 2eff/r^2*(y[3]-(1-4/r)*(y[7]+y[11]))
    g[9]=y[10];   g[10]=W*y[9]  - 2eff/r^2*(y[1]-(1-4/r)*(y[5]+y[9]))
    g[11]=y[12];  g[12]=W*y[11] - 2eff/r^2*(y[3]-(1-4/r)*(y[7]+y[11]))
    g
end

# H1Lorenz consumer convention: arity-5 with a DEAD `l` (verified bit-identical across l in the
# arity gate). One decision, one implementation — delegate to the canonical arity-4 form (AHA).
monopole_ode!(g, rs, y, l, ω) = monopole_ode!(g, rs, y, ω)

"""
    odd_ode!(g, rs, y, l, ω) -> g

Odd-parity (l≥2) system: 2 coupled fields R9,R10 → 8 components (λ = (l+2)(l−1)).
"""
function odd_ode!(g, rs, y, l, ω)
    r, eff, _, _, W = _bg(rs, l, ω)
    lam = (l+2)*(l-1)
    g[1]=y[2];  g[2]=(W + 4eff/r^2*(1-4.5/r))*y[1] - 2eff/r^2*(1-3/r)*y[5]
    g[3]=y[4];  g[4]=(W + 4eff/r^2*(1-4.5/r))*y[3] - 2eff/r^2*(1-3/r)*y[7]
    g[5]=y[6];  g[6]=(W - 2eff/r^2)*y[5] - 2eff/r^2*lam*y[1]
    g[7]=y[8];  g[8]=(W - 2eff/r^2)*y[7] - 2eff/r^2*lam*y[3]
    g
end

"""
    odd_dipole_ode!(g, rs, y, l, ω) -> g

Odd-parity dipole (l=1) special case: 1 field R9 → 4 components.
"""
function odd_dipole_ode!(g, rs, y, l, ω)
    r  = tortoise_schw_inverse(rs)
    eff = schwarzschild_lapse(r); Vl = bls_potential(r, l)
    W = 4Vl + 4eff/r^2*(1-4.5/r) - ω^2
    g[1]=y[2]; g[2]=W*y[1]
    g[3]=y[4]; g[4]=W*y[3]
    g
end

end # module LorenzRadial
