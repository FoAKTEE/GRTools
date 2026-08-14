# Coupled-set descriptors + selection + the γ↔(l,m) mode enumeration — the leaf solve-plan layer of
# the Lorenz field solver (port of H1Lorenz coupled_sets.jl + driver.jl gamma_to_lm/lm_to_gamma/
# mode_list). A `CoupledSet` bundles, per parity sector, the CONSUMED radial RHS, the in/out
# boundary-condition series (already built in Boundary.jl), the R1..R10 sources (Orbit.jl) and the
# gauge callback (Solver.jl). It does NOT re-own the RHS: the `ode!` field POINTS AT
# `MetricReconstruction.Lorenz.LorenzRadial`'s `even_ode!`/…/`odd_dipole_ode!`. No integration here —
# pure wiring. Type-generic {T}: every member is precision-agnostic, the T flows through at call.

using MetricReconstruction.Lorenz.LorenzRadial:
    even_ode!, even_dipole_ode!, even_static_ode!, monopole_ode!, odd_ode!, odd_dipole_ode!

# γ↔(l,m) enumeration primitives are CANONICAL in NumericsCore.modeindex — imported (not re-owned)
# and re-exported through LorenzField for API parity. (NumericsCore already a declared SelfForce dep.)
using NumericsCore: gamma_to_lm, lm_to_gamma, mode_list

"""
    CoupledSet

Per-parity solve-plan descriptor. `kind` ∈ (:monopole,:odd_dipole,:even_dipole,:even,:odd,
:even_static) tags the special handling in the assembler. `num_coupled` is the BLS field count of
the sector, `num_gauge` the algebraically-reconstructed gauge-field count. The callbacks are stored
as function objects (no closures over precision), so a `CoupledSet` is type-generic:

  - `ode!`    — the CONSUMED `LorenzRadial` RHS (`(g,rs,y,l,ω)->g`), referenced not re-owned.
  - `out_bc!` — outer (outgoing) BC series, uniform `(a0s,l,ω,p; rs_out)` signature.
  - `in_bc!`  — inner (near-horizon, ingoing) BC series, `(a0s,l,ω)`.
  - `sources` — the R1..R10 source coefficient functions of the sector (`(l,m,orbit)->Complex`).
  - `gauge!`  — in-place algebraic gauge reconstruction `(inhom,orbit,l,m,ω,num_coupled)->inhom`.
"""
struct CoupledSet
    kind::Symbol
    num_coupled::Int
    num_gauge::Int
    ode!::Function
    out_bc!::Function
    in_bc!::Function
    sources::Vector{Function}
    gauge!::Function
end

# Uniform (a0s,l,ω,p; rs_out) wrappers over the per-sector out-BC series (Boundary.jl). rs_out (the
# analytic far-field seed r*) is forwarded so the integrator can collapse the empty 110000→grid leg
# for BigFloat. The monopole/even-static series carry no `p`; the wrapper drops it.
_out_even(a0s, l, ω, p; rs_out=eltype(a0s)(110000))         = even_out_bcs(a0s, l, ω, p; rs_out=rs_out)
_out_even_dipole(a0s, l, ω, p; rs_out=eltype(a0s)(110000))  = even_dipole_out_bcs(a0s, l, ω, p; rs_out=rs_out)
_out_monopole(a0s, l, ω, p; rs_out=eltype(a0s)(110000))     = monopole_out_bcs(a0s, l, ω; rs_out=rs_out)
_out_odd(a0s, l, ω, p; rs_out=eltype(a0s)(110000))          = odd_out_bcs(a0s, l, ω, p; rs_out=rs_out)
_out_odd_dipole(a0s, l, ω, p; rs_out=eltype(a0s)(110000))   = odd_dipole_out_bcs(a0s, l, ω, p; rs_out=rs_out)
_out_even_static(a0s, l, ω, p; rs_out=eltype(a0s)(110000))  = even_static_out_bcs(a0s, l, ω; rs_out=rs_out)

# Gauge-callback closures wrapping the Solver.jl algebraic reconstructions (forward-referenced — the
# global names resolve at CALL time, after Solver.jl is included). One uniform
# (inhom,orbit,l,m,ω,num_coupled) signature.
_g_even(inhom, orbit, l, m, ω, nc)       = construct_R2_R4_gauge_fields_on_grid!(inhom, orbit, l, ω, nc)
_g_odd(inhom, orbit, l, m, ω, nc)        = construct_R8_gauge_field_on_grid!(inhom, orbit, l, m, ω, nc)
_g_even_static(inhom, orbit, l, m, ω, nc) = construct_R6_R7_R2_R4_gauge_fields_on_grid!(inhom, orbit, l, nc)

# Build the six per-parity descriptors. `ode!` fields are the CONSUMED LorenzRadial RHS (===), the
# sources are the Orbit.jl R*_source functions, the BC series the Boundary.jl wrappers above.
function _make_csets()
    even = CoupledSet(:even, 5, 2, even_ode!, _out_even, even_in_bcs,
        Function[R1_source, R3_source, R5_source, R6_source, R7_source], _g_even)
    odd = CoupledSet(:odd, 2, 1, odd_ode!, _out_odd, odd_in_bcs,
        Function[R9_source, R10_source], _g_odd)
    odd_dipole = CoupledSet(:odd_dipole, 1, 1, odd_dipole_ode!, _out_odd_dipole, odd_dipole_in_bcs,
        Function[R9_source], _g_odd)
    monopole = CoupledSet(:monopole, 3, 1, monopole_ode!, _out_monopole, monopole_in_bcs,
        Function[R1_source, R3_source, R6_source], _g_even)
    even_dipole = CoupledSet(:even_dipole, 4, 2, even_dipole_ode!, _out_even_dipole, even_dipole_in_bcs,
        Function[R1_source, R3_source, R5_source, R6_source], _g_even)
    even_static = CoupledSet(:even_static, 3, 4, even_static_ode!, _out_even_static, even_static_in_bcs,
        Function[R1_source, R3_source, R5_source], _g_even_static)
    (monopole=monopole, odd_dipole=odd_dipole, even_dipole=even_dipole,
     even=even, odd=odd, even_static=even_static)
end

"""Immutable table of the six per-parity `CoupledSet` descriptors (built once at load)."""
const CSETS = _make_csets()

"""
    select_coupled_set(l, m; static=false) -> CoupledSet

Parity-sector descriptor for mode `(l,m)` (port of H1Lorenz `select_correct_coupled_set`).
`(l+m)` even ⇒ even sector (l=0 monopole, l=1 even_dipole, else even); `(l+m)` odd ⇒ odd sector
(l=1 odd_dipole, else odd). With `static=true` and the static-mode condition `(l+m)` even, `m==0`,
`l>=2` (n=0 ⇒ ω=0), the even sector swaps to `even_static`; l=m=0 stays the analytic monopole.
"""
function select_coupled_set(l::Int, m::Int; static::Bool=false)
    if iseven(l + m)
        l == 0 && return CSETS.monopole
        l == 1 && return CSETS.even_dipole
        (static && m == 0) && return CSETS.even_static
        return CSETS.even
    else
        l == 1 && return CSETS.odd_dipole
        return CSETS.odd
    end
end
