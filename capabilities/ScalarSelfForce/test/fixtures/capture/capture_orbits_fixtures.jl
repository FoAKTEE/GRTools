# capture_orbits_fixtures.jl — ONE-SHOT generation recipe for test/fixtures/orbits_circular.jl
# (ESSF node `orbits-adapters`, schedule W2-C). Bakes the Spacetime-direct characterization
# oracle for `ScalarSelfForce.CircularOrbit`: every fixture number is computed by calling the
# EXPORTED Spacetime functions the adapter itself must wire (constants_of_motion,
# kerr_circular_ut_uphi, horizons, kerr_isco), so the baked rows pin BOTH the upstream values
# and the adapter's wiring (argument order, x = +1 equatorial convention, Omega = uphi/ut).
#
# RUN (from the repo root; regenerates the fixture file in place):
#   ~/.juliaup/bin/julia --project=capabilities/ScalarSelfForce \
#       capabilities/ScalarSelfForce/test/fixtures/capture/capture_orbits_fixtures.jl
#
# Grid (node contract): a in {0, ±0.5, ±0.9, ±0.998} × r0 in {isco, 1.001·isco, 4M, 6M, 10M}
# (M = 1; a < 0 = retrograde, motion always +φ — Spacetime sign convention). Sub-ISCO rows at
# retrograde spins are UNSTABLE circular geodesics (still exact solutions; dV_eff/dr = 0 holds).
#
# Literals are printed with `repr` (shortest round-trip Float64), so a re-run on the same
# Spacetime version reproduces the file bit-identically.

using Spacetime
using Printf

const SPINS = [0.0, 0.5, -0.5, 0.9, -0.9, 0.998, -0.998]
radii(a) = (ri = kerr_isco(a); [ri, 1.001 * ri, 4.0, 6.0, 10.0])

rows = String[]
for a in SPINS, r0 in radii(a)
    com = constants_of_motion(a, r0, 0.0, 1.0)          # (E, Lz, Q) — equatorial x = +1
    ut, uphi = kerr_circular_ut_uphi(a, r0)
    Omega = uphi / ut
    rp, rm = horizons(a)
    risco = kerr_isco(a)
    push!(rows,
        "    (a = $(repr(a)), r0 = $(repr(r0)), E = $(repr(com.E)), Lz = $(repr(com.Lz)), " *
        "ut = $(repr(ut)), uphi = $(repr(uphi)), Omega = $(repr(Omega)), " *
        "rp = $(repr(rp)), rm = $(repr(rm)), risco = $(repr(risco))),")
    # diagnostic only (NOT baked): circularity residual dV_eff/dr at r0, normalized by the sum
    # of term magnitudes — informs the test gate, never forced. NOTE the printed V_eff in
    # main.tex:190 carries a LAST-TERM SIGN TYPO (+M(L−aE)²/r³): its derivative does NOT vanish
    # on the paper's own E, L solutions (measured ≈0.5 normalized). The correct circularity
    # condition, derived from the Kerr equatorial radial equation
    # (u^r)² = E²−1 + 2M/r − (L²+a²(1−E²))/r² + 2M(L−aE)²/r³, is
    #   M/r² − (L²+a²(1−E²))/r³ + 3M(L−aE)²/r⁴ = 0.
    t1 = 1 / r0^2
    t2 = -(com.Lz^2 + a^2 * (1 - com.E^2)) / r0^3
    t3 = 3 * (com.Lz - a * com.E)^2 / r0^4
    @printf("a=%+6.3f r0=%9.6f  E=%12.9f  Lz=%+12.9f  |dV'|/Σ|terms| = %.3e\n",
            a, r0, com.E, com.Lz, abs(t1 + t2 + t3) / (abs(t1) + abs(t2) + abs(t3)))
end

fixture_path = joinpath(@__DIR__, "..", "orbits_circular.jl")
open(fixture_path, "w") do io
    print(io, """
# orbits_circular.jl — BAKED Spacetime characterization fixtures for the ESSF node
# `orbits-adapters` (ScalarSelfForce.CircularOrbit, schedule W2-C). DO NOT EDIT BY HAND.
#
# GENERATION RECIPE: see test/fixtures/capture/capture_orbits_fixtures.jl (one-shot capture
# calling ONLY exported Spacetime functions — constants_of_motion(a, r0, 0.0, 1.0) [equatorial
# x = +1], kerr_circular_ut_uphi(a, r0), horizons(a), kerr_isco(a); Omega = uphi/ut). Re-running
# the recipe on the same Spacetime version reproduces this file bit-identically.
#
# CITATION: arXiv:2606.04998 main.tex:176-209 (Sec. II.A — circular equatorial Kerr geodesics:
# E, L via the two circularity conditions; u^t, u^phi closed forms; Omega = 1/(a + sqrt(r0^3/M));
# ISCO condition). Conventions: G = c = 1, M = 1; a < 0 encodes retrograde (motion always +φ).
# Grid: a ∈ {0, ±0.5, ±0.9, ±0.998} × r0 ∈ {isco, 1.001·isco, 4, 6, 10}; sub-ISCO rows at
# retrograde spins are unstable circular geodesics (exact solutions; dV_eff/dr = 0 still holds).
# Characterization tolerance (frozen; may only tighten): rtol 1e-12.
const ORBITS_CIRCULAR_FIXTURES = [
""")
    for row in rows
        println(io, row)
    end
    print(io, "]\n")
end
println("wrote $(abspath(fixture_path)) ($(length(rows)) rows)")
