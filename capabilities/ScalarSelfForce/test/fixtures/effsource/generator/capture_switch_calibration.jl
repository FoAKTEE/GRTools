#!/usr/bin/env julia
# capture_switch_calibration.jl — BAKE the effective-source route-switch calibration table
# (ESSF node `effective-source-series`, W8-B). DO NOT EDIT THE OUTPUT BY HAND.
#
# Regenerate (≈ 11 min single-threaded; needs no network):
#   env JULIA_PKG_OFFLINE=true JULIA_NUM_THREADS=1 julia --project=capabilities/ScalarSelfForce \
#       capabilities/ScalarSelfForce/test/fixtures/effsource/generator/capture_switch_calibration.jl \
#       capabilities/ScalarSelfForce/test/fixtures/effsource/switch_calibration.jl
#
# WHAT IT IS. For each of the 3240 points of the frozen calibration grid
#     6 fixture configs × m ∈ {0,1,2,5,10,20} × λ = 10^(-5 : 0.5 : -1)
#       × 5 ray directions frac ∈ {0, ¼, ½, ¾, 1} × both signs of Δr
# (the point is ϱ = λ exactly, Δr = sgn·frac·λ/√rho2_dr2, Δθ = √(1−frac²)·λ/√rho2_dtheta2),
# it records three numbers:
#     rel_series  = |series  − truth| / |truth|      (Float64 route vs BigFloat(256) truth)
#     rel_direct  = |direct  − truth| / |truth|
#     eps_cancel  = eps(Float64) · effective_source_bundle(...).cancellation
# The TRUTH is `effective_source_m` at BigFloat(256), where the direct assembly's λ^{-2}
# roundoff amplification is ~1e-65 — the same oracle the node's matching-zone gate uses.
#
# WHY IT IS BAKED. The shipped test can then replay the whole calibration against the LIVE
# switch rule (which needs only Float64) and re-assert the three headline numbers — worst
# achieved/optimal factor, worst switched relative error, worst near-crossover route
# disagreement — instead of quoting them from a comment. The values are route-INDEPENDENT
# (they are properties of the two routes, not of the rule), so a change to the switch never
# invalidates them; a change to the series VALUE or to the direct route does, and the
# node's own fixture/identity gates catch that first.
using ScalarSelfForce, Printf
const FIXDIR = joinpath(@__DIR__, "..")
include(joinpath(FIXDIR, "effsource_configs.jl"))

out = ARGS[1]
setprecision(BigFloat, 256)
phi_p = EFFSOURCE_PHIP
const VALS = Float64[]
const NDONE = Ref(0)
for c in EFFSOURCE_CONFIGS, m in (0, 1, 2, 5, 10, 20), e in -5.0:0.5:-1.0,
    frac in (0.0, 0.25, 0.5, 0.75, 1.0), sgn in (1, -1)

    o = CircularOrbit(c.a, c.r0); p = Puncture(o)
    oB = CircularOrbit(big(c.a), big(c.r0)); pB = Puncture(oB)
    L = 10.0^e
    dr = sgn * frac * L / sqrt(o.rho2_dr2)
    dth = sqrt(max(0.0, 1 - frac^2)) * L / sqrt(o.rho2_dtheta2)
    r = o.r0 + dr; z = -sin(dth)
    truth = effective_source_m(pB, m, big(r), big(z); phi_p = big(phi_p))
    ser = effective_source_series_m(p, m, r, z; phi_p = phi_p)
    b = effective_source_bundle(p, m, r, z; phi_p = phi_p)
    push!(VALS, Float64(abs(big(ser) - truth) / abs(truth)))
    push!(VALS, Float64(abs(big(b.source) - truth) / abs(truth)))
    push!(VALS, eps(Float64) * b.cancellation)
    NDONE[] += 1
    NDONE[] % 200 == 0 && @printf("  %d / 3240\n", NDONE[])
end
open(out, "w") do io
    println(io, """
# switch_calibration.jl — BAKED route-switch calibration table (ESSF node
# `effective-source-series`, W8-B). DO NOT EDIT BY HAND. Regenerate with
#   julia --project=capabilities/ScalarSelfForce \\
#     capabilities/ScalarSelfForce/test/fixtures/effsource/generator/capture_switch_calibration.jl \\
#     capabilities/ScalarSelfForce/test/fixtures/effsource/switch_calibration.jl
#
# ORACLE: this package's own `effective_source_m` at BigFloat(256) (relative roundoff ~1e-65
#   there), against its Float64 series and direct routes. No external code or data.
#
# GRID (3240 rows, in exactly this nesting order — the test regenerates the coordinates and
# asserts the count, so the order is part of the contract):
#   for c in EFFSOURCE_CONFIGS, m in (0,1,2,5,10,20), e in -5.0:0.5:-1.0,
#       frac in (0.0,0.25,0.5,0.75,1.0), sgn in (1,-1)
#   lambda = 10^e ; Delta r = sgn*frac*lambda/sqrt(rho2_dr2) ;
#   Delta theta = sqrt(1-frac^2)*lambda/sqrt(rho2_dtheta2) ; r = r0 + Delta r ; z = -sin(Delta theta)
#
# COLUMNS (row k = matrix column M[:, k]):
#   1 rel_series | 2 rel_direct | 3 eps(Float64)*cancellation
const EFFSOURCE_SWITCH_CALIBRATION = reshape(Float64[""")
    for i in 1:3:length(VALS)
        @printf(io, "    %.9e, %.9e, %.9e,\n", VALS[i], VALS[i+1], VALS[i+2])
    end
    println(io, "], 3, $(length(VALS) ÷ 3))")
end
@printf("wrote %s (%d rows, %d bytes)\n", out, length(VALS) ÷ 3, filesize(out))
