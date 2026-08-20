# End-to-end second-order Lorenz-gauge Ricci source assembly + HDF5 writer (legacy
# SecondOrderRicci/assembly.jl). Couples 3 harmonics (Wigner.coupling, evaluated HERE host-side by
# the `r2_term` wrapper — the §5.5/§7.7.1 Kernels-purity hoist: Kernels carries no Wigner machinery)
# into the 10 quadratic r2 combinations on the pure-arithmetic `Kernels.r2_term` host over the h1
# corpus read by ReadCorpus.read_h1_corpus, supplies the analytic h1_mass/h1_spin static fields, and closes
# the SecondOrder submodule. `ZeroPadField`/`read_h1_corpus` come from ReadCorpus, `selection_modes` from
# Modes, `h1_mass`/`h1_spin` from AnalyticFields; the r2 ARITHMETIC stays Kernels-owned (NOT re-owned).

import HDF5
import Kernels

_has_parent_mode(h, dh, ddh, l::Integer, m::Integer) =
    any(i -> haskey(h, (i, l, m)) || haskey(dh, (i, l, m)) || haskey(ddh, (i, l, m)), 1:10)

# Highest multipole present in a read_h1_corpus field corpus (the driver's `h1lmax`).
_h1lmax(h::AbstractDict) = isempty(h) ? -1 : maximum(l for (_, l, _) in keys(h))

# Legacy Ricci.cc command-line bandwidth semantics (L66-75).  A negative dlmax is not an empty
# bandwidth: it first trims that many multipoles from the usable first-order corpus, then selects the
# full remaining bandwidth.  Positive bandwidths are simply capped at the usable corpus maximum.
function _legacy_dlmax(h1lmax::Integer, dlmax::Integer)
    if dlmax < 0
        h1lmax += dlmax
        dlmax = h1lmax
    else
        dlmax = min(dlmax, h1lmax)
    end
    return h1lmax, dlmax
end

"""
    r2_term(n, M, r0, r, f, fp, h1, dh1, ddh1, h2, dh2, ddh2, l3, m3, l1, m1, l2, m2; threaded=false)

The coupling-evaluating HOST WRAPPER over the pure-arithmetic `Kernels.r2_term` GPU host (hoisted
here from Kernels, S4-c: design §5.5/§7.7.1 Kernels purity — SecondOrder owns the Wigner coupling
consumption). Evaluates the `n`-th kernel's 3-harmonic angular coupling array ONCE host-side via
[`coupling`](@ref) (`Wigner.coupling`) following the pure-data contract `Kernels.R2_COUPLING_SPEC[n]`,
then launches `Kernels.r2_term` with the precomputed array. Signature, semantics and output are
IDENTICAL to the pre-hoist `Kernels.r2_term` (and to legacy `SecondOrderRicci.r2_term`).
"""
function r2_term(n::Integer, M, r0, r, f, fp,
                 h1, dh1, ddh1, h2, dh2, ddh2,
                 l3::Integer, m3::Integer, l1::Integer, m1::Integer, l2::Integer, m2::Integer;
                 threaded::Bool = false)
    (1 <= n <= 10) || error("SelfForce.SecondOrder.r2_term: kernel R2_$n out of range (have 1:10)")
    (_has_parent_mode(h1, dh1, ddh1, l1, m1) &&
     _has_parent_mode(h2, dh2, ddh2, l2, m2)) || return zeros(ComplexF64, length(r))
    c = Float64[coupling(l3, m3, s3, l1, m1, s1, l2, m2, s2)
                for (s3, s1, s2) in Kernels.R2_COUPLING_SPEC[n]]
    return Kernels.r2_term(n, M, r0, r, f, fp, h1, dh1, ddh1, h2, dh2, ddh2,
                           l3, m3, l1, m1, l2, m2, c; threaded = threaded)
end

"""
    second_order_ricci(fields, i3, l3, m3; h1lmax, dlmax) -> Vector{ComplexF64}

The `(i3,l3,m3)` component of the self-coupled second-order Lorenz-gauge Ricci source: the `R2_i3`
quadratic kernel summed over all first-order mode pairs `(l1,m1),(l2,m2=m3-m1)` that couple in. Port
of the inner accumulation of the `Ricci` driver (`ref-code/SecondOrderRicci/src/Ricci.cc` L117-159),
self-coupling case `dirA == dirB` so the corpus serves as both kernel arguments.

Normalization contract: for `g(ε) = g₀ + εh₁`, this function returns the second variation
`D²R[h₁,h₁] = 2[ε²]R[g₀ + εh₁]`.  In the no-factorial second-order metric expansion this is the
`2δ²R` Lorenz-equation source term.  `PerturbationEquations` 0.4.0 calls the Taylor coefficient
`[ε²]R = ½D²R` its `d2R` output, so compare this API to **twice** that operator.

`fields` is an [`H1RadialJetCorpus`](@ref) returned by [`read_h1_corpus`](@ref). The `(l1,l2)` summation
runs over
`l3-dlmax ≤ l1,l2 ≤ l3+dlmax` intersected with `0 ≤ l1`, `|l3-l1| ≤ l2 ≤ l3+l1`, and the corpus cap
`l1,l2 ≤ h1lmax`; `m1 ∈ -l1:l1` with `|m3-m1| ≤ l2`. `M = 1`. The radial accumulation is on the
Float64-locked `Kernels.r2_term` host (the GPU r2 host), so the source is `ComplexF64`. A negative
`dlmax` follows the legacy driver sentinel: first set `h1lmax ← h1lmax + dlmax`, then use the full
remaining bandwidth `dlmax ← h1lmax`.
"""
function second_order_ricci(fields, i3::Integer, l3::Integer, m3::Integer;
                            h1lmax::Integer = _h1lmax(fields.h),
                            dlmax::Integer = h1lmax)
    h1lmax, dlmax = _legacy_dlmax(h1lmax, dlmax)
    M = 1.0
    r = fields.r; f = fields.f; fp = fields.fp; r0 = fields.r0
    N = length(r)
    # Kernels.r2_term's _pack does its own zero-padding via haskey, so plain ComplexF64 corpora suffice.
    H = _cf64(fields.h); DH = _cf64(fields.dh); DDH = _cf64(fields.ddh)
    rf = Float64.(r); ff = Float64.(f); fpf = Float64.(fp); r0f = Float64(r0)
    acc = zeros(ComplexF64, N)
    for l1 in max(0, l3 - dlmax):min(h1lmax, l3 + dlmax)
        for l2 in max(abs(l3 - l1), l3 - dlmax):min(min(l3 + l1, h1lmax), l3 + dlmax)
            for m1 in -l1:l1
                abs(m3 - m1) > l2 && continue
                tmp = r2_term(i3, M, r0f, rf, ff, fpf, H, DH, DDH, H, DH, DDH,
                              l3, m3, l1, m1, l2, m3 - m1)
                @. acc += tmp
            end
        end
    end
    return acc
end

# Coerce a Dict{(i,l,m)=>Vector{Complex{T}}} corpus to ComplexF64 for the Float64 r2 host.
_cf64(d::AbstractDict{K,Vector{ComplexF64}}) where {K} = d
_cf64(d::AbstractDict{K,V}) where {K,V} =
    Dict{K,Vector{ComplexF64}}(k => ComplexF64.(v) for (k, v) in d)

"""
    ricci_source(fields; lmax, dlmax) -> Dict{NTuple{3,Int},Vector{ComplexF64}}

Full self-coupled second-order Ricci source: every valid target mode `(i3,l3,m3)` of
`selection_modes(lmax)` mapped to its assembled radial vector. Top-level port of the `Ricci` driver
(`ref-code/SecondOrderRicci/src/Ricci.cc`); `lmax`/`dlmax` are capped at the corpus `h1lmax` exactly
as the driver caps them. For `dlmax < 0`, the driver first lowers the usable first-order cutoff by
`abs(dlmax)`, then uses that cutoff as the bandwidth and caps target `lmax` afterward. Defaults
reproduce the single-directory invocation `dlmax = lmax = h1lmax`. Each returned mode has the
[`second_order_ricci`](@ref) normalization `D²R = 2[ε²]R`.
"""
function ricci_source(fields; lmax::Integer = _h1lmax(fields.h), dlmax::Integer = lmax)
    h1lmax = _h1lmax(fields.h)
    h1lmax, dlmax = _legacy_dlmax(h1lmax, dlmax)
    lmax = min(lmax, h1lmax)
    src = Dict{NTuple{3,Int},Vector{ComplexF64}}()
    for (i3, l3, m3) in selection_modes(lmax)
        src[(i3, l3, m3)] = second_order_ricci(fields, i3, l3, m3; h1lmax = h1lmax, dlmax = dlmax)
    end
    return src
end

# --- static-field (mass / spin) symmetrized coupling: Ricci_M.cc / Ricci_S.cc --------------------

# Split an h1_mass / h1_spin result Dict{(i,l,m)=>(;h,dh,ddh)} into the three (i,l,m)=>Vector dicts
# of the read_h1_corpus layout (ComplexF64), so the R2 kernels see the static field exactly like the dynamic
# one.
function _static_corpus(sfield)
    h = Dict{NTuple{3,Int},Vector{ComplexF64}}()
    dh = Dict{NTuple{3,Int},Vector{ComplexF64}}()
    ddh = Dict{NTuple{3,Int},Vector{ComplexF64}}()
    for (k, v) in sfield
        h[k] = ComplexF64.(v.h); dh[k] = ComplexF64.(v.dh); ddh[k] = ComplexF64.(v.ddh)
    end
    return h, dh, ddh
end

# One target (i3,l3,m3) of the symmetric coupling between a static l1-multipole field (mass: l1=0;
# spin: l1=1) with m1=0 and the dynamic first-order field. Port of the inner loop shared by
# Ricci_M.cc (L103-153) and Ricci_S.cc (L103-153): l2 ∈ [|l3-l1|,l3+l1] ∩ [l3-dlmax,l3+dlmax] ∩
# [0,h1lmax], m1=0, source = R2(static,dynamic) + R2(dynamic,static).
function _static_coupled_mode(M, r0, r, f, fp, SH, SDH, SDDH, DYh, DYdh, DYddh,
                              l1::Integer, dlmax::Integer,
                              i3::Integer, l3::Integer, m3::Integer; h1lmax::Integer)
    acc = zeros(ComplexF64, length(r)); m1 = 0
    for l2 in max(abs(l3 - l1), l3 - dlmax):min(min(l3 + l1, h1lmax), l3 + dlmax)
        abs(m3 - m1) > l2 && continue
        m2 = m3 - m1
        tmp1 = r2_term(i3, M, r0, r, f, fp, SH, SDH, SDDH, DYh, DYdh, DYddh, l3, m3, l1, m1, l2, m2)
        tmp2 = r2_term(i3, M, r0, r, f, fp, DYh, DYdh, DYddh, SH, SDH, SDDH, l3, m3, l2, m2, l1, m1)
        @. acc += tmp1 + tmp2
    end
    return acc
end

function _ricci_static_source(fields, sfield, l1::Integer, dlmax::Integer; lmax::Integer)
    M = 1.0
    h1lmax = _h1lmax(fields.h)
    rf = Float64.(fields.r); ff = Float64.(fields.f); fpf = Float64.(fields.fp); r0f = Float64(fields.r0)
    SH, SDH, SDDH = _static_corpus(sfield)
    DYh = _cf64(fields.h); DYdh = _cf64(fields.dh); DYddh = _cf64(fields.ddh)
    src = Dict{NTuple{3,Int},Vector{ComplexF64}}()
    for (i3, l3, m3) in selection_modes(lmax)
        src[(i3, l3, m3)] = _static_coupled_mode(M, r0f, rf, ff, fpf,
                                                 SH, SDH, SDDH, DYh, DYdh, DYddh,
                                                 l1, dlmax, i3, l3, m3; h1lmax = h1lmax)
    end
    return src
end

"""
    ricci_source_mass(fields; lmax) -> Dict{NTuple{3,Int},Vector{ComplexF64}}

Second-order Ricci source from the symmetric coupling of the first-order MASS perturbation
(`h1_mass`, an l=0 static field) to the dynamic first-order field. Port of `Ricci_M`
(`ref-code/SecondOrderRicci/src/Ricci_M.cc`): fixed `l1=0`, `m1=0`, `dlmax=0`, symmetrized
`R2(mass,dyn)+R2(dyn,mass)` over `selection_modes(lmax)`.
"""
ricci_source_mass(fields; lmax::Integer = _h1lmax(fields.h)) =
    _ricci_static_source(fields, h1_mass(Float64.(fields.r)), 0, 0; lmax = lmax)

"""
    ricci_source_spin(fields; lmax) -> Dict{NTuple{3,Int},Vector{ComplexF64}}

Second-order Ricci source from the symmetric coupling of the first-order SPIN (angular-momentum)
perturbation (`h1_spin`, an l=1, m=0 static field) to the dynamic first-order field. Port of
`Ricci_S` (`ref-code/SecondOrderRicci/src/Ricci_S.cc`): fixed `l1=1`, `m1=0`, `dlmax=1`, symmetrized
`R2(spin,dyn)+R2(dyn,spin)` over `selection_modes(lmax)`.
"""
ricci_source_spin(fields; lmax::Integer = _h1lmax(fields.h)) =
    _ricci_static_source(fields, h1_spin(Float64.(fields.r)), 1, 1; lmax = lmax)

"""
    write_src(path, src; r=nothing) -> path

Write an assembled source `Dict{(i,l,m)=>Vector{Complex}}` to HDF5 at `path` in the upstream driver
layout: one `"src i={i} l={l} m={m}"` dataset per mode (on-disk shape `[N][2]` = real,imag, matching
the C++ writer) plus, if `r` is given, an `"r"` grid dataset. Port of the output block of the `Ricci`
driver (`ref-code/SecondOrderRicci/src/Ricci.cc` L170-197). A Julia `(2,N)` array maps to the
column-reversed on-disk `[N][2]` C layout.
"""
function write_src(path::AbstractString, src::AbstractDict; r = nothing)
    HDF5.h5open(path, "w") do file
        for ((i, l, m), v) in src
            D = Matrix{Float64}(undef, 2, length(v))   # Julia (2,N) ⇒ on-disk [N][2] (C order)
            @inbounds for j in eachindex(v)
                D[1, j] = real(v[j]); D[2, j] = imag(v[j])
            end
            file["src i=$i l=$l m=$m"] = D
        end
        r === nothing || (file["r"] = collect(float.(r)))
    end
    return path
end
