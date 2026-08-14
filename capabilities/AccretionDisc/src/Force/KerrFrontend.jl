# KERR FRONT-END — the a≠0 metric source: MetricReconstruction.Lorenz (Dolan-2024 Kerr Lorenz
# reconstruction) → BLS h̄^i(r) → H1Mode{T} (design §4.10 Force: "a≠0 front-end via
# MetricReconstruction.Lorenz"; §7.3: AccretionDisc reaches the Lorenz machinery ONLY through
# MetricReconstruction — NEVER through SelfForce/Kernels). From DiscEMRIForce/kerr_source.jl.
#
# The a=0 path loads the Lorenz-gauge BLS corpus (`H1Mode` from Loader.jl) and assembles h_μν via
# Assemble.jl. The Kerr (a≠0) path REPLACES THE METRIC SOURCE ONLY: the spin-general
# `Lorenz.to_BLS_kerr(a,r0,l,m,r)` produces the SAME 10 Barack-Lousto-Sago variables h̄¹..h̄¹⁰(r)
# the corpus holds, so the downstream effective-force → 𝒬-source → forced-ODE → F8 pipeline
# (Source/Fluid/Pipeline, all spin-agnostic given h_μν) runs UNCHANGED. Per (l,m), the (l+m)-parity
# selects the physical BLS fields: even ⇒ the 7 even fields [h̄¹,h̄³,h̄⁵,h̄⁶,h̄⁷,h̄²,h̄⁴] (loader
# order), odd ⇒ the 3 odd fields [h̄⁹,h̄¹⁰,h̄⁸] — matching the corpus convention.
#
# PARALLEL (§8): the per-l reconstruction (one banded closure_kerr solve shared across l) and the
# per-grid-point to_BLS_kerr evaluation are the embarrassingly-parallel axes — threaded with
# race-free preallocated slots under the blas_single_thread guard; threaded == serial bit-identical.

import MetricReconstruction
import MetricReconstruction.Lorenz: closure_kerr, smatrices_kerr, reconstruct_kerr, to_BLS_kerr,
                                    _attach_recmeta
import NumericsCore: blas_single_thread, highorder_deriv

# loader field order (Assemble.jl): even col→BLS = [h̄¹,h̄³,h̄⁵,h̄⁶,h̄⁷,h̄²,h̄⁴]; odd = [h̄⁹,h̄¹⁰,h̄⁸].
const _KERR_EVEN_PERM = (1, 3, 5, 6, 7, 2, 4)
const _KERR_ODD_PERM  = (9, 10, 8)

"""
    kerr_mode(a, r_p, l, m; T=Float64, rgrid, recs=nothing, S=nothing, ls=nothing, threaded=true)
        -> H1Mode{T}

Assemble a single Kerr (l,m) metric-perturbation mode into an `H1Mode{T}` (the same container the
a=0 corpus loader returns), from `MetricReconstruction.Lorenz.to_BLS_kerr(a,r_p,l,m,r)` sampled on
the radial grid `rgrid`. The BLS h̄^i(r) are packed into `val_right`/`dr_right` in the loader field
order (parity-aware) so `hmunu`/`force_mode`/`source_F` consume them unchanged; `dr_right` uses
five-point Fornberg differentiation on any strictly increasing grid. Pass shared `recs`/`S`/`ls`
(from one banded `closure_kerr`) to
reuse the reconstruction across modes; `threaded=false` when the CALLER already threads the outer
(l) axis (avoids nested-thread oversubscription).
"""
function kerr_mode(a, r_p, l::Integer, m::Integer; T::Type=Float64, rgrid,
                   recs=nothing, S=nothing, ls=nothing, threaded::Bool=true,
                   bls_evaluator=to_BLS_kerr)
    par = kerr_parity(l, m)
    perm = par === :even ? _KERR_EVEN_PERM : _KERR_ODD_PERM
    nf = length(perm)
    rg = T.(collect(rgrid))
    N = length(rg)
    N >= 2 || throw(ArgumentError("kerr_mode: rgrid needs at least two points"))
    all(rg[i] < rg[i+1] for i in 1:N-1) ||
        throw(ArgumentError("kerr_mode: rgrid must be strictly increasing"))
    # evaluate the 10 BLS fields on the grid (shared recs/S → one closure solve). The grid axis is
    # independent (per-point to_BLS_kerr), threaded with race-free preallocated slots.
    bls = Vector{NTuple{10,Complex{T}}}(undef, N)
    _eval!(i) = (b = bls_evaluator(a, r_p, l, m, rg[i]; recs=recs, S=S, ls=ls);
                 bls[i] = ntuple(j -> Complex{T}(b[j]), 10))
    if threaded
        blas_single_thread() do
            Threads.@threads for i in 1:N; _eval!(i); end
        end
    else
        for i in 1:N; _eval!(i); end
    end
    # Pack val_right in loader field order. Five-point Fornberg weights retain the old fourth-order
    # accuracy on a uniform grid and remain correct on arbitrary strictly increasing grids.
    valR = [Vector{Complex{T}}(undef, N) for _ in 1:nf]
    drR  = [Vector{Complex{T}}(undef, N) for _ in 1:nf]
    for (jf, jbls) in enumerate(perm)
        col = [bls[i][jbls] for i in 1:N]
        valR[jf] .= col
        for i in 1:N
            drR[jf][i] = highorder_deriv(col, rg, i; npts=min(5, N), order=1)
        end
    end
    # left grid mirrors the right grid (the full sampled grid includes r<r_p so force_mode's k±1
    # central FD always has neighbours). r0 = r_p.
    H1Mode{T}(Int(l), Int(m), rg, valR, drR, rg, valR, drR, T(r_p),
              T[], nf, par)
end

"""
    kerr_modes(a, r_p, lrange, m; T=Float64, rgrid, lmax_proj=nothing, return_build=false)
        -> Dict{Int,H1Mode{T}}  [or (modes, build) when return_build]

Assemble all Kerr (l,m) modes for l in `lrange` IN PARALLEL (one banded `closure_kerr` solve at the
projection `lmax_proj` shared across ALL l; the per-l `reconstruct_kerr` — the expensive κ chain —
threaded over the spheroidal-l list with race-free slots). Returns `l => H1Mode{T}`, drop-in for the
a=0 `load_modes`. `rgrid` is the shared radial sampling grid (must bracket the resonances + the
forced-solve span, with ≥2 guard points each end so the k±1/k±2 FD stays in-range).
`return_build=true` additionally returns the shared `(recs, S, ls)` build so callers (e.g. the §8
threaded==serial bit-identity gate) can re-evaluate `kerr_mode` without re-solving the closure.
"""
function kerr_modes(a, r_p, lrange, m::Integer; T::Type=Float64, rgrid, lmax_proj=nothing,
                    return_build::Bool=false)
    ls = collect(lrange)
    # ONE banded closure solve at the max projection lmax, shared across ALL l (recs + S built once).
    lmx = lmax_proj === nothing ? max(2*maximum(ls), maximum(ls) + 8) : lmax_proj
    sol, aux = closure_kerr(a, r_p, m; lmax=lmx, return_aux=true)
    S = smatrices_kerr(a, r_p, m; lmax=lmx)
    auxls = collect(aux.ls)
    # radial span the modes are sampled over → the high-magnitude BigFloat-interp gate in branchK.
    # Pad slightly so FD stencils / jump evals at the grid ends stay inside the interpolant domain.
    rg_lo = Float64(minimum(rgrid)) - 0.5; rg_hi = Float64(maximum(rgrid)) + 0.5
    rspan = (rg_lo, rg_hi)
    recvec = Vector{Any}(undef, length(auxls))
    blas_single_thread() do
        Threads.@threads for j in eachindex(auxls)
            ll = auxls[j]
            recvec[j] = _attach_recmeta(reconstruct_kerr(a, r_p, ll; m=m, aux=aux, sol=sol, rspan=rspan), ll, aux)
        end
    end
    recs = Dict{Int,Any}(auxls[j] => recvec[j] for j in eachindex(auxls))
    # grid assembly: the output-l loop is small (≈|lrange|) so keep it serial and let each kerr_mode
    # thread its (larger) radial grid — all cores busy on the per-point to_BLS_kerr evaluation.
    out = Vector{H1Mode{T}}(undef, length(ls))
    for i in eachindex(ls)
        out[i] = kerr_mode(a, r_p, ls[i], m; T=T, rgrid=rgrid, recs=recs, S=S, ls=auxls, threaded=true)
    end
    modes = Dict(ls[i] => out[i] for i in eachindex(ls))
    return_build ? (modes, (recs=recs, S=S, ls=auxls)) : modes
end
