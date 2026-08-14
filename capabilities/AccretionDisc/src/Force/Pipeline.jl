# PIPELINE — forced master-ODE solve + F8 torque-balance / flux-crosscheck drivers (design §4.10
# Force: "forced master-ODE solve + torque-balance/flux cross-checks").
# From DiscEMRIForce/pipeline.jl (h1lorenz forced_solve / f8_torque_balance / f8_flux_crosscheck / f5_source).
#
#   master ODE:  𝔥'' = 𝔉 − A_m 𝔥' − B_m 𝔥      (A_m, B_m REAL ⇒ exact ang-mom conservation law)
#   integrating factor p(r) = exp(∫A_m dr);   d/dr[p·Im(𝔥*𝔥')] = p·Im(𝔥*𝔉)
#     LHS = ∂_r J̇^r_Adv   (advected-flux divergence of the excited density wave)
#     RHS = ∂T^φ_Mat       (matter torque density deposited by the metric-reconstructed forcing 𝔉)
#   The match ∂_r J̇^r_Adv == ∂T^φ_Mat IS the F8 torque balance (~1e-7).
#
# 𝒬, A_m, B_m_full, lindblad_locations, risco_disc come from ..Modes; the forcing 𝔉 from Source.jl.
# AHA: the legacy in-file `_risco_disc` Bardeen fork is DROPPED — `Modes.risco_disc` (→ the canonical
# Spacetime.kerr_isco) is the single ISCO home (design §4.10 collisions). TYPE-GENERIC {T}.
# PARALLEL (§8): the 𝔉(r) grid precompute is THREADED over grid points (race-free per-slot output)
# under the NumericsCore.blas_single_thread guard (the 𝒬 4×4 LAPACK inv inside a threaded sweep) —
# threaded == serial BIT-IDENTICAL.

using Printf
using ..Modes: RelativisticDisc, disc_OmegaK, A_m, B_m_full, lindblad_locations, risco_disc, csound
import NumericsCore: blas_single_thread

# ---- 𝔉(r) precompute on a grid (THREADED) + linear interpolant -----------------------------------------------
"""
    ForcingInterp{T}

A linear interpolant of the radial forcing 𝔉(r) sampled on grid `rg` with values `Fg`.
"""
struct ForcingInterp{T<:Real}
    rg::Vector{T}; Fg::Vector{Complex{T}}
end

"""
    build_forcing(d, modes, m, r_p, ω0, rg; lmax) -> ForcingInterp{T}

Precompute 𝔉(r) on grid `rg` IN PARALLEL across threads (per-grid-point source assembly, race-free
preallocated slots; BLAS pinned to one thread for the scope — §8 guard), returning a linear
interpolant. Threaded == serial bit-identical (one writer per slot, deterministic single-thread BLAS).
"""
function build_forcing(d::RelativisticDisc, modes::Dict{Int,H1Mode{T}}, m::Integer, r_p, ω0, rg;
                       lmax::Integer) where {T}
    g = modes[m].r_right
    rgT = T.(collect(rg))
    Fg = Vector{Complex{T}}(undef, length(rgT))
    blas_single_thread() do
        Threads.@threads for i in eachindex(rgT)
            k = argmin(abs.(g .- rgT[i]))
            Fg[i] = source_F(d, modes, k, m, r_p, ω0; lmax=lmax)
        end
    end
    ForcingInterp{T}(rgT, Fg)
end

"Evaluate the forcing interpolant at radius `r` (endpoint-clamped linear interpolation)."
function (F::ForcingInterp{T})(r) where {T}
    rT = T(r)
    rT <= F.rg[1] && return F.Fg[1]
    rT >= F.rg[end] && return F.Fg[end]
    i = searchsortedlast(F.rg, rT); t = (rT - F.rg[i])/(F.rg[i+1] - F.rg[i])
    (1 - t)*F.Fg[i] + t*F.Fg[i+1]
end

# ---- forced master-ODE solve (RK4 inward, source-driven, outgoing BC) ----------------------------------------
"""
    forced_solve(d, m, r_p, F, r_start, r_end, N) -> (rs, hs, hps)

Integrate the forced master ODE 𝔥'' = 𝔉 − A_m 𝔥' − B_m 𝔥 by RK4 from `r_start` inward to `r_end`
(N steps), starting from 𝔥=𝔥'=0 (no incoming wave) so the solution is purely source-driven.
Returns the radial grid and the complex 𝔥, 𝔥' trajectories. From DiscEMRIForce/pipeline.jl (h1lorenz rk4).
"""
function forced_solve(d::RelativisticDisc, m::Integer, r_p, F::ForcingInterp{T}, r_start, r_end, N::Integer) where {T}
    Am(r) = A_m(d, r, r_p, m); Bm(r) = B_m_full(d, r, r_p, m)
    rhs(r, y) = Complex{T}[y[2], F(r) - Am(r)*y[2] - Bm(r)*y[1]]
    h = (T(r_end) - T(r_start))/N; r = T(r_start); y = Complex{T}[0, 0]
    rs = Vector{T}(undef, N+1); hs = Vector{Complex{T}}(undef, N+1); hps = Vector{Complex{T}}(undef, N+1)
    rs[1] = r; hs[1] = y[1]; hps[1] = y[2]
    for i in 1:N
        k1 = rhs(r, y); k2 = rhs(r + h/2, y + h/2*k1); k3 = rhs(r + h/2, y + h/2*k2); k4 = rhs(r + h, y + h*k3)
        y = y + h/6*(k1 + 2k2 + 2k3 + k4); r += h
        rs[i+1] = r; hs[i+1] = y[1]; hps[i+1] = y[2]
    end
    rs, hs, hps
end

# =============================================================================================================
# F8 TORQUE BALANCE  ∂_r J̇^r_Adv  ==  ∂T^φ_Mat   (the primary G-D2 gate)
# =============================================================================================================

"""
    f8_torque_balance(d, m, r_p, modes, lmax; rg, N, r_start, r_end, rmin, rmax, verbose)
        -> (mean_rel, max_rel, mean_rel_integral, max_rel_integral, npts, rs, dJ, TRq)

The F8 angular-momentum torque balance: integrate the forced master ODE, form the two independently
computed torque densities ∂_r J̇^r_Adv (advected-flux divergence) and ∂T^φ_Mat (forcing torque), and
report their relative agreement over the physical outer region [rmin,rmax] (differential + integral
forms). G-D2 oracle ~1e-7. From DiscEMRIForce/pipeline.jl.
"""
function f8_torque_balance(d::RelativisticDisc, m::Integer, r_p, modes::Dict{Int,H1Mode{T}}, lmax::Integer;
                           rg=10.2:0.8:16.0, N::Integer=16000, r_start=16.0, r_end=10.4,
                           rmin=10.6, rmax=15.4, verbose::Bool=true) where {T}
    ω0 = m*disc_OmegaK(d, r_p)
    F = build_forcing(d, modes, m, r_p, ω0, rg; lmax=lmax)
    rs, hs, hps = forced_solve(d, m, r_p, F, r_start, r_end, N)
    Am(r) = A_m(d, r, r_p, m); Nn = length(rs)
    # integrating factor p(r)
    ps = zeros(T, Nn); ps[1] = one(T)
    for i in 2:Nn; ps[i] = ps[i-1]*exp(real((Am(rs[i]) + Am(rs[i-1]))/2)*(rs[i] - rs[i-1])); end
    Jw(i)  = ps[i]*imag(conj(hs[i])*hps[i])       # advected flux J̇^r_Adv
    TRq(i) = ps[i]*imag(conj(hs[i])*F(rs[i]))      # torque density ∂T^φ_Mat
    dJ(i)  = (Jw(i+1) - Jw(i-1))/(rs[i+1] - rs[i-1])
    rI, rO = lindblad_locations(d, r_p, m)
    idx = Int[]; for i in 2:Nn-1; (rs[i] < T(rmax) && rs[i] > T(rmin)) && push!(idx, i); end
    if verbose
        @printf("\n=== F8 TORQUE BALANCE  ∂_r J̇^r_Adv == ∂T^φ_Mat  (a=%.2f, m=%d, r_p=%.1f, OLR=%.3f) ===\n",
                Float64(d.a), m, Float64(r_p), Float64(rO))
        @printf("%8s | %16s | %16s | rel.diff\n", "r", "∂_rJ̇_Adv", "∂T^φ_Mat")
    end
    maxrel = zero(T); sr = zero(T); nc = 0
    tinytol = T(1)/T(10)^12
    for i in idx[1:max(1, length(idx)÷10):end]
        rel = abs(dJ(i) - TRq(i))/(abs(TRq(i)) + eps(T))
        if abs(TRq(i)) > tinytol; maxrel = max(maxrel, rel); sr += rel; nc += 1; end
        verbose && @printf("%8.3f | %16.6e | %16.6e | %.2e\n",
                           Float64(rs[i]), Float64(dJ(i)), Float64(TRq(i)), Float64(rel))
    end
    meanrel = nc > 0 ? sr/nc : T(NaN)
    # integral form: J̇_Adv(r) = ∫ from r_start of ∂T  (Jw at i=1 is the BC)
    Jcum = zeros(T, Nn); Jcum[1] = Jw(1)
    for i in 2:Nn; Jcum[i] = Jcum[i-1] + (TRq(i) + TRq(i-1))/2*(rs[i] - rs[i-1]); end
    relint = [abs(Jcum[i] - Jw(i))/(abs(Jw(i)) + eps(T)) for i in idx]
    meanint = sum(relint)/length(relint); maxint = maximum(relint)
    if verbose
        @printf("\nmean rel.diff = %.2e (max %.2e) over %d points with |∂T|>1e-12\n",
                Float64(meanrel), Float64(maxrel), nc)
        @printf("integral form J̇_Adv(r)=∫∂T:  mean rel.diff = %.2e (max %.2e)\n",
                Float64(meanint), Float64(maxint))
    end
    (mean_rel=meanrel, max_rel=maxrel, mean_rel_integral=meanint, max_rel_integral=maxint,
     npts=nc, rs=[rs[i] for i in idx], dJ=[dJ(i) for i in idx], TRq=[TRq(i) for i in idx])
end

"""
    f8_torque_balance_ilr(d, m, r_p, modes, lmax; gap, rstart_off, rend_off, N, verbose) -> (same NamedTuple)

The INNER (ILR-side) F8 torque balance: the SAME operator identity d/dr[p·Im(𝔥*𝔥')]=p·Im(𝔥*𝔉) run
on the inner region [ISCO, r_p) where the inner Lindblad resonance launches the inward-propagating
density wave. Auto-selects the inner sampling grid + forced-solve span from the ISCO and the mode
grid; integrates from just below r_p inward toward the ISCO. The F8 closure is direction/side-
agnostic. (The mode grid must sample r<r_p; the legacy honest-scope note carries: the off-equatorial
source assembly inside the corotation band is slow/fragile at a≠0 — a manual driver, not a CI gate.)
"""
function f8_torque_balance_ilr(d::RelativisticDisc, m::Integer, r_p, modes::Dict{Int,H1Mode{T}}, lmax::Integer;
                               gap=T(2)/T(10), rstart_off=T(3)/T(10), rend_off=T(3)/T(10),
                               N::Integer=16000, verbose::Bool=true) where {T}
    rI, rO = lindblad_locations(d, r_p, m)
    risco = risco_disc(d)
    # inner physical window: a bit inside r_p down to a bit outside the ISCO; within the sampled grid.
    g = modes[m].r_right
    rinner_lo = max(T(risco) + T(rend_off), T(minimum(g)) + T(2)/T(10))
    rinner_hi = T(r_p) - T(rstart_off)
    rmin = rinner_lo + gap; rmax = rinner_hi - gap
    rg = collect(range(rinner_lo, rinner_hi; length=10))
    f8_torque_balance(d, m, r_p, modes, lmax; rg=rg, N=N,
                      r_start=Float64(rinner_hi), r_end=Float64(rinner_lo),
                      rmin=Float64(rmin), rmax=Float64(rmax), verbose=verbose)
end

# =============================================================================================================
# F8 FLUX CROSS-CHECK   J_fluid ∝ J_wave  (the G-D3 gate)
# =============================================================================================================

"""
    f8_flux_crosscheck(d, m, r_p, modes, lmax; rs_sample, rg, N, r_start, r_end, verbose)
        -> (cv_full, cv_wave, ratios_full, ratios_wave)

Cross-check the advected ang-mom flux two independent ways: J_wave = p·Im(𝔥*𝔥') (wave-action form)
vs J_fluid = (ε+p)Re[δu^r δu_φ*] from the FULL fluid reconstruction. The wave-part reconstruction
tracks J_wave to a smooth O(1) WKB factor (small CV); the full reconstruction is forcing-dominated
(large CV). G-D3 oracle ~const. From DiscEMRIForce/pipeline.jl.
"""
function f8_flux_crosscheck(d::RelativisticDisc, m::Integer, r_p, modes::Dict{Int,H1Mode{T}}, lmax::Integer;
                            rs_sample=[11.0,11.6,12.2,12.8,13.6,14.4], rg=10.2:0.8:16.0,
                            N::Integer=12000, r_start=16.0, r_end=10.4, verbose::Bool=true) where {T}
    ω0 = m*disc_OmegaK(d, r_p)
    F = build_forcing(d, modes, m, r_p, ω0, rg; lmax=lmax)
    rs, hs, hps = forced_solve(d, m, r_p, F, r_start, r_end, N)
    Am(r) = A_m(d, r, r_p, m); Nn = length(rs)
    ps = zeros(T, Nn); ps[1] = one(T)
    for i in 2:Nn; ps[i] = ps[i-1]*exp(real((Am(rs[i]) + Am(rs[i-1]))/2)*(rs[i] - rs[i-1])); end
    g = modes[m].r_right
    rf = T[]; rw = T[]
    verbose && @printf("\n=== F8 flux cross-check: J_full / J_wavepart  vs  J_wave = p·Im(𝔥*𝔥') ===\n")
    verbose && @printf("%7s | %13s | %13s %10s | %13s %10s\n",
                       "r", "J_wave", "J_full", "r[full]", "J_wavepart", "r[wave]")
    for rr in rs_sample
        i = argmin(abs.(rs .- T(rr))); k = argmin(abs.(g .- rs[i]))
        Jwave = ps[i]*imag(conj(hs[i])*hps[i])
        Jf = Jadv(d, modes, k, m, r_p, ω0, hs[i], hps[i]; lmax=lmax, wave=false)
        Jwp = Jadv(d, modes, k, m, r_p, ω0, hs[i], hps[i]; lmax=lmax, wave=true)
        push!(rf, Jf/Jwave); push!(rw, Jwp/Jwave)
        verbose && @printf("%7.2f | %13.4e | %13.4e %10.2e | %13.4e %10.2e\n",
                           Float64(rs[i]), Float64(Jwave), Float64(Jf), Float64(Jf/Jwave),
                           Float64(Jwp), Float64(Jwp/Jwave))
    end
    cv(x) = (mu = sum(x)/length(x); sqrt(sum((x .- mu).^2)/length(x))/abs(mu))
    cvf = cv(rf); cvw = cv(rw)
    if verbose
        @printf("\nfull δU      ratio CV = %.2e (forcing-dominated)\n", Float64(cvf))
        @printf("wave-part δU ratio: %.2f→%.2f, CV = %.2e (tracks the wave-action flux)\n",
                Float64(rw[1]), Float64(rw[end]), Float64(cvw))
    end
    (cv_full=cvf, cv_wave=cvw, ratios_full=rf, ratios_wave=rw)
end

# =============================================================================================================
# F5 SOURCE STRUCTURE  |𝔉(r)| peaking at the OLR  (the source-plot driver)
# =============================================================================================================

"""
    f5_source(d, m, r_p, modes, lmax; rg, verbose, highorder, npts, spectral, smooth_lo, smooth_hi, cheb_deg)
        -> (rgrid, absF, r_peak, OLR)

Sample |𝔉(r)| (the radial forcing amplitude) over `rg` (OLR-side). The forcing peaks at the OLR (the
Lindblad resonance that launches the outgoing density wave). `highorder=true` switches to the
wing-precision `source_F_grid` array pass (optionally `spectral`). THREADED over the sample grid
(race-free per-slot; blas_single_thread guard). From DiscEMRIForce/pipeline.jl (h1lorenz f5_source.jl).
"""
function f5_source(d::RelativisticDisc, m::Integer, r_p, modes::Dict{Int,H1Mode{T}}, lmax::Integer;
                   rg=8.2:0.25:15.0, verbose::Bool=true,
                   highorder::Bool=false, npts::Integer=9, spectral::Bool=false,
                   smooth_lo=nothing, smooth_hi=nothing, cheb_deg::Integer=8) where {T}
    ω0 = m*disc_OmegaK(d, r_p); g = modes[m].r_right
    rI, rO = lindblad_locations(d, r_p, m)
    rgrid = T.(collect(rg)); absF = Vector{T}(undef, length(rgrid))
    if highorder
        # WING-PRECISION array pass: one high-order/spectral source evaluation on the WHOLE grid, then sample.
        Fg = source_F_grid(d, modes, m, r_p, ω0; lmax=lmax, npts=npts, spectral=spectral,
                           smooth_lo=smooth_lo, smooth_hi=smooth_hi, cheb_deg=cheb_deg)
        Threads.@threads for i in eachindex(rgrid)
            k = argmin(abs.(g .- rgrid[i]))
            absF[i] = abs(Fg[k])
        end
    else
        blas_single_thread() do
            Threads.@threads for i in eachindex(rgrid)
                k = argmin(abs.(g .- rgrid[i]))
                absF[i] = abs(source_F(d, modes, k, m, r_p, ω0; lmax=lmax))
            end
        end
    end
    ipk = argmax(absF)
    if verbose
        @printf("\n=== F5 source |𝔉(r)| (a=%.2f, m=%d, r_p=%.1f, OLR=%.3f) ===\n",
                Float64(d.a), m, Float64(r_p), Float64(rO))
        @printf("peak |𝔉| at r=%.2f (OLR=%.2f); |𝔉(peak)|=%.3e\n",
                Float64(rgrid[ipk]), Float64(rO), Float64(absF[ipk]))
    end
    (rgrid=rgrid, absF=absF, r_peak=rgrid[ipk], OLR=rO)
end
