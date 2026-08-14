# modesum.jl — time-domain (direct/interpolated) + frequency-domain SPA waveform mode-sum.
#
# Canonical home for `h = Σ_{l,m,n} A_{lmn}·sYlm·e^{-iΦ_{mn}}` (de-dups the FEW-vs-KerrBackwards forks of
# the SAME arithmetic). One work-item per output sample + serial inner mode reduction ⇒ CPU-serial ==
# CPU-threaded bit-identical (design §8). Built tests-first in a subsequent round.
#
# Consolidates: KerrBackwards·summation·{mode_sum,interpolated_mode_sum} + ·fdsummation·{spa_func,
# fd_flip_average,fd_interpolated_mode_sum,_SPA_*}; FEW·summation/{directmodesum,interpolatedmodesum,
# fdinterpolatedmodesum}.
#
# Time-domain kernels (`direct_mode_sum`/`interpolated_mode_sum`) landed first; the FD-SPA variant
# (`fd_interpolated_mode_sum`/`spa_func`/`fd_flip_average`) was hoisted here from the transitional
# `EMRIWaveforms.ModeSum.FEW` local copy (S4-a, §5.5 'Time-domain + FD-SPA + AAK mode-sum live in
# Kernels' / §6.8 single-copy mandate); EMRIWaveforms now delegates through thin facades.

export direct_mode_sum, interpolated_mode_sum
export fd_interpolated_mode_sum, fd_flip_average, spa_func

# Element type of the complex output: float-promote the (possibly BigFloat) real parts of every input so
# the kernel stays type-generic (BigFloat phases/amplitudes flow through; never hard-Float64'd).
_modesum_eltype(teuk, ylms, Phi_phi, Phi_r) =
    complex(float(promote_type(real(eltype(teuk)), real(eltype(ylms)),
                               eltype(Phi_phi), eltype(Phi_r))))

"""
    direct_mode_sum(teuk_modes, ylms, Phi_phi, Phi_r, m_arr, n_arr; threaded=false) -> Vector

Direct time-domain ±m mode-sum waveform (the single canonical primitive consolidating the bit-identical
FEW `DirectModeSum` and KerrBackwards `mode_sum` forks):

    h(tᵢ) = Σⱼ yⱼ Aᵢⱼ e^{-iψᵢⱼ} + Σ_{j: mⱼ>0} y'ⱼ conj(Aᵢⱼ) e^{+iψᵢⱼ},   ψᵢⱼ = mⱼ Φ_φ(tᵢ) + nⱼ Φ_r(tᵢ).

- `teuk_modes`: `nt × K` complex amplitudes `Aᵢⱼ` for the `K` base (`m ≥ 0`) modes.
- `ylms`: length `2K`. First `K` are the base harmonics `yⱼ`; last `K` the `m < 0` partner harmonics `y'ⱼ`
  (the FEW `ylmkeep = [keep_modes, temp2]` convention). Partner entries are consumed only where `mⱼ > 0`.
- `Phi_phi`, `Phi_r`: length `nt` azimuthal/radial phases. `m_arr`, `n_arr`: length `K` base mode indices.

One work-item per output sample `i` + a serial inner mode reduction ⇒ the serial and bounded-threaded
paths are bit-identical (`threaded=true` partitions samples into ≤20 chunks; design §8 bit-equality gate).
"""
function direct_mode_sum(teuk_modes::AbstractMatrix, ylms::AbstractVector,
                         Phi_phi::AbstractVector, Phi_r::AbstractVector,
                         m_arr::AbstractVector{<:Integer}, n_arr::AbstractVector{<:Integer};
                         threaded::Bool = false)
    nt, K = size(teuk_modes)
    length(ylms) == 2K ||
        throw(ArgumentError("ylms must have length 2*size(teuk_modes,2) = $(2K)"))
    (length(m_arr) == K && length(n_arr) == K) ||
        throw(ArgumentError("m_arr and n_arr must have length size(teuk_modes,2) = $K"))
    (length(Phi_phi) == nt && length(Phi_r) == nt) ||
        throw(ArgumentError("Phi_phi and Phi_r must have length size(teuk_modes,1) = $nt"))

    Tc = _modesum_eltype(teuk_modes, ylms, Phi_phi, Phi_r)
    h = Vector{Tc}(undef, nt)
    with_kernel_threads(nt; threaded = threaded) do i
        acc = zero(Tc)
        @inbounds for j in 1:K
            psi = m_arr[j] * Phi_phi[i] + n_arr[j] * Phi_r[i]
            acc += ylms[j] * teuk_modes[i, j] * cis(-psi)
            if m_arr[j] > 0
                acc += ylms[K + j] * conj(teuk_modes[i, j]) * cis(psi)
            end
        end
        @inbounds h[i] = acc
    end
    return h
end

"""
    interpolated_mode_sum(t, teuk_modes, ylms, Phi_phi, Phi_r, m_arr, n_arr;
                          dt=10.0, out_len=nothing, threaded=false) -> Vector

Time-domain interpolated mode-sum: the sparse trajectory `(t, teuk_modes, Phi_φ, Phi_r)` is fit with one
not-a-knot cubic spline (real+imag of every mode and the two phases — the FEW `y_all` stacking via
`NumericsCore.CubicSplineInterpolant`), resampled onto the uniform dense grid `tᵢ = dt·(i-1)`,
`i = 1 … out_len` (`out_len` defaults to `⌊t[end]/dt⌋ + 1`), then reduced by [`direct_mode_sum`]. Because
the spline reproduces cubic data exactly, at nodes that coincide with the dense grid this equals
`direct_mode_sum` on the raw modes. `threaded` is forwarded to the reduction.
"""
function interpolated_mode_sum(t::AbstractVector, teuk_modes::AbstractMatrix,
                               ylms::AbstractVector, Phi_phi::AbstractVector,
                               Phi_r::AbstractVector,
                               m_arr::AbstractVector{<:Integer},
                               n_arr::AbstractVector{<:Integer};
                               dt::Real = 10.0, out_len::Union{Integer,Nothing} = nothing,
                               threaded::Bool = false)
    nt, K = size(teuk_modes)
    length(t) == nt ||
        throw(ArgumentError("t must have length size(teuk_modes,1) = $nt"))
    olen = out_len === nothing ? floor(Int, t[end] / dt) + 1 : Int(out_len)

    Tr = float(promote_type(real(eltype(teuk_modes)), eltype(Phi_phi), eltype(Phi_r), eltype(t)))
    # FEW y_all stacking: rows 1..K real, K+1..2K imag, then Phi_phi, Phi_r.
    Y = Matrix{Tr}(undef, 2K + 2, nt)
    @inbounds for i in 1:nt
        for k in 1:K
            Y[k, i]     = real(teuk_modes[i, k])
            Y[K + k, i] = imag(teuk_modes[i, k])
        end
        Y[2K + 1, i] = Phi_phi[i]
        Y[2K + 2, i] = Phi_r[i]
    end
    spline = CubicSplineInterpolant(collect(Tr, t), Y)

    tm = Matrix{Complex{Tr}}(undef, olen, K)
    pp = Vector{Tr}(undef, olen)
    pr = Vector{Tr}(undef, olen)
    @inbounds for i in 1:olen
        td = dt * (i - 1)
        for k in 1:K
            tm[i, k] = complex(spline(td, k), spline(td, K + k))
        end
        pp[i] = spline(td, 2K + 1)
        pr[i] = spline(td, 2K + 2)
    end

    return direct_mode_sum(tm, ylms, pp, pr, m_arr, n_arr; threaded = threaded)
end

# =====================================================================================================
# FD-SPA — frequency-domain stationary-phase mode-sum (hoisted from EMRIWaveforms.ModeSum.FEW, S4-a).
#
# The SPA primitives (`spa_func`/`_cube_roots`/`_fd_seg_freq_range`, series constants) are moved
# VERBATIM. The mode-sum itself is re-expressed in GATHER form — one work-item per output frequency
# bin + a serial inner (mode, segment, root) reduction — so it satisfies the §8 parallelism contract
# (CPU-serial == CPU-threaded BIT-IDENTICAL). The legacy kernel was a SCATTER (each (mode, segment)
# writes its visited bins AND the mirror ±m partner bins); the gather replicates its per-bin
# floating-point accumulation ORDER exactly (contributions sorted by (mode, segment, source-bin), the
# i==zidx direct/partner root-interleave included), so the output is bit-identical to the pre-hoist
# fixture, not merely close. Float64 tier by construction (the SPA series constants are Float64;
# EMRIWaveforms is the Float64/GPU tier — no BigFloat inheritance requirement here).
# =====================================================================================================

const _SPA_GAMP13 = 2.67893853470774763   # Γ(1/3)
const _SPA_GAMM13 = -4.06235381827920125  # Γ(-1/3)

const _SPA_C1 = (0.5, -0.09375, 0.0050223214285714285714, -0.00012555803571428571429,
    1.8109332074175824176e-6, -1.6977498819539835165e-8, 1.1169407118118312608e-10,
    -5.4396463237589184781e-13, 2.0398673714095944293e-15, -6.0710338434809358015e-18,
    1.4687985105195812423e-20, -2.9454515585285720100e-23, 4.9754249299469121790e-26,
    -7.1760936489618925658e-29)
const _SPA_C2 = (1.0, -0.375, 0.028125, -0.00087890625, 0.000014981356534090909091,
    -1.6051453429383116883e-7, 1.1802539286311115355e-9, -6.3227889033809546546e-12,
    2.5772237377911499951e-14, -8.2603324929203525483e-17, 2.1362928861000911763e-19,
    -4.5517604107246260858e-22, 8.1281435905796894390e-25, -1.2340298973552160079e-27)
const _SPA_CL = (im, 0.069444444444444444444, -0.037133487654320987654im,
    -0.037993059127800640146, 0.057649190412669721333im, 0.11609906402551541102,
    -0.29159139923075051147im, -0.87766696951001691647, 3.0794530301731669934im)

"""
    spa_func(x::Real) -> ComplexF64

Stationary-phase function `i√x e^{-ix} K_{1/3}(-ix)` (legacy `interpolate.cu::SPAFunc`): the
small-`|x|≤7` ascending series and the large-`|x|>7` asymptotic series of the upstream kernel.
"""
function spa_func(x::Real)
    if abs(x) <= 7.0
        xx = complex(float(x))
        pref1 = exp(-2.0 * π * im / 3.0) * xx^(5.0 / 6.0) * _SPA_GAMM13 / 2.0^(1.0 / 3.0)
        pref2 = exp(-π * im / 3.0) * xx^(1.0 / 6.0) * _SPA_GAMP13 / 2.0^(2.0 / 3.0)
        x2 = x * x
        ser1 = _SPA_C1[14]; ser2 = _SPA_C2[14]
        for kk in 13:-1:1
            ser1 = _SPA_C1[kk] + x2 * ser1
            ser2 = _SPA_C2[kk] + x2 * ser2
        end
        return exp(-im * x) * (pref1 * ser1 + pref2 * ser2)
    else
        y = 1.0 / x
        pref = exp(-0.75 * im * π) * sqrt(0.5 * π)
        ser = _SPA_CL[9]
        for kk in 8:-1:1
            ser = _SPA_CL[kk] + y * ser
        end
        return pref * ser
    end
end

"""
    _cube_roots(a, b, c, d) -> Tuple

Finite roots of `a x³ + b x² + c x + d = 0`. Genuine cubics use the legacy Cardano arithmetic;
exactly degenerate leading coefficients dispatch to the quadratic or linear formula. A nonzero
constant has no roots, while the identically-zero polynomial has no isolated roots, so both return
the empty tuple.
"""
function _cube_roots(a, b, c, d)
    if iszero(a)
        if iszero(b)
            return iszero(c) ? () : (complex(-d / c),)
        end
        disc = sqrt(complex(c * c - 4 * b * d))
        nplus = -c + disc
        nminus = -c - disc
        q = (abs(nplus) >= abs(nminus) ? nplus : nminus) / 2
        if iszero(q)
            z = complex(zero(q))
            return (z, z)
        end
        return (q / b, d / q)
    end
    b2 = b * b; b3 = b2 * b; a2 = a * a
    Delta0 = b2 - 3 * a * c
    Delta1 = 2 * b3 - 9 * a * b * c + 27 * a2 * d
    sqrt_disc = sqrt(complex(Delta1 * Delta1 - 4 * Delta0^3))
    C3 = (Delta1 + sqrt_disc) / 2
    # Cardano's conventional `+sqrt` branch cancels exactly for cubics such as
    # x^3 - 1. The conjugate branch represents the same root set without the 0/0.
    iszero(C3) && (C3 = (Delta1 - sqrt_disc) / 2)
    C = C3^(one(a) / 3)
    if iszero(C)
        iszero(Delta0) && iszero(Delta1) ||
            error("_cube_roots: Cardano branch cancellation produced a zero divisor")
        r = complex(-b / (3a))
        return (r, r, r)
    end
    xi = (-one(a) + sqrt(complex(-3 * one(a)))) / 2
    r1 = -(one(a) / (3a)) * (b + C + Delta0 / C)
    r2 = -(one(a) / (3a)) * (b + xi * C + Delta0 / (xi * C))
    r3 = -(one(a) / (3a)) * (b + xi * xi * C + Delta0 / (xi * xi * C))
    return (r1, r2, r3)
end

"""
    _fd_seg_freq_range(a, b, c, d, L) -> (fmin, fmax)

Turnover-aware frequency range of the cubic `F(x) = d + c x + b x² + a x³` over `x ∈ [0, L]`,
including any interior critical point where `F'(x) = 0` (legacy `fdinterp.py` derivative root-solve).
"""
function _fd_seg_freq_range(a, b, c, d, L)
    F(x) = d + c * x + b * x^2 + a * x^3
    fmin = min(F(0.0), F(L)); fmax = max(F(0.0), F(L))
    inner = (2b)^2 - 4 * (3a) * c
    if inner >= 0 && a != 0
        s = sqrt(inner)
        for x in ((-(2b) + s) / (2 * (3a)), (-(2b) - s) / (2 * (3a)))
            if 0.0 < x < L
                fx = F(x); fmin = min(fmin, fx); fmax = max(fmax, fx)
            end
        end
    elseif iszero(a) && !iszero(b)
        x = -c / (2b)
        if 0.0 < x < L
            fx = F(x); fmin = min(fmin, fx); fmax = max(fmax, fx)
        end
    end
    return (fmin, fmax)
end

"""
    fd_flip_average(waveform, frequency) -> (hplus, hcross)

Convert the raw `±f` SPA spectrum into `h₊(f)`/`h₍×₎(f)` via the flip-average symmetry (legacy
`fdinterp.py`): `s = -reverse(waveform)`, then `h₊ = (Re(s+rev s) + i Im(s−rev s))/2`,
`h₍×₎ = (−Im(s+rev s) + i Re(s−rev s))/2`.
"""
function fd_flip_average(waveform::AbstractVector, frequency::AbstractVector)
    fd = -reverse(waveform)
    fr = reverse(fd)
    hp = @. real(fd + fr) / 2 + im * imag(fd - fr) / 2
    hc = @. -imag(fd + fr) / 2 + im * real(fd - fr) / 2
    return (hp, hc)
end

# SPA value of ONE stationary-phase root: given the spline `sp` (rows 1..K re, K+1..2K im, 2K+1 f_φ,
# 2K+2 f_r, 2K+3 Φ_φ, 2K+4 Φ_r), mode `mi` (indices `m`,`n`), segment `seg` spanning [tlo, thi), output
# frequency `f` and root offset `x`, return (true, R·amp·e^{i(2πf t⋆ − ψ)}) if the root lies in the
# segment — expressions verbatim from the legacy kernel (bit-identical contributions).
@inline function _fd_spa_root(sp, K::Int, mi::Int, seg::Int, m, n, tlo, thi, f::Float64, x::Float64)
    tt = tlo + x
    (tt < tlo || tt >= thi) && return (false, zero(ComplexF64))
    x2 = x * x; x3 = x * x2
    @inbounds begin
        R_re = sp.y[mi, seg] + sp.c1[mi, seg] * x + sp.c2[mi, seg] * x2 + sp.c3[mi, seg] * x3
        R_im = sp.y[K+mi, seg] + sp.c1[K+mi, seg] * x + sp.c2[K+mi, seg] * x2 + sp.c3[K+mi, seg] * x3
        Phi_phi_i = sp.y[2K+3, seg] + sp.c1[2K+3, seg] * x + sp.c2[2K+3, seg] * x2 + sp.c3[2K+3, seg] * x3
        Phi_r_i   = sp.y[2K+4, seg] + sp.c1[2K+4, seg] * x + sp.c2[2K+4, seg] * x2 + sp.c3[2K+4, seg] * x3
        fdot_phi = sp.c1[2K+1, seg] + 2 * sp.c2[2K+1, seg] * x + 3 * sp.c3[2K+1, seg] * x2
        fdot_r   = sp.c1[2K+2, seg] + 2 * sp.c2[2K+2, seg] * x + 3 * sp.c3[2K+2, seg] * x2
        fddot_phi = 2 * sp.c2[2K+1, seg] + 6 * sp.c3[2K+1, seg] * x
        fddot_r   = 2 * sp.c2[2K+2, seg] + 6 * sp.c3[2K+2, seg] * x
    end
    phase_term = m * Phi_phi_i + n * Phi_r_i
    fdot  = m * fdot_phi  + n * fdot_r
    fddot = m * fddot_phi + n * fddot_r
    temp_arg = 2.0 * π * fdot^3 / (3.0 * fddot^2)
    amp = -1.0 * fdot / abs(fddot) * 2.0 / sqrt(3.0) *
          spa_func(temp_arg) / sqrt(complex(temp_arg))
    return (true, complex(R_re, R_im) * amp * exp(im * (2.0 * π * f * tt - phase_term)))
end

# One (mode, segment) accumulation pass at output frequency `f`: solve the stationary cubic
# `f_mode(x) = f`, accumulate `w·te` (direct) or `w·conj(te)` (±m partner) for every in-segment root
# in Cardano root order — the exact per-source-bin event the legacy scatter emitted.
@inline function _fd_spa_pass(acc::ComplexF64, sp, K, mi, seg, m, n, tlo, thi,
                              c3::Float64, c2::Float64, c1::Float64, y0::Float64,
                              f::Float64, w, conjugate::Bool)
    roots = _cube_roots(c3, c2, c1, y0 - f)
    for r in roots
        ok, te = _fd_spa_root(sp, K, mi, seg, m, n, tlo, thi, f, real(r))
        ok || continue
        acc += conjugate ? w * conj(te) : w * te
    end
    return acc
end

"""
    fd_interpolated_mode_sum(t, teuk_modes, ylms, Phi_phi, Phi_r, f_phi, f_r, m_arr, n_arr;
                             f_arr, include_minus_m=true, separate_modes=false,
                             threaded=false) -> (hplus, hcross)

Frequency-domain stationary-phase mode-sum (single canonical copy of the FEW
`FDInterpolatedModeSum.sum` / KerrBackwards `fdsummation` kernel, §5.5). One not-a-knot cubic spline
(`NumericsCore.CubicSplineInterpolant`) fits the real/imag mode amplitudes, the mode frequencies
`f_phi, f_r`, and the phases `Φ_φ, Φ_r`; for every output bin of `f_arr`, mode and trajectory segment,
the turnover-aware segment frequency range ([`_fd_seg_freq_range`]) decides membership, the stationary
time solves the cubic `f_mode(t)=f` ([`_cube_roots`]), and the SPA value `R·amp·e^{i(2πft⋆−ψ)}` is
accumulated with the `K_{1/3}` correction ([`spa_func`]).

`include_minus_m` adds the `-m` partner (`conj`, mirror bin about the mandatory exact-`0.0` bin of
`f_arr`). `separate_modes=true` returns the raw per-mode `±f` spectrum matrix (`K × length(f_arr)`;
no flip-average, `include_minus_m` forced off); else `(h₊, h₍×₎)` via [`fd_flip_average`].

GATHER form: one work-item per output frequency bin + serial inner (mode, segment, root) reduction,
replicating the legacy scatter's per-bin accumulation order exactly ⇒ `threaded=true` is
BIT-IDENTICAL to serial (§8 gate) and both are bit-identical to the pre-hoist EMRIWaveforms kernel.
"""
function fd_interpolated_mode_sum(t::AbstractVector, teuk_modes::AbstractMatrix,
                                  ylms::AbstractVector, Phi_phi::AbstractVector,
                                  Phi_r::AbstractVector, f_phi::AbstractVector,
                                  f_r::AbstractVector, m_arr::AbstractVector{<:Integer},
                                  n_arr::AbstractVector{<:Integer};
                                  f_arr::AbstractVector, include_minus_m::Bool = true,
                                  separate_modes::Bool = false, threaded::Bool = false)
    nt, K = size(teuk_modes)
    nseg = nt - 1
    length(t) == nt ||
        throw(ArgumentError("t must have length size(teuk_modes,1) = $nt"))
    length(ylms) == 2K ||
        throw(ArgumentError("ylms must have length 2*size(teuk_modes,2) = $(2K)"))
    (length(m_arr) == K && length(n_arr) == K) ||
        throw(ArgumentError("m_arr and n_arr must have length size(teuk_modes,2) = $K"))
    (length(Phi_phi) == nt && length(Phi_r) == nt && length(f_phi) == nt && length(f_r) == nt) ||
        throw(ArgumentError("Phi_phi, Phi_r, f_phi, f_r must have length size(teuk_modes,1) = $nt"))

    freq = collect(float.(f_arr))
    nf = length(freq)
    zidx = findfirst(==(0.0), freq)
    zidx === nothing && throw(ArgumentError("f_arr must contain exactly one 0.0 bin"))
    f_first = freq[1]
    df = freq[2] - freq[1]

    # single not-a-knot spline: re(1..K), im(K+1..2K), f_phi, f_r, Phi_phi, Phi_r.
    Y = Matrix{Float64}(undef, 2K + 4, nt)
    @views Y[1:K, :]    .= permutedims(real.(teuk_modes))
    @views Y[K+1:2K, :] .= permutedims(imag.(teuk_modes))
    @views Y[2K+1, :]   .= f_phi
    @views Y[2K+2, :]   .= f_r
    @views Y[2K+3, :]   .= Phi_phi
    @views Y[2K+4, :]   .= Phi_r
    sp = CubicSplineInterpolant(collect(Float64, t), Y)

    include_mm = include_minus_m && !separate_modes   # legacy forces the partner off per-mode

    # precompute per (mode, segment): mode-frequency cubic coefficients + visited source-bin range
    # [bl, bh] in the legacy b-index (source bin i ⇔ b = i-1), exactly the scatter's loop bounds.
    fy  = Matrix{Float64}(undef, K, nseg); fc1 = similar(fy); fc2 = similar(fy); fc3 = similar(fy)
    bl  = Matrix{Int}(undef, K, nseg);     bh  = Matrix{Int}(undef, K, nseg)
    @inbounds for mi in 1:K, seg in 1:nseg
        m = m_arr[mi]; n = n_arr[mi]
        L = t[seg+1] - t[seg]
        fy[mi, seg]  = m * sp.y[2K+1, seg]  + n * sp.y[2K+2, seg]
        fc1[mi, seg] = m * sp.c1[2K+1, seg] + n * sp.c1[2K+2, seg]
        fc2[mi, seg] = m * sp.c2[2K+1, seg] + n * sp.c2[2K+2, seg]
        fc3[mi, seg] = m * sp.c3[2K+1, seg] + n * sp.c3[2K+2, seg]
        fmin, fmax = _fd_seg_freq_range(fc3[mi, seg], fc2[mi, seg], fc1[mi, seg], fy[mi, seg], L)
        bl[mi, seg] = trunc(Int, abs((fmin - f_first) / df)) + 1
        bh[mi, seg] = trunc(Int, abs((fmax - f_first) / df))
    end

    if separate_modes
        wfm = zeros(ComplexF64, K, nf)
        with_kernel_threads(nf; threaded = threaded) do i
            f = freq[i]
            @inbounds for mi in 1:K
                m = m_arr[mi]; n = n_arr[mi]; yp = ylms[mi]
                for seg in 1:nseg
                    bl[mi, seg] <= i - 1 <= bh[mi, seg] || continue
                    roots = _cube_roots(fc3[mi, seg], fc2[mi, seg], fc1[mi, seg], fy[mi, seg] - f)
                    for r in roots
                        ok, te = _fd_spa_root(sp, K, mi, seg, m, n, t[seg], t[seg+1], f, real(r))
                        ok && (wfm[mi, i] = yp * te)   # legacy assignment semantics: last root wins
                    end
                end
            end
        end
        return wfm
    end

    wf = zeros(ComplexF64, nf)
    with_kernel_threads(nf; threaded = threaded) do i
        acc = zero(ComplexF64)
        j = 2 * zidx - i    # mirror source bin whose ±m partner lands on bin i (involution about zidx)
        @inbounds for mi in 1:K
            m = m_arr[mi]; n = n_arr[mi]
            yp = ylms[mi]; ym = ylms[K + mi]
            partner = include_mm && m != 0 && 1 <= j <= nf
            for seg in 1:nseg
                dir_in = bl[mi, seg] <= i - 1 <= bh[mi, seg]
                par_in = partner && (bl[mi, seg] <= j - 1 <= bh[mi, seg])
                (dir_in || par_in) || continue
                tlo = t[seg]; thi = t[seg+1]
                c3 = fc3[mi, seg]; c2 = fc2[mi, seg]; c1 = fc1[mi, seg]; y0 = fy[mi, seg]
                if dir_in && par_in && i == j
                    # zero bin: the legacy scatter emitted direct then partner PER ROOT — interleave.
                    roots = _cube_roots(c3, c2, c1, y0 - freq[i])
                    for r in roots
                        ok, te = _fd_spa_root(sp, K, mi, seg, m, n, tlo, thi, freq[i], real(r))
                        ok || continue
                        acc += yp * te
                        acc += ym * conj(te)
                    end
                elseif dir_in && (!par_in || i < j)
                    # direct events precede the partner's source bin (or no partner): b ascending.
                    acc = _fd_spa_pass(acc, sp, K, mi, seg, m, n, tlo, thi, c3, c2, c1, y0,
                                       freq[i], yp, false)
                    par_in && (acc = _fd_spa_pass(acc, sp, K, mi, seg, m, n, tlo, thi, c3, c2, c1,
                                                  y0, freq[j], ym, true))
                else
                    # partner's source bin precedes the direct events (or no direct).
                    acc = _fd_spa_pass(acc, sp, K, mi, seg, m, n, tlo, thi, c3, c2, c1, y0,
                                       freq[j], ym, true)
                    dir_in && (acc = _fd_spa_pass(acc, sp, K, mi, seg, m, n, tlo, thi, c3, c2, c1,
                                                  y0, freq[i], yp, false))
                end
            end
        end
        @inbounds wf[i] = acc
    end
    return fd_flip_average(wf, freq)
end
