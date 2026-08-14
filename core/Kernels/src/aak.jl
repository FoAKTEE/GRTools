# aak.jl — Augmented-Analytic-Kludge harmonic summation.
#
# Canonical home for the AAK kludge mode-sum (`aak_summation`/`_aak_point`/`_aak_rot_coeff`/
# `_aak_lam_fixup!`); one work-item per output sample + serial harmonic reduction (parallel==serial
# bit-identity gate, design §8).
#
# Consolidates: KerrBackwards·summation (AAK path); FEW·summation/aaksummation·{aak_ak_OmegaPhi,
# aak_waveform}. The trajectory-driven `aak_waveform` (run_inspiral/kerr_frequencies) AND the
# Peters–Mathews Bessel builder `ak_bessel_coeffs` live in EMRIWaveforms.AAK (§5.5/§7.7.1 purity
# hoist: Kernels carries NO SpecialFunctions machinery); Kernels owns the pure harmonic-sum kernel
# and receives the Bessel coefficients through the injected `bessel_coeffs` provider. The provider
# CANNOT be a precomputed array: the coefficients depend on the spline-resampled eccentricity at the
# dense output nodes, which only exists inside this kernel — so the coefficient BUILDER is injected
# (host-side dependency inversion), keeping the per-work-item body free of any Bessel dependency.

export aak_summation

"""
    _aak_rot_coeff(iota, theta_S, phi_S, theta_K, phi_K, alpha) -> NTuple{4}

Rotation coefficients `(rot0, rot1, -rot1, rot0)` mapping the kludge wave-frame polarisations to the
detector frame (`gpuAAK.cu::d_RotCoeff`), built from the sky direction `n`, orbital angular momentum `L`,
spin `S`, and the cross products `n×L`, `n×S`.
"""
function _aak_rot_coeff(iota, theta_S, phi_S, theta_K, phi_K, alpha)
    n = (sin(theta_S) * cos(phi_S), sin(theta_S) * sin(phi_S), cos(theta_S))
    S = (sin(theta_K) * cos(phi_K), sin(theta_K) * sin(phi_K), cos(theta_K))
    L = (cos(iota) * sin(theta_K) * cos(phi_K) + sin(iota) * (sin(alpha) * sin(phi_K) - cos(alpha) * cos(theta_K) * cos(phi_K)),
         cos(iota) * sin(theta_K) * sin(phi_K) - sin(iota) * (sin(alpha) * cos(phi_K) + cos(alpha) * cos(theta_K) * sin(phi_K)),
         cos(iota) * cos(theta_K) + sin(iota) * cos(alpha) * sin(theta_K))
    cross(u, v) = (u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1])
    dotp(u, v) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
    vnorm(u) = sqrt(dotp(u, u))
    nxL = cross(n, L); nxS = cross(n, S)
    norm = vnorm(nxL) * vnorm(nxS)
    norm < 1e-6 && (norm = 1e-6)
    cosrot = dotp(nxL, nxS) / norm
    sinrot = (dotp(L, nxS) - dotp(S, nxL)) / norm
    rot0 = 2.0 * cosrot * cosrot - 1.0
    rot1 = cosrot * sinrot
    return (rot0, rot1, -rot1, rot0)
end

"""
    _aak_point(bessel_coeffs, e, Phi, gim, alp, nu, gimdot, OmegaPhi, lam, t;
               M_phys, S_phys, mu, qS, phiS, qK, phiK, dist, nmodes, mich) -> ComplexF64

Single-timestep AAK waveform value `h₊ - i h₍×₎` (or detector `hᴵ - i hᴵᴵ` when `mich`), summing harmonics
`n = 1 … nmodes` from the analytic-kludge basis (`gpuAAK.cu::make_waveform` inner body). Per-harmonic
Peters–Mathews amplitudes `a,b,c` use the five Bessel coefficients supplied by the injected
`bessel_coeffs(n, e) -> (J_{n-2}, J_{n-1}, J_n, J_{n+1}, J_{n+2})` provider (canonical implementation:
`EMRIWaveforms.AAK.ak_bessel_coeffs`; Kernels itself carries no Bessel machinery — §5.5 purity).
"""
function _aak_point(bessel_coeffs, e, Phi, gim, alp, nu, gimdot, OmegaPhi, lam, t;
                    M_phys, S_phys, mu, qS, phiS, qK, phiK, dist, nmodes, mich)
    fill_val = 1e-6
    qS < fill_val && (qS = fill_val); qK < fill_val && (qK = fill_val)
    qS > π - fill_val && (qS = π - fill_val); qK > π - fill_val && (qK = π - fill_val)
    cosqS = cos(qS); sinqS = sin(qS); cosqK = cos(qK); sinqK = sin(qK)
    cosphiK = cos(phiK); sinphiK = sin(phiK)
    halfsqrt3 = sqrt(3.0) / 2.0
    zeta = (mu * MTSUN_SI) / dist / GPCINSEC
    lam > π - fill_val && (lam = π - fill_val); lam < fill_val && (lam = fill_val)
    coslam = cos(lam); sinlam = sin(lam); cosalp = cos(alp); sinalp = sin(alp)
    cosqL = cosqK * coslam + sinqK * sinlam * cosalp
    sinqL = sqrt(1.0 - cosqL * cosqL)
    phiLup = sinqK * sinphiK * coslam - cosphiK * sinlam * sinalp - cosqK * sinphiK * sinlam * cosalp
    phiLdown = sinqK * cosphiK * coslam + sinphiK * sinlam * sinalp - cosqK * cosphiK * sinlam * cosalp
    phiL = atan(phiLup, phiLdown)
    Ldotn = cosqL * cosqS + sinqL * sinqS * cos(phiL - phiS)
    Ldotn2 = Ldotn * Ldotn
    Sdotn = cosqK * cosqS + sinqK * sinqS * cos(phiK - phiS)
    beta = if S_phys == 0.0
        0.0
    else
        betaup = -Sdotn + coslam * Ldotn
        betadown = sinqS * sin(phiK - phiS) * sinlam * cosalp + (cosqK * Sdotn - cosqS) / sinqK * sinlam * sinalp
        atan(betaup, betadown)
    end
    gam = 2.0 * (gim + beta)
    cos2gam = cos(gam); sin2gam = sin(gam)

    local cosorbphs = 0.0
    if mich
        orbphs = 2.0 * π * t / YRSID_SI
        cosorbphs = cos(orbphs - phiS); sinorbphs = sin(orbphs - phiS)
        cosq = 0.5 * cosqS - halfsqrt3 * sinqS * cosorbphs
        phiw = orbphs + atan(halfsqrt3 * cosqS + 0.5 * sinqS * cosorbphs, sinqS * sinorbphs)
        psiup = 0.5 * cosqK - halfsqrt3 * sinqK * cos(orbphs - phiK) - cosq * (cosqK * cosqS + sinqK * sinqS * cos(phiK - phiS))
        psidown = 0.5 * sinqK * sinqS * sin(phiK - phiS) -
                  halfsqrt3 * cos(orbphs) * (cosqK * sinqS * sin(phiS) - cosqS * sinqK * sin(phiK)) -
                  halfsqrt3 * sin(orbphs) * (cosqS * sinqK * cos(phiK) - cosqK * sinqS * cos(phiS))
        psi = atan(psiup, psidown)
        cosq1 = 0.5 * (1.0 + cosq * cosq)
        cos2phi = cos(2.0 * phiw); sin2phi = sin(2.0 * phiw)
        cos2psi = cos(2.0 * psi); sin2psi = sin(2.0 * psi)
        FplusI = cosq1 * cos2phi * cos2psi - cosq * sin2phi * sin2psi
        FcrosI = cosq1 * cos2phi * sin2psi + cosq * sin2phi * cos2psi
        FplusII = cosq1 * sin2phi * cos2psi + cosq * cos2phi * sin2psi
        FcrosII = cosq1 * sin2phi * sin2psi - cosq * cos2phi * cos2psi
    else
        FplusI = 1.0; FcrosI = 0.0; FplusII = 0.0; FcrosII = 1.0
    end

    Amp = (abs(OmegaPhi) * M_phys * MTSUN_SI)^(2.0 / 3.0) * zeta
    rot = _aak_rot_coeff(lam, qS, phiS, qK, phiK, alp)

    hItemp = 0.0; hIItemp = 0.0
    for n in 1:nmodes
        nPhi = if mich
            fn = n * nu + gimdot / π
            Doppler = 2.0 * π * fn * AUsec * sinqS * cosorbphs
            n * Phi + Doppler
        else
            n * Phi
        end
        J0, J1, J2, J3, J4 = bessel_coeffs(n, e)
        a = -n * Amp * (J0 - 2.0 * e * J1 + 2.0 / n * J2 + 2.0 * e * J3 - J4) * cos(nPhi)
        b = -n * Amp * sqrt(1 - e * e) * (J0 - 2.0 * J2 + J4) * sin(nPhi)
        c = 2.0 * Amp * J2 * cos(nPhi)
        Aplus = -(1.0 + Ldotn2) * (a * cos2gam - b * sin2gam) + c * (1 - Ldotn2)
        Acros = 2.0 * Ldotn * (b * cos2gam + a * sin2gam)
        Aplusold = Aplus; Acrosold = Acros
        Aplus = Aplusold * rot[1] + Acrosold * rot[2]
        Acros = Aplusold * rot[3] + Acrosold * rot[4]
        if mich
            hItemp += halfsqrt3 * (FplusI * Aplus + FcrosI * Acros)
            hIItemp += halfsqrt3 * (FplusII * Aplus + FcrosII * Acros)
        else
            hItemp += FplusI * Aplus + FcrosI * Acros
            hIItemp += FplusII * Aplus + FcrosII * Acros
        end
    end
    return complex(hItemp, -hIItemp)
end

"""
    _aak_lam_fixup!(lam; fill_val=1e-6) -> lam

In-place polar-angle fixup: nudge `λ = ι` away from `0, π/2, π` by `fill_val` so the analytic-kludge
geometric projection stays regular (`aakwave.py` AAKSummation.sum `lam` block).
"""
function _aak_lam_fixup!(lam::AbstractVector; fill_val::Real = 1e-6)
    h = π / 2
    @inbounds for i in eachindex(lam)
        l = lam[i]
        if l > π - fill_val
            lam[i] = π - fill_val
        elseif l < fill_val
            lam[i] = fill_val
        elseif abs(l - h) < fill_val
            lam[i] = l > h ? h + fill_val : h - fill_val
        end
    end
    return lam
end

"""
    aak_summation(t, e, Phivec, gimvec, alpvec, nuvec, gimdotvec, OmegaPhi, lam,
                  M, a, mu, qS, phiS, qK, phiK, dist, nmodes;
                  bessel_coeffs, mich=false, dt=10.0, out_len=nothing, threaded=false)
        -> Vector{ComplexF64}

Augmented-analytic-kludge waveform (FEW `AAKSummation.sum` CPU path). The sparse kludge basis
`(e, Φ, γ, α, ν, γ̇, Ω_φ, λ)` is fit with one not-a-knot cubic spline (`NumericsCore.CubicSplineInterpolant`),
resampled onto the dense grid `tᵢ = dt·(i-1)` (`out_len` defaults to `⌊t[end]/dt⌋+1`), and each node
reduced by [`_aak_point`]. `mich=false` returns `h₊ - i h₍×₎`; `mich=true` the LISA `hᴵ - i hᴵᴵ`.
`M`,`a` are MBH mass (M_⊙) and spin; `mu` the secondary mass (M_⊙); `dist` the distance (Gpc).

`bessel_coeffs` is the REQUIRED Peters–Mathews coefficient provider,
`bessel_coeffs(n::Int, e::Float64) -> (J_{n-2}, J_{n-1}, J_n, J_{n+1}, J_{n+2})` at argument `n·e`
with `J_{-1} = -J_1` at `n = 1` (canonical implementation: `EMRIWaveforms.AAK.ak_bessel_coeffs`).
Injected because Kernels carries no SpecialFunctions dep (§5.5 purity: deps ==
{NumericsCore, KernelAbstractions}) and the coefficients depend on the spline-resampled eccentricity
at the dense nodes, so they cannot be precomputed as a bounded array by the caller.

One work-item per dense node + serial harmonic reduction ⇒ the serial and bounded-threaded paths are
bit-identical (`threaded=true`; design §8 bit-equality gate).
"""
function aak_summation(t::AbstractVector, e::AbstractVector, Phivec::AbstractVector,
                       gimvec::AbstractVector, alpvec::AbstractVector, nuvec::AbstractVector,
                       gimdotvec::AbstractVector, OmegaPhi::AbstractVector, lam::AbstractVector,
                       M::Real, a::Real, mu::Real, qS::Real, phiS::Real, qK::Real, phiK::Real,
                       dist::Real, nmodes::Integer;
                       bessel_coeffs, mich::Bool = false, dt::Real = 10.0,
                       out_len::Union{Integer,Nothing} = nothing, threaded::Bool = false)
    nt = length(t)
    olen = out_len === nothing ? floor(Int, t[end] / dt) + 1 : Int(out_len)

    # FEW y_all stacking (aakwave.py): e, Phivec, gimvec, alpvec, nuvec, gimdotvec, OmegaPhi, lam.
    Y = Matrix{Float64}(undef, 8, nt)
    @views Y[1, :] .= e;        @views Y[2, :] .= Phivec
    @views Y[3, :] .= gimvec;   @views Y[4, :] .= alpvec
    @views Y[5, :] .= nuvec;    @views Y[6, :] .= gimdotvec
    @views Y[7, :] .= OmegaPhi; @views Y[8, :] .= lam
    spline = CubicSplineInterpolant(collect(Float64, t), Y)

    h = Vector{ComplexF64}(undef, olen)
    with_kernel_threads(olen; threaded = threaded) do i
        td = dt * (i - 1)
        @inbounds h[i] = _aak_point(bessel_coeffs,
                                    spline(td, 1), spline(td, 2), spline(td, 3), spline(td, 4),
                                    spline(td, 5), spline(td, 6), spline(td, 7), spline(td, 8), td;
                                    M_phys = M, S_phys = a, mu = mu, qS = qS, phiS = phiS,
                                    qK = qK, phiK = phiK, dist = dist, nmodes = nmodes, mich = mich)
    end
    return h
end
