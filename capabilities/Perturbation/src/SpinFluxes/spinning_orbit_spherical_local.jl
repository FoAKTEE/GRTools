# LOCAL spherical-orbit linear-in-spin correction engine + assembly + r/x FD derivatives —
# faithful 1:1 of SpinningOrbit/src/{numeric_correction_spherical,orbit_correction_spherical,
# orbit_derivatives_spherical}.jl (SpinningTrajectory.wl `CorrectionSpherical` L945-1084,
# `KerrSpinOrbitCorrection2Spherical` L1087-1221, `KerrSpinOrbitCorrectionSpherical` L1224-1228),
# rebuilt on Spacetime's geodesic layer. The θ-only (1D-DFT) analog of the eccentric engine in
# spinning_orbit_correction_local.jl: r≡p constant, drdλ=0. Reuses that file's _oc_h1/h2 and
# _oc_Sg**/_oc_dSg** helpers, _d2z_dλ2, and the _chith_hat Darwin map. See the FIREWALL NOTE in
# spinning_secondary_local.jl. Paper arXiv:2506.20726.

# The 22 spherical force functions (.wl L962-983); auto-translated, verbatim from the legacy port.
"`ℛt` t-equation source (.wl L962)."
_sph_Rt(r, z, dz, K, a) = -((2*a^2*dz*r*z*(3*r^2*(2*K + r^2) - 2*a^2*(3*K + 4*r^2)*z^2 + a^4*z^4) + 0*(3*K*r^4 - 9*a^2*r^2*(2*K + r^2)*z^2 + a^4*(3*K + 14*r^2)*z^4 - a^6*z^6))/(sqrt(K)*(r^2 + a^2*z^2)^4))
"`ℛφ` φ-equation source (.wl L963)."
_sph_Rphi(r, z, dz, K, a) = -((a*(-2*dz*r*(a^2 + r^2)*z*(3*r^2*(2*K + r^2) - 2*a^2*(3*K + 4*r^2)*z^2 + a^4*z^4) + 0*(-1 + z^2)*(3*K*r^4 - 9*a^2*r^2*(2*K + r^2)*z^2 + a^4*(3*K + 14*r^2)*z^4 - a^6*z^6)))/(sqrt(K)*(r^2 + a^2*z^2)^4))
"`𝒢θ` r-equation (.wl L964)."
_sph_Gth(r, z, dz, chith, incl, a) = -((2*dz*r*(a^2 + (-2 + r)*r)*sin(incl)*sin(chith))/((-1 + z^2)*(r^2 + a^2*z^2)))
"`ℋ𝔯` r-equation radial-shift metric force (.wl L965)."
_sph_HR(r, z, dz, En, Lz, a) = (1/2)*((1/(a^2 + (-2 + r)*r)^2)*(2*a^2*(-2*a^2*En^2 + 4*a*En*Lz + (-2 + a^2)*Lz^2) + 2*r*(2*a^2*(-(a^2*En^2) - 4*a*En*Lz + Lz^2) + r*(a*(16*a*En^2 + 4*En*Lz - a*Lz^2) - 2*En^2*r*(8 + 2*a^2 - 5*r + r^2)))) + (Lz^2*(a^4 - 4*a^2*r + 2*a^2*r^2 + r^4))/((a^2 + r^2)^2*(-1 + z)) - (Lz^2*(a^4 - 4*a^2*r + 2*a^2*r^2 + r^4))/((a^2 + r^2)^2*(1 + z)) - (8*r^3*(a*(a*En - Lz) + En*r^2)^2)/((a^2 + r^2)*(r^2 + a^2*z^2)^2) + (8*r*(a*(a*En - Lz) + En*r^2)*(a^3*(a*En - Lz) + 3*a^2*En*r^2 + 2*En*r^4))/((a^2 + r^2)^2*(r^2 + a^2*z^2))) + 0^2*((2 - a^2 + (-2 + r)*r)/(a^2 + (-2 + r)*r)^2 + (2*r^2)/(r^2 + a^2*z^2)^2 - 1/(r^2 + a^2*z^2)) + (dz^2*(-(a^2*r^2) + r^4 + a^2*(a^2 + r*(-4 + 3*r))*z^2))/((-1 + z^2)*(r^2 + a^2*z^2)^2)
"`ℋθ` r-equation polar metric force (.wl L966)."
_sph_Hth(r, z, dz, chith, incl, En, Lz, a) = (1/((-1 + z^2)^2*(r^2 + a^2*z^2)^2))*2*r*z*(Lz^2*r^4 + 2*a^2*Lz^2*r^2*z^2 + a^4*Lz^2*z^4 + 2*r^3*(-(a*En) - Lz + a*En*z^2)*(-(a*En) + Lz + a*En*z^2) + 2*a^2*r*(-(a*En) + Lz + a*En*z^2)*(-(a*En) + Lz + (a*En - 2*Lz)*z^2))*sin(incl)*sin(chith) - (2*a^2*0^2*r*z*sin(incl)*sin(chith))/(r^2 + a^2*z^2)^2 - (2*dz^2*r*(a^2 + (-2 + r)*r)*z*(r^2*cot(incl)^2 + a^2*(1 - 2*z^2 + z^4*csc(incl)^2))*csc(chith)*sin(incl))/((-1 + z^2)^2*(r^2 + a^2*z^2)^2)
"`ℐ1θ` r-equation polar inhomogeneity (.wl L967)."
_sph_I1th(r, z, dz, Uth, a) = (2*dz^2*r*(a^2 + (-2 + r)*r)/Uth)/((-1 + z^2)*(r^2 + a^2*z^2))
"`ℐ2` couples δudst (.wl L968)."
_sph_I2(r, z, En, Lz, a) = (1/(a^2 + (-2 + r)*r))*2*(a^3*(a*En - Lz) + r^2*(a*(2*a*En + Lz) + En*(-4 + r)*r)) - (4*r^2*(a*(a*En - Lz) + En*r^2))/(r^2 + a^2*z^2)
"`ℐ3` couples δudsφ (.wl L969)."
_sph_I3(r, z, En, Lz, a) = (2*(Lz*(-2 + r)^2*r^3 - a*En*r^3*(-4 + 3*r)*(-1 + z^2) + a^5*En*z^2*(-1 + z^2) + a^4*Lz*z^2*(1 + (-1 + r)*z^2) + a^2*Lz*r^2*(-1 + (-3 + 2*r)*z^2) - a^3*En*r^2*(-1 + z^4)))/((a^2 + (-2 + r)*r)*(-1 + z^2)*(r^2 + a^2*z^2))
"`𝒥` r-equation source (.wl L970): geo + ℐ2·δudst + ℐ3·δudsφ."
_sph_J(r, z, dut, duphi, En, Lz, K, a) = ((a^2*En - a*Lz + En*r^2)*(-3*K*r^4 + 9*a^2*r^2*(2*K + r^2)*z^2 - a^4*(3*K + 14*r^2)*z^4 + a^6*z^6))/(sqrt(K)*(r^2 + a^2*z^2)^3) + _sph_I2(r, z, En, Lz, a)*dut + _sph_I3(r, z, En, Lz, a)*duphi
"`𝒦θ` norm-equation (.wl L971)."
_sph_Kth(r, z, dz, chith, incl, a) = (2*dz*sin(incl)*sin(chith))/((-1 + z^2)*(r^2 + a^2*z^2))
"`ℳθ` norm-equation polar metric force (.wl L972)."
_sph_Mth(r, z, dz, chith, incl, En, Lz, a) = (2*a^2*0^2*z*sin(chith)*sin(incl))/((a^2 + (-2 + r)*r)*(r^2 + a^2*z^2)^2) - (1/((a^2 + (-2 + r)*r)*(-1 + z^2)^2*(r^2 + a^2*z^2)^2))*2*z*(Lz^2*r^4 + 2*a^2*Lz^2*r^2*z^2 + a^4*Lz^2*z^4 + 2*r^3*(-(a*En) - Lz + a*En*z^2)*(-(a*En) + Lz + a*En*z^2) + 2*a^2*r*(-(a*En) + Lz + a*En*z^2)*(-(a*En) + Lz + (a*En - 2*Lz)*z^2))*sin(chith)*sin(incl) + (2*dz^2*z*(r^2*cot(incl)^2 + a^2*(1 - 2*z^2 + z^4*csc(incl)^2))*sin(incl))/((-1 + z^2)^2*(r^2 + a^2*z^2)^2*sin(chith))
"`𝒩1θ` norm-equation polar inhomogeneity (.wl L973)."
_sph_N1th(r, z, dz, Uth, a) = -((2*dz^2/Uth)/((-1 + z^2)*(r^2 + a^2*z^2)))
"`𝒩2` couples δudst (.wl L974)."
_sph_N2(r, z, En, Lz, a) = 2*En + (4*r*(a*(a*En - Lz) + En*r^2))/((a^2 + (-2 + r)*r)*(r^2 + a^2*z^2))
"`𝒩3` couples δudsφ (.wl L975)."
_sph_N3(r, z, En, Lz, a) = (-2*r*(2*a*En - 2*Lz + Lz*r) - 2*a*(a*Lz - 2*En*r)*z^2)/((a^2 + (-2 + r)*r)*(-1 + z^2)*(r^2 + a^2*z^2))
"`𝒫` norm-equation source (.wl L976): 𝒩2·δudst + 𝒩3·δudsφ."
_sph_P(r, z, dut, duphi, En, Lz, a) = _sph_N2(r, z, En, Lz, a)*dut + _sph_N3(r, z, En, Lz, a)*duphi
"`𝒬θ` θ-equation (.wl L977)."
_sph_Qth(z, chith, incl) = (sin(incl)*sin(chith))/sqrt(1 - z^2)
"`𝒮θ` θ-equation (.wl L978)."
_sph_Sth(r, z, dz, chith, incl, a) = -((2*z*sqrt(1 - z^2)*(r^2*cot(incl)^2 + a^2*(1 - 2*z^2 + z^4*csc(incl)^2))*dz*sin(incl))/((-1 + z^2)^2*(r^2 + a^2*z^2)*sin(chith)))
"`𝒯θ` θ-equation (.wl L979)."
_sph_Tth(r, z, dz, d2z, chith, incl, En, Lz, a) = (1/((a^2 + (-2 + r)*r)*(1 - z^2)^(5/2)*(r^2 + a^2*z^2)^2))*(3*a^6*Lz^2*z^6 + Lz^2*r^6*(1 + 2*z^2) + a^4*Lz^2*r^2*z^4*(7 + 2*z^2) + a^2*Lz^2*r^4*z^2*(5 + 4*z^2) + 2*a^4*r*z^2*(-(-(a*En) + Lz)^2 + 2*(-(a*En) + Lz)^2*z^2 - (a^2*En^2 - 2*a*En*Lz + 4*Lz^2)*z^4) + 2*a^2*r^3*((-(a*En) + Lz)^2 - (5*a*En - 3*Lz)*(a*En - Lz)*z^2 + (7*a^2*En^2 - 10*a*En*Lz - 4*Lz^2)*z^4 + a*En*(-3*a*En + 4*Lz)*z^6) - 2*r^5*(-(a^2*En^2) + Lz^2 + 2*(2*a^2*En^2 + Lz^2)*z^2 + a^2*En^2*z^4*(-5 + 2*z^2)))*sin(incl)*sin(chith) + (a^2*0^2*(a^2*z^2 + r^2*(-1 + 2*z^2))*sin(incl)*sin(chith))/((a^2 + (-2 + r)*r)*sqrt(1 - z^2)*(r^2 + a^2*z^2)^2) - (1/(2*(1 - z^2)^(5/2)*(r^2 + a^2*z^2)^2))*dz^2*csc(incl)*csc(chith)*(2*r^4*(1 + 2*z^2)*cos(incl)^2 + a^2*r^2*(-4*z^6 + 2*z^4*(7 + 2*cos(2*incl)) + z^2*(-5 + 3*cos(2*incl)) + 2*sin(incl)^2) + 2*a^4*z^2*(z^4*(1 + 2*cos(2*incl)) - sin(incl)^2 + 2*z^2*sin(incl)^2)) + (2*z*cos(incl)*cot(incl)*csc(chith)*(dz^2*z + d2z*(-z^2 + sin(incl)^2)))/((1 - z^2)^(3/2)*(-1 + 2*z^2 + cos(2*incl)))
"`𝒰1θ` θ-equation polar inhomogeneity (.wl L980)."
_sph_U1th(r, z, dz, d2z, Uth, a) = -((-2*d2z*(-1 + z^2)*(r^2 + a^2*z^2) + 2*dz^2*z*(r^2 + a^2*(-1 + 2*z^2)))/((1 - z^2)^(3/2)*(r^2 + a^2*z^2)*Uth))
"`𝒰2` couples δudst (.wl L981)."
_sph_U2(r, z, En, Lz, a) = (4*a^2*r*(a^2*En - a*Lz + En*r^2)*z*sqrt(1 - z^2))/((a^2 + (-2 + r)*r)*(r^2 + a^2*z^2))
"`𝒰3` couples δudsφ (.wl L982)."
_sph_U3(r, z, En, Lz, a) = -((2*z*(Lz*(-2 + r)*r^3 + a^4*Lz*z^4 - 2*a^3*En*r*(-1 + z^2)^2 + 2*a^2*Lz*r*(1 + (-2 + r)*z^2)))/((a^2 + (-2 + r)*r)*(1 - z^2)^(3/2)*(r^2 + a^2*z^2)))
"`𝒱` θ-equation source (.wl L983): geo + 𝒰2·δudst + 𝒰3·δudsφ."
_sph_V(r, z, dut, duphi, En, Lz, K, a) = -((2*a*r*z*(3*r^2*(2*K + r^2) - 2*a^2*(3*K + 4*r^2)*z^2 + a^4*z^4)*(Lz + a*En*(-1 + z^2)))/(sqrt(K)*(r^2 + a^2*z^2)^3))/sqrt(1 - z^2) + _sph_U2(r, z, En, Lz, a)*dut + _sph_U3(r, z, En, Lz, a)*duphi

"""
    _correction_spherical(a, p, e, x, kmax) -> NamedTuple

Port of `CorrectionSpherical[orbitGeo,kmax]` (SpinningTrajectory.wl L945-1084): the θ-only numeric
linear-in-spin orbit-correction engine for a spherical Kerr orbit (r≡p, e=0). Samples the 22
spin-curvature forces on the `4kmax` midpoint polar-Mino grid, DFTs them, assembles the M-matrix
(r/θ/norm blocks + a constant-r radial-shift δ𝔯 block), solves M·v=c (least squares). Returns
`(; ΥθS, udt0S, udφ0S, error, kmax, δχθSIm, δ𝔯four, ℛtfourIm, ℛφfourIm)`; `kmax` must be even.
"""
function _correction_spherical(a, p, e, x, kmax::Integer)
    iseven(kmax) || throw(ArgumentError("CorrectionSpherical requires even kmax (got $kmax)"))
    En = geodesic_energy(a, p, e, x); Lz = geodesic_angular_momentum(a, p, e, x)
    Q = carter_constant(a, p, e, x); K = Q + (Lz - a * En)^2
    incl = acos(x)
    Υθ = kerr_frequencies(a, p, e, x; Time="Mino").Uth
    r = float(p)                         # r ≡ p (spherical orbit), drdλ ≡ 0
    stepsθ = 4 * kmax
    pol = geodesic_trajectory_polar(a, p, e, x)

    # --- θ-grid samples at polar Mino phase qz_j = 2π(j-1/2)/stepsθ (.wl L998) ---
    qz = [2pi * (j - 0.5) / stepsθ for j in 1:stepsθ]
    pg = [pol(qz[j]) for j in 1:stepsθ]
    zv = [g.z for g in pg]; dzv = [g.Uz for g in pg]
    d2zv = [_d2z_dλ2(zv[j], En, Lz, Q, a) for j in 1:stepsθ]
    χθv = [(begin
        base = acos(clamp(cos(pg[j].θ) / sin(incl), -1.0, 1.0))
        fl = 2pi * floor(qz[j] / (2pi))
        mod(qz[j], 2pi) < pi ? base + fl : 2pi + fl - base
    end) for j in 1:stepsθ]

    # --- 1D DFT in the polar phase ---
    Ejk = [cis(2pi * k * (j - 0.5) / stepsθ) for j in 1:stepsθ, k in -kmax:kmax]   # stepsθ×(2kmax+1)
    dft(list) = (list' * Ejk)[:] / stepsθ      # length-(2kmax+1) vector, indexed [k+kmax+1]
    idx(v, k) = v[k + kmax + 1]

    # ---------------- t,φ equations: ℛt/ℛφ DFT, then δudst/δudsφ (.wl L998-1001, L960-961) ----------------
    Rtfour = dft([_sph_Rt(r, zv[j], dzv[j], K, a) for j in 1:stepsθ])
    Rφfour = dft([_sph_Rphi(r, zv[j], dzv[j], K, a) for j in 1:stepsθ])
    δu(four, q) = real(sum(k == 0 ? 0.0 + 0im :
                           imag(idx(four, k)) / (-k * Υθ) * cis(-k * q)
                           for k in -kmax:2:kmax))
    δudst = [δu(Rtfour, qz[j]) for j in 1:stepsθ]
    δudsφ = [δu(Rφfour, qz[j]) for j in 1:stepsθ]

    # ---------------- r equation forces (.wl L1003-1016) ----------------
    Gθfour = dft([_sph_Gth(r, zv[j], dzv[j], χθv[j], incl, a)        for j in 1:stepsθ])
    Hθfour = dft([_sph_Hth(r, zv[j], dzv[j], χθv[j], incl, En, Lz, a) for j in 1:stepsθ])
    I1θfour = dft([_sph_I1th(r, zv[j], dzv[j], Υθ, a)                for j in 1:stepsθ])
    I2four = dft([_sph_I2(r, zv[j], En, Lz, a)                       for j in 1:stepsθ])
    I3four = dft([_sph_I3(r, zv[j], En, Lz, a)                       for j in 1:stepsθ])
    Jfour = dft([_sph_J(r, zv[j], δudst[j], δudsφ[j], En, Lz, K, a)  for j in 1:stepsθ])

    # ---------------- θ equation forces (.wl L1018-1031) ----------------
    Qθfour = dft([_sph_Qth(zv[j], χθv[j], incl)                      for j in 1:stepsθ])
    Sθfour = dft([_sph_Sth(r, zv[j], dzv[j], χθv[j], incl, a)        for j in 1:stepsθ])
    Tθfour = dft([_sph_Tth(r, zv[j], dzv[j], d2zv[j], χθv[j], incl, En, Lz, a) for j in 1:stepsθ])
    U1θfour = dft([_sph_U1th(r, zv[j], dzv[j], d2zv[j], Υθ, a)       for j in 1:stepsθ])
    U2four = dft([_sph_U2(r, zv[j], En, Lz, a)                       for j in 1:stepsθ])
    U3four = dft([_sph_U3(r, zv[j], En, Lz, a)                       for j in 1:stepsθ])
    Vfour = dft([_sph_V(r, zv[j], δudst[j], δudsφ[j], En, Lz, K, a)  for j in 1:stepsθ])

    # ---------------- norm equation forces (.wl L1033-1044) ----------------
    Kθfour = dft([_sph_Kth(r, zv[j], dzv[j], χθv[j], incl, a)        for j in 1:stepsθ])
    Mθfour = dft([_sph_Mth(r, zv[j], dzv[j], χθv[j], incl, En, Lz, a) for j in 1:stepsθ])
    N1θfour = dft([_sph_N1th(r, zv[j], dzv[j], Υθ, a)                for j in 1:stepsθ])
    N2four = dft([_sph_N2(r, zv[j], En, Lz, a)                       for j in 1:stepsθ])
    N3four = dft([_sph_N3(r, zv[j], En, Lz, a)                       for j in 1:stepsθ])
    Pfour = dft([_sph_P(r, zv[j], δudst[j], δudsφ[j], En, Lz, a)     for j in 1:stepsθ])

    # ============================ submatrix assembly (.wl L1046-1065) ============================
    keven = collect(-kmax:2:kmax)            # kmax+1 rows
    kodd  = collect(-kmax+1:2:kmax-1)        # kmax-1 rows
    kk    = vcat(-kmax:2:-2, 2:2:kmax)       # nonzero even kk modes (kmax columns)
    HRconst = _sph_HR(r, 0.0, sqrt(Q), En, Lz, a)   # ℋ𝔯[p,0,0,√Qhat] — z=0, dzdλ=√Qhat (.wl L1047)
    band(k, c) = -kmax <= k - c <= kmax

    M1θ = [band(k, c) ? -im*c*Υθ*idx(Gθfour, k - c) + idx(Hθfour, k - c) : 0.0 + 0im for k in keven, c in kk]
    M1𝔯 = [band(k, c) ? (k == c ? 1 : 0)*((-im*c*Υθ)^2 + HRconst) : 0.0 + 0im for k in keven, c in kk]
    M1Υθ = idx.(Ref(I1θfour), keven); M1t = idx.(Ref(I2four), keven); M1φ = idx.(Ref(I3four), keven)
    M1 = hcat(-imag(M1θ), real(M1𝔯), real(M1Υθ), real(M1t), real(M1φ))

    M2θ = [band(k, c) ? -c^2*Υθ^2*idx(Qθfour, k - c) - im*c*Υθ*idx(Sθfour, k - c) + idx(Tθfour, k - c) : 0.0 + 0im for k in kodd, c in kk]
    M2𝔯 = zeros(length(kodd), length(kk))
    M2Υθ = idx.(Ref(U1θfour), kodd); M2t = idx.(Ref(U2four), kodd); M2φ = idx.(Ref(U3four), kodd)
    M2 = hcat(-imag(M2θ), M2𝔯, real(M2Υθ), real(M2t), real(M2φ))

    M3θ = [band(k, c) ? -im*c*Υθ*idx(Kθfour, k - c) + idx(Mθfour, k - c) : 0.0 + 0im for k in keven, c in kk]
    M3𝔯 = zeros(length(keven), length(kk))
    M3Υθ = idx.(Ref(N1θfour), keven); M3t = idx.(Ref(N2four), keven); M3φ = idx.(Ref(N3four), keven)
    M3 = hcat(-imag(M3θ), M3𝔯, real(M3Υθ), real(M3t), real(M3φ))

    M = vcat(M1, M2, M3)

    # ---------------- source vector c (.wl L1065) ----------------
    cJ = real(idx.(Ref(Jfour), keven)); cV = real(idx.(Ref(Vfour), kodd)); cP = real(idx.(Ref(Pfour), keven))
    cvec = -vcat(cJ, cV, cP)

    # ---------------- least squares (Mathematica LeastSquares = min-norm) ----------------
    v = pinv(M) * cvec
    δχθSfourIm = v[1:kmax]
    δ𝔯four = v[kmax+1 : 2*kmax]
    ΥθS = real(v[end-2]); udt0S = real(v[end-1]); udφ0S = real(v[end])
    err = norm(M * v - cvec)

    return (; ΥθS, udt0S, udφ0S, error = err, kmax,
            δχθSIm = δχθSfourIm, δ𝔯four,
            ℛtfourIm = imag.(Rtfour), ℛφfourIm = imag.(Rφfour))
end

"""
    _kerr_spin_orbit_correction2_spherical(a, p, e, x, correction) -> NamedTuple

Port of `KerrSpinOrbitCorrection2Spherical[correction]` (SpinningTrajectory.wl L1087-1221):
assembles the linear-in-spin spherical orbit-correction object from the spherical engine output —
geodesic constants/Mino frequencies, the δχθS sin-series and δ𝔯 cos-series, the ES/LS zero-modes
(θ-only midpoint average), the UtS/UφS-list → 1D θ-DFT → ΔtS/ΔφS deltas, and the frequency
corrections {ΓS,ΥθS,ΥφS}.
"""
function _kerr_spin_orbit_correction2_spherical(a, p, e, x, correction)
    ScriptI = acos(x)
    kmax = correction.kmax
    stepsθ = 4 * kmax

    # --- geodesic constants & Mino frequencies (.wl L1152-1156) ---
    Ehat = geodesic_energy(a, p, e, x)
    Lzhat = geodesic_angular_momentum(a, p, e, x)
    Qhat = carter_constant(a, p, e, x)
    Khat = Qhat + (a * Ehat - Lzhat)^2
    fr = kerr_frequencies(a, p, e, x; Time="Mino")
    Υθhat = fr.Uth; Υφhat = fr.Uph; Γhat = fr.Ut

    pol = geodesic_trajectory_polar(a, p, e, x)
    zhat(wz) = cos(pol(wz).θ)

    # --- correction-engine outputs ---
    δχθSfourIm = correction.δχθSIm           # length kmax
    δ𝔯four = correction.δ𝔯four               # length kmax
    ΥθS = correction.ΥθS
    udt0S = correction.udt0S; udφ0S = correction.udφ0S

    # --- the χθhat Darwin-angle map (.wl L1095) ---
    χθhat(wz) = _chith_hat(a, p, e, x, wz / Υθhat)
    function dχθhat(wz)
        h = 1e-6
        (χθhat(wz + h) - χθhat(wz - h)) / (2h)
    end

    # --- δχθS sin-series and δ𝔯 cos-series (.wl L1096-1099) ---
    δχθS(wz) = 2 * sum(δχθSfourIm[k + kmax ÷ 2] * sin(2k * wz) for k in 1:(kmax ÷ 2))
    dδχθS(wz) = 2 * sum(2k * Υθhat * δχθSfourIm[k + kmax ÷ 2] * cos(2k * wz) for k in 1:(kmax ÷ 2))
    δ𝔯(wz) = 2 * sum(δ𝔯four[k + kmax ÷ 2] * cos(2k * wz) for k in 1:(kmax ÷ 2))
    dδ𝔯(wz) = -2 * sum((2k * Υθhat) * δ𝔯four[k + kmax ÷ 2] * sin(2k * wz) for k in 1:(kmax ÷ 2))

    # --- geodesic polar velocity along the phase variable (.wl L1100); drdλ ≡ 0 (spherical) ---
    Δ(r) = r^2 - 2r + a^2
    Σ(r, z) = r^2 + a^2 * z^2
    dzdλ(wz) = -sign(pi - mod(wz, 2pi)) *
               sqrt(Qhat - zhat(wz)^2 * (a^2 * (1 - Ehat^2) * (1 - zhat(wz)^2) + Lzhat^2 + Qhat))

    # --- spin 4-vector components (.wl L1103-1106), r≡p, drdλ≡0 ---
    St(wz) = -1 / sqrt(Khat) * (-a * (p * dzdλ(wz)) / Σ(p, zhat(wz)))
    Sr(wz) = -1 / sqrt(Khat) * (a * zhat(wz) * ((p^2 + a^2) * Ehat - a * Lzhat) / Δ(p))
    Sz(wz) = -1 / sqrt(Khat) * (a * p * Ehat - p * Lzhat / (1 - zhat(wz)^2))
    Sφ(wz) = -1 / sqrt(Khat) * ((p * (p^2 + a^2) * dzdλ(wz)) / Σ(p, zhat(wz)))

    # --- position corrections rS/zS (.wl L1117-1118) ---
    rS(wz) = δ𝔯(wz)
    zS(wz) = -sin(ScriptI) * δχθS(wz) * sin(χθhat(wz))

    # --- ES, LS scalar zero-modes (θ-only midpoint averages, .wl L1173-1174) ---
    wzg = [(j - 0.5) * 2pi / stepsθ for j in 1:stepsθ]
    ES = -udt0S + sum(_oc_h1(p, zhat(wzg[j]), Khat, a) for j in 1:stepsθ) / stepsθ
    LS = udφ0S - sum(_oc_h2(p, zhat(wzg[j]), Khat, Ehat, Lzhat, a) for j in 1:stepsθ) / stepsθ

    # --- four-velocity correction zero-modes (.wl L1175-1176): NO h3 in udφS (drdλ=0) ---
    udtS(wz) = -ES + _oc_h1(p, zhat(wz), Khat, a)
    udφS(wz) = LS + _oc_h2(p, zhat(wz), Khat, Ehat, Lzhat, a)

    # --- UtS/UφS lists → 1D θ-DFT → Im[ΔtS/ΔφS]four (.wl L1142-1192) ---
    UtS(wz) = (_oc_dSgttdr(p, a) * rS(wz) + _oc_dSgttdz(zhat(wz), a) * zS(wz)) * (-Ehat) +
              _oc_Sgtt(p, zhat(wz), a) * udtS(wz) +
              _oc_dSgtphdr(p, a) * rS(wz) * Lzhat +
              _oc_Sgtph(p, a) * udφS(wz)
    UφS(wz) = (_oc_dSgphphdr(p, a) * rS(wz) + _oc_dSgphphdz(zhat(wz), a) * zS(wz)) * Lzhat +
              _oc_Sgphph(p, zhat(wz), a) * udφS(wz) +
              _oc_dSgtphdr(p, a) * rS(wz) * (-Ehat) +
              _oc_Sgtph(p, a) * udtS(wz)
    Uthat(wz) = _oc_Sgtt(p, zhat(wz), a) * (-Ehat) + _oc_Sgtph(p, a) * Lzhat
    Uφhat(wz) = _oc_Sgphph(p, zhat(wz), a) * Lzhat + _oc_Sgtph(p, a) * (-Ehat)

    UtSlist   = [UtS(wzg[j])   for j in 1:stepsθ]
    UφSlist   = [UφS(wzg[j])   for j in 1:stepsθ]
    Uthatlist = [Uthat(wzg[j]) for j in 1:stepsθ]
    Uφhatlist = [Uφhat(wzg[j]) for j in 1:stepsθ]

    Ejk = [cis(2pi * k * (j - 0.5) / stepsθ) for j in 1:stepsθ, k in -kmax:kmax]   # stepsθ×(2kmax+1)
    dft(list) = (list' * Ejk)[:] / stepsθ      # length-(2kmax+1) vector, indexed [k+kmax+1]
    UtSfour   = dft(UtSlist);   UφSfour   = dft(UφSlist)
    Uthatfour = dft(Uthatlist); Uφhatfour = dft(Uφhatlist)

    # ImΔtSfour[k]/ImΔφSfour[k] (.wl L1191-1192), indexed [k+kmax+1]; k=0 mode → 0
    ImΔtSfour = zeros(Float64, 2 * kmax + 1)
    ImΔφSfour = zeros(Float64, 2 * kmax + 1)
    for k in -kmax:kmax
        k == 0 && continue
        ImΔtSfour[k + kmax + 1] = real(UtSfour[k + kmax + 1]) / (k * Υθhat) -
                                  real(Uthatfour[k + kmax + 1]) / (k * Υθhat)^2 * (k * ΥθS)
        ImΔφSfour[k + kmax + 1] = real(UφSfour[k + kmax + 1]) / (k * Υθhat) -
                                  real(Uφhatfour[k + kmax + 1]) / (k * Υθhat)^2 * (k * ΥθS)
    end
    ΓS  = real(UtSfour[kmax + 1])   # k=0 mode
    ΥφS = real(UφSfour[kmax + 1])

    # --- ΔtS/ΔφS closures (.wl L1121-1122,L1203): the sum is written TWICE in the source —
    #     replicated VERBATIM (a genuine .wl transcription quirk; the single-sum OrbitCorrection
    #     form below is the physical one, consumed by the flux modes). ---
    ΔtS(wz) = sum(2 * ImΔtSfour[k + kmax + 1] * sin(k * wz) for k in 1:kmax) +
              sum(2 * ImΔtSfour[k + kmax + 1] * sin(k * wz) for k in -kmax:kmax)
    ΔφS(wz) = sum(2 * ImΔφSfour[k + kmax + 1] * sin(k * wz) for k in 1:kmax) +
              sum(2 * ImΔφSfour[k + kmax + 1] * sin(k * wz) for k in -kmax:kmax)
    ocΔtS(wz) = 2 * sum(ImΔtSfour[k + kmax + 1] * sin(k * wz) for k in 1:kmax)
    ocΔφS(wz) = 2 * sum(ImΔφSfour[k + kmax + 1] * sin(k * wz) for k in 1:kmax)

    # --- OrbitCorrection: NamedTuple-returning closure (.wl L1123-1141) ---
    function OrbitCorrection(wz)
        δχθSp = δχθS(wz); χθhatp = χθhat(wz)
        return (ΔtS = ocΔtS(wz), rS = rS(wz), zS = zS(wz), ΔφS = ocΔφS(wz),
                UrS = dδ𝔯(wz),
                UzS = -sin(ScriptI) * (dχθhat(wz) * (cos(χθhatp) * δχθSp * Υθhat + sin(χθhatp) * ΥθS) +
                       sin(χθhatp) * dδχθS(wz)))
    end

    FourVelocityCorrection = (
        udt = wz -> udtS(wz),
        ur = wz -> dδ𝔯(wz),
        uz = wz -> -sin(ScriptI) * (dχθhat(wz) * (cos(χθhat(wz)) * δχθS(wz) * Υθhat +
              sin(χθhat(wz)) * ΥθS) + sin(χθhat(wz)) * dδχθS(wz)),
        udφ = wz -> udφS(wz))

    return (; a = float(a), p = float(p), e = float(e), Inclination = x, ScriptI,
            Ehat, Lzhat, Khat,
            TrajectoryCorrection = (wz -> ΔtS(wz), wz -> rS(wz), wz -> zS(wz), wz -> ΔφS(wz)),
            TrajectoryGeo = pol,
            MinoFrequenciesGeo = (Γhat, Υθhat, Υφhat),
            BLFrequenciesGeo = (Υθhat / Γhat, Υφhat / Γhat),
            FourVelocityCorrection,
            Spin4vector = (wz -> St(wz), wz -> Sr(wz), wz -> Sz(wz), wz -> Sφ(wz)),
            TrajectoryDeltas = (Δtz = wz -> pol(wz).Δtz, Δφz = wz -> pol(wz).Δφz),
            MinoFrequenciesCorrection = (ΓS, ΥθS, ΥφS),
            BLFrequenciesCorrection = (ΥθS / Γhat - Υθhat / Γhat^2 * ΓS,
                                       ΥφS / Γhat - Υφhat / Γhat^2 * ΓS),
            ES, LS,
            OrbitCorrection,
            error = correction.error,
            kmax)
end

"""
    KerrSpinOrbitCorrectionSpherical(a, p, e, x, kmax) -> NamedTuple

Port of `KerrSpinOrbitCorrectionSpherical[a,p,e,x,kmax]` (SpinningTrajectory.wl L1224-1228): runs
the spherical correction engine then the `_kerr_spin_orbit_correction2_spherical` assembly.
`kmax` must be even; `e` must be 0 (spherical orbit). Unexported (internal kinematics).
"""
function KerrSpinOrbitCorrectionSpherical(a, p, e, x, kmax::Integer)
    iseven(kmax) || throw(ArgumentError("kmax must be even (got $kmax)"))
    iszero(e) || throw(ArgumentError("KerrSpinOrbitCorrectionSpherical requires a spherical orbit e=0 (got e=$e)"))
    correction = _correction_spherical(a, p, e, x, kmax)
    return _kerr_spin_orbit_correction2_spherical(a, p, e, x, correction)
end

"""
    _orbit_derivatives_spherical(a, p, e, x, kmax; h=1e-5) -> NamedTuple

r/x (≡ ∂p, ∂x) finite-difference derivatives of `KerrSpinOrbitCorrectionSpherical(a,p,e,x,kmax)`,
matching the `orbitDerivatives` Association the .wl `TeukolskySpinModeSphericalCorrection` reads.
Keys: `dEndr,dLzdr,dEndx,dLzdx` (geodesic constants), `dBLFrequenciesdr,dBLFrequenciesdx` (2-vec
geo BL freqs Ωθ,Ωφ), `dMinoFrequenciesdr,dMinoFrequenciesdx` (3-vec geo Mino Γ,Υθ,Υφ), and
`OrbitDerivatives` — a closure `wz -> (; dzdr,dzdx,dUzdr,dUzdx,dΔtdr,dΔφdr,dΔtdx,dΔφdx)` (central
FD of the OrbitCorrection closure over p,x). `e` must be 0. USE THE SAME `h` as the oracle.
"""
function _orbit_derivatives_spherical(a, p, e, x, kmax::Integer; h::Real=1e-5)
    iszero(e) || throw(ArgumentError("_orbit_derivatives_spherical requires a spherical orbit e=0 (got e=$e)"))

    En(pp, xx) = geodesic_energy(a, pp, e, xx)
    Lz(pp, xx) = geodesic_angular_momentum(a, pp, e, xx)
    dEndr = (En(p + h, x) - En(p - h, x)) / (2h)
    dLzdr = (Lz(p + h, x) - Lz(p - h, x)) / (2h)
    dEndx = (En(p, x + h) - En(p, x - h)) / (2h)
    dLzdx = (Lz(p, x + h) - Lz(p, x - h)) / (2h)

    minofreq(pp, xx) = (f = kerr_frequencies(a, pp, e, xx; Time="Mino"); (f.Ut, f.Uth, f.Uph))
    mf(pp, xx) = (mm = minofreq(pp, xx); (Γ = mm[1], BL = (mm[2] / mm[1], mm[3] / mm[1])))
    mfpp = mf(p + h, x); mfpm = mf(p - h, x); mfxp = mf(p, x + h); mfxm = mf(p, x - h)
    mfrp = minofreq(p + h, x); mfrm = minofreq(p - h, x)
    mfxpf = minofreq(p, x + h); mfxmf = minofreq(p, x - h)
    dBLFrequenciesdr = ((mfpp.BL[1] - mfpm.BL[1]) / (2h), (mfpp.BL[2] - mfpm.BL[2]) / (2h))
    dBLFrequenciesdx = ((mfxp.BL[1] - mfxm.BL[1]) / (2h), (mfxp.BL[2] - mfxm.BL[2]) / (2h))
    dMinoFrequenciesdr = ((mfrp[1] - mfrm[1]) / (2h), (mfrp[2] - mfrm[2]) / (2h), (mfrp[3] - mfrm[3]) / (2h))
    dMinoFrequenciesdx = ((mfxpf[1] - mfxmf[1]) / (2h), (mfxpf[2] - mfxmf[2]) / (2h), (mfxpf[3] - mfxmf[3]) / (2h))

    # OrbitCorrection closures at the four shifted points (built ONCE; closures over wz)
    OCpp = KerrSpinOrbitCorrectionSpherical(a, p + h, e, x, kmax).OrbitCorrection
    OCpm = KerrSpinOrbitCorrectionSpherical(a, p - h, e, x, kmax).OrbitCorrection
    OCxp = KerrSpinOrbitCorrectionSpherical(a, p, e, x + h, kmax).OrbitCorrection
    OCxm = KerrSpinOrbitCorrectionSpherical(a, p, e, x - h, kmax).OrbitCorrection

    function OrbitDerivatives(wz)
        cpp = OCpp(wz); cpm = OCpm(wz); cxp = OCxp(wz); cxm = OCxm(wz)
        return (dzdr  = (cpp.zS  - cpm.zS)  / (2h), dzdx  = (cxp.zS  - cxm.zS)  / (2h),
                dUzdr = (cpp.UzS - cpm.UzS) / (2h), dUzdx = (cxp.UzS - cxm.UzS) / (2h),
                dΔtdr = (cpp.ΔtS - cpm.ΔtS) / (2h), dΔφdr = (cpp.ΔφS - cpm.ΔφS) / (2h),
                dΔtdx = (cxp.ΔtS - cxm.ΔtS) / (2h), dΔφdx = (cxp.ΔφS - cxm.ΔφS) / (2h))
    end

    return (; dEndr, dLzdr, dEndx, dLzdx,
            dBLFrequenciesdr, dBLFrequenciesdx,
            dMinoFrequenciesdr, dMinoFrequenciesdx,
            OrbitDerivatives, h = float(h))
end
