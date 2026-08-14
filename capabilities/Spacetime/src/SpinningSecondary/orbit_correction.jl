# Orbit-correction-object assembly + the two public functions it unlocks — consolidates
# SpinningOrbit/src/orbit_correction.jl (port of SpinningTrajectory.wl `KerrSpinOrbitCorrection2`
# L748-906, `KerrSpinOrbitCorrection` L907-911, `KerrSpinOrbitNum` L914-938; renamed
# `_spin_orbit_correction2` / `spin_orbit_correction` / `spin_orbit_numeric` in the GrTool style).
# Given the numeric `Correction` engine output (here the verified `_numeric_correction` NamedTuple),
# assemble the linear-in-spin orbit correction: the geodesic constants/Mino-frequencies, the
# δχrS/δχθS/δ𝔯/δ𝔷 Fourier series, the four-velocity correction pieces, the UtS/UφS → DFT →
# ΔtS/ΔφS trajectory deltas, and the frequency corrections {ΓS,ΥrS,ΥθS,ΥφS}.
# R3 (load-bearing numerics — the headline ΥφS/ΓS come out of this DFT assembly).

# ---- four-velocity / metric helper functions (SpinningTrajectory.wl L753-784); h2 (L754) is the
#      page-long one — MMA→Julia auto-translated (gen_orbcorr_jl.wl, InputForm+head-subst). R3. ----
"`h1` (.wl L753)."
_oc_h1(r, z, K, a) = (r*(-3*a^2*r^2*z^2 + a^4*z^4 + K*(r^2 - 3*a^2*z^2)))/(sqrt(K)*(r^2 + a^2*z^2)^3)
"`h2` (.wl L754) — the large four-velocity-equation source. R3."
_oc_h2(r, z, K, En, Lz, a) = (-(En*Lz*r^6) + a^4*En*Lz*r^2*z^4 + a^2*En*Lz*r^4*(-2 + z^2) - a^6*En*Lz*z^4*(-2 + z^2) + a^7*En^2*z^4*(-1 + z^2) + a*r^3*(Lz^2*r + K*(-1 + z^2) + K*r*(-1 + 2*z^2) + r^3*(z^2 - En^2*(-1 + z^2))) + a^3*(-(En^2*r^4*(-1 + z^2)) + r*z^2*(2*r^3 + 2*K*r*z^2 - 3*K*(-1 + z^2) - 3*r^2*(-1 + z^2))) + a^5*z^4*(K - Lz^2 + r*(-1 + z^2 + r*(2 - z^2 + En^2*(-1 + z^2)))))/(sqrt(K)*(r^2 + a^2*z^2)^3)
"`h3` (.wl L755)."
_oc_h3(r, z, K, a) = (-2*a*r*z)/(sqrt(K)*(r^2 + a^2*z^2)^2)
"`Σgtt` (.wl L777)."
_oc_Sgtt(r, z, a) = -((a^2 + r^2)^2/(a^2 + (-2 + r)*r)) + a^2*(1 - z^2)
"`Σgφφ` (.wl L778)."
_oc_Sgphph(r, z, a) = -(a^2/(a^2 + (-2 + r)*r)) + (1 - z^2)^(-1)
"`Σgtφ` (.wl L779)."
_oc_Sgtph(r, a) = a*(1 - (a^2 + r^2)/(a^2 - 2*r + r^2))
"`dΣgtt/dr` (.wl L780)."
_oc_dSgttdr(r, a) = (2*(a^2 + r^2)*((-1 + r)*(a^2 + r^2) - 2*r*(a^2 - 2*r + r^2)))/(a^2 - 2*r + r^2)^2
"`dΣgtt/dz` (.wl L781)."
_oc_dSgttdz(z, a) = -2*a^2*z
"`dΣgtφ/dr` (.wl L782)."
_oc_dSgtphdr(r, a) = (2*a*((-1 + r)*(a^2 + r^2) - r*(a^2 - 2*r + r^2)))/(a^2 - 2*r + r^2)^2
"`dΣgφφ/dr` (.wl L783)."
_oc_dSgphphdr(r, a) = (2*a^2*(-1 + r))/(a^2 + (-2 + r)*r)^2
"`dΣgφφ/dz` (.wl L784)."
_oc_dSgphphdz(z, a) = (2*z)/(-1 + z^2)^2

"""
    _spin_orbit_correction2(a, p, e, x, correction) -> NamedTuple

Port of `KerrSpinOrbitCorrection2[a,p,e,x,correction]` (SpinningTrajectory.wl L748-906): assembles
the linear-in-spin orbit-correction object from the numeric `Correction` engine output
`correction` (here the `_numeric_correction` NamedTuple). Builds the geodesic constants/Mino
frequencies, the δχ/δ𝔯/δ𝔷 Fourier series, the four-velocity correction pieces, the UtS/UφS-list →
2-torus DFT → ΔtS/ΔφS trajectory-time/φ deltas, and the frequency corrections {ΓS,ΥrS,ΥθS,ΥφS}.

Returns a NamedTuple with keys mirroring the .wl Association: `a, p, e, ScriptI, Ehat, Lzhat, Khat,
MinoFrequenciesGeo=(Γ̂,Υr̂,Υθ̂,Υφ̂), BLFrequenciesGeo, MinoFrequenciesCorrection=(ΓS,ΥrS,ΥθS,ΥφS),
BLFrequenciesCorrection, ES, LS, error, TrajectoryGeo, TrajectoryDeltas, TrajectoryCorrection,
FourVelocityCorrection, Spin4vector, OrbitCorrection`. The function-valued keys are Julia closures of
`(wr,wz)` (the phase variables). R3.
"""
function _spin_orbit_correction2(a, p, e, x, correction)
    ScriptI = acos(x)
    nmax = correction.nmax
    kmax = correction.kmax
    stepsr = 4 * nmax
    stepsθ = 4 * kmax

    # --- geodesic constants & Mino frequencies (.wl L823-828; KerrGeoOrbit "Phases") ---
    Ehat = geodesic_energy(a, p, e, x)
    Lzhat = geodesic_angular_momentum(a, p, e, x)
    Qhat = carter_constant(a, p, e, x)
    Khat = Qhat + (a * Ehat - Lzhat)^2
    fr = kerr_frequencies(a, p, e, x; Time="Mino")
    Υrhat = fr.Ur; Υθhat = fr.Uth; Υφhat = fr.Uph; Γhat = fr.Ut

    # geodesic trajectory closures (phase-variable input, mirror orbit["Trajectory"][[2]]/[[3]])
    rad = geodesic_trajectory_radial(a, p, e, x)
    pol = geodesic_trajectory_polar(a, p, e, x)
    rhat(wr) = rad(wr).r
    zhat(wz) = cos(pol(wz).θ)

    # --- correction-engine outputs (the `_numeric_correction` NamedTuple) ---
    δχrSfourIm = correction.δχrSIm           # length 2nmax
    δχθSfourIm = correction.δχθSIm           # length kmax
    ΥrS = correction.ΥrS; ΥθS = correction.ΥθS
    udt0S = correction.udt0S; udφ0S = correction.udφ0S
    # δ𝔯four/δ𝔷four are flat v-slices; the .wl `Partition[·,kmax]` reshapes them row-major.
    # δ𝔯four → (2nmax+1)×kmax,  δ𝔷four → (2nmax)×kmax  (row-major).
    δ𝔯four = permutedims(reshape(correction.δ𝔯four, kmax, 2 * nmax + 1))
    δ𝔷four = permutedims(reshape(correction.δ𝔷four, kmax, 2 * nmax))

    # --- the χrhat/χθhat Darwin-angle maps (.wl L756-757, wr/wz-phase based with JacobiSN).
    #     Equivalent (verified) to the Mino-time _chir_hat/_chith_hat with λ = w/Υ. ---
    χrhat(wr) = _chir_hat(a, p, e, x, wr / Υrhat)
    χθhat(wz) = _chith_hat(a, p, e, x, wz / Υθhat)
    # derivative of the Darwin-angle map wrt the phase: dχ/dw = (dχ/dλ)/Υ. The .wl uses χrhat'[wr];
    # χr(w) is piecewise-smooth so a central finite difference of the phase-map is robust and exact
    # enough (the maps are analytic away from the turning points where Sin[χ]=0 kills the term).
    function dχrhat(wr)
        h = 1e-6
        (χrhat(wr + h) - χrhat(wr - h)) / (2h)
    end
    function dχθhat(wz)
        h = 1e-6
        (χθhat(wz + h) - χθhat(wz - h)) / (2h)
    end

    # --- δχrS/δχθS and δ𝔯/δ𝔷 Fourier series (.wl L758-765).
    #     NB: the .wl indexes δχrSfourIm[[n+nmax]] (k=1..nmax → SECOND half of the 2nmax vector) and
    #     δχθSfourIm[[k+kmax/2]] (k=1..kmax/2 → second half of the kmax vector). ---
    δχrS(wr) = 2 * sum(δχrSfourIm[n + nmax] * sin(n * wr) for n in 1:nmax)
    dδχrS(wr) = 2 * sum(n * Υrhat * δχrSfourIm[n + nmax] * cos(n * wr) for n in 1:nmax)
    δχθS(wz) = 2 * sum(δχθSfourIm[k + kmax ÷ 2] * sin(2k * wz) for k in 1:(kmax ÷ 2))
    dδχθS(wz) = 2 * sum(2k * Υθhat * δχθSfourIm[k + kmax ÷ 2] * cos(2k * wz) for k in 1:(kmax ÷ 2))
    # δ𝔯[wr,wz] = 2 Σ_{n=-nmax..nmax, k=1..kmax/2} δ𝔯four[n+nmax+1, k+kmax/2] cos(n wr + 2k wz)  (.wl L762)
    δ𝔯(wr, wz) = 2 * sum(δ𝔯four[n + nmax + 1, k + kmax ÷ 2] * cos(n * wr + 2k * wz)
                          for n in -nmax:nmax, k in 1:(kmax ÷ 2))
    dδ𝔯(wr, wz) = -2 * sum((n * Υrhat + 2k * Υθhat) * δ𝔯four[n + nmax + 1, k + kmax ÷ 2] *
                            sin(n * wr + 2k * wz) for n in -nmax:nmax, k in 1:(kmax ÷ 2))
    # δ𝔷[wr,wz] = 2 Σ_{n=1..nmax, k=-kmax/2..kmax/2-1} δ𝔷four[n+nmax, k+kmax/2+1] cos(n wr +(2k+1) wz) (.wl L764)
    δ𝔷(wr, wz) = 2 * sum(δ𝔷four[n + nmax, k + kmax ÷ 2 + 1] * cos(n * wr + (2k + 1) * wz)
                          for n in 1:nmax, k in (-(kmax ÷ 2)):(kmax ÷ 2 - 1))
    dδ𝔷(wr, wz) = -2 * sum((n * Υrhat + (2k + 1) * Υθhat) * δ𝔷four[n + nmax, k + kmax ÷ 2 + 1] *
                            sin(n * wr + (2k + 1) * wz) for n in 1:nmax, k in (-(kmax ÷ 2)):(kmax ÷ 2 - 1))

    # --- geodesic velocities along the phase variable (.wl L766-768) ---
    Δ(r) = r^2 - 2r + a^2
    Σ(r, z) = r^2 + a^2 * z^2
    drdλ(wr) = sign(pi - mod(wr, 2pi)) *
               sqrt((Ehat * (rhat(wr)^2 + a^2) - a * Lzhat)^2 - Δ(rhat(wr)) * (rhat(wr)^2 + Khat))
    dzdλ(wz) = -sign(pi - mod(wz, 2pi)) *
               sqrt(Qhat - zhat(wz)^2 * (a^2 * (1 - Ehat^2) * (1 - zhat(wz)^2) + Lzhat^2 + Qhat))

    # --- spin 4-vector components (.wl L771-774) ---
    St(wr, wz) = -1 / sqrt(Khat) * (-a * (rhat(wr) * dzdλ(wz) + zhat(wz) * drdλ(wr)) / Σ(rhat(wr), zhat(wz)))
    Sr(wr, wz) = -1 / sqrt(Khat) * (a * zhat(wz) * ((rhat(wr)^2 + a^2) * Ehat - a * Lzhat) / Δ(rhat(wr)))
    Sz(wr, wz) = -1 / sqrt(Khat) * (a * rhat(wr) * Ehat - rhat(wr) * Lzhat / (1 - zhat(wz)^2))
    Sφ(wr, wz) = -1 / sqrt(Khat) *
                 ((a^2 * zhat(wz) * (1 - zhat(wz)^2) * drdλ(wr) + rhat(wr) * (rhat(wr)^2 + a^2) * dzdλ(wz)) /
                  Σ(rhat(wr), zhat(wz)))

    # --- position corrections rS/zS (.wl L785-786) ---
    rS(wr, wz) = e * p * δχrS(wr) * sin(χrhat(wr)) / (1 + e * cos(χrhat(wr)))^2 + δ𝔯(wr, wz)
    zS(wr, wz) = -sin(ScriptI) * δχθS(wz) * sin(χθhat(wz)) + δ𝔷(wr, wz)

    # --- ES, LS scalar zero-modes (midpoint Mino-grid averages, .wl L855-856) ---
    # grid points wr_i = (i-1/2)·2π/stepsr, wz_j = (j-1/2)·2π/stepsθ (phase variables on [0,2π])
    wrg = [(i - 0.5) * 2pi / stepsr for i in 1:stepsr]
    wzg = [(j - 0.5) * 2pi / stepsθ for j in 1:stepsθ]
    h1grid = sum(_oc_h1(rhat(wrg[i]), zhat(wzg[j]), Khat, a) for i in 1:stepsr, j in 1:stepsθ)
    ES = -udt0S + h1grid / (stepsr * stepsθ)
    h2grid = sum(_oc_h2(rhat(wrg[i]), zhat(wzg[j]), Khat, Ehat, Lzhat, a) +
                 _oc_h3(rhat(wrg[i]), zhat(wzg[j]), Khat, a) * drdλ(wrg[i]) * dzdλ(wzg[j])
                 for i in 1:stepsr, j in 1:stepsθ)
    LS = udφ0S - h2grid / (stepsr * stepsθ)

    # --- four-velocity correction zero-modes (.wl L857-858) ---
    udtS(wr, wz) = -ES + _oc_h1(rhat(wr), zhat(wz), Khat, a)
    udφS(wr, wz) = LS + (_oc_h2(rhat(wr), zhat(wz), Khat, Ehat, Lzhat, a) +
                         _oc_h3(rhat(wr), zhat(wz), Khat, a) * drdλ(wr) * dzdλ(wz))

    # --- UtS/UφS lists → 2-torus DFT → Im[ΔtS/ΔφS]four (.wl L818-875) ---
    UtS(wr, wz) = (_oc_dSgttdr(rhat(wr), a) * rS(wr, wz) + _oc_dSgttdz(zhat(wz), a) * zS(wr, wz)) * (-Ehat) +
                  _oc_Sgtt(rhat(wr), zhat(wz), a) * udtS(wr, wz) +
                  _oc_dSgtphdr(rhat(wr), a) * rS(wr, wz) * Lzhat +
                  _oc_Sgtph(rhat(wr), a) * udφS(wr, wz)
    UφS(wr, wz) = (_oc_dSgphphdr(rhat(wr), a) * rS(wr, wz) + _oc_dSgphphdz(zhat(wz), a) * zS(wr, wz)) * Lzhat +
                  _oc_Sgphph(rhat(wr), zhat(wz), a) * udφS(wr, wz) +
                  _oc_dSgtphdr(rhat(wr), a) * rS(wr, wz) * (-Ehat) +
                  _oc_Sgtph(rhat(wr), a) * udtS(wr, wz)
    Uthat(wr, wz) = _oc_Sgtt(rhat(wr), zhat(wz), a) * (-Ehat) + _oc_Sgtph(rhat(wr), a) * Lzhat
    Uφhat(wr, wz) = _oc_Sgphph(rhat(wr), zhat(wz), a) * Lzhat + _oc_Sgtph(rhat(wr), a) * (-Ehat)

    UtSlist   = [UtS(wrg[i], wzg[j])   for i in 1:stepsr, j in 1:stepsθ]
    UφSlist   = [UφS(wrg[i], wzg[j])   for i in 1:stepsr, j in 1:stepsθ]
    Uthatlist = [Uthat(wrg[i], wzg[j]) for i in 1:stepsr, j in 1:stepsθ]
    Uφhatlist = [Uφhat(wrg[i], wzg[j]) for i in 1:stepsr, j in 1:stepsθ]

    Eni = [cis(2pi * n * (i - 0.5) / stepsr) for n in -nmax:nmax, i in 1:stepsr]
    Ejk = [cis(2pi * k * (j - 0.5) / stepsθ) for j in 1:stepsθ, k in -kmax:kmax]
    dft(list) = Eni * list * Ejk / (stepsr * stepsθ)
    UtSfour   = dft(UtSlist);   UφSfour   = dft(UφSlist)
    Uthatfour = dft(Uthatlist); Uφhatfour = dft(Uφhatlist)

    # ImΔtSfour[n,k] / ImΔφSfour[n,k] (.wl L874-875), indexed [n+nmax+1, k+kmax+1]
    ImΔtSfour  = zeros(Float64, 2 * nmax + 1, 2 * kmax + 1)
    ImΔφSfour  = zeros(Float64, 2 * nmax + 1, 2 * kmax + 1)
    for n in -nmax:nmax, k in -kmax:kmax
        (n == 0 && k == 0) && continue
        den = n * Υrhat + k * Υθhat
        ImΔtSfour[n + nmax + 1, k + kmax + 1] =
            real(UtSfour[n + nmax + 1, k + kmax + 1]) / den -
            real(Uthatfour[n + nmax + 1, k + kmax + 1]) / den^2 * (n * ΥrS + k * ΥθS)
        ImΔφSfour[n + nmax + 1, k + kmax + 1] =
            real(UφSfour[n + nmax + 1, k + kmax + 1]) / den -
            real(Uφhatfour[n + nmax + 1, k + kmax + 1]) / den^2 * (n * ΥrS + k * ΥθS)
    end
    ΓS  = real(UtSfour[nmax + 1, kmax + 1])
    ΥφS = real(UφSfour[nmax + 1, kmax + 1])

    # --- ΔtS/ΔφS trajectory-time/φ deltas (.wl L789-790) ---
    ΔtS(wr, wz) =
        sum(2 * ImΔtSfour[nmax + 1, k + kmax + 1] * sin(k * wz) for k in 1:kmax) +
        sum(2 * ImΔtSfour[n + nmax + 1, k + kmax + 1] * sin(n * wr + k * wz)
            for n in 1:nmax, k in -kmax:kmax)
    ΔφS(wr, wz) =
        sum(2 * ImΔφSfour[nmax + 1, k + kmax + 1] * sin(k * wz) for k in 1:kmax) +
        sum(2 * ImΔφSfour[n + nmax + 1, k + kmax + 1] * sin(n * wr + k * wz)
            for n in 1:nmax, k in -kmax:kmax)

    # --- OrbitCorrection: assoc-returning closure (.wl L791-817, returns ΔtS/rS/zS/ΔφS/UrS/UzS) ---
    function OrbitCorrection(wr, wz)
        δχrSp = δχrS(wr); χrhatp = χrhat(wr)
        δχθSp = δχθS(wz); χθhatp = χθhat(wz)
        return (ΔtS = ΔtS(wr, wz), rS = rS(wr, wz), zS = zS(wr, wz), ΔφS = ΔφS(wr, wz),
                UrS = p * e * (dχrhat(wr) * (cos(χrhatp) * δχrSp * Υrhat + sin(χrhatp) * ΥrS +
                       2 * e * sin(χrhatp)^2 * δχrSp * Υrhat / (1 + e * cos(χrhatp))) +
                       sin(χrhatp) * dδχrS(wr)) / (1 + e * cos(χrhatp))^2 + dδ𝔯(wr, wz),
                UzS = -sin(ScriptI) * (dχθhat(wz) * (cos(χθhatp) * δχθSp * Υθhat + sin(χθhatp) * ΥθS) +
                       sin(χθhatp) * dδχθS(wz)) + dδ𝔷(wr, wz))
    end

    # four-velocity correction tuple (udtS, UrS-component, UzS-component, udφS) (.wl L890)
    FourVelocityCorrection = (
        udt = (wr, wz) -> udtS(wr, wz),
        ur = (wr, wz) -> p * e * (dχrhat(wr) * (cos(χrhat(wr)) * δχrS(wr) * Υrhat + sin(χrhat(wr)) * ΥrS +
              2 * e * sin(χrhat(wr))^2 * δχrS(wr) * Υrhat / (1 + e * cos(χrhat(wr)))) +
              sin(χrhat(wr)) * dδχrS(wr)) / (1 + e * cos(χrhat(wr)))^2 + dδ𝔯(wr, wz),
        uz = (wr, wz) -> -sin(ScriptI) * (dχθhat(wz) * (cos(χθhat(wz)) * δχθS(wz) * Υθhat +
              sin(χθhat(wz)) * ΥθS) + sin(χθhat(wz)) * dδχθS(wz)) + dδ𝔷(wr, wz),
        udφ = (wr, wz) -> udφS(wr, wz))

    return (; a = float(a), p = float(p), e = float(e), ScriptI,
            Ehat, Lzhat, Khat,
            TrajectoryCorrection = ((wr, wz) -> ΔtS(wr, wz), (wr, wz) -> rS(wr, wz),
                                    (wr, wz) -> zS(wr, wz), (wr, wz) -> ΔφS(wr, wz)),
            TrajectoryGeo = (rad, pol),
            MinoFrequenciesGeo = (Γhat, Υrhat, Υθhat, Υφhat),
            BLFrequenciesGeo = (Υrhat / Γhat, Υθhat / Γhat, Υφhat / Γhat),
            FourVelocityCorrection,
            Spin4vector = ((wr, wz) -> St(wr, wz), (wr, wz) -> Sr(wr, wz),
                           (wr, wz) -> Sz(wr, wz), (wr, wz) -> Sφ(wr, wz)),
            TrajectoryDeltas = (Δtr = (wr) -> rad(wr).Δtr, Δφr = (wr) -> rad(wr).Δφr,
                                Δtz = (wz) -> pol(wz).Δtz, Δφz = (wz) -> pol(wz).Δφz),
            MinoFrequenciesCorrection = (ΓS, ΥrS, ΥθS, ΥφS),
            BLFrequenciesCorrection = (ΥrS / Γhat - Υrhat / Γhat^2 * ΓS,
                                       ΥθS / Γhat - Υθhat / Γhat^2 * ΓS,
                                       ΥφS / Γhat - Υφhat / Γhat^2 * ΓS),
            ES, LS,
            δχrS, δχθS, δr = δ𝔯, δz = δ𝔷,
            OrbitCorrection,
            error = correction.error,
            nmax, kmax)
end

"""
    spin_orbit_correction(a, p, e, x, nmax, kmax) -> NamedTuple

Port of `KerrSpinOrbitCorrection[a,p,e,x,nmax,kmax]` (SpinningTrajectory.wl L907-911): runs the
numeric `Correction` engine (`_numeric_correction`) then the `_spin_orbit_correction2` assembly.
Returns the orbit-correction NamedTuple (see `_spin_orbit_correction2`). `kmax` must be even. R3.
"""
function spin_orbit_correction(a, p, e, x, nmax::Integer, kmax::Integer)
    iseven(kmax) || throw(ArgumentError("kmax must be even (got $kmax)"))
    nmax > 0 || throw(ArgumentError("nmax must be a positive integer (got $nmax)"))
    correction = _numeric_correction(a, p, e, x, nmax, kmax)
    return _spin_orbit_correction2(a, p, e, x, correction)
end

"""
    spin_orbit_numeric(orbitCorrection, spar) -> NamedTuple

Port of `KerrSpinOrbitNum[orbitCorrection,spar]` (SpinningTrajectory.wl L914-938): the thin
frequency assembly — `{Γ,Υr,Υθ,Υφ} = MinoFrequenciesGeo + spar·MinoFrequenciesCorrection`, then the
BL frequencies `Υ/Γ`. Carries through `a,p,e,ScriptI,s,Ehat,Lzhat,Khat,TrajectoryGeo,
TrajectoryDeltas,MinoFrequencies(Geo/Correction),ES,LS,OrbitCorrection`. `orbitCorrection` is the
NamedTuple returned by `spin_orbit_correction`. R3.
"""
function spin_orbit_numeric(orbitCorrection, spar)
    mfg = orbitCorrection.MinoFrequenciesGeo
    mfc = orbitCorrection.MinoFrequenciesCorrection
    Γ  = mfg[1] + spar * mfc[1]
    Υr = mfg[2] + spar * mfc[2]
    Υθ = mfg[3] + spar * mfc[3]
    Υφ = mfg[4] + spar * mfc[4]
    return (; a = orbitCorrection.a, p = orbitCorrection.p, e = orbitCorrection.e,
            ScriptI = orbitCorrection.ScriptI, s = spar,
            Ehat = orbitCorrection.Ehat, Lzhat = orbitCorrection.Lzhat, Khat = orbitCorrection.Khat,
            BLFrequencies = (Υr / Γ, Υθ / Γ, Υφ / Γ),
            TrajectoryGeo = orbitCorrection.TrajectoryGeo,
            TrajectoryDeltas = orbitCorrection.TrajectoryDeltas,
            MinoFrequencies = (Γ, Υr, Υθ, Υφ),
            MinoFrequenciesGeo = mfg,
            MinoFrequenciesCorrection = mfc,
            ES = orbitCorrection.ES, LS = orbitCorrection.LS,
            OrbitCorrection = orbitCorrection.OrbitCorrection)
end
