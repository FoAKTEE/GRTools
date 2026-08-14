# Numeric linear-in-spin orbit-correction engine (spin_orbit_correction via the 2-torus DFT) —
# consolidates SpinningOrbit/src/numeric_correction.jl (port of SpinningTrajectory.wl `Correction`,
# L507-746). The geodesic spin-curvature force is Fourier-decomposed over the (λr,λθ) 2-torus, then
# the orbit correction is solved linearly: the 47 algebraic force/coefficient terms (design §4.1's
# "~60 spin-curvature force terms" together with the 4 worldline deviations and the 11 _oc_*
# four-velocity/metric helpers), the Darwin-angle maps, the DFT machinery and the M-matrix
# least-squares solve. Σ=r²+a²z², Δ=r²−2r+a², K=K̃ (Carter). The algebraic force terms are pure
# arithmetic ⇒ type-generic {T}; the DFT/solve engine is Float64-tier (LAPACK pinv), faithful to
# the oracle. Serial (no threaded LAPACK region ⇒ no §8 BLAS-guard needed). GENERIC (inclined)
# bound orbits; the oracle evaluates it at x²<1 — at x²=1 the Darwin θ-angle map degenerates
# (sin(acos(x))=0), exactly as in the oracle. R3 (load-bearing numerics).

"`𝓡t` source for the BL-time spin shift δu_{ds,t} (SpinningTrajectory.wl L528). R3."
function _Rt_force(r, z, dr, dz, K, a)
    Σ = r^2 + a^2 * z^2
    num = 2 * a^2 * dz * r * z * (3 * r^2 * (2K + r^2) - 2 * a^2 * (3K + 4 * r^2) * z^2 + a^4 * z^4) +
          dr * (3K * r^4 - 9 * a^2 * r^2 * (2K + r^2) * z^2 + a^4 * (3K + 14 * r^2) * z^4 - a^6 * z^6)
    return -num / (sqrt(K) * Σ^4)
end

"`𝓡φ` source for the azimuthal spin shift δu_{ds,φ} (SpinningTrajectory.wl L529). R3."
function _Rphi_force(r, z, dr, dz, K, a)
    Σ = r^2 + a^2 * z^2
    num = a * (-2 * dz * r * (a^2 + r^2) * z * (3 * r^2 * (2K + r^2) - 2 * a^2 * (3K + 4 * r^2) * z^2 + a^4 * z^4) +
               dr * (-1 + z^2) * (3K * r^4 - 9 * a^2 * r^2 * (2K + r^2) * z^2 + a^4 * (3K + 14 * r^2) * z^4 - a^6 * z^6))
    return -num / (sqrt(K) * Σ^4)
end

"`𝓘2` coefficient (couples δu_{ds,t} into the energy/Carter correction; .wl L543). R3."
function _I2_coef(r, z, En, Lz, K, a)
    return 2 * (a^3 * (a * En - Lz) + r^2 * (a * (2 * a * En + Lz) + En * (-4 + r) * r)) / (a^2 + (-2 + r) * r) -
           (4 * r^2 * (a * (a * En - Lz) + En * r^2)) / (r^2 + a^2 * z^2)
end

"`𝓘3` coefficient (couples δu_{ds,φ}; .wl L544). R3."
function _I3_coef(r, z, En, Lz, K, a)
    return (2 * (Lz * (-2 + r)^2 * r^3 - a * En * r^3 * (-4 + 3r) * (-1 + z^2) + a^5 * En * z^2 * (-1 + z^2) +
                 a^4 * Lz * z^2 * (1 + (-1 + r) * z^2) + a^2 * Lz * r^2 * (-1 + (-3 + 2r) * z^2) -
                 a^3 * En * r^2 * (-1 + z^4))) / ((a^2 + (-2 + r) * r) * (-1 + z^2) * (r^2 + a^2 * z^2))
end

"`𝓙` energy/Carter source (.wl L545): geodesic part + 𝓘2·δu_t + 𝓘3·δu_φ."
function _J_source(r, z, dut, duphi, En, Lz, K, a)
    geo = ((a^2 * En - a * Lz + En * r^2) *
           (-3K * r^4 + 9 * a^2 * r^2 * (2K + r^2) * z^2 - a^4 * (3K + 14 * r^2) * z^4 + a^6 * z^6)) /
          (sqrt(K) * (r^2 + a^2 * z^2)^3)
    return geo + _I2_coef(r, z, En, Lz, K, a) * dut + _I3_coef(r, z, En, Lz, K, a) * duphi
end

# --- geodesic 2nd derivatives wrt Mino λ (SpinningTrajectory.wl L520-521) ---
"`d²r/dλ²` along the geodesic (.wl L520). R3."
_d2r_dλ2(r, En, Lz, K, a) =
    0.5 * (-((-2 + 2r) * (K + r^2)) - 2r * (a^2 - 2r + r^2) + 4 * En * r * (-a * Lz + En * (a^2 + r^2)))
"`d²z/dλ²` along the geodesic (z=cosθ; .wl L521, uses Carter Q). R3."
_d2z_dλ2(z, En, Lz, Q, a) =
    (2 * a^2 * (1 - En^2) * z^3 - 2z * (Lz^2 + Q + a^2 * (1 - En^2) * (1 - z^2))) / 2

# --- the simpler r/θ/norm force components (SpinningTrajectory.wl L533,535,541-542,547,549) ---
"`𝓖ℛ` radial-metric force (.wl L533). R3."
_GR_force(r, z, dr, a) = (-2 * (-1 + r) / (a^2 + (-2 + r) * r) - 2r / (r^2 + a^2 * z^2)) * dr
"`𝓖ℨ` polar-metric force (.wl L535). R3."
_GZ_force(r, z, dz, a) = (2r * (a^2 + (-2 + r) * r) * dz) / ((-1 + z^2) * (r^2 + a^2 * z^2))
"`𝓘1r` radial inhomogeneity (.wl L541, Υr the geodesic radial Mino freq). R3."
_I1r_coef(r, z, dr, d2r, Υr, a) =
    (2 * (r^2 * (dr^2 * (3 - 2r) + d2r * (-2 + r) * r) + a^4 * d2r * z^2 +
          a^2 * (d2r * r * (r + (-2 + r) * z^2) - dr^2 * (r + (-1 + r) * z^2)))) /
    ((a^2 + (-2 + r) * r) * (r^2 + a^2 * z^2) * Υr)
"`𝓘1θ` polar inhomogeneity (.wl L542, Υθ the geodesic polar Mino freq). R3."
_I1θ_coef(r, z, dz, Υθ, a) =
    (2 * dz^2 * r * (a^2 + (-2 + r) * r) / Υθ) / ((-1 + z^2) * (r^2 + a^2 * z^2))
"`𝓚ℛ` radial norm-force (.wl L547). R3."
_KR_force(r, z, dr, a) = (2 * dr) / ((a^2 + (-2 + r) * r) * (r^2 + a^2 * z^2))
"`𝓚ℨ` polar norm-force (.wl L549). R3."
_KZ_force(r, z, dz, a) = -((2 * dz) / ((-1 + z^2) * (r^2 + a^2 * z^2)))

# --- 2nd-order metric forces 𝓗ℛ (.wl L538), 𝓗ℨ (.wl L540) — MMA→Julia auto-translated
#     (gen_numcorr_jl.wl, InputForm+head-subst); verbatim-oracle verified. R3. ---
"`𝓗ℛ` 2nd-order radial metric force (.wl L538). R3."
_HR_force(r, z, dr, dz, En, Lz, a) = (dz^2*(-(a^2*r^2) + r^4 + a^2*(a^2 + r*(-4 + 3*r))*z^2))/((-1 + z^2)*(r^2 + a^2*z^2)^2) + dr^2*((2 - a^2 + (-2 + r)*r)/(a^2 + (-2 + r)*r)^2 + (2*r^2)/(r^2 + a^2*z^2)^2 - (r^2 + a^2*z^2)^(-1)) + ((2*a^2*(-2*a^2*En^2 + 4*a*En*Lz + (-2 + a^2)*Lz^2) + 2*r*(2*a^2*(-(a^2*En^2) - 4*a*En*Lz + Lz^2) + r*(a*(16*a*En^2 + 4*En*Lz - a*Lz^2) - 2*En^2*r*(8 + 2*a^2 - 5*r + r^2))))/(a^2 + (-2 + r)*r)^2 + (Lz^2*(a^4 - 4*a^2*r + 2*a^2*r^2 + r^4))/((a^2 + r^2)^2*(-1 + z)) - (Lz^2*(a^4 - 4*a^2*r + 2*a^2*r^2 + r^4))/((a^2 + r^2)^2*(1 + z)) - (8*r^3*(a*(a*En - Lz) + En*r^2)^2)/((a^2 + r^2)*(r^2 + a^2*z^2)^2) + (8*r*(a*(a*En - Lz) + En*r^2)*(a^3*(a*En - Lz) + 3*a^2*En*r^2 + 2*En*r^4))/((a^2 + r^2)^2*(r^2 + a^2*z^2)))/2

"`𝓗ℨ` 2nd-order polar metric force (.wl L540). R3."
_HZ_force(r, z, dr, dz, En, Lz, a) = (2*a^2*dr^2*r*z)/(r^2 + a^2*z^2)^2 - (2*dz^2*r*(a^2 + (-2 + r)*r)*z*(r^2 + a^2*(-1 + 2*z^2)))/((-1 + z^2)^2*(r^2 + a^2*z^2)^2) + (2*r*z*(-(Lz^2*r^4) - 2*a^2*Lz^2*r^2*z^2 - a^4*Lz^2*z^4 - 2*r^3*(-(a*En) - Lz + a*En*z^2)*(-(a*En) + Lz + a*En*z^2) - 2*a^2*r*(-(a*En) + Lz + a*En*z^2)*(-(a*En) + Lz + (a*En - 2*Lz)*z^2)))/((-1 + z^2)^2*(r^2 + a^2*z^2)^2)

# --- energy/Lz/Carter-correction coefficient layer (𝓝/𝓤/𝓠/𝓜ℨ + 𝓥/𝓟 sources); auto-translated. R3. ---
"`𝓝1r` (.wl L554)."
_N1r_coef(r, z, dr, Ur, a) = (2*dr^2)/((a^2 + (-2 + r)*r)*Ur*(r^2 + a^2*z^2))
"`𝓝1θ` (.wl L555)."
_N1th_coef(r, z, dz, Uth, a) = (-2*dz^2)/(Uth*(-1 + z^2)*(r^2 + a^2*z^2))
"`𝓝2` (.wl L556)."
_N2_coef(r, z, En, Lz, a) = 2*En + (4*r*(a*(a*En - Lz) + En*r^2))/((a^2 + (-2 + r)*r)*(r^2 + a^2*z^2))
"`𝓝3` (.wl L557)."
_N3_coef(r, z, En, Lz, a) = (-2*r*(2*a*En - 2*Lz + Lz*r) - 2*a*(a*Lz - 2*En*r)*z^2)/((a^2 + (-2 + r)*r)*(-1 + z^2)*(r^2 + a^2*z^2))
"`𝓤2` (.wl L574)."
_U2_coef(r, z, En, Lz, a) = (4*a^2*r*(a^2*En - a*Lz + En*r^2)*z*sqrt(1 - z^2))/((a^2 + (-2 + r)*r)*(r^2 + a^2*z^2))
"`𝓤3` (.wl L575)."
_U3_coef(r, z, En, Lz, a) = (-2*z*(Lz*(-2 + r)*r^3 + a^4*Lz*z^4 - 2*a^3*En*r*(-1 + z^2)^2 + 2*a^2*Lz*r*(1 + (-2 + r)*z^2)))/((a^2 + (-2 + r)*r)*(1 - z^2)^(3/2)*(r^2 + a^2*z^2))
"`𝓠ℨ` (.wl L560)."
_QZ_force(z) = -(1/sqrt(1 - z^2))
"`𝓜ℨ` polar Carter force (.wl L553)."
_MZ_force(r, z, dr, dz, En, Lz, a) = (-2*a^2*dr^2*z)/((a^2 + (-2 + r)*r)*(r^2 + a^2*z^2)^2) + (2*dz^2*z*(r^2 + a^2*(-1 + 2*z^2)))/((-1 + z^2)^2*(r^2 + a^2*z^2)^2) + (2*z*(Lz^2*r^4 + 2*a^2*Lz^2*r^2*z^2 + a^4*Lz^2*z^4 + 2*r^3*(-(a*En) - Lz + a*En*z^2)*(-(a*En) + Lz + a*En*z^2) + 2*a^2*r*(-(a*En) + Lz + a*En*z^2)*(-(a*En) + Lz + (a*En - 2*Lz)*z^2)))/((a^2 + (-2 + r)*r)*(-1 + z^2)^2*(r^2 + a^2*z^2)^2)
"`𝓜𝔯` norm-equation 2nd-order radial metric force (.wl L551); MMA→Julia auto-translated (gen_numcorr_MR.wl). R3."
_MR_force(r, z, dr, dz, En, Lz, a) = (2*dz^2*r)/((-1 + z^2)*(r^2 + a^2*z^2)^2) + (2*dr^2*(-(r*(a^2 + r*(-3 + 2*r))) - a^2*(-1 + r)*z^2))/((a^2 + (-2 + r)*r)^2*(r^2 + a^2*z^2)^2) + (2*(Lz^2*r^5 + a^4*Lz^2*r*z^4 + En^2*r^6*(-1 + z^2) - a^4*(-(a*En) + Lz)^2*z^2*(-1 + z^2) - r^4*(2*(a*En - 2*Lz)*(a*En - Lz) + a*En*z^2*(-3*a*En + 6*Lz + a*En*z^2)) + a^2*r^2*(-(-(a*En) + Lz)^2 - (a*En + Lz)*z^2*(-3*a*En + 3*Lz + 2*a*En*z^2)) + 2*r^3*(2*(-(a*En) + Lz)^2 + a*z^2*(-4*a*En^2 + 4*En*Lz + a*Lz^2 + 2*a*En^2*z^2))))/((a^2 + (-2 + r)*r)^2*(-1 + z^2)*(r^2 + a^2*z^2)^2)
"`𝓥`geo (.wl L576 geodesic term)."
_Vgeo(r, z, En, Lz, K, a) = (-2*a*r*z*(3*r^2*(2*K + r^2) - 2*a^2*(3*K + 4*r^2)*z^2 + a^4*z^4)*(Lz + a*En*(-1 + z^2)))/(sqrt(K)*sqrt(1 - z^2)*(r^2 + a^2*z^2)^3)
"`𝓥` Lz/Carter source (.wl L576): geo + 𝓤2·δu_t + 𝓤3·δu_φ."
_V_source(r, z, dut, duphi, En, Lz, K, a) = _Vgeo(r, z, En, Lz, K, a) + _U2_coef(r, z, En, Lz, a)*dut + _U3_coef(r, z, En, Lz, a)*duphi
"`𝓟` source (.wl L558): 𝓝2·δu_t + 𝓝3·δu_φ."
_P_source(r, z, dut, duphi, En, Lz, a) = _N2_coef(r, z, En, Lz, a)*dut + _N3_coef(r, z, En, Lz, a)*duphi

# --- χ/sin(I)-dependent r/θ inhomogeneous forces (Darwin-angle-driven); auto-translated. R3. ---
"`𝓕r` (.wl L530)."
_Fr_force(chir, e, p) = (e*p*sin(chir))/(1 + e*cos(chir))^2
"`𝓖r` (.wl L532)."
_Gr_force(r, z, dr, chir, a, e, p) = (2*dr*e*r*(r^3*(p*(-a^2 + p) - (a^2*(-1 + e^2) + p^2)*r + (-1 + e^2 + p)*r^2) + a^2*(a^2*p^2 + r*(-(p*(3*a^2 + p)) + r*(-2*a^2*(-1 + e^2) + 4*p - r*(3 - 3*e^2 + p + (-1 + e^2)*r))))*z^2)*sin(chir))/(p*(a^2 + (-2 + r)*r)*(p^2 + r*(-2*p + r - e^2*r))*(r^2 + a^2*z^2))
"`𝓖θ` (.wl L534)."
_Gth_force(r, z, dz, chith, incl, a) = (-2*dz*r*(a^2 + (-2 + r)*r)*sin(chith)*sin(incl))/((-1 + z^2)*(r^2 + a^2*z^2))
"`𝓚r` (.wl L546)."
_Kr_force(r, z, dr, chir, a, e, p) = (2*dr*e*r^2*sin(chir))/(p*(a^2 + (-2 + r)*r)*(r^2 + a^2*z^2))
"`𝓚θ` (.wl L548)."
_Kth_force(r, z, dz, chith, incl, a) = (2*dz*sin(chith)*sin(incl))/((-1 + z^2)*(r^2 + a^2*z^2))
"`𝓠θ` (.wl L559)."
_Qth_force(z, chith, incl) = (sin(chith)*sin(incl))/sqrt(1 - z^2)
"`𝓢r` (.wl L561)."
_Sr_force(r, z, dr, chir, a, e, p) = (2*a^2*dr*e*r^2*z*sqrt(1 - z^2)*sin(chir))/(p*(a^2 + (-2 + r)*r)*(r^2 + a^2*z^2))
"`𝓢ℛ` (.wl L562)."
_SR_force(r, z, dr, a) = (2*a^2*dr*z*sqrt(1 - z^2))/((a^2 + (-2 + r)*r)*(r^2 + a^2*z^2))
"`𝓢θ` (.wl L563)."
_Sth_force(r, z, dz, chith, incl, a) = (-2*dz*z*sqrt(1 - z^2)*csc(chith)*(r^2*cot(incl)^2 + a^2*(1 - 2*z^2 + z^4*csc(incl)^2))*sin(incl))/((-1 + z^2)^2*(r^2 + a^2*z^2))
"`𝓢ℨ` (.wl L564)."
_SZ_force(r, z, dz, a) = (-2*dz*z*(r^2 + a^2*(-1 + 2*z^2)))/((1 - z^2)^(3/2)*(r^2 + a^2*z^2))
"`𝓤1r` (.wl L572)."
_U1r_coef(r, z, dr, Ur, a) = (2*a^2*dr^2*z*sqrt(1 - z^2))/((a^2 + (-2 + r)*r)*Ur*(r^2 + a^2*z^2))

# --- heaviest 2nd-order χ/sin(I)-dependent forces (auto-translated from the .wl SOURCE, gen_numcorr_jl5.wl). R3. ---
"`𝓗r` (.wl L536-537)."
_Hr_force(r, z, dr, dz, d2r, chir, En, Lz, a, e, p) = ((d2r*(-p^2 + 3*p*r + 2*(-1 + e^2)*r^2) + dr^2*(-2*(-1 + e^2 + p) - (-1 + e^2)*r - (e*p^2)/(2*(p + (-1 + e)*r)) + (e*p^2)/(2*p - 2*(1 + e)*r) + (2*(-1 + a^2)*(2*a^2*(-1 + e^2 + p) + ((-4 + a^2)*(-1 + e^2) - 4*p + p^2)*r))/(a^2 + (-2 + r)*r)^2 + (((-2 + a^2)*(-1 + e^2) - 2*p + p^2)*(-2 + r))/(a^2 + (-2 + r)*r) + (2*r^3*(p + (-1 + e)*r)*(-p + (1 + e)*r))/(r^2 + a^2*z^2)^2 + (r*(3*p^2 - 8*p*r - 5*(-1 + e^2)*r^2))/(r^2 + a^2*z^2)))*csc(chir))/(e*p*r) + (dz^2*e*r^2*(-(a^2*r^2) + r^4 + a^2*(a^2 + r*(-4 + 3*r))*z^2)*sin(chir))/(p*(-1 + z^2)*(r^2 + a^2*z^2)^2) + (e*r^2*(Lz^2*(-2 + r)^2*r^6 + 4*a*En*Lz*r^6*(-1 + z^2) + 2*En^2*(-3 + r)*r^8*(-1 + z^2) - 4*a^5*En*Lz*r^2*z^2*(-10 + 8*r - z^2)*(-1 + z^2) - 4*a^3*En*Lz*r^3*(-1 + z^2)*(r*(-1 + 2*r) + 2*(-4 + r)*(-1 + r)*z^2) + a^2*r^3*(4*En^2*r^3*(-1 + z^2)*(r + (-3 + r)*z^2) + Lz^2*(r*(2 + 3*(-2 + r)*r) + (-16 + r*(22 + (-10 + r)*r))*z^2)) + 4*a^7*En*Lz*z^2*(-z^2 + z^4 - 2*r*(-1 + z^4)) + a^6*z^2*(-4*En^2*r^2*(-1 + z^2)*(5 - 3*r + (-4 + r)*z^2) + Lz^2*(-4*r + (2 + r*(-6 + 5*r))*z^2 - (2 + (-2 + r)*r)*z^4)) + a^8*z^2*(Lz^2*z^4 - 2*En^2*(-1 + z^2)*(z^2 + r*(-2 + z^2))) + a^4*r^2*(-2*En^2*r*(-1 + z^2)*(-((-1 + r)*r) - 2*(4 + r*(-4 + 3*r))*z^2 + (8 + (-5 + r)*r)*z^4) + Lz^2*z^2*(20 - r*(20 + r*(-7 + z^2)))))*sin(chir))/(p*(a^2 + (-2 + r)*r)^2*(-1 + z^2)*(r^2 + a^2*z^2)^2)
"`𝓗θ` (.wl L539)."
_Hth_force(r, z, dr, dz, chith, incl, En, Lz, a) = (-2*dz^2*r*(a^2 + (-2 + r)*r)*z*csc(chith)*(r^2*cot(incl)^2 + a^2*(1 - 2*z^2 + z^4*csc(incl)^2))*sin(incl))/((-1 + z^2)^2*(r^2 + a^2*z^2)^2) - (2*a^2*dr^2*r*z*sin(chith)*sin(incl))/(r^2 + a^2*z^2)^2 + (2*r*z*(Lz^2*r^4 + 2*a^2*Lz^2*r^2*z^2 + a^4*Lz^2*z^4 + 2*r^3*(-(a*En) - Lz + a*En*z^2)*(-(a*En) + Lz + a*En*z^2) + 2*a^2*r*(-(a*En) + Lz + a*En*z^2)*(-(a*En) + Lz + (a*En - 2*Lz)*z^2))*sin(chith)*sin(incl))/((-1 + z^2)^2*(r^2 + a^2*z^2)^2)
"`𝓜r` (.wl L550)."
_Mr_force(r, z, dr, dz, chir, En, Lz, a, e, p) = (2*dr^2*(r^3*((a^2 - p)*p + r*(a^2*(-1 + e^2) + p^2 - (-1 + e^2 + p)*r)) + a^2*(-(a^2*p^2) + r*(p*(3*a^2 + p) + r*(2*a^2*(-1 + e^2) - 4*p + r*(3 - 3*e^2 + p + (-1 + e^2)*r))))*z^2)*csc(chir))/(e*p*r*(a^2 + (-2 + r)*r)^2*(r^2 + a^2*z^2)^2) + (2*dz^2*e*r^3*sin(chir))/(p*(-1 + z^2)*(r^2 + a^2*z^2)^2) + (2*e*r^2*(Lz^2*r^5 + a^4*Lz^2*r*z^4 + En^2*r^6*(-1 + z^2) - a^4*(-(a*En) + Lz)^2*z^2*(-1 + z^2) - r^4*(2*(a*En - 2*Lz)*(a*En - Lz) + a*En*z^2*(-3*a*En + 6*Lz + a*En*z^2)) + a^2*r^2*(-(-(a*En) + Lz)^2 - (a*En + Lz)*z^2*(-3*a*En + 3*Lz + 2*a*En*z^2)) + 2*r^3*(2*(-(a*En) + Lz)^2 + a*z^2*(-4*a*En^2 + 4*En*Lz + a*Lz^2 + 2*a*En^2*z^2)))*sin(chir))/(p*(a^2 + (-2 + r)*r)^2*(-1 + z^2)*(r^2 + a^2*z^2)^2)
"`𝓜θ` (.wl L552)."
_Mth_force(r, z, dr, dz, chith, incl, En, Lz, a) = (2*dz^2*z*csc(chith)*(r^2*cot(incl)^2 + a^2*(1 - 2*z^2 + z^4*csc(incl)^2))*sin(incl))/((-1 + z^2)^2*(r^2 + a^2*z^2)^2) + (2*a^2*dr^2*z*sin(chith)*sin(incl))/((a^2 + (-2 + r)*r)*(r^2 + a^2*z^2)^2) - (2*z*(Lz^2*r^4 + 2*a^2*Lz^2*r^2*z^2 + a^4*Lz^2*z^4 + 2*r^3*(-(a*En) - Lz + a*En*z^2)*(-(a*En) + Lz + a*En*z^2) + 2*a^2*r*(-(a*En) + Lz + a*En*z^2)*(-(a*En) + Lz + (a*En - 2*Lz)*z^2))*sin(chith)*sin(incl))/((a^2 + (-2 + r)*r)*(-1 + z^2)^2*(r^2 + a^2*z^2)^2)
"`𝓣r` (.wl L565-566)."
_Tr_force(r, z, dr, dz, chir, En, Lz, a, e, p) = (-2*a^2*dz^2*e*r^3*z*sin(chir))/(p*sqrt(1 - z^2)*(r^2 + a^2*z^2)^2) + (2*a^2*dr^2*e*r*z*sqrt(1 - z^2)*(r^3*(p*(-a^2 + p) - (a^2*(-1 + e^2) + p^2)*r + (-1 + e^2 + p)*r^2) + a^2*(a^2*p^2 + r*(-(p*(3*a^2 + p)) + r*(-2*a^2*(-1 + e^2) + 4*p - r*(3 - 3*e^2 + p + (-1 + e^2)*r))))*z^2)*sin(chir))/((a^2 + (-2 + r)*r)^2*(p^3 + p*r*(-2*p + r - e^2*r))*(r^2 + a^2*z^2)^2) + (2*a^2*e*r^2*z*(-(Lz^2*r^5) - a^4*Lz^2*r*z^4 - En^2*r^6*(-1 + z^2) + a^4*(-(a*En) + Lz)^2*z^2*(-1 + z^2) + r^4*(2*(a*En - 2*Lz)*(a*En - Lz) + a*En*z^2*(-3*a*En + 6*Lz + a*En*z^2)) + a^2*r^2*((-(a*En) + Lz)^2 + (a*En + Lz)*z^2*(-3*a*En + 3*Lz + 2*a*En*z^2)) + 2*r^3*(-2*(-(a*En) + Lz)^2 - a*z^2*(-4*a*En^2 + 4*En*Lz + a*Lz^2 + 2*a*En^2*z^2)))*sin(chir))/(p*(a^2 + (-2 + r)*r)^2*sqrt(1 - z^2)*(r^2 + a^2*z^2)^2)
"`𝓣ℛ` (.wl L567)."
_TR_force(r, z, dr, dz, En, Lz, a) = (-2*a^2*dz^2*r*z)/(sqrt(1 - z^2)*(r^2 + a^2*z^2)^2) - (2*a^2*dr^2*z*sqrt(1 - z^2)*(r*(a^2 + r*(-3 + 2*r)) + a^2*(-1 + r)*z^2))/((a^2 + (-2 + r)*r)^2*(r^2 + a^2*z^2)^2) + (2*a^2*z*(-(Lz^2*r^5) - a^4*Lz^2*r*z^4 - En^2*r^6*(-1 + z^2) + a^4*(-(a*En) + Lz)^2*z^2*(-1 + z^2) + r^4*(2*(a*En - 2*Lz)*(a*En - Lz) + a*En*z^2*(-3*a*En + 6*Lz + a*En*z^2)) + a^2*r^2*((-(a*En) + Lz)^2 + (a*En + Lz)*z^2*(-3*a*En + 3*Lz + 2*a*En*z^2)) + 2*r^3*(-2*(-(a*En) + Lz)^2 - a*z^2*(-4*a*En^2 + 4*En*Lz + a*Lz^2 + 2*a*En^2*z^2))))/((a^2 + (-2 + r)*r)^2*sqrt(1 - z^2)*(r^2 + a^2*z^2)^2)
"`𝓣θ` (.wl L568-569)."
_Tth_force(r, z, dr, dz, d2z, chith, incl, En, Lz, a) = (a^2*dr^2*(a^2*z^2 + r^2*(-1 + 2*z^2))*sin(chith)*sin(incl))/((a^2 + (-2 + r)*r)*sqrt(1 - z^2)*(r^2 + a^2*z^2)^2) + ((3*a^6*Lz^2*z^6 + Lz^2*r^6*(1 + 2*z^2) + a^4*Lz^2*r^2*z^4*(7 + 2*z^2) + a^2*Lz^2*r^4*z^2*(5 + 4*z^2) + 2*a^4*r*z^2*(-(-(a*En) + Lz)^2 + 2*(-(a*En) + Lz)^2*z^2 - (a^2*En^2 - 2*a*En*Lz + 4*Lz^2)*z^4) + 2*a^2*r^3*((-(a*En) + Lz)^2 - (5*a*En - 3*Lz)*(a*En - Lz)*z^2 + (7*a^2*En^2 - 10*a*En*Lz - 4*Lz^2)*z^4 + a*En*(-3*a*En + 4*Lz)*z^6) - 2*r^5*(-(a^2*En^2) + Lz^2 + 2*(2*a^2*En^2 + Lz^2)*z^2 + a^2*En^2*z^4*(-5 + 2*z^2)))*sin(chith)*sin(incl))/((a^2 + (-2 + r)*r)*(1 - z^2)^(5/2)*(r^2 + a^2*z^2)^2) + (2*z*cos(incl)*cot(incl)*csc(chith)*(dz^2*z + d2z*(-z^2 + sin(incl)^2)))/((1 - z^2)^(3/2)*(-1 + 2*z^2 + cos(2*incl))) - (dz^2*csc(chith)*csc(incl)*(2*r^4*(1 + 2*z^2)*cos(incl)^2 + a^2*r^2*(-4*z^6 + 2*z^4*(7 + 2*cos(2*incl)) + z^2*(-5 + 3*cos(2*incl)) + 2*sin(incl)^2) + 2*a^4*z^2*(z^4*(1 + 2*cos(2*incl)) - sin(incl)^2 + 2*z^2*sin(incl)^2)))/(2*(1 - z^2)^(5/2)*(r^2 + a^2*z^2)^2)
"`𝓣ℨ` (.wl L570-571)."
_TZ_force(r, z, dr, dz, d2z, incl, En, Lz, a) = (a^2*dr^2*(r^2 - (a^2 + 2*r^2)*z^2))/((a^2 + (-2 + r)*r)*sqrt(1 - z^2)*(r^2 + a^2*z^2)^2) + (-3*a^6*Lz^2*z^6 - Lz^2*r^6*(1 + 2*z^2) - a^4*Lz^2*r^2*z^4*(7 + 2*z^2) - a^2*Lz^2*r^4*z^2*(5 + 4*z^2) + 2*a^4*r*z^2*((-(a*En) + Lz)^2 - 2*(-(a*En) + Lz)^2*z^2 + (a^2*En^2 - 2*a*En*Lz + 4*Lz^2)*z^4) + 2*a^2*r^3*(-(-(a*En) + Lz)^2 + (5*a*En - 3*Lz)*(a*En - Lz)*z^2 + (-7*a^2*En^2 + 10*a*En*Lz + 4*Lz^2)*z^4 + a*En*(3*a*En - 4*Lz)*z^6) + 2*r^5*(-(a^2*En^2) + Lz^2 + 2*(2*a^2*En^2 + Lz^2)*z^2 + a^2*En^2*z^4*(-5 + 2*z^2)))/((a^2 + (-2 + r)*r)*(1 - z^2)^(5/2)*(r^2 + a^2*z^2)^2) + (z*(-d2z + z*(-dz^2 + d2z*z)*csc(incl)^2))/((1 - z^2)^(3/2)*(1 - z^2*csc(incl)^2)) + (dz^2*(2*a^4*z^2*(-3*z^6 + z^4*(3 - 2*cos(2*incl)) + z^2*(-2 + cos(2*incl)) + sin(incl)^2) + 2*r^4*(-(z^2*(1 + z^2 + cos(2*incl))) + sin(incl)^2) - a^2*r^2*(4*z^6 + z^2*(-5 + 3*cos(2*incl)) + z^4*(6 + 4*cos(2*incl)) + 2*sin(incl)^2)))/((1 - z^2)^(5/2)*(r^2 + a^2*z^2)^2*(-1 + 2*z^2 + cos(2*incl)))
"`𝓤1θ` (.wl L573)."
_U1th_coef(r, z, dz, d2z, chith, En, Lz, a, Uth) = -((-2*d2z*(-1 + z^2)*(r^2 + a^2*z^2) + 2*dz^2*z*(r^2 + a^2*(-1 + 2*z^2)))/(Uth*(1 - z^2)^(3/2)*(r^2 + a^2*z^2)))

# --- Mino-time → Darwin-angle reparametrisation χr(λr), χθ(λθ) (SpinningTrajectory.wl L524-525):
#     the angle the χ-dependent forces are evaluated at on the 2-torus grid. R3. ---
"`χrhat(λr)` — radial Darwin angle from the geodesic r(λr) (.wl L524). R3."
function _chir_hat(a, p, e, x, λr)
    Υr = kerr_frequencies(a, p, e, x; Time="Mino").Ur; Λr = 2pi / Υr
    r = geodesic_trajectory_radial(a, p, e, x)(Υr * λr).r
    base = acos(clamp((p / r - 1) / e, -1.0, 1.0))
    fl = 2pi * floor(λr / Λr)
    return mod(λr, Λr) < Λr / 2 ? base + fl : 2pi + fl - base
end
"`χθhat(λθ)` — polar Darwin angle from the geodesic θ(λθ) (.wl L525). R3."
function _chith_hat(a, p, e, x, λθ)
    Υθ = kerr_frequencies(a, p, e, x; Time="Mino").Uth; Λθ = 2pi / Υθ
    θ = geodesic_trajectory_polar(a, p, e, x)(Υθ * λθ).θ
    base = acos(clamp(cos(θ) / sin(acos(x)), -1.0, 1.0))
    fl = 2pi * floor(λθ / Λθ)
    return mod(λθ, Λθ) < Λθ / 2 ? base + fl : 2pi + fl - base
end

"""
    _numeric_force_dft(a, p, e, x, nmax, kmax) -> (Rtfour, Rphifour)

The 2-torus DFT machinery of the numeric correction engine (SpinningTrajectory.wl L592-600): build
the `𝓡t`/`𝓡φ` spin-curvature force on the `stepsr×stepsθ=(4nmax)×(4kmax)` midpoint Mino grid, then
the 2D DFT `Exp_{ni}·list·Exp_{jk}/(stepsr·stepsθ)` → the `(n,k)` Fourier-coefficient matrices. The
same machinery DFTs every force in the engine. R3.
"""
function _numeric_force_dft(a, p, e, x, nmax::Integer, kmax::Integer)
    En = geodesic_energy(a, p, e, x); Lz = geodesic_angular_momentum(a, p, e, x)
    Q = carter_constant(a, p, e, x); K = Q + (Lz - a * En)^2
    fr = kerr_frequencies(a, p, e, x; Time="Mino"); Υr = fr.Ur; Υθ = fr.Uth
    Λr = 2pi / Υr; Λθ = 2pi / Υθ
    stepsr = 4 * nmax; stepsθ = 4 * kmax
    rad = geodesic_trajectory_radial(a, p, e, x); pol = geodesic_trajectory_polar(a, p, e, x)
    # geodesic grid samples (midpoint rule): rhat,drdλ at λr_i; zhat,dzdλ at λθ_j
    rg = [rad(Υr * (i - 0.5) * Λr / stepsr) for i in 1:stepsr]
    zg = [pol(Υθ * (j - 0.5) * Λθ / stepsθ) for j in 1:stepsθ]
    Rt = [_Rt_force(rg[i].r, zg[j].z, rg[i].Ur, zg[j].Uz, K, a) for i in 1:stepsr, j in 1:stepsθ]
    Rφ = [_Rphi_force(rg[i].r, zg[j].z, rg[i].Ur, zg[j].Uz, K, a) for i in 1:stepsr, j in 1:stepsθ]
    Eni = [cis(2pi * n * (i - 0.5) / stepsr) for n in -nmax:nmax, i in 1:stepsr]   # (2nmax+1)×stepsr
    Ejk = [cis(2pi * k * (j - 0.5) / stepsθ) for j in 1:stepsθ, k in -kmax:kmax]   # stepsθ×(2kmax+1)
    return (Eni * Rt * Ejk / (stepsr * stepsθ), Eni * Rφ * Ejk / (stepsr * stepsθ))
end

"""
    _delta_u_st(four, Υr, Υθ, nmax, kmax, λr, λθ) -> Real

The spin-shift potential δu_{ds,t/φ}(λr,λθ) (SpinningTrajectory.wl L526-527): the periodic solution
of the t/φ spin-curvature equation, `Re Σ_{n,k≠0,even k} Im[𝓡four_{nk}]/(−kΥθ−nΥr)·e^{−i(kΥθλθ+nΥrλr)}`
from the force Fourier coefficients `four` (= Rtfour or Rφfour). Aligned spin ⇒ even-k sum. R3.
"""
function _delta_u_st(four, Υr, Υθ, nmax::Integer, kmax::Integer, λr, λθ)
    s = zero(ComplexF64)
    for n in -nmax:nmax, k in -kmax:2:kmax
        (n == 0 && k == 0) && continue
        s += imag(four[n + nmax + 1, k + kmax + 1]) / (-k * Υθ - n * Υr) *
             cis(-k * Υθ * λθ - n * Υr * λr)
    end
    return real(s)
end

# =============================================================================
# Final stage: the M-matrix least-squares assembly (SpinningTrajectory.wl
# `Correction[orbitGeo,nmax,kmax]`, L596-745). Every spin-curvature force is
# DFT'd over the (λr,λθ) 2-torus, the convolution-structured submatrices M1/M2/M3
# (r/θ/norm equations) are built from those Fourier coefficients, and the linear
# system M·v = c is solved (least squares) for the frequency/four-velocity
# corrections ΥrS,ΥθS,udt0S,udφ0S. R3 (load-bearing numerics).
# =============================================================================

"""
    _numeric_correction(a, p, e, x, nmax, kmax) -> NamedTuple

Port of `Correction[orbitGeo,nmax,kmax]` (SpinningTrajectory.wl L596-745): the numeric
linear-in-spin orbit-correction engine. Builds the Fourier-coefficient matrix of every
spin-curvature force on the `(4nmax)×(4kmax)` midpoint Mino grid, assembles the M-matrix
least-squares system, and solves for the linear-in-spin Mino-frequency corrections and the
four-velocity zero-modes.

Returns `(; ΥrS, ΥθS, udt0S, udφ0S, error, nmax, kmax)` where `error = ‖M·v − c‖` (the
least-squares residual). The δχ/δ𝔯/δ𝔷 mode coefficients are also returned (δχrSIm, δχθSIm,
δ𝔯four, δ𝔷four) for completeness. R3.
"""
function _numeric_correction(a, p, e, x, nmax::Integer, kmax::Integer)
    En = geodesic_energy(a, p, e, x); Lz = geodesic_angular_momentum(a, p, e, x)
    Q = carter_constant(a, p, e, x); K = Q + (Lz - a * En)^2
    incl = acos(x)
    fr = kerr_frequencies(a, p, e, x; Time="Mino"); Υr = fr.Ur; Υθ = fr.Uth
    Λr = 2pi / Υr; Λθ = 2pi / Υθ
    stepsr = 4 * nmax; stepsθ = 4 * kmax
    rad = geodesic_trajectory_radial(a, p, e, x); pol = geodesic_trajectory_polar(a, p, e, x)

    # --- geodesic grid samples on the midpoint Mino grid (mirror _numeric_force_dft) ---
    # λr_i = (i-1/2)Λr/stepsr  (i=1..stepsr);  λθ_j = (j-1/2)Λθ/stepsθ  (j=1..stepsθ)
    λrg = [(i - 0.5) * Λr / stepsr for i in 1:stepsr]
    λθg = [(j - 0.5) * Λθ / stepsθ for j in 1:stepsθ]
    rg = [rad(Υr * λrg[i]) for i in 1:stepsr]   # .r, .Ur (=drdλ)
    zg = [pol(Υθ * λθg[j]) for j in 1:stepsθ]   # .z, .Uz (=dzdλ)
    rv = [g.r for g in rg]; drv = [g.Ur for g in rg]
    zv = [g.z for g in zg]; dzv = [g.Uz for g in zg]
    d2rv = [_d2r_dλ2(rv[i], En, Lz, K, a) for i in 1:stepsr]
    d2zv = [_d2z_dλ2(zv[j], En, Lz, Q, a) for j in 1:stepsθ]
    χrv = [_chir_hat(a, p, e, x, λrg[i]) for i in 1:stepsr]
    χθv = [_chith_hat(a, p, e, x, λθg[j]) for j in 1:stepsθ]

    # --- DFT matrices and helper (mirror _numeric_force_dft exactly) ---
    Eni = [cis(2pi * n * (i - 0.5) / stepsr) for n in -nmax:nmax, i in 1:stepsr]
    Ejk = [cis(2pi * k * (j - 0.5) / stepsθ) for j in 1:stepsθ, k in -kmax:kmax]
    dft(list) = Eni * list * Ejk / (stepsr * stepsθ)   # (2nmax+1)×(2kmax+1), indexed [n+nmax+1,k+kmax+1]

    # ---------------- t,φ equations: source for δu_{ds,t/φ} (L597-600) ----------------
    Rtlist = [_Rt_force(rv[i], zv[j], drv[i], dzv[j], K, a) for i in 1:stepsr, j in 1:stepsθ]
    Rφlist = [_Rphi_force(rv[i], zv[j], drv[i], dzv[j], K, a) for i in 1:stepsr, j in 1:stepsθ]
    Rtfour = dft(Rtlist); Rφfour = dft(Rφlist)

    # δu_{ds,t}(λr,λθ), δu_{ds,φ}(λr,λθ) on the grid (L526-527)
    δudst = [_delta_u_st(Rtfour, Υr, Υθ, nmax, kmax, λrg[i], λθg[j]) for i in 1:stepsr, j in 1:stepsθ]
    δudsφ = [_delta_u_st(Rφfour, Υr, Υθ, nmax, kmax, λrg[i], λθg[j]) for i in 1:stepsr, j in 1:stepsθ]

    # ---------------- r equation forces (L602-616) ----------------
    Frlist = [_Fr_force(χrv[i], e, p)                                                   for i in 1:stepsr, j in 1:stepsθ]
    F𝔯list = [1.0                                                                       for i in 1:stepsr, j in 1:stepsθ]  # 𝓕𝔯func[]:=1 (.wl L531)
    Grlist = [_Gr_force(rv[i], zv[j], drv[i], χrv[i], a, e, p)                          for i in 1:stepsr, j in 1:stepsθ]
    G𝔯list = [_GR_force(rv[i], zv[j], drv[i], a)                                        for i in 1:stepsr, j in 1:stepsθ]
    Gθlist = [_Gth_force(rv[i], zv[j], dzv[j], χθv[j], incl, a)                         for i in 1:stepsr, j in 1:stepsθ]
    G𝔷list = [_GZ_force(rv[i], zv[j], dzv[j], a)                                        for i in 1:stepsr, j in 1:stepsθ]
    Hrlist = [_Hr_force(rv[i], zv[j], drv[i], dzv[j], d2rv[i], χrv[i], En, Lz, a, e, p) for i in 1:stepsr, j in 1:stepsθ]
    H𝔯list = [_HR_force(rv[i], zv[j], drv[i], dzv[j], En, Lz, a)                        for i in 1:stepsr, j in 1:stepsθ]
    Hθlist = [_Hth_force(rv[i], zv[j], drv[i], dzv[j], χθv[j], incl, En, Lz, a)         for i in 1:stepsr, j in 1:stepsθ]
    H𝔷list = [_HZ_force(rv[i], zv[j], drv[i], dzv[j], En, Lz, a)                        for i in 1:stepsr, j in 1:stepsθ]
    I1rlist = [_I1r_coef(rv[i], zv[j], drv[i], d2rv[i], Υr, a)                          for i in 1:stepsr, j in 1:stepsθ]
    I1θlist = [_I1θ_coef(rv[i], zv[j], dzv[j], Υθ, a)                                   for i in 1:stepsr, j in 1:stepsθ]
    I2list = [_I2_coef(rv[i], zv[j], En, Lz, K, a)                                      for i in 1:stepsr, j in 1:stepsθ]
    I3list = [_I3_coef(rv[i], zv[j], En, Lz, K, a)                                      for i in 1:stepsr, j in 1:stepsθ]
    Jlist = [_J_source(rv[i], zv[j], δudst[i, j], δudsφ[i, j], En, Lz, K, a)            for i in 1:stepsr, j in 1:stepsθ]

    Frfour = dft(Frlist); F𝔯four = dft(F𝔯list)
    Grfour = dft(Grlist); G𝔯four = dft(G𝔯list); Gθfour = dft(Gθlist); G𝔷four = dft(G𝔷list)
    Hrfour = dft(Hrlist); H𝔯four = dft(H𝔯list); Hθfour = dft(Hθlist); H𝔷four = dft(H𝔷list)
    I1rfour = dft(I1rlist); I1θfour = dft(I1θlist); I2four = dft(I2list); I3four = dft(I3list)
    Jfour = dft(Jlist)

    # ---------------- θ equation forces (L633-647) ----------------
    Qθlist = [_Qth_force(zv[j], χθv[j], incl)                                           for i in 1:stepsr, j in 1:stepsθ]
    Q𝔷list = [_QZ_force(zv[j])                                                          for i in 1:stepsr, j in 1:stepsθ]
    Srlist = [_Sr_force(rv[i], zv[j], drv[i], χrv[i], a, e, p)                          for i in 1:stepsr, j in 1:stepsθ]
    S𝔯list = [_SR_force(rv[i], zv[j], drv[i], a)                                        for i in 1:stepsr, j in 1:stepsθ]
    Sθlist = [_Sth_force(rv[i], zv[j], dzv[j], χθv[j], incl, a)                         for i in 1:stepsr, j in 1:stepsθ]
    S𝔷list = [_SZ_force(rv[i], zv[j], dzv[j], a)                                        for i in 1:stepsr, j in 1:stepsθ]
    Trlist = [_Tr_force(rv[i], zv[j], drv[i], dzv[j], χrv[i], En, Lz, a, e, p)          for i in 1:stepsr, j in 1:stepsθ]
    T𝔯list = [_TR_force(rv[i], zv[j], drv[i], dzv[j], En, Lz, a)                        for i in 1:stepsr, j in 1:stepsθ]
    Tθlist = [_Tth_force(rv[i], zv[j], drv[i], dzv[j], d2zv[j], χθv[j], incl, En, Lz, a) for i in 1:stepsr, j in 1:stepsθ]
    T𝔷list = [_TZ_force(rv[i], zv[j], drv[i], dzv[j], d2zv[j], incl, En, Lz, a)         for i in 1:stepsr, j in 1:stepsθ]
    U1rlist = [_U1r_coef(rv[i], zv[j], drv[i], Υr, a)                                   for i in 1:stepsr, j in 1:stepsθ]
    U1θlist = [_U1th_coef(rv[i], zv[j], dzv[j], d2zv[j], χθv[j], En, Lz, a, Υθ)         for i in 1:stepsr, j in 1:stepsθ]
    U2list = [_U2_coef(rv[i], zv[j], En, Lz, a)                                         for i in 1:stepsr, j in 1:stepsθ]
    U3list = [_U3_coef(rv[i], zv[j], En, Lz, a)                                         for i in 1:stepsr, j in 1:stepsθ]
    Vlist = [_V_source(rv[i], zv[j], δudst[i, j], δudsφ[i, j], En, Lz, K, a)            for i in 1:stepsr, j in 1:stepsθ]

    Qθfour = dft(Qθlist); Q𝔷four = dft(Q𝔷list)
    Srfour = dft(Srlist); S𝔯four = dft(S𝔯list); Sθfour = dft(Sθlist); S𝔷four = dft(S𝔷list)
    Trfour = dft(Trlist); T𝔯four = dft(T𝔯list); Tθfour = dft(Tθlist); T𝔷four = dft(T𝔷list)
    U1rfour = dft(U1rlist); U1θfour = dft(U1θlist); U2four = dft(U2list); U3four = dft(U3list)
    Vfour = dft(Vlist)

    # ---------------- norm equation forces (L664-676) ----------------
    Krlist = [_Kr_force(rv[i], zv[j], drv[i], χrv[i], a, e, p)                          for i in 1:stepsr, j in 1:stepsθ]
    K𝔯list = [_KR_force(rv[i], zv[j], drv[i], a)                                        for i in 1:stepsr, j in 1:stepsθ]
    Kθlist = [_Kth_force(rv[i], zv[j], dzv[j], χθv[j], incl, a)                         for i in 1:stepsr, j in 1:stepsθ]
    K𝔷list = [_KZ_force(rv[i], zv[j], dzv[j], a)                                        for i in 1:stepsr, j in 1:stepsθ]
    Mrlist = [_Mr_force(rv[i], zv[j], drv[i], dzv[j], χrv[i], En, Lz, a, e, p)          for i in 1:stepsr, j in 1:stepsθ]
    M𝔯list = [_MR_force(rv[i], zv[j], drv[i], dzv[j], En, Lz, a)                        for i in 1:stepsr, j in 1:stepsθ]
    Mθlist = [_Mth_force(rv[i], zv[j], drv[i], dzv[j], χθv[j], incl, En, Lz, a)         for i in 1:stepsr, j in 1:stepsθ]
    M𝔷list = [_MZ_force(rv[i], zv[j], drv[i], dzv[j], En, Lz, a)                        for i in 1:stepsr, j in 1:stepsθ]
    N1rlist = [_N1r_coef(rv[i], zv[j], drv[i], Υr, a)                                   for i in 1:stepsr, j in 1:stepsθ]
    N1θlist = [_N1th_coef(rv[i], zv[j], dzv[j], Υθ, a)                                  for i in 1:stepsr, j in 1:stepsθ]
    N2list = [_N2_coef(rv[i], zv[j], En, Lz, a)                                         for i in 1:stepsr, j in 1:stepsθ]
    N3list = [_N3_coef(rv[i], zv[j], En, Lz, a)                                         for i in 1:stepsr, j in 1:stepsθ]
    Plist = [_P_source(rv[i], zv[j], δudst[i, j], δudsφ[i, j], En, Lz, a)               for i in 1:stepsr, j in 1:stepsθ]

    Krfour = dft(Krlist); K𝔯four = dft(K𝔯list); Kθfour = dft(Kθlist); K𝔷four = dft(K𝔷list)
    Mrfour = dft(Mrlist); M𝔯four = dft(M𝔯list); Mθfour = dft(Mθlist); M𝔷four = dft(M𝔷list)
    N1rfour = dft(N1rlist); N1θfour = dft(N1θlist); N2four = dft(N2list); N3four = dft(N3list)
    Pfour = dft(Plist)

    # ============================ submatrix assembly ============================
    # band-conditioned convolution submatrices, exactly mirroring L691-717.
    # Rows iterate (n,k) with n outer (−nmax..nmax) and k inner (the equation's k-parity).
    # Columns of the convolution blocks iterate the (nn) or (nn,kk) sum indices.

    keven = collect(-kmax:2:kmax)            # k modes for r/norm equations  (kmax+1 values)
    kodd  = collect(-kmax+1:2:kmax-1)        # k modes for θ equation        (kmax-1 values)
    nn1   = vcat(-nmax:-1, 1:nmax)           # nn over ±1..±nmax             (2nmax values)
    nnall = collect(-nmax:nmax)              # nn over −nmax..nmax           (2nmax+1 values)
    kk2   = vcat(-kmax:2:-2, 2:2:kmax)       # kk even ±2..±kmax             (kmax values)
    kkodd = collect(-kmax+1:2:kmax-1)        # kk odd −kmax+1..kmax−1        (kmax values)

    idx(m, n, k) = m[n + nmax + 1, k + kmax + 1]   # Fourier-matrix [n,k] accessor

    # ---- block builders. nrows over (n,krange); columns over the listed sum index ----
    # single-inner-index block (Flatten[...,{{1,2},{3}}]): cols = nn  OR  kk
    function block3_nn(f, krange, nns)        # f(n,k,nn) with band on n
        rows = [(n, k) for n in -nmax:nmax for k in krange]
        [f(nk[1], nk[2], nn) for nk in rows, nn in nns]
    end
    function block3_kk(f, krange, kks)        # f(n,k,kk) with band on k
        rows = [(n, k) for n in -nmax:nmax for k in krange]
        [f(nk[1], nk[2], kk) for nk in rows, kk in kks]
    end
    # double-inner-index block (Flatten[...,{{1,2},{3,4}}]): cols = (nn,kk), nn outer kk inner
    function block4(f, krange, nns, kks)      # f(n,k,nn,kk)
        rows = [(n, k) for n in -nmax:nmax for k in krange]
        cols = [(nn, kk) for nn in nns for kk in kks]
        [f(nk[1], nk[2], c[1], c[2]) for nk in rows, c in cols]
    end
    band_n(n, nn) = -nmax <= n - nn <= nmax
    band_k(k, kk) = -kmax <= k - kk <= kmax

    # column-vector blocks (Flatten[mat[[All, sel]], {{1,2}}]) → length-(rows) vector,
    # n outer, k inner over the selected k columns. Returns a column vector matching the
    # (n,krange) row ordering of the convolution blocks above.
    function colvec(mat, krange)
        [idx(mat, n, k) for n in -nmax:nmax for k in krange]
    end

    # ---------------- M1 : r equation (k even) ----------------
    M1r = block3_nn((n, k, nn) -> band_n(n, nn) ?
            -nn^2 * Υr^2 * idx(Frfour, n - nn, k) - im * nn * Υr * idx(Grfour, n - nn, k) + idx(Hrfour, n - nn, k) : 0.0,
            keven, nn1)
    M1θ = block3_kk((n, k, kk) -> band_k(k, kk) ?
            -im * kk * Υθ * idx(Gθfour, n, k - kk) + idx(Hθfour, n, k - kk) : 0.0,
            keven, kk2)
    M1𝔯 = block4((n, k, nn, kk) -> (band_n(n, nn) && band_k(k, kk)) ?
            (-im * nn * Υr - im * kk * Υθ)^2 * idx(F𝔯four, n - nn, k - kk) -
            im * (nn * Υr + kk * Υθ) * idx(G𝔯four, n - nn, k - kk) + idx(H𝔯four, n - nn, k - kk) : 0.0,
            keven, nnall, kk2)
    M1𝔷 = block4((n, k, nn, kk) -> (band_n(n, nn) && band_k(k, kk)) ?
            -im * (nn * Υr + kk * Υθ) * idx(G𝔷four, n - nn, k - kk) + idx(H𝔷four, n - nn, k - kk) : 0.0,
            keven, nn1, kkodd)
    M1Υr = colvec(I1rfour, keven); M1Υθ = colvec(I1θfour, keven)
    M1t = colvec(I2four, keven);   M1φ = colvec(I3four, keven)
    M1 = hcat(-imag(M1r), -imag(M1θ), real(M1𝔯), real(M1𝔷),
              real(M1Υr), real(M1Υθ), real(M1t), real(M1φ))

    # ---------------- M2 : θ equation (k odd) ----------------
    M2r = block3_nn((n, k, nn) -> band_n(n, nn) ?
            -im * nn * Υr * idx(Srfour, n - nn, k) + idx(Trfour, n - nn, k) : 0.0,
            kodd, nn1)
    M2θ = block3_kk((n, k, kk) -> band_k(k, kk) ?
            -kk^2 * Υθ^2 * idx(Qθfour, n, k - kk) - im * kk * Υθ * idx(Sθfour, n, k - kk) + idx(Tθfour, n, k - kk) : 0.0,
            kodd, kk2)
    M2𝔯 = block4((n, k, nn, kk) -> (band_n(n, nn) && band_k(k, kk)) ?
            -im * (nn * Υr + kk * Υθ) * idx(S𝔯four, n - nn, k - kk) + idx(T𝔯four, n - nn, k - kk) : 0.0,
            kodd, nnall, kk2)
    M2𝔷 = block4((n, k, nn, kk) -> (band_n(n, nn) && band_k(k, kk)) ?
            (-im * nn * Υr - im * kk * Υθ)^2 * idx(Q𝔷four, n - nn, k - kk) -
            im * (nn * Υr + kk * Υθ) * idx(S𝔷four, n - nn, k - kk) + idx(T𝔷four, n - nn, k - kk) : 0.0,
            kodd, nn1, kkodd)
    M2Υr = colvec(U1rfour, kodd); M2Υθ = colvec(U1θfour, kodd)
    M2t = colvec(U2four, kodd);   M2φ = colvec(U3four, kodd)
    M2 = hcat(-imag(M2r), -imag(M2θ), real(M2𝔯), real(M2𝔷),
              real(M2Υr), real(M2Υθ), real(M2t), real(M2φ))

    # ---------------- M3 : norm equation (k even) ----------------
    M3r = block3_nn((n, k, nn) -> band_n(n, nn) ?
            -im * nn * Υr * idx(Krfour, n - nn, k) + idx(Mrfour, n - nn, k) : 0.0,
            keven, nn1)
    M3θ = block3_kk((n, k, kk) -> band_k(k, kk) ?
            -im * kk * Υθ * idx(Kθfour, n, k - kk) + idx(Mθfour, n, k - kk) : 0.0,
            keven, kk2)
    M3𝔯 = block4((n, k, nn, kk) -> (band_n(n, nn) && band_k(k, kk)) ?
            -im * (nn * Υr + kk * Υθ) * idx(K𝔯four, n - nn, k - kk) + idx(M𝔯four, n - nn, k - kk) : 0.0,
            keven, nnall, kk2)
    M3𝔷 = block4((n, k, nn, kk) -> (band_n(n, nn) && band_k(k, kk)) ?
            -im * (nn * Υr + kk * Υθ) * idx(K𝔷four, n - nn, k - kk) + idx(M𝔷four, n - nn, k - kk) : 0.0,
            keven, nn1, kkodd)
    M3Υr = colvec(N1rfour, keven); M3Υθ = colvec(N1θfour, keven)
    M3t = colvec(N2four, keven);   M3φ = colvec(N3four, keven)
    M3 = hcat(-imag(M3r), -imag(M3θ), real(M3𝔯), real(M3𝔷),
              real(M3Υr), real(M3Υθ), real(M3t), real(M3φ))

    M = vcat(M1, M2, M3)

    # ---------------- source vector c (L719) ----------------
    cJ = [idx(Jfour, n, k) for n in -nmax:nmax for k in keven]   # even k
    cV = [idx(Vfour, n, k) for n in -nmax:nmax for k in kodd]    # odd  k
    cP = [idx(Pfour, n, k) for n in -nmax:nmax for k in keven]   # even k
    c = real(-vcat(cJ, cV, cP))

    # ---------------- least squares (Mathematica LeastSquares = min-norm solution) ----------------
    v = pinv(M) * c

    nu = 2 * nmax + kmax + kmax * (2 * nmax + 1) + 2 * nmax * kmax + 4   # # unknowns
    δχrSIm = v[1:2*nmax]
    δχθSIm = v[2*nmax+1 : 2*nmax+kmax]
    δ𝔯four = v[2*nmax+kmax+1 : 2*nmax+kmax+kmax*(2*nmax+1)]
    δ𝔷four = v[2*nmax+kmax+kmax*(2*nmax+1)+1 : 2*nmax+kmax+kmax*(2*nmax+1)+2*nmax*kmax]
    ΥrS = v[end-3]; ΥθS = v[end-2]; udt0S = v[end-1]; udφ0S = v[end]
    err = norm(M * v - c)

    return (; ΥrS, ΥθS, udt0S, udφ0S, error = err, nmax, kmax,
            δχrSIm, δχθSIm, δ𝔯four, δ𝔷four)
end
