# Closed-form NP/GHP source kernels Sourceace{l,n,m,m†}2± and Sourcenorm2±
# (Sourceincoordinates.nb "Analytical formulas", Polcar & Witzany 2025). Verbatim port of the
# saved notebook output cells (xAct derivation not re-run): equatorial Schwarzschild θ=π/2, r=R,
# with the higher RadHZ derivatives already eliminated via the radial Teukolsky equation, so each
# kernel is a rational function of (MB,R,ω,mi,li), the Hertz potential value/first-derivative
# (`rhz`,`drhz` for the +ω mode; `rhzb`,`drhzb` for the −ω/conjugate mode), and the spin-weighted
# spherical harmonics `yv(s,m) = {}_sY_{li,m}(π/2)`.
#
# `yv` is any callable `(s,m) -> {}_sY_{li,m}(π/2)`, supplied by the caller (SourceAmplitudes wires
# in Harmonics) — these kernels therefore carry no Harmonics dependency. The +kernels read
# `yv(s,mi)`; the −kernels read `yv(s,-mi)`. Snorm2± are the LHS-normalisation kernels (continuity
# check); the source-amplitude pipeline itself sets Snorm≡0 (see source_amplitudes / upstream
# Sourcesamplitudes). All arithmetic is type-generic so BigFloat inputs are preserved.

function _source_kernel_l_plus(MB, R, omega, mi, li, rhz, drhz, yv)
    return (rhz*(2*im*li*(1 + li)*sqrt(-2 + li + li^2)*sqrt(MB)*(2*MB - R)*R^(3/2)*yv(-1, mi) + sqrt(li*(1 + li))*sqrt(-2 + li + li^2)*R^3*(2 - 4*im*MB*omega + 3*im*omega*R)*yv(0, mi) + 8*im*sqrt(-2 + li + li^2)*MB^(3/2)*R^(3/2)*yv(1, mi) - 12*sqrt(-2 + li + li^2)*MB^(3/2)*omega*R^(5/2)*yv(1, mi) - 4*sqrt(-2 + li + li^2)*sqrt(MB)*omega*R^(7/2)*yv(1, mi) - 4*im*sqrt(-2 + li + li^2)*sqrt(MB)*omega^2*R^(9/2)*yv(1, mi) + 16*MB^3*yv(2, mi) - 24*MB^2*R*yv(2, mi) + 4*li*MB^2*R*yv(2, mi) + 4*li^2*MB^2*R*yv(2, mi) - 8*im*MB^2*omega*R^2*yv(2, mi) - 4*im*MB*omega*R^3*yv(2, mi) + im*li*MB*omega*R^3*yv(2, mi) + im*li^2*MB*omega*R^3*yv(2, mi) + 2*MB*omega^2*R^4*yv(2, mi)) - drhz*(2*MB - R)*R*(sqrt(li*(1 + li))*sqrt(-2 + li + li^2)*R^2*yv(0, mi) + 4*im*sqrt(-2 + li + li^2)*sqrt(MB)*sqrt(R)*(MB + im*omega*R^2)*yv(1, mi) + MB*(-4*MB + R*(-2 + li + li^2 - 2*im*omega*R))*yv(2, mi)))/(16*R*(-3*MB + R))
end

function _source_kernel_l_minus(MB, R, omega, mi, li, rhzb, drhzb, yv)
    return ((-1)^mi*(drhzb*(2*MB - R)*R*(MB*(4*MB + R*(2 - li - li^2 - 2*im*omega*R))*yv(-2, -mi) - sqrt(-2 + li + li^2)*sqrt(R)*(4*sqrt(MB)*(im*MB + omega*R^2)*yv(-1, -mi) + sqrt(li*(1 + li))*R^(3/2)*yv(0, -mi))) + rhzb*(MB*(16*MB^2 + 4*MB*R*(-6 + li + li^2 + 2*im*omega*R) - im*omega*R^3*(-4 + li + li^2 + 2*im*omega*R))*yv(-2, -mi) + sqrt(-2 + li + li^2)*R^(3/2)*(4*sqrt(MB)*(MB*(2*im + 3*omega*R) + omega*R^2*(1 - im*omega*R))*yv(-1, -mi) + sqrt(li*(1 + li))*R^(3/2)*(2 + 4*im*MB*omega - 3*im*omega*R)*yv(0, -mi) + 2*im*li*(1 + li)*sqrt(MB)*(2*MB - R)*yv(1, -mi)))))/(16*R*(-3*MB + R))
end

function _source_kernel_n_plus(MB, R, omega, mi, li, rhz, drhz, yv)
    return (drhz*(2*MB - R)*R*(sqrt(li*(1 + li))*sqrt(-2 + li + li^2)*R*(-2*MB + R)^2*yv(0, mi) + 4*sqrt(-2 + li + li^2)*sqrt(MB)*(2*MB - R)*sqrt(R)*(-(im*MB) + omega*R^2)*yv(1, mi) + MB*(8*MB^2 + 2*(-4 + li + li^2)*MB*R + R^2*(2 - li - li^2 - 2*im*omega*R + 4*omega^2*R^2))*yv(2, mi)) + rhz*(-(im*sqrt(li*(1 + li))*sqrt(-2 + li + li^2)*R^2*(-2*MB + R)^2*(-2*im + omega*R)*yv(0, mi)) - 2*im*sqrt(-2 + li + li^2)*sqrt(MB)*R^(3/2)*(4*MB^2*(-2 + li + li^2 - 3*im*omega*R) + 2*MB*R*(2 - 2*li - 2*li^2 + im*omega*R + 2*omega^2*R^2) + R^2*(li + li^2 - 2*omega*R*(-im + omega*R)))*yv(1, mi) + MB*(32*MB^3 + 8*MB^2*R*(-8 + li + li^2 - 4*im*omega*R) + R^3*(-8 - 4*im*omega*R - 6*omega^2*R^2 - 4*im*omega^3*R^3 + li*(4 + 3*im*omega*R) + li^2*(4 + 3*im*omega*R)) - 2*MB*R^2*(li*(6 + 3*im*omega*R) + li^2*(6 + 3*im*omega*R) + 4*(-5 - 4*im*omega*R + omega^2*R^2)))*yv(2, mi)))/(32*R^2*(-3*MB + R))
end

function _source_kernel_n_minus(MB, R, omega, mi, li, rhzb, drhzb, yv)
    return ((-1)^mi*(drhzb*(2*MB - R)*R*(MB*(8*MB^2 + 2*(-4 + li + li^2)*MB*R + R^2*(2 - li - li^2 + 2*im*omega*R + 4*omega^2*R^2))*yv(-2, -mi) + sqrt(-2 + li + li^2)*sqrt(R)*(-2*MB + R)*(4*sqrt(MB)*(im*MB + omega*R^2)*yv(-1, -mi) + sqrt(li*(1 + li))*sqrt(R)*(-2*MB + R)*yv(0, -mi))) + rhzb*(MB*(32*MB^3 + 8*MB^2*R*(-8 + li + li^2 + 4*im*omega*R) + R^3*(-8 + 4*im*omega*R - 6*omega^2*R^2 + 4*im*omega^3*R^3 + li*(4 - 3*im*omega*R) + li^2*(4 - 3*im*omega*R)) - 2*MB*R^2*(li*(6 - 3*im*omega*R) + li^2*(6 - 3*im*omega*R) + 4*(-5 + 4*im*omega*R + omega^2*R^2)))*yv(-2, -mi) + im*sqrt(-2 + li + li^2)*R^(3/2)*(-2*MB + R)*(2*sqrt(MB)*(2*MB*(-2 + li + li^2 + 3*im*omega*R) + R*(-li - li^2 + 2*omega*R*(im + omega*R)))*yv(-1, -mi) + sqrt(li*(1 + li))*sqrt(R)*(-2*MB + R)*(2*im + omega*R)*yv(0, -mi)))))/(32*R^2*(-3*MB + R))
end

function _source_kernel_m_plus(MB, R, omega, mi, li, rhz, drhz, yv)
    return (2*drhz*(2*MB - R)*R*(im*sqrt(li*(1 + li))*sqrt(-2 + li + li^2)*sqrt(MB)*(2*MB - R)*sqrt(R)*yv(0, mi) + 2*sqrt(-2 + li + li^2)*(MB^2 - im*MB*omega*R^2 + im*omega*R^3)*yv(1, mi) + sqrt(MB)*(MB + im*omega*R^2)*(-4*omega*R^(3/2)*yv(2, mi) + sqrt(-6 + li + li^2)*sqrt(MB)*yv(3, mi))) + rhz*(2*sqrt(li*(1 + li))*sqrt(-2 + li + li^2)*sqrt(MB)*(2*MB - R)*R^(3/2)*(-2*im + omega*R)*yv(0, mi) + sqrt(-2 + li + li^2)*(16*MB^3 + 12*MB^2*R*(-2 - im*omega*R) + 2*MB*R^2*(2 + li + li^2 + 6*im*omega*R - 2*omega^2*R^2) + R^3*(-li - li^2 + 4*omega*R*(-2*im + omega*R)))*yv(1, mi) + sqrt(MB)*(8*MB^2 + 2*MB*R*(-6 + li + li^2 - 3*im*omega*R) + R^2*(2 - li - li^2 - 2*im*omega*R + 2*omega^2*R^2))*(-4*omega*R^(3/2)*yv(2, mi) + sqrt(-6 + li + li^2)*sqrt(MB)*yv(3, mi))))/(16*sqrt(2)*R*(-3*MB + R))
end

function _source_kernel_m_minus(MB, R, omega, mi, li, rhzb, drhzb, yv)
    return -1/16*((-1)^mi*sqrt(-2 + li + li^2)*(2*drhzb*sqrt(MB)*(2*MB - R)*R*(sqrt(MB)*(MB - im*omega*R^2)*yv(-1, -mi) - im*sqrt(li*(1 + li))*sqrt(R)*(-2*MB + R)*yv(0, -mi)) + rhzb*(MB*(8*MB^2 + 2*MB*R*(-6 + li + li^2 + 3*im*omega*R) + R^2*(2 - li - li^2 + 2*im*omega*R + 2*omega^2*R^2))*yv(-1, -mi) + R*(-2*MB + R)*(2*sqrt(li*(1 + li))*sqrt(MB)*sqrt(R)*(2*im + omega*R)*yv(0, -mi) - li*(1 + li)*(2*MB - R)*yv(1, -mi)))))/(sqrt(2)*R*(-3*MB + R))
end

function _source_kernel_mdag_plus(MB, R, omega, mi, li, rhz, drhz, yv)
    return (sqrt(-2 + li + li^2)*(2*drhz*sqrt(MB)*(2*MB - R)*R*(-(im*sqrt(li*(1 + li))*sqrt(R)*(-2*MB + R)*yv(0, mi)) + sqrt(MB)*(MB + im*omega*R^2)*yv(1, mi)) + rhz*(li*(1 + li)*R*(-2*MB + R)^2*yv(-1, mi) + 2*sqrt(li*(1 + li))*sqrt(MB)*(2*MB - R)*R^(3/2)*(-2*im + omega*R)*yv(0, mi) + MB*(8*MB^2 + 2*MB*R*(-6 + li + li^2 - 3*im*omega*R) + R^2*(2 - li - li^2 - 2*im*omega*R + 2*omega^2*R^2))*yv(1, mi))))/(16*sqrt(2)*R*(-3*MB + R))
end

function _source_kernel_mdag_minus(MB, R, omega, mi, li, rhzb, drhzb, yv)
    return -1/16*((-1)^mi*(2*drhzb*(2*MB - R)*R*(sqrt(-6 + li + li^2)*MB*(MB - im*omega*R^2)*yv(-3, -mi) + 4*sqrt(MB)*omega*R^(3/2)*(MB - im*omega*R^2)*yv(-2, -mi) + sqrt(-2 + li + li^2)*(2*(MB^2 + im*MB*omega*R^2 - im*omega*R^3)*yv(-1, -mi) + im*sqrt(li*(1 + li))*sqrt(MB)*(2*MB - R)*sqrt(R)*yv(0, -mi))) + rhzb*(sqrt(-6 + li + li^2)*MB*(8*MB^2 + 2*MB*R*(-6 + li + li^2 + 3*im*omega*R) + R^2*(2 - li - li^2 + 2*im*omega*R + 2*omega^2*R^2))*yv(-3, -mi) + 4*sqrt(MB)*omega*R^(3/2)*(8*MB^2 + 2*MB*R*(-6 + li + li^2 + 3*im*omega*R) + R^2*(2 - li - li^2 + 2*im*omega*R + 2*omega^2*R^2))*yv(-2, -mi) + sqrt(-2 + li + li^2)*((16*MB^3 + 12*im*MB^2*R*(2*im + omega*R) + 2*MB*R^2*(2 + li + li^2 - 6*im*omega*R - 2*omega^2*R^2) + R^3*(-li - li^2 + 4*omega*R*(2*im + omega*R)))*yv(-1, -mi) - 2*sqrt(li*(1 + li))*sqrt(MB)*(2*MB - R)*R^(3/2)*(2*im + omega*R)*yv(0, -mi)))))/(sqrt(2)*R*(-3*MB + R))
end

function _source_kernel_norm_plus(MB, R, omega, mi, li, rhz, drhz, yv)
    return (2*drhz*sqrt(MB)*(2*MB - R)*R*(-(im*sqrt(-2 + li + li^2)*sqrt(R)*(-2*MB + R)*yv(1, mi)) + sqrt(MB)*(MB + im*omega*R^2)*yv(2, mi)) + rhz*(sqrt(li*(1 + li))*sqrt(-2 + li + li^2)*R*(-2*MB + R)^2*yv(0, mi) + 2*sqrt(-2 + li + li^2)*sqrt(MB)*(2*MB - R)*R^(3/2)*(-2*im + omega*R)*yv(1, mi) + MB*(8*MB^2 + 2*MB*R*(-6 + li + li^2 - 3*im*omega*R) + R^2*(2 - li - li^2 - 2*im*omega*R + 2*omega^2*R^2))*yv(2, mi)))/(24*MB - 8*R)
end

function _source_kernel_norm_minus(MB, R, omega, mi, li, rhzb, drhzb, yv)
    return ((-1)^mi*(2*drhzb*sqrt(MB)*(2*MB - R)*R*(sqrt(MB)*(MB - im*omega*R^2)*yv(-2, -mi) - im*sqrt(-2 + li + li^2)*sqrt(R)*(-2*MB + R)*yv(-1, -mi)) + rhzb*(MB*(8*MB^2 + 2*MB*R*(-6 + li + li^2 + 3*im*omega*R) + R^2*(2 - li - li^2 + 2*im*omega*R + 2*omega^2*R^2))*yv(-2, -mi) + sqrt(-2 + li + li^2)*R*(-2*MB + R)*(2*sqrt(MB)*sqrt(R)*(2*im + omega*R)*yv(-1, -mi) + sqrt(li*(1 + li))*(-2*MB + R)*yv(0, -mi)))))/(24*MB - 8*R)
end
