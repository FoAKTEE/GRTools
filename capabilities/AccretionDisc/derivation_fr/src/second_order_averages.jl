"""
Canonical FR24 v1 (23)–(32) / FD25 shared AMF and torque certificate.
At fixed r, Sigma is real, amplitudes are arbitrary complex, Phi_ml real,
m,m'≥1, l,l' integers, n=Omega_p=n_p>0. Use physical real fields and
the disc torque sign -r∮deltaSigma*dphi(Phi). No field equation or
assumption of AMF conservation is used. pi_exact is mathematical π.
"""
function certify_fr24_second_order_averages()
    c = certify_fourier_orthogonality()
    @variables ar ai br bi sr si Phi r Sigma m n pi_exact C S delta_m delta_l t
    A, B, density = ar+im*ai, br+im*bi, sr+im*si
    p = real_product_coefficients(A,B)
    # C=cos((l-l')nt), S=sin((l-l')nt). The surviving pm term has
    # exp(-i(l-l')nt), mp its conjugate. delta_m is δ[m,m'].
    phi_integral = 2pi_exact*delta_m*(p.pm*(C-im*S)+p.mp*(C+im*S))
    flux = r^2*Sigma*phi_integral
    printed_flux = pi_exact*r^2*Sigma*delta_m*(real(A*conj(B))*C+imag(A*conj(B))*S)
    push!(c, certify_symbolic("FR24 (23)-(26) instantaneous AMF from real fields",
        "2407.20398v1::eq:AMF-def; 2407.20398v1::eq26-instantaneous-AMF",
        flux-printed_flux))
    push!(c, certify_symbolic("FR24 azimuthal orthogonality kills unequal m", "2407.20398v1::eq24-eq26",
        Symbolics.substitute(flux,Dict(delta_m=>0))))

    theta = m*t # derivative with respect to this azimuth placeholder
    push!(c, certify_symbolic("FR24 potential azimuth derivative amplitude is i m Phi", "2407.20398v1::eq27-torque-definition",
        Symbolics.derivative(Phi*cos(theta),t)+m*Phi*sin(theta)))
    q = real_product_coefficients(density,im*m*Phi)
    torque = -r*2pi_exact*delta_m*(q.pm*(C-im*S)+q.mp*(C+im*S))
    printed_torque = -pi_exact*r*m*Phi*delta_m*(si*C-sr*S)
    push!(c, certify_symbolic("FR24 (27)-(28) instantaneous torque from real fields",
        "2407.20398v1::eq:time_dependent_dTdr", torque-printed_torque))
    push!(c, certify_symbolic("FR24 torque azimuthal orthogonality kills unequal m", "2407.20398v1::eq:time_dependent_dTdr",
        Symbolics.substitute(torque,Dict(delta_m=>0))))

    # Time averaging uses the same character lemma in M=nt; the explicit
    # normalization cancels the Jacobian. Thus <C>=delta_l and <S>=0.
    push!(c, certify_symbolic("FR24 orbit-average measure has unit constant mode", "2407.20398v1::eq29-eq30",
        n/(2pi_exact)*(2pi_exact/n)-1))
    time_average = Dict(C=>delta_l,S=>0)
    avg_flux = Symbolics.substitute(flux,time_average)
    avg_torque = Symbolics.substitute(torque,time_average)
    push!(c, certify_symbolic("FR24 (29)-(30) only equal (m,l) survive",
        "2407.20398v1::eq29-eq30; 2506.19917v3::eq:orbit_avg_amf_dtdr",
        [avg_flux-pi_exact*r^2*Sigma*delta_m*delta_l*real(A*conj(B)),
         avg_torque+pi_exact*r*m*Phi*delta_m*delta_l*si]))
    diagonal = Dict(delta_m=>1,delta_l=>1)
    modal_flux = Symbolics.substitute(avg_flux,diagonal)
    modal_torque = Symbolics.substitute(avg_torque,diagonal)
    # Attach EVERY source-qualified label owned here, rather than treating
    # a duplicate label as an unqualified assertion of shared physics.
    for paper in ("2407.20398v1","2506.19917v3")
        push!(c, certify_symbolic("$paper AMF definition and real-field expansion",
            "$paper::eq:AMF-def", flux-printed_flux))
        push!(c, certify_symbolic("$paper F_J,ml = pi r^2 Sigma Re[ur uphi*]",
            "$paper::eq:F_jml", modal_flux-pi_exact*r^2*Sigma*real(A*conj(B))))
        push!(c, certify_symbolic("$paper dT_ml/dr = -pi r m Phi Im[dSigma]",
            "$paper::eq:dTdr_jml", modal_torque+pi_exact*r*m*Phi*imag(density)))
    end
    push!(c, certify_symbolic("FD25 torque definition shares the disc sign", "2506.19917v3::eq:dTdr_def",
        torque-printed_torque))
    push!(c, certify_symbolic("FR real amplitude / 2 gives 4pi complex-pair convention", "FR-to-complex-reality-pair normalization",
        4pi_exact*r^2*Sigma*real((A/2)*conj(B/2))-modal_flux))
    println("ASSUMPTIONS m,m'>=1; l,l' in Z; fixed-r axisymmetric real Sigma; real Phi_ml. ",
        "FR24 Omega_p and FD25 n_p both map to n; FR22 lm maps to ml without changing indices.")
    println("AVERAGING phi integration keeps cross-l terms at fixed m; the normalized ",
        "physical-time average kills l!=l'. Finite sums extend to L2-convergent fields. ",
        "Reality pairs are (m,l),(-m,-l); the one-sided Re convention gives pi, not 4pi.")
    println("SOURCE CONVENTION FR24 (26): both cosine and sine terms carry pi*r^2*Sigma ",
        "and both sums. The result certifies wave AMF and disc torque, not equality of their radial gradients.")
    c
end
