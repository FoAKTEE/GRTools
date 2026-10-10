"""
Coefficients of the four characters in Re[A exp(iθ)] Re[B exp(iθ')].
The keys give the signs of θ and θ'. Keeping all four terms makes reality
pairing and its factors explicit. Used for cosine, sine, AMF and torque.
"""
function real_product_coefficients(A, B)
    (pp=A*B/4, pm=A*conj(B)/4, mp=conj(A)*B/4, mm=conj(A)*conj(B)/4)
end

"""
Common analytic orthogonality lemma, with exact symbolic algebraic checks.
For an integer k≠0, the endpoint values sin(2πk)=sin(0)=0 and
cos(2πk)=cos(0)=1 follow from periodicity. Check the primitives and their
endpoint differences separately; the k=0 branch is integrated separately.
Thus ∮exp(ikx)dx=2π δ[k,0]. No finite mode scan or numerical quadrature is
used as a proof. pi_exact denotes mathematical π, kept symbolic throughout.
"""
function certify_fourier_orthogonality()
    @variables x k pi_exact ar ai br bi z w
    c = Certificate[]
    primitives = [sin(k*x)/k, -cos(k*x)/k]
    push!(c, certify_symbolic("Fourier nonzero-character primitive", "integer-character lemma",
        [Symbolics.derivative(primitives[1],x)-cos(k*x),
         Symbolics.derivative(primitives[2],x)-sin(k*x)]))
    endpoint = Dict(sin(k*x)=>0, cos(k*x)=>1)
    # At x=0 the two primitive values are exactly 0 and -1/k.
    push!(c, certify_symbolic("Fourier nonzero-character endpoint cancellation", "integer-character lemma",
        [Symbolics.substitute(primitives[1],endpoint),
         Symbolics.substitute(primitives[2],endpoint)+1/k]))
    push!(c, certify_symbolic("Fourier zero-character integral is 2pi", "integer-character lemma",
        Symbolics.substitute(x,Dict(x=>2pi_exact))-Symbolics.substitute(x,Dict(x=>0))-2pi_exact))
    A, B = ar+im*ai, br+im*bi
    p = real_product_coefficients(A,B)
    product = (A*z+conj(A)/z)*(B*w+conj(B)/w)/4
    push!(c, certify_symbolic("Fourier four reality-paired product characters", "integer-character lemma",
        product-(p.pp*z*w+p.pm*z/w+p.mp*w/z+p.mm/(z*w))))
    # The difference characters are conjugates, even for complex responses.
    push!(c, certify_symbolic("Fourier surviving pair gives one half Re[A B*]", "integer-character lemma",
        p.pm+p.mp-real(A*conj(B))/2))
    println("ANALYTIC LEMMA integer k; k!=0 uses the checked primitive and periodic ",
        "endpoints; k=0 gives 2pi. Apply independently in phi and M. For m,m'>=1, ",
        "m+m'>0 kills pp/mm; pm/mp survive iff m=m'. In time they survive iff l=l'.")
    c
end

"""
FR22 v2 (9)–(12), with E-ECC-2 correcting the physical-time inverse (22).
Inputs: ECC::fr22-eq-kepler, pericentre at phi=t=0, n>0, 0≤e<1,
real direct potential and time-reversal-even softening. The source's sum is
the m≥1 projection: it omits the entire m=0 sector, which need not vanish.
All identities hold for finite sums. L² Fourier completeness and Cauchy–
Schwarz extend projections and integrated products to convergent L² fields;
pointwise convergence at a singular unsoftened source is not asserted.
"""
function certify_fr22_fourier()
    c = certify_fourier_orthogonality()
    @variables r rp eps2 Gmu cp sp cpp spp phi t m l n z potential amplitude
    @variables pi_exact delta_m delta_l ar ai Omega_p n_p
    ref = "2207.14637v2::eq:cosine_expansion"
    # Kepler gives rp(-t)=rp(t), cos(phi_p) even, sin(phi_p) odd.
    # cp/sp denote cos/sin(phi); cpp/spp denote cos/sin(phi_p).
    separation2 = r^2+rp^2-2r*rp*(cp*cpp+sp*spp)+eps2
    reversal = Dict(sp=>-sp,spp=>-spp)
    push!(c, certify_symbolic("FR22 direct potential under (phi,t)->(-phi,-t)", ref,
        Symbolics.substitute(-Gmu/sqrt(separation2),reversal)+Gmu/sqrt(separation2)))
    theta = m*phi-l*n*t
    push!(c, certify_symbolic("FR22 Fourier phase reverses sign", ref,
        Symbolics.substitute(theta,Dict(phi=>-phi,t=>-t))+theta))
    # z=exp(iθ). The torus reversal preserves dphi dM and the potential,
    # while the sine kernel changes sign, so I_s=-I_s and I_s=0.
    sine = (z-1/z)/(2im)
    cosine = (z+1/z)/2
    push!(c, certify_symbolic("FR22 odd sine projection cancels on the torus", ref,
        potential*(sine+Symbolics.substitute(sine,Dict(z=>1/z)))))
    push!(c, certify_symbolic("FR22 even cosine projection is retained", ref,
        Symbolics.substitute(cosine,Dict(z=>1/z))-cosine))
    push!(c, certify_symbolic("FR22 real coefficient reality pair equals cosine", ref,
        amplitude*z/2+amplitude/(2z)-amplitude*cosine))

    cc = real_product_coefficients(Num(1),Num(1))
    cs = real_product_coefficients(Num(1),-im*Num(1)) # Re[-i exp(iθ)]=sin θ
    # delta_m=δ[m,m'], delta_l=δ[l,l']; pp/mm vanish since m+m'>0.
    inner = (2pi_exact)^2/n*(cc.pm+cc.mp)*delta_m*delta_l
    push!(c, certify_symbolic("FR22 double cosine orthogonality: 2pi^2/n", "2207.14637v2::eq:fourier_expansion",
        inner-2pi_exact^2/n*delta_m*delta_l))
    push!(c, certify_symbolic("FR22 cosine-sine orthogonality", "2207.14637v2::eq:fourier_expansion",
        (2pi_exact)^2/n*(cs.pm+cs.mp)*delta_m*delta_l))
    push!(c, certify_symbolic("FR22 corrected inverse has unit modal gain", "2207.14637v2::eq:phi_mode_integral; E-ECC-2",
        n/(2pi_exact^2)*inner-delta_m*delta_l))
    push!(c, certify_symbolic("E-ECC-2 printed inverse has gain 1/n", "2207.14637v2::eq:phi_mode_integral; E-ECC-2",
        inner/(2pi_exact^2)-delta_m*delta_l/n))
    @variables M time_unit length_unit
    push!(c, certify_symbolic("FR22 physical-time to mean-anomaly measure", "2207.14637v2::eq:phi_mode_integral; E-ECC-2",
        [Symbolics.derivative(n*t,t)-n, n*Symbolics.derivative(M/n,M)-1]))
    push!(c, certify_symbolic("E-ECC-2 printed units have extra time; corrected units", "2207.14637v2::eq:phi_mode_integral; E-ECC-2",
        [(length_unit^2/time_unit^2*time_unit)/(length_unit^2/time_unit^2)-time_unit,
         (1/time_unit)*(length_unit^2/time_unit^2*time_unit)-length_unit^2/time_unit^2]))

    speed = l*n/m
    push!(c, certify_symbolic("FR22 pattern speed and total frequency", "2207.14637v2::eq:olm",
        [m*(phi-speed*t)-theta, Symbolics.derivative(theta,t)+l*n,
         Symbolics.derivative(theta,t)+speed*Symbolics.derivative(theta,phi)]))
    # Differentiate the actual real modal ansatz, retaining complex response.
    A = ar+im*ai
    response = ar*cos(theta)-ai*sin(theta)
    dphi_amplitude, dt_amplitude = im*m*A, -im*l*n*A
    push!(c, certify_symbolic("FR22 modal ansatz: azimuth and time derivatives", "2207.14637v2::fr22:eq12-modal-ansatz",
        [Symbolics.derivative(response,phi)-(real(dphi_amplitude)*cos(theta)-imag(dphi_amplitude)*sin(theta)),
         Symbolics.derivative(response,t)-(real(dt_amplitude)*cos(theta)-imag(dt_amplitude)*sin(theta))]))
    # Each source uses the same m>=1,l in Z convention; only the frequency
    # name and the order of subscripts change (FR22 lm -> FR24/FD25 ml).
    for (paper,frequency,phase) in (("2407.20398v1",Omega_p,m*phi-l*Omega_p*t),
                                  ("2506.19917v3",n_p,m*(phi-l*n_p/m*t)))
        push!(c, certify_symbolic("$paper shared potential expansion crosswalk",
            "$paper::eqn:potential_expansion",
            amplitude*(Symbolics.substitute(phase,Dict(frequency=>n))-theta)))
        push!(c, certify_symbolic("$paper shared modal response crosswalk",
            "$paper::eqn:modal_repsonse",
            Symbolics.substitute(ar*cos(phase)-ai*sin(phase),Dict(frequency=>n))-response))
    end
    println("ASSUMPTIONS FR22 parity is joint (phi,t) reversal, not phi reversal at fixed t; ",
        "the Kepler certificate fixes pericentre and proves the orbit parities. Positive ",
        "separation radicand fixes the square-root branch. Softening must be even in t.")
    println("SCOPE m>=1, l in Z (including l=0 and negative l); m=0 is omitted, not zero. ",
        "Reality pairing is (m,l)<->(-m,-l). No independent negative-m modes are summed.")
    println("ERRATUM E-ECC-2: FR22 (22), 3_numerical_procedure.tex:30, needs n/(2pi^2) ",
        "for physical dt. Printed 1/(2pi^2) returns amplitude/n and has an extra time unit; ",
        "only n=1 units hide the error. The M=n*t inverse has 1/(2pi^2).")
    c
end
