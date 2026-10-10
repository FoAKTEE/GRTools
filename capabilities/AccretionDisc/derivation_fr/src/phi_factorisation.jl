"""
Derive, rather than attribute to the paper, the factorisation of corrected
FR22 (22)–(23). With A=m*phi_p(M)-l*M and dimensional softening epsilon,
Phi_lm = -G*Mp/(2pi) int dM cos(A)*b_(1/2)^m(r/rp;epsilon/rp)/rp.

The analytic integral steps use translation invariance of a full-period
integral, odd-function cancellation, and the definition of b; the symbolic
certificates check their algebra and normalization exactly. Positive rp and
epsilon make the kernel smooth and bounded, justify Fubini/differentiation
under the integral and fix the square-root branch. Epsilon is independent
of the field radius r. No random identity testing is needed.
"""
function certify_fr22_phi_factorisation()
    c = certify_fourier_orthogonality()
    @variables r rp epsilon cp sp t phi_p psi m l M n Gmu pi_exact z w b delta_lm
    ref = "2207.14637v2::eq:phi_mode_integral; E-ECC-2"
    # A rational unit-circle chart for phi_p proves the geometric shift.
    C, S = (1-t^2)/(1+t^2), 2t/(1+t^2)
    shifted_cos, shifted_sin = C*cp-S*sp, S*cp+C*sp
    push!(c, certify_symbolic("FR22 azimuth shift preserves separation",ref,
        shifted_cos*C+shifted_sin*S-cp))
    push!(c, certify_symbolic("FR22 shifted azimuth Jacobian and modal phase",ref,
        [Symbolics.derivative(phi_p+psi,psi)-1,
         m*(phi_p+psi)-l*M-(m*psi+(m*phi_p-l*M))]))
    D = r^2+rp^2-2r*rp*cp+epsilon^2
    Q = 1+(r/rp)^2-2(r/rp)*cp+(epsilon/rp)^2
    push!(c, certify_symbolic("FR22 positive-root kernel scaling D=rp^2*Q",ref,D-rp^2*Q))
    @variables u
    kernel_chart = Symbolics.substitute(D,Dict(cp=>(1-u^2)/(1+u^2)))
    push!(c, certify_symbolic("FR22 Plummer kernel is even under psi reversal",ref,
        Symbolics.substitute(kernel_chart^(-1//2),Dict(u=>-u))-kernel_chart^(-1//2)))
    # z=exp(i*m*psi), w=exp(i*A). This is the exact addition formula.
    Cpsi, Spsi = (z+1/z)/2, (z-1/z)/(2im)
    CA, SA = (w+1/w)/2, (w-1/w)/(2im)
    shifted = (z*w+1/(z*w))/2
    push!(c, certify_symbolic("FR22 cosine addition in shifted azimuth",ref,
        shifted-(Cpsi*CA-Spsi*SA)))
    @variables kernel
    push!(c, certify_symbolic("FR22 even Plummer kernel cancels sine projection",ref,
        kernel*(Spsi+Symbolics.substitute(Spsi,Dict(z=>1/z)))))
    push!(c, certify_symbolic("FR22 even-paired integrand retains cos(A) cos(m*psi)",ref,
        kernel*(shifted+Symbolics.substitute(shifted,Dict(z=>1/z)))/2-kernel*Cpsi*CA))
    push!(c, certify_symbolic("FR22 corrected double integral factorises with -GMp/(2pi)",ref,
        -Gmu*n/(2pi_exact^2)*(1/n)*(pi_exact*b/rp)*CA+Gmu/(2pi_exact)*b/rp*CA))
    push!(c, certify_symbolic("FR22 circular orbit selects l=m by character orthogonality",ref,
        -Gmu/(2pi_exact)*(b/rp)*(2pi_exact*delta_lm)+Gmu*b/rp*delta_lm))
    @variables alpha eta B0 B1 B2 alpha0
    # An arbitrary twice-differentiable b has this second jet at alpha0.
    # rp,epsilon fixed under d/dr; eta=epsilon/rp therefore stays fixed.
    jet = (B0+B1*(r/rp-alpha0)+B2*(r/rp-alpha0)^2/2)/rp
    push!(c, certify_symbolic("FR22 radial derivatives have b'/rp^2 and b''/rp^3",ref,
        [Symbolics.substitute(Symbolics.derivative(jet,r),Dict(r=>rp*alpha0))-B1/rp^2,
         Symbolics.derivative(Symbolics.derivative(jet,r),r)-B2/rp^3,
         Symbolics.derivative(epsilon/rp,r)]))
    q = 1+alpha^2-2alpha*cp+eta^2
    k = q^(-1//2)
    # The differentiated integrands used by the runtime Laplace kernels.
    first = (cp-alpha)*q^(-3//2)
    second = (1//2+2*(alpha^2-2alpha*cp)+3//2*(2cp^2-1)-eta^2)*q^(-5//2)
    push!(c, certify_symbolic("runtime Laplace first derivative integrand", "1811.09628::eq:laplace_coef",
        Symbolics.derivative(k,alpha)-first))
    push!(c, certify_symbolic("runtime Laplace second derivative integrand", "1811.09628::eq:laplace_coef",
        Symbolics.derivative(first,alpha)-second))
    println("DERIVATION phi=phi_p+psi has Jacobian 1. Integer m makes the shifted ",
        "integrand 2pi-periodic, so translate its interval to [-pi,pi]. The Plummer ",
        "kernel is even in psi; its sine integral is zero and its cosine integral ",
        "is pi*b/rp by MR19 (13), the continuum definition of runtime laplace_b_soft.")
    println("SCOPE a,rp,n>0, 0<=e<1, integer m>=1 and l in Z; epsilon>0 fixed with respect ",
        "to r (it may depend on orbital phase). sqrt(D)=rp*sqrt(Q) by positive roots. ",
        "The rational angle chart extends at its omitted point by continuity. ",
        "E-ECC-2 supplies n/(2pi^2), dt=dM/n. No indirect term is included.")
    println("CIRCULAR rp=a, phi_p=M and constant epsilon give -GMp/a*b*delta_lm ",
        "by the checked integer-character lemma: MR19 (12), runtime Phi_m with indirect=0 ",
        "in a=GMp=1 units. The runtime midpoint sum approximates the continuum b; ",
        "its finite-quadrature agreement is numerical evidence, not an exact certificate.")
    c
end
