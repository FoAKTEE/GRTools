"""
Pressureless FR22 (26), eq:linblad_resonances (source spelling retained).
Let x=sqrt(r/a)>0, n>0 and Omega=kappa=n/x^3. All algebra is rational
in x; the unique positive cube root gives r/a=((m+s)/l)^(2/3).
The unsquared resonance equation fixes the branch and excludes spurious
positive radii obtained by taking fractional powers of negative ratios.
"""
function certify_fr22_lindblad()
    @variables x m l n Omega
    c = Certificate[]
    ref = "2207.14637v2::eq:linblad_resonances"
    wt = l*n-m*n/x^3
    for s in (-1,1)
        push!(c, certify_symbolic("FR22 Lindblad sign $s: unsquared radius equation",ref,
            x^3/n*(wt-s*n/x^3)-(l*x^3-(m+s))))
        push!(c, certify_symbolic("FR22 Lindblad sign $s: positive-branch frequency",ref,
            l*n-m*(l*n/(m+s))-s*l*n/(m+s)))
    end
    push!(c, certify_symbolic("FR22 corotation: x^3=m/l",ref,
        [x^3/n*wt-(l*x^3-m), l*n-m*(l*n/m)]))
    # l>0: these positive gaps in x^3 preserve order in r=a*x^2.
    push!(c, certify_symbolic("E-ECC-8 r_minus < r_corotation < r_plus",ref,
        [m/l-(m-1)/l-1/l, (m+1)/l-m/l-1/l]))
    push!(c, certify_symbolic("FR22 m=1,l>0: inner residual is l*n, never zero",ref,
        Symbolics.substitute(wt+n/x^3,Dict(m=>1))-l*n))
    push!(c, certify_symbolic("FR22 m=1: outer radius has x^3=2/l",ref,
        Symbolics.substitute(x^3/n*(wt-n/x^3),Dict(m=>1))-(l*x^3-2)))
    @variables L
    push!(c, certify_symbolic("FR22 l=-L<=0: no finite corotation",ref,
        Symbolics.substitute(wt,Dict(l=>-L))+L*n+m*n/x^3))
    for s in (-1,1)
        push!(c, certify_symbolic("FR22 l=-L<0: sign $s has no positive radius",ref,
            Symbolics.substitute(x^3/n*(wt-s*n/x^3),Dict(l=>-L))+(L*x^3+m+s)))
    end
    push!(c, certify_symbolic("FR22 (m,l)=(1,0): apsidal resonance at every radius",ref,
        Symbolics.substitute(wt+n/x^3,Dict(m=>1,l=>0))))
    push!(c, certify_symbolic("FR22 l=0: only the m=1 inner branch can vanish",ref,
        [Symbolics.substitute(wt+n/x^3,Dict(l=>0))+(m-1)*n/x^3,
         Symbolics.substitute(wt-n/x^3,Dict(l=>0))+(m+1)*n/x^3]))
    println("ASSUMPTIONS a,n,x>0; x=sqrt(r/a); pressureless Omega=kappa=n/x^3; integer m>=1. ",
        "For l>0 and m>1: r_minus/a=((m-1)/l)^(2/3), rc/a=(m/l)^(2/3), ",
        "r_plus/a=((m+1)/l)^(2/3). Positive cube roots are unique and monotone. ",
        "The frequency-form certificates require m+s!=0; m=1,s=-1 is checked separately.")
    println("E-ECC-8: minus is INNER (wt=-kappa); plus is OUTER (wt=+kappa). ",
        "The prose following FR22 (26) interchanges these labels. m=1 has no finite inner ",
        "resonance for l>0. For l=-L<0, L*x^3+m+s>0 excludes both branches; ",
        "(1,0) is everywhere apsidal, not a localized Lindblad radius. No finite ",
        "corotation exists for l<=0; m>1,l=0 has no Lindblad radius.")
    c
end
