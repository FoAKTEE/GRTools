"""
First variation of inviscid polar continuity and Euler about ur=0, uphi=r*Omega.
The local radial coordinate x carries arbitrary first jets; lambda counts the
perturbation. z is the formal character exp(i*m*(phi-omega_lm*t)), so its
derivatives implement the already certified Fourier ansatz. Linear real
operators extend complex-linearly, hence arbitrary complex amplitudes are allowed.
"""
function fr22_linearised_residuals()
    @variables x lambda z r Sigma Sigma_r P P_r Omega Omega_r Phi0_r
    @variables ds ds_r dp dp_r ur ur_r uphi uphi_r Phi Phi_r m omega_lm kappa2
    radius = r+x
    density = Sigma+Sigma_r*x+lambda*(ds+ds_r*x)*z
    pressure = P+P_r*x+lambda*(dp+dp_r*x)*z
    vr = lambda*(ur+ur_r*x)*z
    vp = radius*(Omega+Omega_r*x)+lambda*(uphi+uphi_r*x)*z
    potential = Phi0_r*x+lambda*(Phi+Phi_r*x)*z
    dr(f) = Symbolics.derivative(f,x)
    dphi(f) = im*m*z*Symbolics.derivative(f,z)
    dt(f) = -im*m*omega_lm*z*Symbolics.derivative(f,z)
    nonlinear = [dt(density)+dr(radius*density*vr)/radius+dphi(density*vp)/radius,
        dt(vr)+vr*dr(vr)+vp*dphi(vr)/radius-vp^2/radius+dr(pressure)/density+dr(potential),
        dt(vp)+vr*dr(vp)+vp*dphi(vp)/radius+vr*vp/radius+dphi(pressure)/(radius*density)+dphi(potential)/radius]
    base = Dict(x=>0,lambda=>0,z=>1)
    # Differentiate scalar real/imaginary parts: derivative(complex(...),lambda)
    # otherwise leaves an unevaluated complex-constructor differential.
    derived = [Symbolics.substitute(Symbolics.derivative(real(f),lambda),base)+
        im*Symbolics.substitute(Symbolics.derivative(imag(f),lambda),base) for f in nonlinear]
    equilibrium = [Symbolics.substitute(Symbolics.substitute(f,base),
        Dict(Phi0_r=>r*Omega^2-P_r/Sigma)) for f in nonlinear]
    wt = m*(omega_lm-Omega)
    # Independent transcription: FR22 2_framework.tex:123,128,133, RHS moved left.
    printed = [-im*wt*ds+(Sigma*ur+r*Sigma_r*ur+r*Sigma*ur_r)/r+im*m*Sigma*uphi/r,
        -im*wt*ur-2Omega*uphi+dp_r/Sigma-P_r*ds/Sigma^2+Phi_r,
        -im*wt*uphi+kappa2*ur/(2Omega)+im*m/r*(dp/Sigma+Phi)]
    epicycle = Dict(kappa2=>2Omega/r*(2r*Omega+r^2*Omega_r))
    (; derived, printed, equilibrium, epicycle,
       symbols=(; r,Sigma,Sigma_r,P_r,Omega,Omega_r,ds,dp,dp_r,ur,ur_r,uphi,Phi,Phi_r,m,omega_lm,kappa2))
end

"""
FR22 (13)–(15), MR20 (2)–(4), and the shared FR24/FD25 equations.
Sources: FR22 Sections/2_framework.tex:120–137; MR20 src/main.tex:71–77;
FR24 Sections/2_linear_framework.tex:105–116; FD25 Sections/2_linear_theory.tex:39–49.
Assume r,Sigma,Omega>0, smooth axisymmetric equilibrium, no viscosity or
disc self-gravity. No thermodynamic closure is used to derive these three
equations: arbitrary radial temperature gradients and delta P are retained.
"""
function certify_fr22_linearised()
    system = fr22_linearised_residuals()
    (; r,Sigma,Sigma_r,P_r,Omega,ds,dp,dp_r,ur,ur_r,uphi,Phi,Phi_r,m,omega_lm,kappa2) = system.symbols
    c = Certificate[]
    push!(c, certify_symbolic("FR22 stationary background satisfies polar fluid equations",
        "2207.14637v2::eq:framework:orbital_frequency",system.equilibrium))
    labels = ("eq:density_pert_eqn","eq:ur_eqn","eq:uphi_eqn")
    for (j,label) in enumerate(labels)
        residual = system.derived[j]-Symbolics.substitute(system.printed[j],system.epicycle)
        push!(c, certify_symbolic("FR22 first variation of polar equation $j",
            "2207.14637v2::$label",residual))
    end
    @variables omega_p Phi_m Phi_m_r l n omega_ml
    # MR20 main.tex:71–75, transcribed in its own source notation.
    mr20 = [-im*m*(omega_p-Omega)*ds+(Sigma*ur+r*Sigma_r*ur+r*Sigma*ur_r)/r+im*m*Sigma*uphi/r,
        -im*m*(omega_p-Omega)*ur-2Omega*uphi+dp_r/Sigma-P_r*ds/Sigma^2+Phi_m_r,
        -im*m*(omega_p-Omega)*uphi+kappa2*ur/(2Omega)+im*m/r*(dp/Sigma+Phi_m)]
    crosswalk = Dict(omega_p=>omega_lm,Phi_m=>Phi,Phi_m_r=>Phi_r)
    for j in 1:3
        push!(c, certify_symbolic("MR20 equation $(j+1) equals FR22 equation $(j+12)",
            "1911.01428::eq:pert$j; 2207.14637v2::$(labels[j])",
            Symbolics.substitute(mr20[j],crosswalk)-system.printed[j]))
    end
    push!(c, certify_symbolic("FR22 Doppler frequency m(omega_lm-Omega)=l*n-m*Omega",
        "2207.14637v2::eq:olm",m*(l*n/m-Omega)-(l*n-m*Omega)))
    # FR24 and FD25 print the same three operators (source locations above),
    # changing lm -> ml only. Subscripts designate mode indices, not a swap of values.
    for paper in ("2407.20398v1","2506.19917v3"), (j,label) in enumerate(labels)
        shared = Symbolics.substitute(system.printed[j],Dict(omega_lm=>omega_ml))
        push!(c, certify_symbolic("$paper shared fluid equation $j", "$paper::$label",
            Symbolics.substitute(shared,Dict(omega_ml=>omega_lm, kappa2=>system.epicycle[kappa2]))-system.derived[j]))
    end
    println("ASSUMPTIONS r,Sigma,Omega>0; smooth axisymmetric inviscid 2D equilibrium; ",
        "u_r=0, u_phi=r*Omega, Phi0'=r*Omega^2-P'/Sigma. First variation only; ",
        "no EoS closure or q=0 restriction. Fourier frequency is l*n, not the pattern speed.")
    c
end
