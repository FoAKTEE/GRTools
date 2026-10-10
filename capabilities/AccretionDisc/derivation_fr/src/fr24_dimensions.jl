"""
Dimensional audit of every numbered FR24 v1 formula and the dimension-bearing
inline definitions used by them, including the appendix normalisations.
Each term is replaced by its M,L,T (and temperature) monomial. Addition,
subtraction, real/imaginary parts, sums and averages preserve dimensions;
each additive term is checked separately. Numerical coefficients and p,q,h,e,
b,m,l are dimensionless. d/dr divides by L and dt integration multiplies by T.
Log derivatives mean X'/X (equivalently log(X/Xref)), not log of a unit.
This audits units only, not coefficients, signs, averaging or fitting accuracy.
Source anchors are the begin-equation/align lines in the immutable v1 mirror.
"""
function certify_fr24_dimensions()
    @variables mass length time temperature
    U = Num(1)
    G = length^3/(mass*time^2)
    sigma = mass/length^2
    cs = length/time
    pressure = mass/time^2
    omega = 1/time
    potential = cs^2
    detuning = omega^2
    flux = mass*length^2/time^2
    angular = mass*length^2/time
    invLT = 1/length
    c = Certificate[]
    # A single verifier adapter for the >=55 source equations: no dimensional
    # arithmetic implementation, just the same R3 symbolic zero certificate.
    dim(name, ref, lhs, terms...) = push!(c, certify_symbolic(
        "FR24 dimensions: "*name, ref, Num[lhs-term for term in terms]))

    dim("direct potential", "2_linear_framework.tex:17", potential, G*mass/length)
    dim("surface density", "2_linear_framework.tex:26", sigma, sigma*U)
    dim("sound speed", "2_linear_framework.tex:33", cs, cs*U)
    dim("EoS", "2_linear_framework.tex:40", pressure, cs^2*sigma)
    # Tdisc ∝ P/Sigma contains the implicit mean-particle-mass/k_B factor.
    dim("temperature proportionality", "2_linear_framework.tex:47", temperature,
        (temperature/cs^2)*pressure/sigma)
    dim("aspect ratio", "2_linear_framework.tex:54", U, (cs/omega)/length)
    dim("Omega squared", "2_linear_framework.tex:63", omega^2,
        G*mass/length^3, (pressure/length)/(length*sigma), omega^2*U)
    dim("kappa squared", "2_linear_framework.tex:70", omega^2,
        (omega/length)*(length^2*omega/length), omega^2*U)
    dim("potential expansion", "2_linear_framework.tex:83", potential, potential*U)
    dim("pattern speed", "2_linear_framework.tex:90", omega, U*omega)
    dim("modal response phase", "2_linear_framework.tex:97", U, omega*time)
    dim("continuity", "2_linear_framework.tex:104", sigma/time,
        omega*sigma, (length*sigma*cs/length)/length, sigma*cs/length)
    dim("radial Euler", "2_linear_framework.tex:104", cs/time,
        omega*cs, pressure/(sigma*length), pressure*sigma/(length*sigma^2), potential/length)
    dim("azimuthal Euler", "2_linear_framework.tex:104", cs/time,
        omega*cs, omega^2/omega*cs, pressure/(length*sigma), potential/length)
    dim("master equation all terms", "2_linear_framework.tex:117", potential/length^2,
        potential/length^2, invLT*potential/length,
        omega/(length*omega)*invLT*potential,
        invLT^2*potential, potential/length^2, detuning/cs^2*potential,
        potential/length^2, invLT*potential/length,
        omega/(length*omega)*invLT*potential, potential/length^2)
    dim("thermal scale", "2_linear_framework.tex:128", invLT, (cs^2/length)/cs^2)
    dim("detuning", "2_linear_framework.tex:135", detuning, omega^2, omega^2)
    dim("density reconstruction", "2_linear_framework.tex:142", sigma, sigma*potential/cs^2)
    dim("radial velocity reconstruction", "2_linear_framework.tex:142", cs,
        omega*potential/(length*detuning), omega*potential/(length*detuning),
        omega*invLT*potential/detuning)
    dim("azimuthal velocity reconstruction", "2_linear_framework.tex:142", cs,
        omega^2*potential/(omega*length*detuning), omega*potential/(length*detuning),
        omega^2*invLT*potential/(omega*detuning))
    dim("Lindblad radii", "2_linear_framework.tex:165", U, length/length)
    # Printed eq:phi_mode_integral is NOT dimensionally correct in physical t.
    # Certify the extra T explicitly; PL00 eq (21), paplar99.tex:580–583,
    # supplies omega/(2pi²(1+delta_m0)); here m>=1, so delta_m0=0.
    dim("printed Fourier integral has extra time", "2_linear_framework.tex:196",
        time, (potential*time)/potential, (G*mass*time/length)/potential)
    dim("Fourier integral corrected by Omega_p", "2_linear_framework.tex:196",
        potential, omega*potential*time, omega*G*mass*time/length)
    println("SOURCE DISCREPANCY FR24 eq:phi_mode_integral lacks Omega_p for dimensional t; ",
            "PL00 (21) includes omega. In M=Omega_p*t or Omega_p=1 units the discrepancy is hidden.")
    dim("AMF definition", "2_linear_framework.tex:222", flux, length^2*sigma*cs^2)
    dim("radial velocity expansion", "2_linear_framework.tex:229", cs, cs*U)
    dim("azimuthal velocity expansion", "2_linear_framework.tex:232", cs, cs*U)
    # The common pi*r²*Sigma and sums apply to both sine/cosine terms.
    dim("time dependent AMF", "2_linear_framework.tex:238", flux,
        length^2*sigma*cs^2, length^2*sigma*cs^2)
    dim("torque density definition", "2_linear_framework.tex:245", flux/length,
        length*sigma*potential)
    dim("time dependent torque density", "2_linear_framework.tex:251", flux/length,
        length*potential*sigma, length*potential*sigma)
    dim("orbit averaged AMF", "2_linear_framework.tex:260", flux, flux, omega*flux*time)
    dim("orbit averaged torque density", "2_linear_framework.tex:260", flux/length,
        flux/length, omega*(flux/length)*time)
    dim("modal AMF", "2_linear_framework.tex:268", flux, length^2*sigma*cs^2)
    dim("modal torque density", "2_linear_framework.tex:268", flux/length, length*potential*sigma)
    dim("one sided torque", "2_linear_framework.tex:277", flux, (flux/length)*length)
    dim("net modal torque", "2_linear_framework.tex:286", flux,
        (flux/length)*length, flux, flux)
    dim("corotation radius", "2_linear_framework.tex:294", length, length*U)
    dim("corotation AMF jump", "2_linear_framework.tex:300", flux, flux, flux)
    dim("corotation torque jump", "2_linear_framework.tex:300", flux, flux, flux)
    dim("corotation torque", "2_linear_framework.tex:307", flux, flux, flux)
    dim("Lindblad torque", "2_linear_framework.tex:313", flux, flux, flux)

    dim("F_J0", "3_circular_benchmark.tex:33", flux, sigma*length^4*omega^2*U)
    dim("vortensity", "3_circular_benchmark.tex:42", omega/sigma,
        (length^2*omega/length)/(length*sigma), omega/sigma*U)
    dim("WKB torque", "3_circular_benchmark.tex:73", flux,
        sigma/omega^2*(length*potential/length)^2,
        sigma/omega^2*(omega*potential/omega)^2)
    dim("WKB pressure correction", "3_circular_benchmark.tex:73", U, cs^2/(length*omega)^2)
    dim("net torque fit", "3_circular_benchmark.tex:100", flux, U*flux)
    dim("Lindblad torque fit", "3_circular_benchmark.tex:100", flux, U*flux)
    # Compare squares of the root expression; positive unit scales fix the root.
    dim("planet angular momentum squared", "4_eccentric_torques.tex:107", angular^2,
        mass^2*G*mass*length, (mass*omega*length^2)^2)
    dim("orbital evolution logarithmic rates", "4_eccentric_torques.tex:114", omega,
        (angular/time)/angular, (length/time)/length, U/time)
    dim("inverse timescales", "4_eccentric_torques.tex:121", omega,
        (angular/time)/angular, (length/time)/length, U/time)
    dim("timescale relation", "4_eccentric_torques.tex:128", omega, omega, U*omega)
    dim("angular momentum rate", "4_eccentric_torques.tex:135", omega, flux/angular)
    dim("semimajor axis rate", "4_eccentric_torques.tex:141", omega,
        omega*flux/(mass*length^2*omega^2))
    dim("tau0 inverse", "4_eccentric_torques.tex:158", omega, flux/(mass*omega*length^2))
    dim("PL00 angular momentum fit", "5_discussion.tex:61", omega, U*omega)
    dim("PL00 eccentricity fit", "5_discussion.tex:66", omega, U*omega)
    dim("CN06 angular momentum fit", "5_discussion.tex:79", omega, U*omega)
    dim("CN06 eccentricity fit", "5_discussion.tex:86", omega, U*omega)

    dim("inline Kepler and scale height", "inline:section2.1", omega^2, G*mass/length^3)
    dim("inline Doppler and resonance frequencies", "inline:section2.2-resonances", omega,
        U*omega, omega)
    dim("inline softening and corotation width", "inline:section2.3", length, U*cs/omega, U*length)
    dim("inline azimuths and Fourier orders", "inline:section2.3-angles", U, U, length/(cs/omega))
    dim("inline potential eccentricity scaling", "inline:section2.3-small-e", potential, potential*U)
    dim("inline thermal mass", "inline:section2.2", mass, U^3*mass)
    dim("inline pseudo enthalpy", "inline:section2.2", potential, pressure/sigma)
    dim("inline Jacobi energy", "inline:section4.4", flux, G*mass^2/length, omega*angular)
    dim("inline energy and momentum transfer", "inline:section4.4", flux/time, omega*flux)
    dim("inline angular momentum transfer", "inline:section4.4-torque", angular/time, flux)
    dim("inline torque and AMF gradients", "inline:sections3-4", flux/length, flux/length)
    dim("inline scaled e and Mach number", "inline:section4", U, U/U, cs/cs)
    dim("inline orbital period and time", "inline:section5", time, 1/omega)
    dim("inline physical disc and planet parameters", "inline:section5-units", sigma, mass/length^2)
    dim("figure and appendix normalised quantities", "inline:figures-and-appendices", U,
        flux/flux, (flux/length)/(flux/length), omega/omega,
        length/length, sigma/sigma, time*omega, mass/mass)
    println("AUDIT CONVENTIONS temperature proportionality includes particle-mass/k_B; ",
            "c_{s^2} in the printed master/reconstruction equations means c_s^2; ",
            "the AMF sine/cosine terms share their prefactor and sums.")
    println("AUDIT CONVENTIONS inline zero conditions apply separately to each quantity; ",
            "dimensional example values carry their adjacent AU, g cm^-2, year, or solar/Earth mass units; ",
            "power-law proportionalities include a dimensional reference amplitude. All remaining ",
            "parameter/order comparisons involve dimensionless e,h,p,q,b,m,l,Q_e or their ratios.")
    c
end
