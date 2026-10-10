"""
C-NORM owns amplitudes, phases, units and reality-pair counting only.
C-HARM (ECC::st-worldline-rates) supplies omega=m*Omega_phi+k*Omega_r
and the coordinate-time projector. No orbit or projector is derived here.
"""
module EccentricConventions
using Symbolics
using ..R3Derivation: Certificate, certify_symbolic, certify_bool

# Reuse the admitted Newtonian character algebra and the canonical mode adapter.
include("../../derivation_fr/src/fourier.jl")
include("../../src/ForcingModes.jl")

"One representative of (m,k),(-m,-k); (0,0) is excluded from dynamic feedback."
independent_dynamic_pair(m::Integer, k::Integer) = m > 0 || (m == 0 && k > 0)

"Azimuth integral/time-average weight for one independent nonconstant pair."
function reality_pair_weight(m::Integer, k::Integer, ::Type{T}=Float64) where {T<:AbstractFloat}
    independent_dynamic_pair(m,k) || throw(DomainError((m,k), "need m>0 or m=0,k>0"))
    4*T(pi)
end

"""
Azimuth integral and normalized radial-phase average of two real paired modes.
Characters are exp(i(m*alpha-k*beta)); orthogonality is the admitted FR22 lemma.
The integer conditions below are exact Kronecker deltas, not mode truncations.
For m=m'=0,k,k'>0 the sum characters vanish by beta averaging even though
azimuth integration cannot remove them. A static coefficient is a single real
coefficient, not the doubled A*exp(iθ)+conj(A)*exp(-iθ) representation here.
"""
function paired_product_integral(A, B, m, k, mp, kp, pi_exact)
    terms = real_product_coefficients(2A,2B)
    sum_zero = m+mp == 0 && k+kp == 0
    difference_zero = m-mp == 0 && k-kp == 0
    2pi_exact*((terms.pp+terms.mm)*sum_zero + (terms.pm+terms.mp)*difference_zero)
end

function convention_certificates(; verbose::Bool=true)
    C = certify_fourier_orthogonality()
    @variables ar ai br bi pi_exact phi t w m k l n Phi G M c Sigma cs2
    @variables omega_phi omega_r r mu href a sigmap
    A, B = ar+im*ai, br+im*bi
    phase = m*phi-w*t
    field = ar*cos(phase)-ai*sin(phase)
    da, dt = im*m*A, -im*w*A
    push!(C, certify_symbolic("C-NORM phase exp(i m phi-i omega t)", "G26 (45); shared-forcing-api",
        [Symbolics.derivative(field,phi)-(real(da)*cos(phase)-imag(da)*sin(phase)),
         Symbolics.derivative(field,t)-(real(dt)*cos(phase)-imag(dt)*sin(phase))]; verbose))
    # C-HARM supplies this frequency; only its Newtonian notation is mapped.
    inherited_phase = m*phi-(m*omega_phi+k*omega_r)*t
    push!(C, certify_symbolic("C-NORM consumed C-HARM phase -> FR l=m+k", "C-HARM input; FR22 (9)",
        Symbolics.substitute(inherited_phase,Dict(omega_phi=>n,omega_r=>n,k=>l-m))-
        (m*phi-l*n*t); verbose))
    push!(C, certify_symbolic("C-NORM real amplitudes are twice complex coefficients", "FR24 pi / G26 4pi",
        4pi_exact*real((A/2)*conj(B/2))-pi_exact*real(A*conj(B)); verbose))
    phi_hat = Phi/2
    h_hat = -2phi_hat/c^2
    push!(C, certify_symbolic("C-NORM weak-field h_tt=-Phi_lm in c=1 units", "G26 (52); SW newtonian.jl",
        [h_hat+Phi/c^2, Symbolics.substitute(h_hat,Dict(c=>1))+Phi]; verbose))
    # Consumes FR22's inverse: physical time dt versus M_anom=n*t.
    gain = (2pi_exact)*(2pi_exact/n)/2
    push!(C, certify_symbolic("C-NORM physical-time inverse n/(2pi^2)", "ECC::fr22-eq-fourier",
        [n*gain/(2pi_exact^2)-1, gain/(2pi_exact^2)-1/n]; verbose))
    push!(C, certify_symbolic("C-NORM torque pi reconciles with metric half amplitude", "ECC::fr24-second-order-averages",
        4pi_exact*real((A/2)*conj(im*m*h_hat))/2 +
        pi_exact*real(A*conj(im*m*Phi))/c^2; verbose))

    # Symbolic amplitude identity for each integer-character sector. The proof
    # for arbitrary indices is the two delta conditions in the helper above.
    for (mm,kk) in ((1,0),(2,-3),(1,-1),(0,1),(0,3))
        val = paired_product_integral(A,B,mm,kk,mm,kk,pi_exact)
        push!(C, certify_symbolic("C-NORM 4pi dynamic pair ($mm,$kk)", "G26 (73) two-index extension",
            val-4pi_exact*real(A*conj(B)); verbose))
    end
    cross = [paired_product_integral(A,B,0,1,0,2,pi_exact),
             paired_product_integral(A,B,2,-1,2,1,pi_exact),
             paired_product_integral(A,B,1,2,2,1,pi_exact)]
    push!(C, certify_symbolic("C-NORM cross modes vanish after both averages", "integer-character lemma", cross; verbose))
    push!(C, certify_bool("C-NORM independent half lattice excludes static (0,0)", "C-NORM",
        all(independent_dynamic_pair(mm,kk) != independent_dynamic_pair(-mm,-kk)
            for mm in -5:5 for kk in -5:5 if (mm,kk)!=(0,0)) &&
        !independent_dynamic_pair(0,0),
        "sign partition proof: m!=0 picks positive m; m=0,k!=0 picks positive k; finite regression included"; verbose))

    # Restore G,c,M: L0=GM/c², t0=L0/c, Sigma0=M/L0².
    L0 = G*M/c^2
    t0 = L0/c
    sigma0 = M/L0^2
    push!(C, certify_symbolic("C-NORM G=c=M=1 phase and LI units", "DESIGN v3 units",
        [w*t-(w*t0)*(t/t0), (cs2/c^2)*(Sigma*c^2)-Sigma*cs2,
         (G*(Sigma*c^2)/c^4)*L0-Sigma/sigma0]; verbose))
    FJphys = sigmap*a^4*n^2*(mu/M)^2/href^3
    FJcode = (sigmap/sigma0)*(a/L0)^4*(n*t0)^2*(mu/M)^2/href^3
    push!(C, certify_symbolic("C-NORM torque unit Mc^2 and inverse-time unit 1/t0", "FR24 F_J0, tau_0",
        [FJphys-M*c^2*FJcode,
         FJphys/(mu*n*a^2)-(FJcode/((mu/M)*(n*t0)*(a/L0)^2))/t0]; verbose))
    C
end
end
