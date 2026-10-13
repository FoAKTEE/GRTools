# Conserved-baryon angular-momentum flux, with OPTION A perturbations.
# Stress origin: derivation/src/diagnostics.jl:163–172, emitted by li_emission.jl.
# No adiabatic U/V/W tensor enters this file.

_diagnostic_stress(::Val{LocallyIsothermalOptionA},c,r,y,geo;hmodes,landau) =
    linear_stress(c,r,y;hmodes,landau).dT

"""
    reynolds_flux(mode; measure, enthalpy, dur, du_phi)

Connected quadratic Reynolds contribution to the physical angular-momentum
current: `4π sqrt(-g) (epsilon+p) Re(delta u^r conj(delta u_phi))`.
`du_phi` is COVARIANT and includes `g_phi,a delta u^a + h_phi,a u^a`; with a Landau
shift the caller passes its partner continued in omega_d (`li_flux_decomposition`).
Inputs are half amplitudes, as in C-NORM. The definition is independent of the
thermodynamic closure; `li_flux_decomposition` substitutes OPTION A solutions.

Newtonian limit: measure=r, enthalpy=Sigma, dur=du_r_FR/2 and
du_phi=r*du_phi_FR/2, giving FR24 v1 eq:F_jml exactly:
`π r² Sigma Re(du_r_FR conj(du_phi_FR))` (2_linear_framework.tex:268–287).
Its source-free LI transport law is F_J/c_T² constant, not F_J constant.
"""
function reynolds_flux(mode::HarmonicMode;measure,enthalpy,dur,du_phi)
    reality_pair_weight(mode)*measure*enthalpy*real(dur*conj(du_phi))
end

"""
    li_flux_decomposition(context,r,y;hmodes,landau=true)

Expand sqrt(-g) T^r_phi = sqrt(-g) (epsilon+p) u^r u_phi about u^r=0.
The terms quadratic in the first-order fields are

    quadratic = sqrt(-g) < (epsilon+p) du^r du_phi
                         + u_phi [d(epsilon+p)+(epsilon+p) dlog(sqrt(-g))] du^r >.

Subtract the transport of BACKGROUND angular momentum by the conserved baryon
current `h u_phi Mdot`, where `h=(epsilon+p)/rho` (owner decision 2026-10-13,
gr-li-o-amf-definition; certificate: derivation/scripts/certify_li_amf_definition.jl):

    mean_transport = sqrt(-g) u_phi <[d epsilon+(epsilon+p) dlog(sqrt(-g))] du^r>.

Relation (19) gives `h d rho = d epsilon + xi^r (epsilon' - h rho')`;
the displacement term averages to zero against `du^r`. Thus
`wave=quadratic-mean_transport` is (72a): the connected Reynolds product PLUS
`sqrt(-g) u_phi <dp du^r>`. Evaluate it directly to avoid subtractive cancellation.
`dlog(sqrt(-g))=g^ab h_ab/2`.
Any unknown second-order mean radial velocity contributes the same
sqrt(-g)(epsilon+p) u_phi <u2^r> to both terms and cancels. The reported
`quadratic`/`mean_transport` contain only the known first-order bilinears;
they do not solve the second-order mean flow. Pressure gives no mixed r-phi
Kronecker term. Metric/velocity cross terms survive in du_phi.

This uses the same flux expression as the adiabatic branch, with OPTION A
perturbations and zero bath four-force. Finite-GR wave/mean balance
does not imply a conserved wave flux or identify its edge jump as total torque.
With the Landau shift on, the conjugated du_phi is continued in omega_d as in the
adiabatic and frozen circular AMF (one convention for every closure); du^r and the
mean-transport factor carry no omega_d.
"""
function li_flux_decomposition(c::ModeContext,r,y;hmodes,landau::Bool=true)
    c.closure==LocallyIsothermalOptionA || throw(ArgumentError("LI flux requires OPTION A context"))
    inp=master_inputs(c,r;hmodes,landau,y)
    f=_context_reconstruction(Val(c.closure),inp); partner=_continued_partner(c,inp)
    geo=_diagnostic_geometry(c,r,hmodes)
    enthalpy=eps_profile(c.disc,r)+pressure(c.disc,r)
    u_phi=sum(geo.g[3,:] .* geo.u)
    dlogmeasure=sum(geo.gi .* geo.metric)/2
    wave=reynolds_flux(c.mode;measure=geo.measure,enthalpy,dur=f.dur,du_phi=partner.du_phi)+
        reality_pair_weight(c.mode)*geo.measure*u_phi*real(f.dp*conj(f.dur))
    mean_transport=reality_pair_weight(c.mode)*geo.measure*u_phi*
        real((f.deps+enthalpy*dlogmeasure)*conj(f.dur))
    (;wave,mean_transport,quadratic=wave+mean_transport)
end

_wave_flux(::Val{LocallyIsothermalOptionA},c,r,y;hmodes,landau) =
    li_flux_decomposition(c,r,y;hmodes,landau).wave

"Wave/mean bilinears at each solved radius; closure is checked by the point method."
li_flux_decomposition(sol::EccentricModeSolution) =
    [li_flux_decomposition(sol.context,r,view(sol.y,:,i);hmodes=sol.source(r),landau=sol.params.landau)
     for (i,r) in enumerate(sol.r)]
