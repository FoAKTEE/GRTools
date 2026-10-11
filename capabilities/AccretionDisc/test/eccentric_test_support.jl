using Test
using AccretionDisc.SpiralWaves
using AccretionDisc.ForcingModes
using Spacetime: EquatorialOrbit, orbit_frequency, kerr_equatorial_potentials,
                 kerr_equatorial_potentials_dr

"A fixed circular sample, also used to capture the pre-context bit-pattern golden."
function eccentric_legacy_sample()
    bg = DiscBackground(χ=0.4, aEMRI=20.0, h0=0.06, 𝔭=1.5, 𝔮=0.75)
    values = Float64[]
    for r in (9.0, 17.0, 23.0, 61.0)
        append!(values, (eps_profile(bg,r), csound(bg,r), pressure(bg,r), gamma1(bg,r),
                        values_tuple(rotation_profile(bg,r))..., dOmega_dr(bg,r)))
        h = singular_field_mode(bg,3,r; N=128)
        C,S = master_matrices(bg,3,r;hmodes=h,rc=20.0)
        for z in (vec(C)...,S...,delta_u_phi(bg,3,r,(0.3+0.1im,-0.2+0.4im);hmodes=h,rc=20.0))
            append!(values, (real(z),imag(z)))
        end
    end
    reinterpret(UInt64,values)
end

values_tuple(x::NamedTuple) = Tuple(x)

"Complex covariant source jets; odd radial components are deliberately nonzero."
function eccentric_test_source(r)
    z = complex(one(r),one(r)/3) / r^2
    (; htt=z, htt_r=-2z/r, htphi=z*r/7, htphi_r=-z/7,
       hphiphi=z*r^2/13, hphiphi_r=zero(z), hrr=z/2, hzz=-z/5,
       htr=z*(1+2im)/11, hrphi=z*r*(2-im)/17)
end
