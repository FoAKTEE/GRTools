using GrTools.AccretionDisc.ForcingModes: SofteningPolicy

"""
    FR22Setup(T=Float64; softening=nothing)

The printed FR22 v2 setup, in units `a = n = 1`. `values` also retains printed
definitions and convergence targets; storing a target does not verify it.
Repeated paper statements refer to the same values in this object.

`softening` is a shared ForcingModes.SofteningPolicy choice (the landed API's
name for the proposed SofteningRule), or `nothing` for the unresolved `0.3H`
reading. No forcing calculation may interpret `nothing` as a selected policy.
SemimajorScaleHeight, InstantaneousScaleHeight and LocalScaleHeight remain
distinct choices. This setup does not decide which scale height FR22 intended.
"""
struct FR22Setup{S<:Union{Nothing,SofteningPolicy},V<:NamedTuple}
    softening::S
    values::V
end

function FR22Setup(::Type{T}=Float64; softening::Union{Nothing,SofteningPolicy}=nothing) where {T<:AbstractFloat}
    eccentricities = (T(1)/10, T(1)/4, T(1)/2)
    phases = (zero(T), T(1)/4, T(1)/2, T(3)/4)
    mmax, dlmax = 150, 8
    values = (
        p=zero(T), q=zero(T), h_p=T(1)/10, epsilon=T(3)/10,
        n_phi=1024, n_t=1024, r_min=T(1)/10, r_max=T(5), n_r=1000,
        n_r_grid="logarithmically spaced", interp="cubic", n_refine=2, n_phase_corr=1,
        m_max=mmax, lm_max=dlmax, e_list=eccentricities, e_list_results=eccentricities,
        snapshots_text=phases, snapshots_caption=phases[1:3],
        wavelet_interval=(1, 800), wavelet_orbits=5,
        zz22_mass=(T(3)/10^6, T(3)/10^3), convergence_phase=zero(T),
        conv_m_fixed_lm=(dlmax, 10), conv_lm_fixed_m=(mmax, 1), conv_adopted=(mmax, dlmax),
        conv_m_claim=(80, 2, 10), conv_m_range=(80, 100), conv_lm_claim=4,
        conv_lm_visual=(4, dlmax), M_th_def=3, H_def=raw"H = c_{s,\text{iso}}/\Omega_K",
        colour_range=(-T(2), T(2)), column_titles=phases)
    FR22Setup(softening, values)
end
