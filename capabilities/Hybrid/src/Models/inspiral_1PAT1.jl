# 1PAT1 inspiral forcing (`InspiralModels/1PAT1.m`): first post-adiabatic (1PA) evolution of a
# circular Schwarzschild EMRI. Adds the 1PA forcing `F1[r0]` on top of the 0PA `F0[r0]`:
# `F1 = prefactor · (F2_I(r0) + 2(1-3M/r0)^{3/2}/(1-6M/r0) · F1_flux(r0) · EFLx'(M/r0))`, where
# `F2_I` is the 2SF energy flux (`Re[Exp[interp[Log r0]]]`, order-8) and `EFLx'` is the
# first-law binding-energy derivative over the redshift interpolant. The state `(Ω, φ)` evolves
# with `dΩ/dt = (ν F0 + ν² F1)/M`; integration stops at `r0 ≤ max(6.25, 6.05+4.5ν-1.5ν²)`.

"""
    Forcing1PAT1

1PAT1 forcing model: the 1SF flux interpolants (`F_H`, `F_I`, order-8), the 2SF flux interpolant
(`F2`, complex `{Log r0, Log flux}`, order-8) and the redshift interpolant (`z`, order-8). Build
with [`forcing_1PAT1`](@ref). Evaluate the forcing with [`forcing_F0`](@ref) (0PA) and
[`forcing_F1`](@ref) (1PA).
"""
struct Forcing1PAT1
    FEH::Interp1{Float64,Float64}            # horizon 1SF flux
    FEI::Interp1{Float64,Float64}            # infinity 1SF flux
    F2::Interp1{Float64,ComplexF64}          # 2SF flux, {Log r0, Log flux}
    z::Interp1{Float64,Float64}              # redshift
end

"""
    forcing_1PAT1(sf_datadir) -> Forcing1PAT1

Load the 1SF flux, 2SF flux and redshift tables under `sf_datadir`
(`InspiralModels/sf_data`) and build the order-8 interpolants used by the 1PA forcing.
"""
function forcing_1PAT1(sf_datadir::AbstractString)
    d1 = joinpath(sf_datadir, "1SF_Flux", "Schwarz_Circ")
    d2 = joinpath(sf_datadir, "2SF_Flux", "Schwarz_Circ")
    dz = joinpath(sf_datadir, "1SF_Local_Invariants", "Schwarz_Circ")
    FEH = interp1_order8(read_wl_table(joinpath(d1, "1SFCircShwarzDotEHorizon.m")))
    FEI = interp1_order8(read_wl_table(joinpath(d1, "1SFCircShwarzDotEInf.m")))
    F2 = interp1_order8(read_wl_table_complex(joinpath(d2, "2SFCircShwarzDotEInf.m")))
    z = interp1_order8(read_wl_table(joinpath(dz, "1SFCircShwarzRedshfit.m")))
    return Forcing1PAT1(FEH, FEI, F2, z)
end

"""
    flux_1SF(model::Forcing1PAT1, r0) -> Float64

Total 1SF energy flux `F1(r0) = F_I(r0) + F_H(r0)` (`1PAT1.m`).
"""
flux_1SF(m::Forcing1PAT1, r0::Real) = m.FEI(r0) + m.FEH(r0)

"""
    flux_1SF_deriv(model::Forcing1PAT1, r0) -> Float64

`d/dr0` of the total 1SF energy flux, `F1'(r0) = F_I'(r0) + F_H'(r0)` (the order-8 interpolant
derivatives), needed by the `δm` forcing channel (`1PAT1ea.m` `F1δM`).
"""
flux_1SF_deriv(m::Forcing1PAT1, r0::Real) = derivative(m.FEI, r0) + derivative(m.FEH, r0)

"""
    flux_2SF(model::Forcing1PAT1, r0) -> Float64

2SF energy flux `F2_I(r0) = Re[Exp[interp(Log r0)]]` (`1PAT1.m`), from the complex
`{Log r0, Log flux}` table.
"""
flux_2SF(m::Forcing1PAT1, r0::Real) = real(exp(m.F2(log(r0))))

"""
    forcing_F0(model::Forcing1PAT1, r0) -> Float64

0PA forcing `F0[r0]` (identical closed form to 0PASchwarz; `1PAT1.m`).
"""
forcing_F0(m::Forcing1PAT1, r0::Real) = _forcing_prefactor(r0) * flux_1SF(m, r0)

"""
    forcing_F1(model::Forcing1PAT1, r0) -> Float64

1PA forcing `F1[r0] = prefactor·(F2_I(r0) + 2(1-3M/r0)^{3/2}/(1-6M/r0)·F1_flux(r0)·EFLx'(M/r0))`
(M = 1; `1PAT1.m`). Reuses [`Hybrid.BindingEnergy.E_FLx_deriv`](@ref) for `EFLx'`.
"""
function forcing_F1(m::Forcing1PAT1, r0::Real)
    pref = _forcing_prefactor(r0)
    coupling = (2 * (1 - 3 / r0)^(3 / 2)) / (1 - 6 / r0)
    return pref * (flux_2SF(m, r0) + coupling * flux_1SF(m, r0) * E_FLx_deriv(m.z, 1 / r0))
end

"""
    omega_dot(model::Forcing1PAT1, r0, ν, M) -> Float64

`dΩ/dt = (ν F0(r0) + ν² F1(r0))/M` (`1PAT1.m`).
"""
omega_dot(m::Forcing1PAT1, r0::Real, ν::Real, M::Real) =
    (ν * forcing_F0(m, r0) + ν^2 * forcing_F1(m, r0)) / M

"""
    stop_r0(model::Forcing1PAT1, ν) -> Float64

Terminating radius `max(6.25, 6.05 + 4.5ν - 1.5ν²)` (`1PAT1.m`); ν-dependent.
"""
stop_r0(::Forcing1PAT1, ν::Real) = max(6.25, 6.05 + 4.5 * ν - 1.5 * ν^2)
