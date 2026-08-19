# SpiralWaves/newtonian.jl — the NEWTONIAN-LIMIT pipeline of arXiv:2608.17003, i.e. the denominator
# T_Newt of Fig. 5 (`fig:TGR_vs_TNewt`, M:927-932: "ratio of the relativistic torque to the Newtonian
# torque ... as a function of the spin of the SMBH").
#
# It is the SAME numerical scheme (implicit midpoint (64) on the App.-C grid, outgoing-wave rows
# (66)-(68), the (75) split) driven by the NEWTONIAN master equations (50a)-(50g) on a NEWTONIAN
# background with a NEWTONIAN (Plummer-softened) forcing potential. Every piece is the M/r → 0 limit
# of the corresponding relativistic piece and is gated against it in the tests, so this file
# introduces no independent physics — only the limit:
#
#   background   Ω_K = r^{−3/2}, Ω•,N = a•^{−3/2} (13 with χ = 0), Σ = ε₀(r/a•)^{−𝔭} (12a),
#                c_s = h₀a•Ω•,N(r/a•)^{−𝔮/2} (12c), p = c_s²Σ (12b), Γ₁ = 1 + c_s² (22), and the
#                rotation law from the NEWTONIAN radial force balance rΩ² = M/r² + p′/Σ, i.e.
#                Ω² = 1/r³ − (𝔭+𝔮)c_s²/r², which is exactly the M/r → 0 limit of (11)
#                (ν′ → 1/r², ψ → ln r, ω̃ → 0, v → rΩ);
#   epicycle     κ² = 4Ω² + 2rΩΩ′ (the fluid rotation law, `_kappa2_fluid` — the same combination the
#                relativistic C reduces to; = Ω² exactly for 𝔭 + 𝔮 = 0, the Keplerian κ = Ω);
#   forcing      Φ(r,Ψ) = −m•/√(ε_soft² + r² + a•² − 2ra•cosΨ), the Plummer-softened (59) Newtonian
#                potential of the secondary — the M/r → 0 limit of the singular field (52), whose
#                s0² (B1) tends to the flat coordinate distance (P_r, P_t → 1) and whose h_tt → 2m•/s0
#                gives Φ = −h_tt/2. Its cosine modes (60) are computed by the SAME spectrally accurate
#                periodic-trapezoid projection as `source.jl` (`_kernel_integrals` on a
#                `SingularFieldCoeffs` whose quadratic form is the flat one), never a second integrator;
#   diagnostics  F_J (72a) with √−g → r, u_φ → r²Ω, ε + p → Σ and the COVARIANT δû_φ = r × (50a)
#                (ERRATUM E-M1), and the torque density (72b), whose Newtonian limit is the
#                Goldreich–Tremaine form √−g⟨𝒯⟩ = −r⟨δΣ ∂_φΦ⟩ = −4π m r Φ̂ Im[δΣ̂] (the h_tt piece of
#                ½δ_ξT^{αβ}∂_φh_αβ; the ξ-independent piece ∝ ΣΦ∂_φΦ integrates to zero over φ).
#
# χ IS IGNORED — deliberately. `T_Newt` is a property of the disc (𝔭, 𝔮, h₀, a•, b), not of the hole:
# every function here uses Ω_K = r^{−3/2} and the χ = 0 fiducial inner edge r_in = r_ISCO(0) = 6M of
# (C1), so `newtonian_total_torque(bg, ...)` returns the same bits for every `bg.χ`. That is what
# makes Fig. 5 a measurement of relativity: the whole χ-dependence of T_GR/T_Newt is in the numerator
# (VERIFIED: the measured T_GR(χ=+0.9)/T_GR(χ=−0.9) at a• = 12M, h₀ = 0.05 is 1.578 / 2.422 / 0.838
# for (𝔭,𝔮) = (0,0) / (1.5,0) / (−1.5,0), against 1.617 / 2.67 / 0.843 read off the published Fig. 5
# ratios — a χ-DEPENDENT denominator cannot produce that).
#
# GRID [PRELIMINARY, recorded]: the Newtonian solve runs on `build_grid` of the χ = 0 background with
# the NEWTONIAN corotation radius passed explicitly (`rc`), so the sinh block (C4)-(C5) is centred on
# the Newtonian pole. Node PLACEMENT and the (C2)-(C3) domain trimming still use the χ = 0
# relativistic WKB profile of grid.jl; that is a discretisation choice, not a model choice (the
# wavenumbers differ by the O(M/r) epicyclic and redshift factors), and it is gated by the resolution
# study in the tests, where T changes by < 1e-3 between (nL,nR,nC) = (30,30,1500) and (120,120,12000).
# Type-generic `{T}`; BigFloat inherits. Additive (EGD::sw-newtonian-limit).

# ---------------------------------------------------------------------------------------------
# Newtonian background (the M/r → 0 limit of background.jl)
# ---------------------------------------------------------------------------------------------

"""
    newtonian_Omega_emri(bg) -> Ω•,N

Newtonian orbital frequency of the secondary, Ω•,N = a•^{−3/2} (M = 1) — (13) with χ = 0, i.e. the
M/a• → 0 limit of [`Omega_emri`](@ref). The forcing frequency of mode m is ω = mΩ•,N (46).
"""
newtonian_Omega_emri(bg::DiscBackground{T}) where {T} = bg.aEMRI^(-3 // 2)

"""
    newtonian_csound(bg, r) -> c_s

Sound speed (12c) evaluated with the NEWTONIAN Ω•,N: c_s = h₀ a• Ω•,N (r/a•)^{−𝔮/2} = h₀ a•^{−1/2}
(r/a•)^{−𝔮/2}. Identical to [`csound`](@ref) at χ = 0; χ-independent by construction.
"""
newtonian_csound(bg::DiscBackground, r) =
    bg.h0 * bg.aEMRI * newtonian_Omega_emri(bg) * (r / bg.aEMRI)^(-bg.𝔮 / 2)

"Height-integrated pressure (12b) of the Newtonian disc: p = c_s²Σ with [`newtonian_csound`](@ref)."
newtonian_pressure(bg::DiscBackground, r) = newtonian_csound(bg, r)^2 * eps_profile(bg, r)

"Adiabatic index (22) of the Newtonian disc: Γ₁ = (Σ+p)c_s²/p = 1 + c_s² ([`newtonian_csound`](@ref))."
newtonian_gamma1(bg::DiscBackground, r) = 1 + newtonian_csound(bg, r)^2

"""
    newtonian_rotation(bg, r) -> (Ω, dΩ, κ2)

Equilibrium rotation of the NEWTONIAN disc from the radial force balance rΩ² = M/r² + p′/Σ with the
power laws (12): p′/Σ = −(𝔭+𝔮)c_s²/r, so

    Ω² = 1/r³ − (𝔭+𝔮) c_s²/r²,   (Ω²)′ = −3/r⁴ + (𝔭+𝔮)(𝔮+2) c_s²/r³,   Ω′ = (Ω²)′/(2Ω),

and κ² = 4Ω² + 2rΩΩ′ = 1/r³ + (𝔭+𝔮)(𝔮−2)c_s²/r² ([`_kappa2_fluid`](@ref) — one implementation; the
closed form is asserted in the tests). Reduces to the Keplerian Ω = r^{−3/2}, κ = Ω for 𝔭 + 𝔮 = 0.
This is the M/r → 0 limit of [`rotation_profile`](@ref): the relativistic pressure correction carries
(1/ν′ − 1/ψ′) = r(r−3) where the Newtonian one carries r, so the two agree to O(M/r) × O(h₀²).
"""
function newtonian_rotation(bg::DiscBackground{T}, r) where {T}
    S = promote_type(T, typeof(r))
    rr = convert(S, r)
    cs2 = newtonian_csound(bg, rr)^2
    pq = bg.𝔭 + bg.𝔮
    Ω2 = 1 / rr^3 - pq * cs2 / rr^2
    Ω2 > 0 || throw(DomainError(rr, "newtonian_rotation: Ω² = $Ω2 ≤ 0 (pressure gradient too strong)"))
    Ω = sqrt(Ω2)
    dΩ2 = -3 / rr^4 + pq * (bg.𝔮 + 2) * cs2 / rr^3
    dΩ = dΩ2 / (2 * Ω)
    (Ω = Ω, dΩ = dΩ, κ2 = _kappa2_fluid(Ω, rr, dΩ))
end

"""
    newtonian_omega_d(bg, m, r) -> ω_d

Doppler-shifted forcing frequency (51) of the Newtonian disc, ω_d = m[Ω•,N − Ω_N(r)]
([`newtonian_rotation`](@ref)); the M/r → 0 limit of [`omega_d`](@ref).
"""
newtonian_omega_d(bg::DiscBackground, m::Integer, r) =
    m * (newtonian_Omega_emri(bg) - newtonian_rotation(bg, r).Ω)

"""
    newtonian_corotation_radius(bg) -> r_c

Corotation radius (61) of the Newtonian disc: the root of Ω_N(r) = Ω•,N on [a•/2, 2a•]
(`NumericsCore.brent`; Ω_N is monotonically decreasing). Exactly a• for 𝔭 + 𝔮 = 0 and slightly
inside/outside it for a pressure-supported/pressure-drained disc — the same statement as
[`corotation_radius`](@ref), whose relativistic root differs from this one at O(h₀²M/r), which is
comparable to the width σ_wc = 10⁻³r_c of the corotation block (C5); hence the Newtonian solve builds
its grid with THIS radius (`build_grid(...; rc)`). Independent of m.
"""
function newtonian_corotation_radius(bg::DiscBackground{T}) where {T}
    f(r) = newtonian_Omega_emri(bg) - newtonian_rotation(bg, r).Ω
    brent(f, bg.aEMRI / 2, 2 * bg.aEMRI)
end

# Landau pole displacement (63) with the Newtonian Ω_K(r_c) = r_c^{−3/2}.
_newtonian_landau_shift(bg::DiscBackground, r, rc) =
    bg.εp * rc^(-3 // 2) * exp(-(r - rc)^2 / (2 * bg.σε^2))

function _newtonian_omega_d_shifted(bg::DiscBackground, m::Integer, r, landau::Bool, rc)
    ωd = newtonian_omega_d(bg, m, r)
    landau ? Complex(ωd, _newtonian_landau_shift(bg, r, rc)) : Complex(ωd, zero(ωd))
end

# ---------------------------------------------------------------------------------------------
# the Newtonian (Plummer-softened) forcing potential and its cosine modes (60)
# ---------------------------------------------------------------------------------------------

"""
    newtonian_potential_coeffs(bg) -> SingularFieldCoeffs

The worldline data of the NEWTONIAN forcing potential in the SAME container the relativistic kernel
uses ([`SingularFieldCoeffs`](@ref)), so that its quadrature (`_kernel_integrals`, `_trapezoid_nodes`,
`_strip_halfwidth`) is reused verbatim: the flat squared separation

    s_N²(r,Ψ) = ε_soft² + (r cosΨ − a•)² + r² sin²Ψ = ε_soft² + a•² + r² − 2 a• r cosΨ

is the P_r = P_t = 1 member of the family, i.e. (A0, A2, B1, C2) = (ε_soft² + a•², 1, −2a•, 0), which
is exactly what (B1) tends to as M/a• → 0 (P_r = a•²/Δ_a → 1, P_t = A2 − C2 → 1) — asserted in the
tests. `Ktt = 2` reproduces the Newtonian limit h_tt = 2m•/s0 of (52), so that Φ = −h_tt/2 = −m•/s_N.
"""
function newtonian_potential_coeffs(bg::DiscBackground{T}) where {T}
    a = bg.aEMRI
    ε = softening_length(bg)
    ΩE = newtonian_Omega_emri(bg)
    SingularFieldCoeffs{T}(a, zero(T), ε,
                           ε^2 + a^2, one(T), -2 * a, zero(T),      # A0, A2, B1, C2
                           one(T), one(T),                          # P_r, P_t
                           T(2), zero(T), zero(T), zero(T), zero(T), # K_tt = 2, K_tφ = K_φφ = K_rr = K_zz = 0
                           one(T), ΩE, -one(T), a^2 * ΩE)           # u^t, u^φ, v_t, v_φ (Newtonian)
end

"""
    newtonian_potential_mode(bg, m, r; N = nothing, mEMRI = 1) -> (Φ, Φ_r, N)

The cosine Fourier mode (60) of the Newtonian forcing potential and its analytic r-derivative,

    Φ̂_m(r) = −m• I_m(r),   Φ̂_m′(r) = −m• I_m′(r),   I_m = (1/2π)∫₀^{2π} cos(mΨ)/s_N(r,Ψ) dΨ,

by the periodic trapezoid rule of `source.jl` (`_kernel_integrals`; `N` defaults to
`_trapezoid_nodes`, which sizes the rule from the analyticity strip of 1/s_N and therefore from the
softening b h₀). REAL by construction (the kernel is even in Ψ), and negative for m ≥ 0.
"""
function newtonian_potential_mode(bg::DiscBackground{T}, m::Integer, r;
                                  N::Union{Nothing,Integer} = nothing, mEMRI = 1) where {T}
    S = promote_type(T, typeof(r), typeof(mEMRI))
    _newtonian_potential_mode(newtonian_potential_coeffs(bg), m, convert(S, r), N, convert(S, mEMRI))
end

function _newtonian_potential_mode(K::SingularFieldCoeffs, m::Integer, r::S,
                                   N::Union{Nothing,Integer}, mE::S) where {S<:Real}
    m ≥ 0 || throw(DomainError(m, "newtonian_potential_mode: need m ≥ 0"))
    NN = N === nothing ? _trapezoid_nodes(K, m, r) : Int(N)
    I, Ir = _kernel_integrals(K, m, r, NN)
    (Φ = -mE * I, Φ_r = -mE * Ir, N = NN)
end

"""
    newtonian_potential_modes(bg, m, rgrid; threads = false, N = nothing, mEMRI = 1) -> Vector

[`newtonian_potential_mode`](@ref) at every radius of `rgrid` (the counterpart of
[`source_modes`](@ref); `threads = true` distributes the radii, bitwise identically).
"""
function newtonian_potential_modes(bg::DiscBackground{T}, m::Integer, rgrid::AbstractVector;
                                   threads::Bool = false, N::Union{Nothing,Integer} = nothing,
                                   mEMRI = 1) where {T}
    S = promote_type(T, eltype(rgrid), typeof(mEMRI))
    K = newtonian_potential_coeffs(bg)
    mE = convert(S, mEMRI)
    proto = _newtonian_potential_mode(K, m, convert(S, first(rgrid)), 2, mE)
    out = Vector{typeof(proto)}(undef, length(rgrid))
    if threads
        Threads.@threads for i in eachindex(rgrid)
            @inbounds out[i] = _newtonian_potential_mode(K, m, convert(S, rgrid[i]), N, mE)
        end
    else
        for i in eachindex(rgrid)
            @inbounds out[i] = _newtonian_potential_mode(K, m, convert(S, rgrid[i]), N, mE)
        end
    end
    out
end

# ---------------------------------------------------------------------------------------------
# the Newtonian master system (50) on the Newtonian background
# ---------------------------------------------------------------------------------------------

"""
    newtonian_master_matrices(bg, m, r; hmode, landau = true, rc) -> (C, S)

The master system (50b)-(50g) at radius `r` on the NEWTONIAN background of
[`newtonian_rotation`](@ref), with `hmode` = the potential mode NamedTuple `(Φ, Φ_r, …)` of
[`newtonian_potential_mode`](@ref) and the Landau-regularised ω_d (62)-(63) about `rc`
([`newtonian_corotation_radius`](@ref)). The (50) algebra itself is `_newtonian_entries` of
`master.jl` — the SAME code path as [`master_matrices_newtonian`](@ref), which evaluates it on the
relativistic background instead. Same return layout as [`master_matrices`](@ref).
"""
function newtonian_master_matrices(bg::DiscBackground{T}, m::Integer, r; hmode,
                                   landau::Bool = true, rc = newtonian_corotation_radius(bg)) where {T}
    S = promote_type(T, typeof(r))
    rr = convert(S, r)
    rot = newtonian_rotation(bg, rr)
    Σ = eps_profile(bg, rr); p = newtonian_pressure(bg, rr)
    _newtonian_entries(m, rr, rot.Ω, rot.κ2, Σ, -bg.𝔭 * Σ / rr, p, -(bg.𝔭 + bg.𝔮) * p / rr,
                       newtonian_gamma1(bg, rr), _newtonian_omega_d_shifted(bg, m, rr, landau, rc),
                       convert(S, hmode.Φ), convert(S, hmode.Φ_r))
end

"""
    newtonian_delta_u_phi(bg, m, r, y, Φ; landau = true, rc) -> δû_φ

The COVARIANT azimuthal velocity perturbation of the Newtonian disc: δû_φ = r × (50a) =
m(δÛ + Φ̂)/ω_d + r(−iδû^r)κ²/(2ω_dΩ) in the master variables y = (δÛ, −iδû^r). ERRATUM E-M1: (50a) as
printed is the ORTHONORMAL component; the AMF (72a) consumes the covariant one (with u_φ = r²Ω), which
is r times it — exactly what the relativistic [`delta_u_phi`](@ref) reduces to.
"""
function newtonian_delta_u_phi(bg::DiscBackground{T}, m::Integer, r, y, Φ;
                               landau::Bool = true, rc = newtonian_corotation_radius(bg)) where {T}
    S = promote_type(T, typeof(r))
    rr = convert(S, r)
    rot = newtonian_rotation(bg, rr)
    ωd = _newtonian_omega_d_shifted(bg, m, rr, landau, rc)
    rr * _newtonian_delta_u_phi(m, rr, rot.Ω, rot.κ2, ωd,
                                (convert(Complex{S}, y[1]), convert(Complex{S}, y[2])), convert(S, Φ))
end

"""
    newtonian_delta_sigma(bg, m, r, y; landau = true, rc) -> δΣ̂

Eulerian surface-density perturbation of the Newtonian disc, from the adiabatic relation (20)-(21)
Δp = Γ₁(p/Σ)ΔΣ with Δ = δ + ξ^r∂_r and the Lagrangian displacement ξ̂^r = −y₂/ω_d (from
δû^r = −iω_dξ̂^r, (16)):

    δΣ̂ = (Σ/(Γ₁p))(δp̂ + ξ̂^r p′) − ξ̂^r Σ′,     δp̂ = Σ y₁   ((48) with ε + p → Σ).

Consumed by [`newtonian_torque_density_mode`](@ref); the AMF/torque balance (71) is what gates it.
"""
function newtonian_delta_sigma(bg::DiscBackground{T}, m::Integer, r, y;
                               landau::Bool = true, rc = newtonian_corotation_radius(bg)) where {T}
    S = promote_type(T, typeof(r))
    rr = convert(S, r)
    Σ = eps_profile(bg, rr); p = newtonian_pressure(bg, rr); Γ1 = newtonian_gamma1(bg, rr)
    dΣ = -bg.𝔭 * Σ / rr; dp = -(bg.𝔭 + bg.𝔮) * p / rr
    ωd = _newtonian_omega_d_shifted(bg, m, rr, landau, rc)
    y1 = convert(Complex{S}, y[1]); y2 = convert(Complex{S}, y[2])
    ξr = -y2 / ωd
    (Σ / (Γ1 * p)) * (Σ * y1 + ξr * dp) - ξr * dΣ
end

# ---------------------------------------------------------------------------------------------
# the Newtonian solve: the SAME scheme (64) + (66)-(68) as solver.jl
# ---------------------------------------------------------------------------------------------

"""
    newtonian_grid(bg, m; kwargs...) -> SpiralGrid

The App.-C grid for the Newtonian solve: [`build_grid`](@ref) of the χ = 0 background with the same
(𝔭, 𝔮, h₀, a•, b, ε_p, σ_ε, ε₀) and the NEWTONIAN corotation radius
([`newtonian_corotation_radius`](@ref)) passed as `rc`, so the sinh block (C4)-(C5) is centred on the
Newtonian pole and the fiducial inner edge (C1) is r_ISCO(0) = 6M for every `bg.χ`. `kwargs` are
forwarded to [`build_grid`](@ref) (`nL`, `nR`, `nC`, `n_anc`, `extent`, …). See the file header for
the recorded [PRELIMINARY] choice of keeping the χ = 0 relativistic WKB profile for node placement.
"""
function newtonian_grid(bg::DiscBackground{T}, m::Integer; kwargs...) where {T}
    bg0 = DiscBackground{T}(; χ = zero(T), aEMRI = bg.aEMRI, h0 = bg.h0, 𝔭 = bg.𝔭, 𝔮 = bg.𝔮,
                            b = bg.b, εp = bg.εp, σε = bg.σε, ε0 = bg.ε0)
    build_grid(bg0, m; rc = newtonian_corotation_radius(bg), kwargs...)
end

"""
    newtonian_solve_mode(bg, m; grid = newtonian_grid(bg, m), landau = true, bc = :outgoing,
                         linsolve = :auto, mEMRI = 1, hmodes_N = nothing, threads = false)
        -> ModeSolution

Solve the NEWTONIAN master system (50) for azimuthal number `m` with the same implicit-midpoint
scheme (64) and outgoing-wave boundary rows (66)-(68) as [`solve_mode`](@ref) — the matrices at the
cell midpoints come from [`newtonian_master_matrices`](@ref) and the boundary branch from
[`outgoing_branch`](@ref) evaluated with the Newtonian ω_d. Returns the same [`ModeSolution`](@ref)
container (so `cumulative_torque`, `corotation_split`, `m_cut` … apply unchanged); `params.rc` is the
NEWTONIAN corotation radius, which the Newtonian diagnostics below read back.
"""
function newtonian_solve_mode(bg::DiscBackground{T}, m::Integer;
                              grid::SpiralGrid{T} = newtonian_grid(bg, m),
                              landau::Bool = true, bc::Symbol = :outgoing, linsolve::Symbol = :auto,
                              mEMRI = 1, hmodes_N::Union{Nothing,Integer} = nothing,
                              threads::Bool = false) where {T}
    m ≥ 1 || throw(ArgumentError("newtonian_solve_mode: m ≥ 1 required"))
    t0 = time()
    r = grid.r
    N = length(r)
    rc = grid.r_c
    rmid = [(r[i] + r[i+1]) / 2 for i in 1:(N - 1)]
    hmid = newtonian_potential_modes(bg, m, rmid; threads = threads, N = hmodes_N, mEMRI = mEMRI)
    Cmid = Vector{Matrix{Complex{T}}}(undef, N - 1)
    Smid = Vector{Vector{Complex{T}}}(undef, N - 1)
    for i in 1:(N - 1)
        Cmid[i], Smid[i] = newtonian_master_matrices(bg, m, rmid[i]; hmode = hmid[i], landau = landau, rc = rc)
    end
    Cin, _ = newtonian_master_matrices(bg, m, r[1];
                                       hmode = newtonian_potential_mode(bg, m, r[1]; N = hmodes_N, mEMRI = mEMRI),
                                       landau = landau, rc = rc)
    Cout, _ = newtonian_master_matrices(bg, m, r[N];
                                        hmode = newtonian_potential_mode(bg, m, r[N]; N = hmodes_N, mEMRI = mEMRI),
                                        landau = landau, rc = rc)
    k_in, Y_in, W_in = outgoing_branch(Cin, newtonian_omega_d(bg, m, r[1]); outward = false, bc = bc)
    k_out, Y_out, W_out = outgoing_branch(Cout, newtonian_omega_d(bg, m, r[N]); outward = true, bc = bc)
    y = solve_midpoint(r, Cmid, Smid, W_in, zero(Complex{T}), W_out, zero(Complex{T}); linsolve = linsolve)
    ra, rr = midpoint_residual(r, Cmid, Smid, y)
    bcres = abs(W_in[1] * y[1, 1] + W_in[2] * y[2, 1]) + abs(W_out[1] * y[1, N] + W_out[2] * y[2, N])
    params = (landau = landau, bc = bc, linsolve = linsolve, mEMRI = convert(T, mEMRI),
              hmodes_N = hmodes_N === nothing ? nothing : Int(hmodes_N), threads = threads, rc = rc)
    ModeSolution{T}(Int(m), collect(r), y, k_in, k_out, Y_in, Y_out, ra, rr, bcres, params, time() - t0)
end

# ---------------------------------------------------------------------------------------------
# the Newtonian diagnostics: (72a), (72b) in the M/r → 0 limit
# ---------------------------------------------------------------------------------------------

"""
    newtonian_amf_mode(bg, m, r, y, Φ; landau = true, rc) -> F_J,m(r)

Angular-momentum flux (72a) of one Newtonian mode: with √−g → r, u_φ → r²Ω, ε + p → Σ and
δp̂ = Σ δÛ,

    F_J,m = 4π r Re[ r²Ω (Σ y₁) (i y₂)* + Σ (i y₂) δû_φ* ],   δû_φ = [`newtonian_delta_u_phi`](@ref),

using the single-mode Fourier weight 4π of (73) (m ≥ 1). Real by construction.
"""
function newtonian_amf_mode(bg::DiscBackground{T}, m::Integer, r, y, Φ;
                            landau::Bool = true, rc = newtonian_corotation_radius(bg)) where {T}
    S = promote_type(T, typeof(r))
    rr = convert(S, r)
    rot = newtonian_rotation(bg, rr)
    Σ = eps_profile(bg, rr)
    y1 = convert(Complex{S}, y[1]); y2 = convert(Complex{S}, y[2])
    δur = im * y2                                             # δû^r = i·(−iδû^r)
    δuφ = newtonian_delta_u_phi(bg, m, rr, (y1, y2), Φ; landau = landau, rc = rc)
    4 * S(π) * rr * (rr^2 * rot.Ω * real(Σ * y1 * conj(δur)) + Σ * real(δur * conj(δuφ)))
end

"""
    newtonian_torque_density_mode(bg, m, r, y, Φ; landau = true, rc) -> √−g⟨𝒯⟩_m(r)

Torque density (72b) of one Newtonian mode — the Goldreich–Tremaine form

    √−g⟨𝒯⟩_m = −r ⟨δΣ ∂_φΦ⟩_m = −4π m r Φ̂ Im[δΣ̂],

the M/r → 0 limit of ½δ_ξT^{αβ}∂_φh_αβ (only h_tt = −2Φ survives, T^{tt} → Σ, and the ξ-independent
piece of (40)-(42) contributes ∝ ∮ΣΦ∂_φΦ dφ = 0). Φ̂ is REAL ([`newtonian_potential_mode`](@ref)), so
only Im δΣ̂ — the part of the density response that lags the forcing — carries torque. Sign convention:
the torque ON THE DISC (M:898), as in [`torque_density_mode`](@ref).
"""
function newtonian_torque_density_mode(bg::DiscBackground{T}, m::Integer, r, y, Φ;
                                       landau::Bool = true, rc = newtonian_corotation_radius(bg)) where {T}
    S = promote_type(T, typeof(r))
    rr = convert(S, r)
    δΣ = newtonian_delta_sigma(bg, m, rr, y; landau = landau, rc = rc)
    -4 * S(π) * m * rr * convert(S, Φ) * imag(δΣ)
end

# per-node potential modes with the settings the solve used (mirrors `_node_modes` of diagnostics.jl)
_newtonian_node_modes(bg::DiscBackground, sol::ModeSolution, rs) =
    newtonian_potential_modes(bg, sol.m, rs; threads = sol.params.threads, N = sol.params.hmodes_N,
                              mEMRI = sol.params.mEMRI)

"""
    newtonian_amf(bg, sol; hmodes = <potential modes at the nodes>) -> Vector{T}
    newtonian_amf(bg, sol, idx; hmodes = <potential modes at `sol.r[idx]`>) -> Vector{T}

[`newtonian_amf_mode`](@ref) at every node of a Newtonian [`ModeSolution`](@ref), or only at the node
indices `idx` (the (75) reduction needs four of them) — the counterpart of [`amf`](@ref).
"""
function newtonian_amf(bg::DiscBackground{T}, sol::ModeSolution{T};
                       hmodes::AbstractVector = _newtonian_node_modes(bg, sol, sol.r)) where {T}
    p = sol.params
    [newtonian_amf_mode(bg, sol.m, sol.r[i], (sol.y[1, i], sol.y[2, i]), hmodes[i].Φ;
                        landau = p.landau, rc = p.rc) for i in eachindex(sol.r)]
end

function newtonian_amf(bg::DiscBackground{T}, sol::ModeSolution{T}, idx::AbstractVector{<:Integer};
                       hmodes::AbstractVector = _newtonian_node_modes(bg, sol, sol.r[idx])) where {T}
    p = sol.params
    [newtonian_amf_mode(bg, sol.m, sol.r[i], (sol.y[1, i], sol.y[2, i]), hmodes[j].Φ;
                        landau = p.landau, rc = p.rc) for (j, i) in enumerate(idx)]
end

"""
    newtonian_torque_density(bg, sol; hmodes = <potential modes at the nodes>) -> Vector{T}

[`newtonian_torque_density_mode`](@ref) at every node of a Newtonian [`ModeSolution`](@ref) — the
counterpart of [`torque_density`](@ref); feed it to `cumulative_torque` / `corotation_split`, which
are model-independent (they consume `sol.r` and the arrays only).
"""
function newtonian_torque_density(bg::DiscBackground{T}, sol::ModeSolution{T};
                                  hmodes::AbstractVector = _newtonian_node_modes(bg, sol, sol.r)) where {T}
    p = sol.params
    [newtonian_torque_density_mode(bg, sol.m, sol.r[i], (sol.y[1, i], sol.y[2, i]), hmodes[i].Φ;
                                   landau = p.landau, rc = p.rc) for i in eachindex(sol.r)]
end

"""
    newtonian_balance_residual(bg, sol; flux = newtonian_amf(bg, sol)) -> (r, lhs, rhs, rel)

The balance law (71) dF_J/dr = d𝔗/dr for a Newtonian mode, cell by cell in the discretisation's own
midpoint form — the Newtonian counterpart of [`balance_residual`](@ref) and the gate that ties
[`newtonian_amf_mode`](@ref), [`newtonian_delta_sigma`](@ref) and
[`newtonian_torque_density_mode`](@ref) together (a wrong δΣ̂, δû_φ or Fourier weight breaks it).
"""
function newtonian_balance_residual(bg::DiscBackground{T}, sol::ModeSolution{T};
                                    flux::AbstractVector = newtonian_amf(bg, sol)) where {T}
    p = sol.params
    r = sol.r
    n = length(r) - 1
    rm = [(r[i] + r[i+1]) / 2 for i in 1:n]
    hm = _newtonian_node_modes(bg, sol, rm)
    lhs = Vector{T}(undef, n); rhs = Vector{T}(undef, n)
    @inbounds for i in 1:n
        lhs[i] = (flux[i+1] - flux[i]) / (r[i+1] - r[i])
        ymid = ((sol.y[1, i] + sol.y[1, i+1]) / 2, (sol.y[2, i] + sol.y[2, i+1]) / 2)
        rhs[i] = newtonian_torque_density_mode(bg, sol.m, rm[i], ymid, hm[i].Φ; landau = p.landau, rc = p.rc)
    end
    scale = maximum(abs, rhs)
    (rm, lhs, rhs, iszero(scale) ? zeros(T, n) : abs.(lhs .- rhs) ./ scale)
end
