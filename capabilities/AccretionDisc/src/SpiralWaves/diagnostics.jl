# SpiralWaves/diagnostics.jl — the angular-momentum-flux / torque diagnostics of arXiv:2608.17003
# Sec. II.E (main.tex:784-851): for one solved azimuthal mode,
#
#   F_J(r)  = ⟨√−g j^r_φ⟩ = √−g [u_φ⟨δp δu^r⟩ + (ε+p)⟨δu^r δu_φ⟩]                            (72a)
#   𝔗(r)    = ∫ dr' √−g ⟨𝒯⟩ = ½∫ dr' √−g ⟨δ_ξT^{αβ} ∂_φ h_αβ⟩                                (72b)
#   T       = 𝔗(r_out) − 𝔗(r_in)                                                              (74)
#   T_L     = F_J(r_out) − F_J(r_in)                                                          (75)
#
# with ⟨·⟩ = ∫₀^{2π}dφ evaluated by the Fourier identity (73) ⟨XY⟩ = 2πX̂₀Ŷ₀ + 4πΣ_{m≥1}Re[X̂_mŶ_m*] —
# an INTEGRAL over φ, not an average, so ONE mode m ≥ 1 contributes 4πRe[X̂_mŶ_m*] (and m = 0 would
# contribute 2πX̂₀Ŷ₀; `solve_mode` requires m ≥ 1, the forced problem has no m = 0 response).
#
# (72a) is evaluated through the identity (69) j^r_J = δpδu^r u_φ + (ε+p)δu^rδu_φ, which the CAS node
# EGD::deriv-diagnostics certified ON-SHELL and φ-integrated (certificates [K69], [K70]; it is NOT a
# pointwise off-shell identity — ERRATUM E-M3, so it may only be applied to solutions of the master
# system, which is exactly what a `ModeSolution` is). The per-mode integrands
# `_r2_amf` / `_r2_torque_density` are GENERATED from the route-R2 derivation (r2_generated.jl, see
# master.jl); √−g = e^{ν+ψ+μ} for the line element (4).
#
# The AMF–torque balance (71) dF_J/dr = d𝔗/dr holds on-shell wherever the Landau term is negligible
# ([K71]); it is the paper's Fig. 3 and the strongest end-to-end gate of the whole pipeline
# ([`balance_residual`](@ref)). Near r_c the +iε shift makes the system non-Hermitian and deposits
# angular momentum, which is the corotation torque of the (75) split.
# Type-generic `{T}`; BigFloat inherits. Additive (EGD::sw-diagnostics).

# `lindblad_torque` is the EXISTING Modes generic (the Newtonian Lin-Shu/Goldreich-Tremaine resonance
# torque of arXiv:2601.19123); the spiral-wave torque below is the same CONCEPT in a different model,
# so it is added as a METHOD on that generic (one name = one concept), selected by its ModeSolution
# argument — exactly as background.jl does for disc_OmegaK / kappar2 / csound / eden0.
import ..Modes: lindblad_torque

"""
    sqrt_minus_g(bg::DiscBackground, r) -> √−g

Metric volume factor of the equatorial line element (4): √−g = e^{ν+ψ+μ} (g_zz = 1). Equals `r`
identically for the Kerr (A1) potentials — asserted in the tests, not assumed here.
"""
function sqrt_minus_g(bg::DiscBackground{T}, r) where {T}
    pot = kerr_equatorial_potentials(convert(promote_type(T, typeof(r)), bg.χ), r)
    exp(pot.ν + pot.ψ + pot.μ)
end

# Fourier weight of (73) for a single mode: 4π for m ≥ 1, 2π for m = 0.
_fourier_weight(::Type{T}, m::Integer) where {T} = m == 0 ? 2 * T(π) : 4 * T(π)

"""
    amf_mode(bg, m, r, y; hmodes, landau = true, rc = nothing) -> F_J,m(r)

Contribution of the single azimuthal mode `m` to the radial angular-momentum flux (72a) at radius `r`,
for the master state `y = (δÛ, −iδû^r)`:

    F_J,m(r) = √−g(r) · w_m · Re[ u_φ δp̂ δû^{r*} + (ε+p) δû^r δû_φ* ],   w_m = 4π (m ≥ 1), 2π (m = 0),

with the bracket generated from route R2 (`_r2_amf`) and the covariant δû_φ of [`delta_u_phi`](@ref)
(ERRATUM E-M1: the COVARIANT component, i.e. r×(50a) in the Newtonian limit, with u_φ = r²Ω). Real by
construction.
"""
function amf_mode(bg::DiscBackground, m::Integer, r, y; hmodes, landau::Bool = true, rc = nothing)
    inp, ω̃, vp = r2_context(bg, m, r; hmodes = hmodes, landau = landau, rc = rc)
    T = typeof(inp.v)
    sqrt_minus_g(bg, r) * _fourier_weight(T, m) * real(_r2_amf(inp, ω̃, vp, y[1], y[2]))
end

"""
    torque_density_mode(bg, m, r, y; hmodes, landau = true, rc = nothing) -> √−g⟨𝒯⟩_m(r)

Contribution of the single azimuthal mode `m` to the torque density of (72b) at radius `r`:

    √−g⟨𝒯⟩_m = √−g · w_m · Re[ ½ Σ_{αβ} δ_ξT̂^{αβ} (i m ĥ_αβ)^* ],

the integrand of 𝔗(r). Generated from route R2 (`_r2_torque_density`, built on δ_ξT^{αβ} = δT^{αβ} −
δ_{ξ=0}T^{αβ} of (40)-(42) in the master variables). Real by construction. The paper's sign convention
is the torque ON THE DISC due to the EMRI (M:898).
"""
function torque_density_mode(bg::DiscBackground, m::Integer, r, y; hmodes, landau::Bool = true, rc = nothing)
    inp, ω̃, vp = r2_context(bg, m, r; hmodes = hmodes, landau = landau, rc = rc)
    T = typeof(inp.v)
    sqrt_minus_g(bg, r) * _fourier_weight(T, m) * real(_r2_torque_density(inp, ω̃, vp, y[1], y[2]))
end

# per-node source modes with the settings the solve used
_node_modes(bg::DiscBackground, sol::ModeSolution, rs) =
    source_modes(bg, sol.m, rs; threads = sol.params.threads, N = sol.params.hmodes_N,
                 mEMRI = sol.params.mEMRI)

"""
    amf(bg::DiscBackground, sol::ModeSolution; hmodes = <source modes at the nodes>) -> Vector{T}
    amf(bg::DiscBackground, sol::ModeSolution, idx; hmodes = <source modes at `sol.r[idx]`>) -> Vector{T}

The angular-momentum flux F_J(r) (72a) of a solved mode at every grid node ([`amf_mode`](@ref)), or
only at the node indices `idx`. Pass `hmodes` (the singular-field Fourier modes at the corresponding
radii, [`source_modes`](@ref)) to share ONE projection with [`torque_density`](@ref) — the projection
(60) is the dominant cost of the diagnostics and is otherwise repeated per call. The indexed method
exists because the (75) reduction needs F_J at FOUR nodes only (the two domain edges and the two edges
of the corotation window), never on the whole grid; both methods evaluate exactly the same expression,
so they agree bitwise.
"""
function amf(bg::DiscBackground{T}, sol::ModeSolution{T};
             hmodes::AbstractVector = _node_modes(bg, sol, sol.r)) where {T}
    p = sol.params
    [amf_mode(bg, sol.m, sol.r[i], (sol.y[1, i], sol.y[2, i]); hmodes = hmodes[i], landau = p.landau, rc = p.rc)
     for i in eachindex(sol.r)]
end

function amf(bg::DiscBackground{T}, sol::ModeSolution{T}, idx::AbstractVector{<:Integer};
             hmodes::AbstractVector = _node_modes(bg, sol, sol.r[idx])) where {T}
    p = sol.params
    [amf_mode(bg, sol.m, sol.r[i], (sol.y[1, i], sol.y[2, i]); hmodes = hmodes[j], landau = p.landau, rc = p.rc)
     for (j, i) in enumerate(idx)]
end

"""
    torque_density(bg::DiscBackground, sol::ModeSolution; hmodes = <source modes at the nodes>) -> Vector{T}

The torque density √−g⟨𝒯⟩(r) (the integrand of 72b) of a solved mode at every grid node
([`torque_density_mode`](@ref)); `hmodes` as in [`amf`](@ref).
"""
function torque_density(bg::DiscBackground{T}, sol::ModeSolution{T};
                        hmodes::AbstractVector = _node_modes(bg, sol, sol.r)) where {T}
    p = sol.params
    [torque_density_mode(bg, sol.m, sol.r[i], (sol.y[1, i], sol.y[2, i]); hmodes = hmodes[i], landau = p.landau, rc = p.rc)
     for i in eachindex(sol.r)]
end

"""
    cumulative_torque(bg, sol; density = torque_density(bg, sol)) -> Vector{T}

𝔗(r) = ∫ dr' √−g⟨𝒯⟩ (72b) accumulated along the grid by the trapezoid rule (second order, the order
of the solver itself) with 𝔗(r_in) = 0. The paper writes the lower limit as r = a•; only DIFFERENCES
of 𝔗 enter (74)-(75), so the constant is immaterial — pass `density` to reuse an already computed
integrand.
"""
function cumulative_torque(bg::DiscBackground{T}, sol::ModeSolution{T};
                           density::AbstractVector = torque_density(bg, sol)) where {T}
    r = sol.r
    out = Vector{T}(undef, length(r))
    out[1] = zero(T)
    @inbounds for i in 2:length(r)
        out[i] = out[i-1] + (r[i] - r[i-1]) * (density[i] + density[i-1]) / 2
    end
    out
end

"""
    total_torque(bg, sol; density = ...) -> T

Total torque on the grid (74): T = 𝔗(r_out) − 𝔗(r_in) — the full torque the EMRI exerts on the disc
in this mode, corotation deposition included.
"""
total_torque(bg::DiscBackground{T}, sol::ModeSolution{T}; kwargs...) where {T} =
    last(cumulative_torque(bg, sol; kwargs...))

"""
    lindblad_torque(bg, sol; flux = amf(bg, sol)) -> T

Lindblad (wave) torque (75): T_L = F_J(r_out) − F_J(r_in) — the total torque minus the angular
momentum deposited at corotation, T_L = T − (Δ_cT − Δ_cF_J), by the balance law (71) integrated on
either side of r_c. Equal to the angular momentum carried away by the spiral waves through the two
domain edges.

NOTE the name clash with the Newtonian `Modes.lindblad_torque` of arXiv:2601.19123 (a different
model): inside `SpiralWaves` this method is selected by its `ModeSolution` argument.
"""
lindblad_torque(bg::DiscBackground{T}, sol::ModeSolution{T};
                flux::AbstractVector = amf(bg, sol)) where {T} = flux[end] - flux[1]

"""
    torque_units(bg::DiscBackground; q = 1) -> (F_J0, F_J1)

The Goldreich–Tremaine normalisations (76): F_{J,0} = ε₀a•⁴Ω•²q²/h₀³ and F_{J,1} = F_{J,0}/a• — named
accessor over the [`F_J0`](@ref) / [`F_J1`](@ref) of background.jl (no second definition). Torques are
quoted in units of F_{J,0}h₀ in the paper's figures.
"""
torque_units(bg::DiscBackground; q = 1) = (F_J0(bg; q = q), F_J1(bg; q = q))

"""
    m_cut(Tm::AbstractVector, frac) -> m_K

The cutoff mode number m_K of M:988-1000: the SMALLEST m such that the cumulative sum of the per-mode
torques `Tm` (`Tm[i]` = the torque of mode m = i) reaches `frac` of the total, e.g.
`m_cut(Tm, 0.95)` = m₉₅ and `m_cut(Tm, 0.45)` = m₄₅. Returns `length(Tm)` if the fraction is never
reached (an under-resolved mode sum — the caller should extend `Tm`).
"""
function m_cut(Tm::AbstractVector, frac::Real)
    0 < frac ≤ 1 || throw(ArgumentError("m_cut: need 0 < frac ≤ 1, got $frac"))
    tot = sum(Tm)
    iszero(tot) && throw(ArgumentError("m_cut: the per-mode torques sum to zero"))
    target = frac * tot
    acc = zero(tot)
    for (i, t) in enumerate(Tm)
        acc += t
        (tot > 0 ? acc ≥ target : acc ≤ target) && return i
    end
    length(Tm)
end

"""
    corotation_split(bg, sol; halfwidth = 3bg.σε, flux = amf(bg, sol), density = torque_density(bg, sol))
        -> (T, T_L, T_c, Δ𝔗, ΔF, window)

The corotation/Lindblad decomposition (75) of a solved mode: the total torque T (74), the Lindblad
(wave) torque T_L = F_J(r_out) − F_J(r_in), the corotation torque T_c = Δ_c𝔗 − Δ_cF_J (the angular
momentum DEPOSITED near corotation), and the two jumps across the corotation neighbourhood
`window` = the index range of |r − r_c| ≤ `halfwidth`.

The paper does not fix the width of that neighbourhood (only differences matter for the second
equality of (75)); the default 3σ_ε is the Landau profile's own scale (63). FINDING (this node,
measured at χ = 0, a• = 15, h0 = 0.05, 𝔭 = 𝔮 = 0, m = 3): the deposition is NOT confined to the sinh
corotation BLOCK of (C4)-(C5) — that block spans only r_c ± 10⁻³r_c, and ≈ 0.7 % of the jump happens
just outside it, which is the accuracy floor of the split when the block edges are used as the window.
With `halfwidth = 3σ_ε` the identity T − T_c = T_L closes at the truncation error of the scheme and
converges with resolution.
"""
function corotation_split(bg::DiscBackground{T}, sol::ModeSolution{T};
                          halfwidth = 3 * bg.σε,
                          flux::AbstractVector = amf(bg, sol),
                          density::AbstractVector = torque_density(bg, sol)) where {T}
    𝔗 = cumulative_torque(bg, sol; density = density)
    win = corotation_window(sol, halfwidth)
    i1, i2 = first(win), last(win)
    Δ𝔗 = 𝔗[i2] - 𝔗[i1]
    ΔF = flux[i2] - flux[i1]
    (T = 𝔗[end] - 𝔗[1], T_L = flux[end] - flux[1], T_c = Δ𝔗 - ΔF, Δ𝔗 = Δ𝔗, ΔF = ΔF, window = i1:i2)
end

"""
    corotation_window(sol::ModeSolution, halfwidth) -> UnitRange

Index range of the corotation neighbourhood |r − r_c| ≤ `halfwidth` on the solution's own grid
(`r_c` = `sol.params.rc`, the radius the Landau shift was centred on). The ONE definition of that
window, shared by [`corotation_split`](@ref) and by the survey reduction (which evaluates the AMF at
its two edges only). Throws when no node lies inside.
"""
function corotation_window(sol::ModeSolution, halfwidth)
    rc = sol.params.rc
    idx = findall(r -> abs(r - rc) ≤ halfwidth, sol.r)
    isempty(idx) && throw(ArgumentError("corotation_window: no grid node within $halfwidth of r_c = $rc"))
    first(idx):last(idx)
end

"""
    balance_residual(bg, sol; flux = amf(bg, sol)) -> (r, lhs, rhs, rel)

The AMF–torque balance law (71) dF_J/dr = d𝔗/dr evaluated CELL BY CELL in the discretisation's own
(second-order, midpoint) form — the numerical content of the paper's Fig. 3:

  * `lhs[i] = (F_J(r_{i+1}) − F_J(r_i))/Δr_i`,
  * `rhs[i] = √−g⟨𝒯⟩` at the cell midpoint r_{i+1/2}, evaluated on the midpoint state
    (y_i + y_{i+1})/2 — the very state the midpoint scheme (64) advances,
  * `rel[i] = |lhs − rhs| / max|rhs|` (normalised by the global scale of the torque density, so that
    the zero crossings of 𝒯 do not produce spurious "relative" errors),
  * `r` = the cell midpoints.

Both sides are exact functions of the SAME discrete solution, so `rel` measures the truncation error
of the scheme (it converges at second order), never a modelling assumption: the identity itself is an
exact on-shell symbolic result ([K71] of EGD::deriv-diagnostics).
"""
function balance_residual(bg::DiscBackground{T}, sol::ModeSolution{T};
                          flux::AbstractVector = amf(bg, sol)) where {T}
    p = sol.params
    r = sol.r
    n = length(r) - 1
    rm = [(r[i] + r[i+1]) / 2 for i in 1:n]
    hm = _node_modes(bg, sol, rm)
    lhs = Vector{T}(undef, n); rhs = Vector{T}(undef, n)
    @inbounds for i in 1:n
        lhs[i] = (flux[i+1] - flux[i]) / (r[i+1] - r[i])
        ymid = ((sol.y[1, i] + sol.y[1, i+1]) / 2, (sol.y[2, i] + sol.y[2, i+1]) / 2)
        rhs[i] = torque_density_mode(bg, sol.m, rm[i], ymid; hmodes = hm[i], landau = p.landau, rc = p.rc)
    end
    scale = maximum(abs, rhs)
    (rm, lhs, rhs, iszero(scale) ? zeros(T, n) : abs.(lhs .- rhs) ./ scale)
end
