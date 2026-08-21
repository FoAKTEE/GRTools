# stage1.jl — Julia re-certification of the background relations of arXiv:2608.17003 Sec. II.A.
#
# Route R3 stage 1.  Every identity is checked in BOTH modalities where that is affordable, and the
# table records which one is the admission evidence for each.  Nothing here consumes the Mathematica
# certificates of EGD::deriv-background — this is an independent Julia execution of the same algebra.

"""
    metric_compatibility() -> ∇_λ g_{μν}

∇_λ g_{μν} = ∂_λ g_{μν} − Γ^σ_{λμ}g_{σν} − Γ^σ_{λν}g_{μσ}; must vanish for the Levi-Civita connection.
"""
function metric_compatibility()
    [∂bg(METRIC[μ, ν], λ) - sum(GAMMA[σ, λ, μ] * METRIC[σ, ν] + GAMMA[σ, λ, ν] * METRIC[μ, σ]
                                for σ in 1:4)
     for λ in 1:4, μ in 1:4, ν in 1:4]
end

"""
    euler_kinematic_identity() -> a_α − [u^t u_φ ∂_αΩ − ∂_α log u^t]

Eq. (7), first equality — a pure kinematic identity of the rotating equilibrium (5),(8).
∂_α log u^t is written as (∂_α u^t)/u^t (rational, no logarithm).
"""
euler_kinematic_identity() =
    [ACC_LOW[α] - (UT * ULOW[3] * ∂bg(Ω_expr, α) - ∂bg(UT, α) / UT) for α in 1:4]

"""
    hydrostatic_from_euler() -> p′ solved from ℰ^r = 0, minus (ε+p)·(11)

Eq. (7) second equality ⇒ Eq. (11): the r-component of the background Euler equation (3c) is linear
in p′; solving it must return exactly the printed hydrostatic balance.
"""
function hydrostatic_from_euler()
    E_r = EULER_UP_BG[IR]
    c1 = Symbolics.derivative(E_r, p_r)              # ℰ^r = c1 p′ + c0 (linear in p′)
    c0 = Symbolics.substitute(E_r, Dict(p_r => Num(0)))
    (-c0 / c1) - (eps + p) * HYDRO_RHS
end

"""
    stage1_certificates(; npts=50, seed=2608, verbose=true) -> Vector{Certificate}

The stage-1 certificate suite: metric (4) and its inverse, √−g, the (8)–(10) decomposition, the
derivation table, ∇g = 0, and Eqs. (6), (7), (3b), (11).
"""
function stage1_certificates(; npts::Int = 50, seed::Int = 2608, verbose::Bool = true)
    C = Certificate[]
    push!(C, certify_symbolic("metric-inverse  g·g⁻¹ − 1", "(4)",
                              vec(METRIC * METRIC_INV - Num[i == j for i in 1:4, j in 1:4]);
                              verbose))
    push!(C, certify_symbolic("volume element  det g + (e^{ν+ψ+μ})²", "(4)",
                              Symbolics.det(METRIC) + SQRT_MINUS_G^2; verbose))
    push!(C, certify_symbolic("metric compatibility  ∇_λ g_{μν}", "(4)",
                              vec(metric_compatibility()); verbose))
    push!(C, certify_symbolic("normalisation  g_{μν}u^μu^ν + 1", "(8)", NORM_BG; verbose))
    push!(C, certify_symbolic("u_t, u_φ", "(9)",
                              [ULOW[IT] + W * (enu + epsi * omt * v), ULOW[IPHI] - epsi * v * W];
                              verbose))
    push!(C, certify_symbolic("ZAMO speed  v = (Ω−ω̃)e^{ψ−ν}", "(10)",
                              (Ω_expr - omt) * epsi / enu - v; verbose))
    push!(C, certify_symbolic("derivation table  ddr(Ω) = Ω′", "(10)",
                              ddr(Ω_expr) - Ω_r_expr; verbose))
    push!(C, certify_symbolic("expansion-free flow  ∇_αu^α", "(6)", DIV_U; verbose))
    push!(C, certify_exact("expansion-free flow  ∇_αu^α", "(6)", DIV_U; npts, seed, verbose))
    push!(C, certify_symbolic("acceleration identity  a_α − u^tu_φ∂_αΩ + ∂_α log u^t", "(7)",
                              euler_kinematic_identity(); verbose))
    push!(C, certify_exact("acceleration identity  a_α − u^tu_φ∂_αΩ + ∂_α log u^t", "(7)",
                           euler_kinematic_identity(); npts, seed, verbose))
    push!(C, certify_symbolic("background energy eq.", "(3b)", ENERGY_BG; verbose))
    push!(C, certify_symbolic("hydrostatic balance  p′ from ℰ^r", "(11)",
                              hydrostatic_from_euler(); verbose))
    push!(C, certify_exact("hydrostatic balance  p′ from ℰ^r", "(11)",
                           hydrostatic_from_euler(); npts, seed, verbose))
    push!(C, certify_symbolic("background Euler ℰ^μ on-shell", "(3c),(11)",
                              onshell.(EULER_UP_BG); verbose))
    push!(C, certify_exact("background Euler ℰ^μ on-shell", "(3c),(11)",
                           onshell.(EULER_UP_BG); npts, seed, onshell = true, verbose))
    C
end
