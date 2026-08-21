# background.jl — the equilibrium disc of arXiv:2608.17003 Sec. II.A, symbolically.
#
# Line element (4)  ds² = −e^{2ν}dt² + e^{2ψ}(dφ − ω̃ dt)² + e^{2μ}dr² + dz²  in (t,r,φ,z), with the
# GENERIC potentials of symbols.jl.  Rotating equilibrium (5) u^α = u^t(t^α + Ω φ^α), the (8)–(10)
# decomposition and the hydrostatic balance (11).  Nothing Kerr-specific appears: (A1a–e) enters only
# when the emitted code is evaluated, through `Spacetime.kerr_equatorial_potentials`.

# ---------------------------------------------------------------------------------------------
# metric (4) and its inverse
# ---------------------------------------------------------------------------------------------

"""
    METRIC

g_{μν} of the line element (4) in (t,r,φ,z).  g_tt = −e^{2ν} + e^{2ψ}ω̃², g_tφ = −e^{2ψ}ω̃,
g_φφ = e^{2ψ}, g_rr = e^{2μ}, g_zz = 1.
"""
const METRIC = Num[-enu^2+epsi^2*omt^2  0        -epsi^2*omt  0
                   0                    emu^2     0           0
                   -epsi^2*omt          0         epsi^2      0
                   0                    0         0           1]

"""
    METRIC_INV

g^{μν}, in closed form (det of the (t,φ) block is −e^{2ν+2ψ}): g^{tt} = −e^{−2ν},
g^{tφ} = −ω̃e^{−2ν}, g^{φφ} = e^{−2ψ} − ω̃²e^{−2ν}, g^{rr} = e^{−2μ}, g^{zz} = 1.
Certified against `METRIC` by `certify_background` (g·g⁻¹ = 1).
"""
const METRIC_INV = Num[-1/enu^2        0         -omt/enu^2                0
                       0               1/emu^2    0                        0
                       -omt/enu^2      0          1/epsi^2 - omt^2/enu^2   0
                       0               0          0                        1]

"√(−g) = e^{ν+ψ+μ} (the volume element of (4))."
const SQRT_MINUS_G = enu * epsi * emu

# ---------------------------------------------------------------------------------------------
# equilibrium four-velocity, (5), (8)–(10)
# ---------------------------------------------------------------------------------------------

"u^t = e^{−ν}/√(1−v²) = W/e^ν  (8)."
const UT = W / enu

"u^α = u^t(1, 0, Ω, 0)  (5), with Ω = ω̃ + v e^{ν−ψ} (10)."
const UUP = Num[UT, 0, UT * Ω_expr, 0]

"u_α = g_{αμ}u^μ."
const ULOW = lower(METRIC, UUP)

"q^{μν} = u^μu^ν + g^{μν}  (2)."
const QUP = Num[UUP[μ] * UUP[ν] + METRIC_INV[μ, ν] for μ in 1:4, ν in 1:4]

"∂_α g_{μν} of the background metric (only α = r survives)."
const DMETRIC = metric_partials(METRIC, ∂bg)

"Γ^λ_{μν} of the background metric (4) — from the generic `christoffel_from_metric`."
const GAMMA = christoffel_from_metric(METRIC_INV, DMETRIC)

# ---------------------------------------------------------------------------------------------
# hydrostatic balance (11)
# ---------------------------------------------------------------------------------------------

"""
    HYDRO_RHS

Right-hand side of (11): p′/(ε+p) = −(ν′ − v²ψ′ + e^{ψ−ν}vω̃′)/(1−v²), written with W² = 1/(1−v²)
so that it stays a rational function of the alphabet.
"""
const HYDRO_RHS = -W^2 * (nu_r - v^2 * psi_r + epsi * v * omt_r / enu)

"""
    HYDRO_SUB

The ON-SHELL substitution p′ → (ε+p)·`HYDRO_RHS`.  Applying it puts an expression on the background
Euler equation (11); route R2 used the same device (its "raw" vs "on-shell" forms).
"""
const HYDRO_SUB = Dict(p_r => (eps + p) * HYDRO_RHS)

"Put `x` on-shell w.r.t. (11)."
onshell(x::Num) = Symbolics.substitute(x, HYDRO_SUB)
onshell(x::Complex{Num}) = onshell(real(x)) + im * onshell(imag(x))
onshell(x::AbstractArray) = map(onshell, x)

# ---------------------------------------------------------------------------------------------
# background equations of motion (3b), (3c) — the SAME functions the perturbation uses
# ---------------------------------------------------------------------------------------------

"a^μ = u^ν∇_νu^μ of the equilibrium."
const ACC_UP = four_acceleration_up(GAMMA, UUP, ∂bg)
"a_α = g_{αμ}a^μ."
const ACC_LOW = lower(METRIC, ACC_UP)
"∇_αu^α of the equilibrium — must vanish identically, Eq. (6)."
const DIV_U = divergence_u(GAMMA, UUP, ∂bg)
"ℰ^μ (3c) of the equilibrium."
const EULER_UP_BG = euler_up(METRIC_INV, GAMMA, UUP, p, eps, ∂bg)
"Energy equation (3b) of the equilibrium — vanishes because ∇_αu^α = 0 and u^r = 0."
const ENERGY_BG = energy_eq(GAMMA, UUP, p, eps, ∂bg)
"g_{μν}u^μu^ν + 1 — vanishes by (8)."
const NORM_BG = norm_constraint(METRIC, UUP)
