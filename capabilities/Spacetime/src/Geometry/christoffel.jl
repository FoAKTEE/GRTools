# Geometry — Kerr Christoffel symbols Γ^λ_{μν} (BL, analytic ∂g) and the a=0/M-only
# Schwarzschild connection + its r/θ derivatives, rehomed from Punctures/worldline_self_field
# (§4.5 Rev-2: geometry is Spacetime physics, not NumericsCore). The full-Kerr and Schwarzschild
# methods coexist by ARITY: `christoffel(a, r, θ; coords=:BL)` (3 positional) vs the M-only
# `christoffel(r, θ; M=1)` (2 positional + keyword mass). Type-generic {T}.

# ---- analytic ∂_r g, ∂_θ g for the BL Kerr metric (M=1) --------------------------------------
function _dg_BL(a::T, r, θ) where {T}
    sθ = sin(θ); cθ = cos(θ); s2 = sθ^2
    Σ = r^2 + a^2*cθ^2; Δ = r^2 - 2r + a^2
    Σr = 2r;            Σθ = -2*a^2*sθ*cθ
    Δr = 2r - 2
    gr = zeros(T, 4, 4); gθ = zeros(T, 4, 4)
    gr[1,1] = 2*(Σ - r*Σr)/Σ^2
    gθ[1,1] = 2*(-r*Σθ)/Σ^2
    gr[1,4] = gr[4,1] = -2*a*s2*(Σ - r*Σr)/Σ^2
    gθ[1,4] = gθ[4,1] = -2*a*r*(2*sθ*cθ*Σ - s2*Σθ)/Σ^2
    gr[2,2] = (Σr*Δ - Σ*Δr)/Δ^2
    gθ[2,2] = Σθ/Δ
    gr[3,3] = Σr
    gθ[3,3] = Σθ
    A  = r^2 + a^2 + 2*a^2*r*s2/Σ
    Ar = 2r + 2*a^2*(s2*Σ - r*s2*Σr)/Σ^2
    Aθ = 2*a^2*r*(2*sθ*cθ*Σ - s2*Σθ)/Σ^2
    gr[4,4] = Ar*s2
    gθ[4,4] = Aθ*s2 + A*2*sθ*cθ
    (gr, gθ)
end

"""
    christoffel(a, r, θ; coords=:BL) -> 4×4×4 Γ[λ,μ,ν]

Kerr Christoffel symbols Γ^λ_{μν} (M = 1), λ the upper index, order (t,r,θ,φ). Computed from
the analytic ∂g (∂_t = ∂_φ = 0). Only the Boyer–Lindquist chart `coords=:BL` is provided.
"""
function christoffel(a::Real, r::Real, θ::Real; coords::Symbol = :BL)
    coords === :BL || throw(ArgumentError("christoffel: only coords=:BL implemented (got $(repr(coords)))"))
    aa, rr, tt = _geometry_promote(a, r, θ)
    T = typeof(aa)
    gi = kerr_metric_inv(aa, rr, tt; coords = :BL)
    gr, gθ = _dg_BL(aa, rr, tt)
    dg = zeros(T, 4, 4, 4)   # dg[α,μ,ν] = ∂_α g_{μν}
    dg[2,:,:] = gr; dg[3,:,:] = gθ
    Γ = zeros(T, 4, 4, 4)
    @inbounds for λ in 1:4, μ in 1:4, ν in 1:4
        acc = zero(T)
        for σ in 1:4
            gi[λ,σ] == 0 && continue
            acc += gi[λ,σ]*(dg[μ,σ,ν] + dg[ν,σ,μ] - dg[σ,μ,ν])
        end
        Γ[λ,μ,ν] = acc/2
    end
    Γ
end

# ---- Schwarzschild (a=0) M-only connection + r/θ derivatives (from Punctures) -----------------
"""
    christoffel(r, θ; M=1) -> 4×4×4 Γ[λ,μ,ν]

Schwarzschild Christoffel symbols Γ^λ_{μν} for a hole of mass `M` (closed form), λ the upper
index, order (t,r,θ,φ). The a→0 limit of the Kerr method (with explicit mass M).
"""
function christoffel(r, θ; M = one(float(typeof(r))))
    T = float(promote_type(typeof(r), typeof(θ), typeof(M)))
    Γ = zeros(T, 4, 4, 4)
    f = 1 - 2M/r; s = sin(θ); c = cos(θ)
    set(d, a, b, v) = (Γ[d,a,b] = v; Γ[d,b,a] = v)
    set(1, 1, 2,  M/(r^2*f))          # Γ^t_{tr}
    set(2, 1, 1,  M*f/r^2)            # Γ^r_{tt}
    set(2, 2, 2, -M/(r^2*f))          # Γ^r_{rr}
    set(2, 3, 3, -(r - 2M))           # Γ^r_{θθ}
    set(2, 4, 4, -(r - 2M)*s^2)       # Γ^r_{φφ}
    set(3, 2, 3,  1/r)                # Γ^θ_{rθ}
    set(3, 4, 4, -s*c)                # Γ^θ_{φφ}
    set(4, 2, 4,  1/r)                # Γ^φ_{rφ}
    set(4, 3, 4,  c/s)                # Γ^φ_{θφ}
    Γ
end

"""
    christoffel_dr(r, θ; M=1) -> 4×4×4 ∂_r Γ[λ,μ,ν]

Radial derivative of the Schwarzschild Christoffel symbols (closed form).
"""
function christoffel_dr(r, θ; M = one(float(typeof(r))))
    T = float(promote_type(typeof(r), typeof(θ), typeof(M)))
    G = zeros(T, 4, 4, 4)
    s = sin(θ); den = r^2 - 2M*r
    set(d, a, b, v) = (G[d,a,b] = v; G[d,b,a] = v)
    set(1, 1, 2, -M*2*(r - M)/den^2)      # ∂_r Γ^t_{tr}
    set(2, 1, 1,  M*(6M - 2r)/r^4)        # ∂_r Γ^r_{tt}
    set(2, 2, 2,  M*2*(r - M)/den^2)      # ∂_r Γ^r_{rr}
    set(2, 3, 3, -one(T))                 # ∂_r Γ^r_{θθ}
    set(2, 4, 4, -s^2)                    # ∂_r Γ^r_{φφ}
    set(3, 2, 3, -1/r^2)                  # ∂_r Γ^θ_{rθ}
    set(4, 2, 4, -1/r^2)                  # ∂_r Γ^φ_{rφ}
    G
end

"""
    christoffel_dtheta(r, θ; M=1) -> 4×4×4 ∂_θ Γ[λ,μ,ν]

Polar derivative of the Schwarzschild Christoffel symbols (closed form).
"""
function christoffel_dtheta(r, θ; M = one(float(typeof(r))))
    T = float(promote_type(typeof(r), typeof(θ), typeof(M)))
    G = zeros(T, 4, 4, 4)
    s = sin(θ); c = cos(θ)
    set(d, a, b, v) = (G[d,a,b] = v; G[d,b,a] = v)
    set(2, 4, 4, -(r - 2M)*2*s*c)         # ∂_θ Γ^r_{φφ}
    set(3, 4, 4, -(c^2 - s^2))            # ∂_θ Γ^θ_{φφ}
    set(4, 3, 4, -1/s^2)                  # ∂_θ Γ^φ_{θφ}
    G
end
