# =============================================================================================================
# Assembly — the CCK mode-sum surface (legacy MetricReconstruction/api.jl mode-sum half): sum coordinate h_μν
# over an (l,m) mode list, and evaluate the linearized-Einstein residual over an (r,θ) grid. Both PARALLEL:
# Threads.@threads over the flattened work (per-mode / per-point output slot; the reduction summation stays in
# fixed order k=1..n so parallel==serial bit-identical). ω may be a scalar or a function ω(l,m). Type-generic {T}.
# =============================================================================================================

import Base.Threads

export h_coord_modesum, einstein_grid

"""
    h_coord_modesum(gauge, coords, a, ω, modes, r, θ; whichweyl=:unit, t=0, φ=0, bc=:Up) -> 4×4 Complex

Sum of the coordinate-frame h_μν over a list of (l,m) modes (PARALLEL over modes with per-mode slots).
`modes` is a vector of (l,m) tuples; `ω` may be a scalar or a function `ω(l,m)`. Each per-mode amplitude is
`weyl_amplitudes(whichweyl, ω_lm, l, m, a)`. The per-mode slots are summed in fixed order (parallel==serial).
"""
function h_coord_modesum(gauge::Symbol, coords::Symbol, a, ω, modes::AbstractVector, r, θ;
                         whichweyl::Symbol=:unit, t=0, φ=0, bc::Symbol=:Up)
    n = length(modes)
    parts = Vector{Any}(undef, n)
    Threads.@threads for k in 1:n
        l, m = modes[k]
        ωlm = ω isa Function ? ω(l, m) : ω
        w = weyl_amplitudes(whichweyl, ωlm, l, m, a)
        parts[k] = h_coordinate(gauge, coords, a, ωlm, l, m, r, θ, t, φ, w; bc=bc)
    end
    acc = zeros(eltype(parts[1]), 4, 4)
    for k in 1:n; acc .+= parts[k]; end
    acc
end

"""
    einstein_grid(gauge, a, ω, l, m, pts; weyl=nothing, bc=:Up) -> Vector{Matrix}

Linearized-Einstein residual G⁽¹⁾_μν[h] at each (r,θ) sample (PARALLEL over the flattened grid; per-point slot).
`pts` is a vector of (r,θ) tuples; t=φ=0. For a vacuum CCK reconstruction each matrix ≈ 0 (G-M4). Since the
linearized-Einstein operator is LINEAR, a vanishing per-mode residual on the grid certifies the assembled mode
sum is vacuum at each grid point.
"""
function einstein_grid(gauge::Symbol, a, ω, l::Integer, m::Integer, pts::AbstractVector; weyl=nothing, bc::Symbol=:Up)
    w = weyl === nothing ? weyl_amplitudes(:unit, ω, l, m, a) : weyl
    n = length(pts)
    out = Vector{Any}(undef, n)
    Threads.@threads for k in 1:n
        r, θ = pts[k]
        out[k] = lin_einstein_residual(gauge, a, ω, l, m, r, θ, zero(a), zero(a), w; bc=bc)
    end
    out
end
