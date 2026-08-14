# =============================================================================================================
# API — the user-facing single-mode public surface (legacy MetricReconstruction/api.jl single-mode half) plus the
# unified `reconstruct`/`psi` entry points tying the CCK reconstruction together. Thin wrappers over the green
# primitives (PlusMinus.component_value / Projection.tetrad_matrix / Projection.h_coordinate /
# WeylScalars.weyl_scalar); each defaults `weyl` to `weyl_amplitudes(:unit,…)`. `reconstruct` consolidates the
# single-mode coordinate h into one user entry (coords as a keyword — AHA); `psi` dispatches :psi0/:psi4 over the
# regenerated Weyl scalar. Type-generic {T} (BigFloat-inherited).
# =============================================================================================================

export metric_tetrad, h_tetrad, h_coord, reconstruct, psi

"""
    metric_tetrad(gauge, comp, a, ω, l, m, r, θ; weyl=weyl_amplitudes(:unit,…), t=0, φ=0, bc=:Up) -> Complex

Single tetrad-frame component h_comp(t,r,θ,φ) for one (s=±2) (l,m) mode. `gauge` ∈ :IRG,:ORG; `comp` a nonzero
component symbol (IRG: :nn,:nm,:nmbar,:mm,:mbarmbar; ORG: :ll,:lm,:lmbar,:mm,:mbarmbar).
"""
function metric_tetrad(gauge::Symbol, comp::Symbol, a, ω, l::Integer, m::Integer, r, θ;
                       weyl=nothing, t=0, φ=0, bc::Symbol=:Up)
    w = weyl === nothing ? weyl_amplitudes(:unit, ω, l, m, a) : weyl
    component_value(gauge, comp, a, ω, l, m, r, θ, t, φ, w; bc=bc)
end

"""
    h_tetrad(gauge, a, ω, l, m, r, θ; weyl, t=0, φ=0, bc=:Up) -> 4×4 Complex

Full tetrad-frame matrix h_ab (index order l,n,m,m̄) for one (l,m) mode.
"""
function h_tetrad(gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ; weyl=nothing, t=0, φ=0, bc::Symbol=:Up)
    w = weyl === nothing ? weyl_amplitudes(:unit, ω, l, m, a) : weyl
    tetrad_matrix(gauge, a, ω, l, m, r, θ, t, φ, w; bc=bc)
end

"""
    h_coord(gauge, coords, a, ω, l, m, r, θ; weyl, t=0, φ=0, bc=:Up) -> 4×4 Complex

Coordinate-frame metric perturbation h_μν (coords ∈ :BL,:Ingoing,:Outgoing) for one (l,m) mode.
"""
function h_coord(gauge::Symbol, coords::Symbol, a, ω, l::Integer, m::Integer, r, θ;
                 weyl=nothing, t=0, φ=0, bc::Symbol=:Up)
    w = weyl === nothing ? weyl_amplitudes(:unit, ω, l, m, a) : weyl
    h_coordinate(gauge, coords, a, ω, l, m, r, θ, t, φ, w; bc=bc)
end

"""
    reconstruct(gauge, a, ω, l, m, r, θ; coords=:BL, weyl, t=0, φ=0, bc=:Up) -> 4×4 Complex h_μν

Unified public single-mode reconstruction entry: the coordinate-frame metric perturbation h_μν for one CCK
(gauge ∈ :IRG,:ORG) (l,m) mode, with `coords` (∈ :BL,:Ingoing,:Outgoing) as a keyword. Consolidates `h_coord`
into the user-facing surface; the default `weyl` is `weyl_amplitudes(:unit,…)`.
"""
function reconstruct(gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ;
                     coords::Symbol=:BL, weyl=nothing, t=0, φ=0, bc::Symbol=:Up)
    h_coord(gauge, coords, a, ω, l, m, r, θ; weyl=weyl, t=t, φ=φ, bc=bc)
end

"""
    psi(which, gauge, a, ω, l, m, r, θ; weyl, t=0, φ=0, bc=:Up) -> Complex

Unified public entry for the gauge-invariant regenerated Weyl scalar of the reconstructed metric. `which` ∈
(:psi0, :psi4). Dispatches `WeylScalars.weyl_scalar`; the default `weyl` is `weyl_amplitudes(:unit,…)`.
"""
function psi(which::Symbol, gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ;
             weyl=nothing, t=0, φ=0, bc::Symbol=:Up)
    weyl_scalar(which, gauge, a, ω, l, m, r, θ; weyl=weyl, t=t, φ=φ, bc=bc)
end
