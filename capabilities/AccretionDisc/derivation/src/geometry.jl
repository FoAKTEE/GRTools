# geometry.jl — the ONE generic differential-geometry helper this project needs.
#
# REUSE NOTE (binding, REUSE.md item 1): `Spacetime.christoffel` has exactly two methods, both
# Kerr-specific and both keyed on (a, r, θ) / (r, θ; M) — `christoffel(a,r,θ;coords)` and
# `christoffel(r,θ;M)` (capabilities/Spacetime/src/Geometry/christoffel.jl:36,63).  There is NO
# metric-argument method anywhere in the suite (checked across capabilities/*), so the generic
# connection below is written ONCE, here, and is a PROMOTION CANDIDATE for
# `Spacetime.christoffel(ginv, dg)`; the Kerr-specific ones are never duplicated or forked.
#
# Everything is element-type generic: it is used with `Num` (background) and with the λ-jets of
# jets.jl (perturbed metric), which is why no `AbstractFloat` appears in any signature.

"""
    christoffel_from_metric(ginv, dg) -> Γ[λ,μ,ν]

Levi-Civita connection Γ^λ_{μν} = ½ g^{λσ}(∂_μ g_{σν} + ∂_ν g_{σμ} − ∂_σ g_{μν}) from the INVERSE
metric `ginv[λ,σ] = g^{λσ}` and the metric partials `dg[α][μ,ν] = ∂_α g_{μν}` (a length-n vector of
n×n matrices).  Element-type generic (`Num`, `Jet`, floats).  Symmetric in (μ,ν) by construction.
"""
function christoffel_from_metric(ginv::AbstractMatrix{T}, dg::AbstractVector) where {T}
    n = size(ginv, 1)
    Γ = Array{T}(undef, n, n, n)
    for λ in 1:n, μ in 1:n, ν in μ:n
        acc = zero(T)
        for σ in 1:n
            gi = ginv[λ, σ]
            iszero_entry(gi) && continue
            acc = acc + gi * (dg[μ][σ, ν] + dg[ν][σ, μ] - dg[σ][μ, ν])
        end
        Γ[λ, μ, ν] = acc / 2
        Γ[λ, ν, μ] = Γ[λ, μ, ν]
    end
    Γ
end

"Cheap structural zero test used to skip terms (never a numerical tolerance)."
iszero_entry(x::Num) = isequal(x, Num(0))
iszero_entry(x::Complex{Num}) = isequal(real(x), Num(0)) && isequal(imag(x), Num(0))
iszero_entry(x) = false

"""
    metric_partials(g, ∂) -> dg

`dg[α][μ,ν] = ∂_α g_{μν}` for the partial-derivative operator `∂(x, α)`.
"""
metric_partials(g::AbstractMatrix, ∂) = [[∂(g[μ, ν], α) for μ in axes(g, 1), ν in axes(g, 2)]
                                         for α in axes(g, 1)]
