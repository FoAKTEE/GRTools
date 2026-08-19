# reproduce/arXiv-2608.17003/_spiral_panels.jl — PRESENTATION ONLY for the two spiral-wake figures
# of arXiv:2608.17003 (Fig. 1 `spiral_grid.png`, Fig. 2 `spirals_two_rows.png`).
#
# Zero physics: the δU maps come from `AccretionDisc.SpiralWaves.synthesize_cartesian` in the drivers
# (fig01.jl, fig02.jl); this file only knows how the paper DRAWS them, following sections 0–2 of the
# VLM style manual `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md` (emri-disk repo):
#
#   * diverging `RdBu` (NOT reversed: δU < 0 red, δU > 0 blue), symmetric limits ±0.77, colorbar
#     `extend='both'` with 9 ticks at multiples of 0.154 (−0.616 … 0.616);
#   * square panels, black box frame, ticks INWARD on all four sides, minor ticks on, no grid;
#   * black filled horizon disc r₊ and a black dash-dot circle at r_ISCO, with an unframed
#     `r_ISCO` legend entry in the top-left panel only;
#   * parameter stamp inside the axes at the top-left, CM serif with inline math.
#
# It is shared by fig01.jl and fig02.jl so the two figures cannot drift apart (the manual's
# cross-figure consistency rule 2: the δU colour encoding is semantic and fixed).
module SpiralPanels

using Plots
using Plots: cgrad, Shape, RGB, mm

export RDBU, CBLIM, CBTICKS, wake_panel, colorbar_panel, panel_stamp!

"Diverging colormap of the δU maps (manual section 0, palette 6): matplotlib/ColorBrewer `RdBu`."
const RDBU = cgrad(:RdBu)

"Symmetric colour limits of Figs. 1–2 (manual section 1: ticks are ±0.154·k, k ≤ 4 ⇒ vmax = 5×0.154)."
const CBLIM = 0.77

"The 9 labelled colorbar ticks of Figs. 1–2: −0.616, −0.462, … , +0.616."
const CBTICKS = [-0.616, -0.462, -0.308, -0.154, 0.0, 0.154, 0.308, 0.462, 0.616]

const CIRCLE_θ = range(0, 2π; length = 721)

"""
    wake_panel(xs, ys, Z; lim, r_isco, r_plus, xticks, yticks, xlabel, ylabel, title,
               stamp, stamp_box, isco_legend) -> Plots.Plot

One δU panel of Fig. 1 / Fig. 2: `Z[j, i]` is the field at `(xs[i], ys[j])` (already transposed for
`heatmap`), `NaN` outside the disc (drawn white). Draws the RdBu map at fixed ±[`CBLIM`](@ref), the
filled black horizon disc of radius `r_plus`, the black dash-dot ISCO circle, the parameter stamp and
(when `isco_legend`) the unframed `r_ISCO` legend entry.
"""
function wake_panel(xs, ys, Z; lim = 40.0, r_isco = nothing, r_plus = nothing,
                    xticks = nothing, yticks = nothing, xlabel = "", ylabel = "", title = "",
                    stamp = "", stamp_box::Bool = false, isco_legend::Bool = false,
                    labelsize = 10, ticksize = 8, stampsize = 8)
    p = heatmap(xs, ys, Z; c = RDBU, clims = (-CBLIM, CBLIM), colorbar = false,
                aspect_ratio = 1, xlims = (-lim, lim), ylims = (-lim, lim),
                framestyle = :box, grid = false, tick_direction = :in, minorticks = 4,
                widen = false, background_color_inside = :white,
                xticks = xticks === nothing ? :auto : xticks,
                yticks = yticks === nothing ? :auto : yticks,
                xlabel = xlabel, ylabel = ylabel, title = title,
                guidefontsize = labelsize, tickfontsize = ticksize, titlefontsize = labelsize,
                legend = isco_legend ? :bottomleft : false, legendfontsize = ticksize,
                foreground_color_legend = nothing, background_color_legend = nothing)
    if r_plus !== nothing && r_plus > 0
        plot!(p, Shape(r_plus .* cos.(CIRCLE_θ), r_plus .* sin.(CIRCLE_θ));
              color = :black, linecolor = :black, linewidth = 0, label = "")
    end
    if r_isco !== nothing
        plot!(p, r_isco .* cos.(CIRCLE_θ), r_isco .* sin.(CIRCLE_θ);
              color = :black, linestyle = :dashdot, linewidth = 1.1,
              label = isco_legend ? "\$r_\\mathrm{ISCO}\$" : "")
    end
    isempty(stamp) || panel_stamp!(p, stamp, lim; box = stamp_box, size = stampsize)
    p
end

"""
    panel_stamp!(p, s, lim; box = false, box_width = 0.66, size = 8)

Parameter stamp of Figs. 1–2: plain text inside the axes at the top-left (manual section 0, "panel
annotations"), optionally on a white background box (Fig. 2's left column masks the field with one;
`box_width` is its width as a fraction of the panel width, ≈ 0.66 in the published figure).
"""
function panel_stamp!(p, s::AbstractString, lim; box::Bool = false, box_width = 0.66,
                      size::Integer = 8)
    x0, y0 = -0.94 * lim, 0.93 * lim
    if box
        w, h = box_width * 2 * lim, 0.10 * lim
        plot!(p, Shape([x0, x0 + w, x0 + w, x0], [y0 - h, y0 - h, y0 + 0.35h, y0 + 0.35h]);
              color = :white, linecolor = :white, linewidth = 0, label = "")
    end
    annotate!(p, x0, y0, text(s, :left, :top, size, "Computer Modern"))
    p
end

"""
    colorbar_panel(; label, lim = CBLIM, ticks = CBTICKS, ticksize = 8, labelsize = 10) -> Plots.Plot

The single tall RdBu colorbar shared by all six panels (manual section 0, "colorbars"): a vertical
gradient with the 9 labelled ticks on the RIGHT, the rotated label to the right of them, and the
TRIANGULAR `extend='both'` caps at either end. The frame is drawn by hand (rectangle + the two
caps), since a `Plots` box frame cannot be pointed.
"""
function colorbar_panel(; label::AbstractString, lim = CBLIM, ticks = CBTICKS,
                        ticksize = 8, labelsize = 10)
    z = collect(range(-lim, lim; length = 513))
    cap = 0.075 * lim                                     # height of the triangular extend caps
    p = heatmap([0.0, 1.0], z, repeat(z, 1, 2); c = RDBU, clims = (-lim, lim), colorbar = false,
                framestyle = :none, grid = false, legend = false, widen = false,
                xlims = (0.0, 1.0), ylims = (-lim - 1.35cap, lim + 1.35cap))
    # triangular caps in the extreme colours
    plot!(p, Shape([0.0, 1.0, 0.5], [lim, lim, lim + cap]);
          color = RDBU[1.0], linecolor = RDBU[1.0], linewidth = 0, label = "")
    plot!(p, Shape([0.0, 1.0, 0.5], [-lim, -lim, -lim - cap]);
          color = RDBU[0.0], linecolor = RDBU[0.0], linewidth = 0, label = "")
    # hand-drawn outline: rectangle sides + the two caps
    plot!(p, [0.0, 0.0, 0.5, 1.0, 1.0, 0.5, 0.0], [-lim, lim, lim + cap, lim, -lim, -lim - cap, -lim];
          color = :black, linewidth = 0.8, label = "")
    # ticks (short inward marks) + labels, on the right edge
    for t in ticks
        plot!(p, [0.82, 1.0], [t, t]; color = :black, linewidth = 0.7, label = "")
        annotate!(p, 1.25, t, text(string(t), :left, :vcenter, ticksize, "Computer Modern"))
    end
    annotate!(p, 3.55, 0.0, text(label, :center, :vcenter, labelsize, "Computer Modern"; rotation = 90))
    plot!(p; xlims = (0.0, 4.2))
    p
end

end # module SpiralPanels
