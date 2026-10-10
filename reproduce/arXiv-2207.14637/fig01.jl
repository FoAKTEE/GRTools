# reproduce/arXiv-2207.14637/fig01.jl — ECC::fr22-fig1.
#
# arXiv:2207.14637 Fig. 1 (`fig:cartoon`, src/Sections/2_framework.tex:12-17; reference drawing
# src/Figures/cartoon.eps, a 368 x 365 pt cairo vector cartoon on page 2): the geometry of the
# problem. The star M_* sits at the origin of the polar grid, the planet's eccentric orbit is the
# solid black line, the planet is at (r_p, phi_p) and an arbitrary disc point is (r, phi).
#
# SCHEMATIC, NO DATA. The acceptance is the same visual structure as the published cartoon under
# the vlm-figure-repro loop (style manual `_style/manuals/arXiv-2207.14637.md` section 1).
#
# ZERO physics here. The orbit is `Modes.planet_orbit` on a `Modes.EccentricOrbit`
# (ECC::fr22-kepler-orbit, FR22 (24)-(25), pericentre at M = 0 on the +x axis). Its eccentricity and
# the units a = n = 1 come from `FR22Setup` (ECC::fr22-setup-parameters): e is the middle printed
# eccentricity, e_list[2] = 0.25. The star is therefore at the focus, a*e to the right of the
# ellipse centre. The published cartoon is not to scale: its drawn ellipse has b/a = 0.90 with the
# star 0.19 a from the centre, which no Kepler ellipse satisfies (b/a = 0.90 needs e = 0.44, a focal
# offset of 0.44 a). The Kepler ellipse here has b/a = 0.97; this orbit aspect is the one structural
# difference from the cartoon found by the VLM review, left unresolved pending an owner decision
# (emri-disk obligation fr22-o-fig1-kepler-vs-cartoon). e = 0.5, the other non-circular FR22 value,
# puts the apocentre (1.5 a) outside the dashed boundary.
#
# Every other number below is PRESENTATION only (canvas, colours, line widths, label offsets, the
# schematic placement of the planet and of the disc point); none of them is an FR22 setting. The
# pink radial gradient and the greys are colour samples of the reference raster (cartoon.eps at
# 150 dpi). The golden registers the plotted geometry (orbit samples, peri/apocentre, focal offset,
# planet and disc-point coordinates, the two angles); the PNG beside it is the VLM artefact.
# `--check` renders to a scratch path and leaves the tracked PNG untouched.

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

const PAPER = "2207.14637"
const ID = "fr22_fig01"

include(joinpath(@__DIR__, "_fr22_setup.jl"))
using GrTools.AccretionDisc.Modes: EccentricOrbit, planet_orbit

const SETUP = FR22Setup()

# ───────────────────────────── presentation constants (drawing only) ─────────────────────────────
const STYLE = (
    planet_orbit_fraction = 1 // 12,   # planet drawn T_p/12 after pericentre (schematic placement)
    point_azimuth_over_pi = 7 // 10,   # azimuth of the arbitrary disc point (r, phi), in pi
    r_out = 1.4,                       # dashed disc boundary, in a; the disc point lies on it
    arc_phip = 0.367,                  # radius of the phi_p arc (and wedge), centred on the star
    # the cartoon's phi arc is NOT centred on the star: an arc of the circle with this centre
    # (offset from the star) and radius, both in a, from the +x axis to the r line
    arc_phi_centre = (-0.124, -0.065), arc_phi_radius = 0.581,
    n_orbit = 721, n_register = 73,    # orbit samples drawn / registered (every 10th)
    n_arc = 181, n_grad = 401,         # arc samples; gradient raster per axis
    canvas = (560, 560),               # px; the published drawing is 368 x 365 pt
    pad = (0.005, 0.085),              # canvas beyond the boundary: left/bottom, right/top (arrow tips)
    star_ms = 15, planet_ms = 8,       # marker sizes (cartoon radii 0.071 a and 0.035 a)
    lw_orbit = 3.1, lw_arc = 3.0, lw_axis = 2.6, lw_boundary = 2.4,
    n_dash = 131, dash_fraction = 0.72,  # dashes round the boundary and their filled fraction
    dot_ms = 1.3, dot_spacing = 0.038,  # square dots of the dotted radii (size, spacing in a)
    # swept-back axis arrowheads: tip beyond the boundary, length, half-width, notch depth (in a)
    arrow_tip = 0.08, arrow_len = 0.18, arrow_halfwidth = 0.067, arrow_notch = 0.05,
    label_size = 22,
    # label centres as measured on the cartoon, in a: offsets from the star, the planet and the
    # disc point; the two angle labels in polar form (radius, fraction of the angle)
    label_star = (-0.15, -0.15), label_planet = (0.161, 0.046), label_rp = (-0.203, -0.010),
    label_point = (0.17, -0.02), label_phi = (0.40, 0.80), label_phip = (0.235, 0.40),
)
# radial pink gradient sampled on the cartoon raster along the 225-degree ray: (r / r_out, RGB)
const PINK = [(0.0, (212, 126, 158)), (0.054, (215, 132, 163)), (0.162, (228, 151, 180)),
              (0.269, (231, 164, 189)), (0.378, (236, 179, 200)), (0.485, (239, 191, 209)),
              (0.593, (242, 205, 219)), (0.701, (245, 218, 229)), (0.809, (249, 232, 239)),
              (0.917, (252, 245, 247)), (0.971, (255, 252, 253)), (1.0, (255, 255, 255))]
const AXIS_GREY = (128, 128, 128)
const WEDGE_GREY = (90, 90, 90)

# ───────────────────────────── the geometry (package calls only) ─────────────────────────────────
v = SETUP.values
const T = typeof(v.h_p)
const ORBIT = EccentricOrbit(one(T), v.e_list[2], one(T))   # FR22 units a = n = 1, e from FR22Setup
const TWO_PI = 2T(pi)

M_grid = range(zero(T), TWO_PI; length = STYLE.n_orbit)
orbit = [planet_orbit(M, ORBIT) for M in M_grid]
xo = [o.r_p * cos(o.phi_p) for o in orbit]
yo = [o.r_p * sin(o.phi_p) for o in orbit]
peri = planet_orbit(zero(T), ORBIT)
apo = planet_orbit(T(pi), ORBIT)
x_peri, x_apo = peri.r_p * cos(peri.phi_p), apo.r_p * cos(apo.phi_p)
centre_x = (x_peri + x_apo) / 2
planet = planet_orbit(TWO_PI * T(STYLE.planet_orbit_fraction), ORBIT)
xp, yp = planet.r_p * cos(planet.phi_p), planet.r_p * sin(planet.phi_p)
phi_pt = T(pi) * T(STYLE.point_azimuth_over_pi)
r_pt = T(STYLE.r_out) * ORBIT.a
xq, yq = r_pt * cos(phi_pt), r_pt * sin(phi_pt)

println(@sprintf("[fr22_fig01] e = %.4f (FR22Setup e_list[2]), a = %.1f, n = %.1f", ORBIT.e, ORBIT.a, ORBIT.n))
println(@sprintf("[fr22_fig01] pericentre x = %.6f, apocentre x = %.6f, ellipse centre x = %.6f (focal offset a e)",
                 x_peri, x_apo, centre_x))
println(@sprintf("[fr22_fig01] planet at M = 2pi/12: r_p = %.6f, phi_p = %.4f rad (%.2f deg)",
                 planet.r_p, planet.phi_p, rad2deg(planet.phi_p)))
println(@sprintf("[fr22_fig01] disc point: r = %.3f, phi = %.4f rad (%.1f deg)", r_pt, phi_pt, rad2deg(phi_pt)))

stride = (STYLE.n_orbit - 1) ÷ (STYLE.n_register - 1)
reg = 1:stride:STYLE.n_orbit
@repro_result "e" ORBIT.e
@repro_result "a" ORBIT.a
@repro_result "orbit_x" xo[reg]
@repro_result "orbit_y" yo[reg]
@repro_result "r_pericentre" peri.r_p
@repro_result "r_apocentre" apo.r_p
@repro_result "orbit_r_min_sampled" minimum(o.r_p for o in orbit)
@repro_result "orbit_r_max_sampled" maximum(o.r_p for o in orbit)
@repro_result "orbit_y_max_sampled" maximum(yo)
@repro_result "ellipse_centre_x" centre_x
@repro_result "planet_mean_anomaly" TWO_PI * T(STYLE.planet_orbit_fraction)
@repro_result "planet_r_p" planet.r_p
@repro_result "planet_phi_p" planet.phi_p
@repro_result "planet_xy" [xp, yp]
@repro_result "point_r" r_pt
@repro_result "point_phi" phi_pt
@repro_result "point_xy" [xq, yq]

# ───────────────────────────── the figure (presentation only) ────────────────────────────────────
using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
gr()

rgb(c) = RGB(c[1] / 255, c[2] / 255, c[3] / 255)
rgba(c, alpha) = RGBA(c[1] / 255, c[2] / 255, c[3] / 255, alpha)
const PINK_GRAD = cgrad([rgb(c) for (_, c) in PINK], [s for (s, _) in PINK])

arc(r, θ1, θ2) = (θ = range(θ1, θ2; length = STYLE.n_arc); (r .* cos.(θ), r .* sin.(θ)))
wedge(r, θ1, θ2) = (c = arc(r, θ1, θ2); Shape(vcat(0.0, c[1], 0.0), vcat(0.0, c[2], 0.0)))

"Distance from the star along azimuth `θ` to the circle (centre `c`, radius `ρ`) enclosing the star."
ray_hit(c, ρ, θ) = (uc = c[1] * cos(θ) + c[2] * sin(θ); uc + sqrt(uc^2 - c[1]^2 - c[2]^2 + ρ^2))

"The cartoon's phi arc: the circle (centre `c`, radius `ρ`) between its hits on the rays at 0 and `θ`."
function offset_arc(c, ρ, θ)
    ends = [(t = ray_hit(c, ρ, φ); (t * cos(φ), t * sin(φ))) for φ in (0.0, θ)]
    α1, α2 = (atan(p[2] - c[2], p[1] - c[1]) for p in ends)
    α = range(α1, α2; length = STYLE.n_arc)
    (c[1] .+ ρ .* cos.(α), c[2] .+ ρ .* sin.(α))
end

"Swept-back (stealth) arrowhead at `tip` along the unit vector `u`: tip, wing, notch, wing."
function arrowhead(tip, u, len, w, notch)
    base = tip .- len .* u
    nrm = (-u[2], u[1])
    inner = tip .- (len - notch) .* u
    Shape([tip[1], base[1] + w * nrm[1], inner[1], base[1] - w * nrm[1]],
          [tip[2], base[2] + w * nrm[2], inner[2], base[2] - w * nrm[2]])
end

"Dotted radius from the star to `tip`, drawn as evenly spaced square dots (the cartoon's dot style)."
function dotted!(plt, tip)
    n = max(2, round(Int, hypot(tip...) / STYLE.dot_spacing))
    t = range(0.0, 1.0; length = n)
    scatter!(plt, tip[1] .* t, tip[2] .* t; marker = :rect, ms = STYLE.dot_ms, color = :black,
             markerstrokewidth = 0)
end

label(x, y, s) = annotate!(x, y, text(s, STYLE.label_size, :black, "Computer Modern"))

function draw(png)
    S = STYLE
    R = S.r_out
    lo, hi = -R - S.pad[1], R + S.pad[2]
    plt = plot(; size = S.canvas, aspect_ratio = :equal, xlims = (lo, hi), ylims = (lo, hi),
               framestyle = :none, grid = false, legend = false, ticks = false,
               background_color = :white, margin = -2Plots.mm)
    # pink disc: radial gradient, darkest at the star, white at the dashed boundary
    g = range(-R, R; length = S.n_grad)
    Z = [hypot(x, y) <= R ? hypot(x, y) / R : NaN for y in g, x in g]
    heatmap!(plt, g, g, Z; c = PINK_GRAD, clims = (0, 1), colorbar = false)
    # translucent grey x and y axes, swept-back arrowheads beyond the boundary
    # (the cartoon draws the y axis and its head darker than the x axis: line / head opacities)
    tip = R + S.arrow_tip
    stem = tip - S.arrow_len + S.arrow_notch
    for (line, t, u, a_line, a_head) in (([-R, stem], (tip, 0.0), (1.0, 0.0), 0.5, 0.45),
                                         ([-R, stem], (0.0, tip), (0.0, 1.0), 0.65, 0.8))
        xs, ys = u[1] == 1 ? (line, [0, 0]) : ([0, 0], line)
        plot!(plt, xs, ys; color = rgba(AXIS_GREY, a_line), lw = S.lw_axis)
        plot!(plt, arrowhead(t, u, S.arrow_len, S.arrow_halfwidth, S.arrow_notch); seriestype = :shape,
              fillcolor = rgba(AXIS_GREY, a_head), linecolor = rgba(AXIS_GREY, a_head), linewidth = 0)
    end
    # dashed disc boundary: short dashes drawn as one NaN-separated series
    Δ = 2pi / S.n_dash
    dash = [arc(R, k * Δ, (k + S.dash_fraction) * Δ) for k in 0:S.n_dash-1]
    bx = reduce(vcat, [vcat(d[1][[1, end]], NaN) for d in dash])
    by = reduce(vcat, [vcat(d[2][[1, end]], NaN) for d in dash])
    plot!(plt, bx, by; color = RGB(0.49, 0.49, 0.49), lw = S.lw_boundary)
    # translucent grey wedges of phi (light) and phi_p (overlap darker), bounded by red arcs
    wg = rgba(WEDGE_GREY, 0.22)
    ox, oy = offset_arc(S.arc_phi_centre, S.arc_phi_radius, phi_pt)
    plot!(plt, Shape(vcat(0.0, ox, 0.0), vcat(0.0, oy, 0.0)); seriestype = :shape, fillcolor = wg, linewidth = 0)
    plot!(plt, wedge(S.arc_phip, 0.0, planet.phi_p); seriestype = :shape, fillcolor = wg, linewidth = 0)
    plot!(plt, ox, oy; color = :red, lw = S.lw_arc)
    plot!(plt, arc(S.arc_phip, 0.0, planet.phi_p)...; color = :red, lw = S.lw_arc)
    # the orbit (package geometry)
    plot!(plt, xo, yo; color = :black, lw = S.lw_orbit)
    # dotted radii to the planet and to the disc point (arrowhead at the point)
    dotted!(plt, (xp, yp))
    u = (cos(phi_pt), sin(phi_pt))
    dotted!(plt, (xq - 0.08u[1], yq - 0.08u[2]))
    plot!(plt, arrowhead((xq, yq), u, 0.09, 0.035, 0.03); seriestype = :shape, fillcolor = :black,
          linecolor = :black, linewidth = 0)
    # star and planet
    scatter!(plt, [0.0], [0.0]; ms = S.star_ms, color = :black, markerstrokewidth = 0)
    scatter!(plt, [xp], [yp]; ms = S.planet_ms, color = :black, markerstrokewidth = 0)
    # labels (CM italics, tight subscripts, placed as in the cartoon)
    label(S.label_star..., "\$M_{\\!\\!*}\$")
    label(xp + S.label_planet[1], yp + S.label_planet[2], "\$M_{\\!\\!p}\$")
    label(xp + S.label_rp[1], yp + S.label_rp[2], "\$r_{\\!\\!p}\$")
    label(xq + S.label_point[1], yq + S.label_point[2], "\$r\$")
    label(S.label_phi[1] * cos(S.label_phi[2] * phi_pt), S.label_phi[1] * sin(S.label_phi[2] * phi_pt), "\$\\phi\$")
    label(S.label_phip[1] * cos(S.label_phip[2] * planet.phi_p),
          S.label_phip[1] * sin(S.label_phip[2] * planet.phi_p), "\$\\phi_{\\!\\!p}\$")
    return save_figure!(plt, png; id = ID)
end

const GOLD_PNG = joinpath(@__DIR__, "expected", ID * ".png")
publish = opts.mode === :regen
out_png = publish ? GOLD_PNG : joinpath(mktempdir(), ID * ".png")
saved = draw(out_png)
publish || println("[", ID, "] --check: rendered to a scratch path, ",
                   relpath(GOLD_PNG, ReproHarness.REPRO_ROOT), " untouched")
@repro_result "png_within_budget" saved.ok

exit(ReproHarness.conclude!(opts; paper = PAPER, id = ID, script = @__FILE__, rtol = 1e-12))
