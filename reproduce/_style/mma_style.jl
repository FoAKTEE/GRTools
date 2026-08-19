# Shared Mathematica-matching plot style — for EXACT figure reproduction (figure_reproduction_criterion).
# MMA_RAINBOW is the EXACT ColorData["Rainbow"] colormap used by the papers, extracted verbatim via
#   wolframscript -code 'Do[c=ColorData["Rainbow"][t];Print[...],{t,0,1,0.0625}]'
# (presentation attribute only — the physics data is computed by the Julia port). Font set to Computer Modern /
# serif to match the papers' LaTeX-rendered titles & axis labels.
using Plots
using Plots: RGB, cgrad

const _MMA_RAINBOW_RGB = [
    (0.471412, 0.108766, 0.527016), (0.311060, 0.117580, 0.664469),
    (0.250728, 0.225386, 0.769152), (0.244080, 0.361242, 0.816084),
    (0.266122, 0.486664, 0.802529), (0.305919, 0.585575, 0.739666),
    (0.360480, 0.655759, 0.645692), (0.429842, 0.701849, 0.540321),
    (0.513417, 0.729920, 0.440682), (0.607651, 0.743718, 0.358588),
    (0.705038, 0.742591, 0.299167), (0.794549, 0.721158, 0.260829),
    (0.863512, 0.670771, 0.236564), (0.901014, 0.582826, 0.216542),
    (0.902853, 0.453964, 0.192014), (0.878107, 0.293208, 0.160481),
    (0.857359, 0.131106, 0.132128),
]
const MMA_RAINBOW = cgrad([RGB(r, g, b) for (r, g, b) in _MMA_RAINBOW_RGB])

# Exact ColorData["TemperatureMap"] — the diverging blue→white→yellow→red map used by the FSl{1,2,3}0k (l,k) figures.
const _MMA_TEMPERATUREMAP_RGB = [
    (0.178927, 0.305394, 0.933501), (0.276291, 0.407730, 0.939046),
    (0.381032, 0.504453, 0.945500), (0.500578, 0.605431, 0.953827),
    (0.642359, 0.720535, 0.964988), (0.775578, 0.824607, 0.978266),
    (0.877842, 0.905431, 0.988211), (0.947822, 0.960607, 0.973208),
    (0.984192, 0.987731, 0.911643), (0.992509, 0.991171, 0.773801),
    (0.993893, 0.989345, 0.576615), (0.983368, 0.955559, 0.389888),
    (0.955963, 0.863115, 0.283425), (0.917161, 0.709278, 0.252204),
    (0.881316, 0.553966, 0.222680), (0.848134, 0.370981, 0.193726),
    (0.817319, 0.134127, 0.164218),
]
const MMA_TEMPERATUREMAP = cgrad([RGB(r, g, b) for (r, g, b) in _MMA_TEMPERATUREMAP_RGB])

# Exact Mathematica default line colors ColorData[97][1,2] — the blue & gold used by the line/marker figures.
const MMA_BLUE  = RGB(0.368417, 0.506779, 0.709798)
const MMA_GOLD  = RGB(0.880722, 0.611041, 0.142051)
const MMA_GREEN = RGB(0.560181, 0.691569, 0.194885)   # ColorData[97][3] — the "Analytical" curve
const MMA_RED   = RGB(0.922526, 0.385626, 0.209179)   # ColorData[97][4]
const MMA_PURPLE = RGB(0.528488, 0.470624, 0.701351)  # ColorData[97][5]
const MMA_BROWN = RGB(0.772079, 0.431554, 0.102387)   # ColorData[97][6]

# Mathematica-like defaults: serif/CM font, boxed frame, no grid (matrix plots have no gridlines).
function mma_defaults!()
    default(fontfamily = "Computer Modern", grid = false, framestyle = :box,
            titlefontsize = 14, guidefontsize = 13, tickfontsize = 11,
            colorbar_titlefontsize = 12, foreground_color_legend = nothing)
    return nothing
end
