# FR22 Figs 2–4 row 3: computed ZZ22 tracks. The package owns the wavelet physics.
# Usage: julia --project=reproduce reproduce/arXiv-2207.14637/_zz22_wavelets.jl <outdir>
module ZZ22Wavelets
using Printf
using GrTools.AccretionDisc.Modes: EccentricOrbit, zz22_track
include("_fr22_setup.jl")

const SETUP = FR22Setup()  # softening is unresolved and unused by the kinematic wavelets

"FR22 snapshots in the fifth orbit; emission from pericentre t=0 (landed vector-path inference)."
function track(e, phase, dir; setup=SETUP)
    v = setup.values
    T = typeof(v.h_p)
    period = 2T(pi)
    zz22_track(EccentricOrbit(one(T), T(e), one(T)),
               (v.wavelet_orbits-1+T(phase))*period, dir;
               cs=v.h_p, dt=period*v.wavelet_interval[1]/v.wavelet_interval[2])
end

function main(outdir; setup=SETUP)
    mkpath(outdir)
    for e in setup.values.e_list, phase in setup.values.snapshots_text
        path = joinpath(outdir, @sprintf("zz22_e%g_T%.2f.csv", e, phase))
        open(path, "w") do io
            println(io, "family,k,log10r,y")
            for (fam, dir) in (("inner", -1), ("outer", +1)), (k, x, y) in track(e, phase, dir; setup)
                @printf(io, "%s,%d,%.6f,%.6f\n", fam, k, x, y)
            end
        end
        println(path)
    end
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    ZZ22Wavelets.main(only(ARGS))
end
