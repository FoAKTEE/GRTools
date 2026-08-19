# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/npghp_moveD/TeukEvolution_run/minrun/build_fixture.jl (untracked one-shot, mtime 2026-07-04)
#               sha256(original) = 29585e66895f6a49e0d58e6bdeb67e33e1b95689ee1e3bde944cefca90e7c6f5
# fixture(s)  : td_hphc_upstream.jl
# oracle      : the out_min/ CSVs produced by run.jl above (same TeukEvolution.jl @ 3c8e47f grid)
# run         : julia --project=<TeukEvolution env> test/fixtures/capture/build_fixture.jl
#
# Paths rewritten to the ENV-overridable GRTOOLS_TEUKEVOLUTION_ROOT / GRTOOLS_REFCODE_ROOT
# (defaults reproduce the capture host).  Body otherwise VERBATIM from the original one-shot.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const TEUK_ROOT = get(ENV, "GRTOOLS_TEUKEVOLUTION_ROOT",
                      joinpath(get(ENV, "GRTOOLS_REFCODE_ROOT", "/data/haiyangw/claude/emri-disk/ref-code"),
                               "TeukEvolution.jl"))
# ----------------------------------------------------------------------------------------------

# Build a plain-text characterization fixture from the out_min CSVs.
cd(TEUK_ROOT)
using Printf
include(joinpath(TEUK_ROOT, "src") * "/Sphere.jl");  import .Sphere

# --- reproduce grid exactly as TeukEvolution.launch does ---
nx=48; ny=24; cl=1.0; cfl=0.5; bhs=0.7; bhm=1.0
minr = bhm*(1 + sqrt(1 - (bhs/bhm)^2))
maxR = 1/minr
dr   = maxR/(nx-1)
dt   = min(cfl*dr*bhm^2, 6/ny^2)
Rv   = [dr*i for i=0:nx-1]
Yv   = Sphere.Y_vals(ny)              # GL roots, y = -cosθ
rvals = [Rv[i]==0 ? Inf : cl^2/Rv[i] for i=1:nx]   # BL areal radius

# --- CSV reader: line -> (t, values[i,j]) ---
function read_csv(path)
    rows = Vector{Tuple{Float64,Matrix{Float64}}}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        f = split(line, ',')
        t = parse(Float64, f[1]); nxr = parse(Int,f[2]); nyr = parse(Int,f[3])
        @assert nxr==nx && nyr==ny
        V = Array{Float64}(undef,nx,ny)
        for i=1:nx, j=1:ny
            V[i,j] = parse(Float64, f[3 + (i-1)*ny + j])
        end
        push!(rows,(t,V))
    end
    return rows
end

fre = read_csv("out_min/lin_f_re_2.csv")
fim = read_csv("out_min/lin_f_im_2.csv")
pre = read_csv("out_min/lin_p_re_2.csv")
pim = read_csv("out_min/lin_p_im_2.csv")
times = [r[1] for r in fre]

# --- subsample ---
isub = [1,12,24,36,48]
jsub = [6,12,18,24]

io = open("../../te_evolution_fixture.txt","w")
pr(s...) = (println(io, s...); )
pr("# ============================================================================")
pr("# Move-D characterization oracle: TeukEvolution.jl linear psi4 evolution")
pr("# ============================================================================")
pr("# PROVENANCE")
pr("#   ref-code repo   : ref-code/TeukEvolution.jl (Ripley, JLRipley314/TeukEvolution.jl)")
pr("#   git rev         : 3c8e47f4202a92778e4d21bd0bfc56839288a1b1 (\"update local changes to remote\")")
pr("#   run driver      : scratchpad/npghp_moveD/TeukEvolution_run/minrun/run.jl")
pr("#   julia version   : ", VERSION)
pr("#   resolved deps   : FastGaussQuadrature v0.5.1, HDF5 v0.16.16, Jacobi v0.5.1, Polynomials v2.0.25")
pr("#   threads         : 1")
pr("# CONFIG")
pr("#   runtype=linear_field  id_kind=gaussian  psi_spin=-2  precision=Float64")
pr("#   nx=$nx  ny=$ny  nt=600  ts=200   cl=$cl  cfl=$cfl  bhm=$bhm  bhs(a)=$bhs")
pr("#   m_vals=[2]  id_l_ang=[2]  id_ru=[4.0]  id_rl=[2.5]  id_width=[1.5]  ingoing=true  id_amp=[[0.1,0.0]]")
pr("# DERIVED GRID")
pr("#   minr(r+)=$(minr)   maxR=1/r+=$(maxR)   dr=$(dr)   dt=min(cfl*dr, 6/ny^2)=$(dt)")
pr("#   Rv[i]=dr*(i-1),  i=1..nx  (R=0 is scri+, R=maxR is horizon);  r_BL[i]=cl^2/Rv[i]")
pr("#   Yv[j] = Gauss-Legendre roots on [-1,1], y=-cos(theta)")
pr("# GRID VALUES AT SUBSAMPLE INDICES")
for i in isub
    pr("#   i=$i  R=$(Rv[i])  r_BL=$(rvals[i])")
end
for j in jsub
    pr("#   j=$j  y=$(Yv[j])  cos(theta)=$(-Yv[j])")
end
pr("# SNAPSHOT TIMES (t/bhm): ", join(times, ", "))
pr("# ----------------------------------------------------------------------------")
pr("# COLUMNS: field  time_index  t  i  j  value_re  value_im")
pr("#   field in {lin_f, lin_p};  time_index 1..4 -> t above")
pr("# ============================================================================")

for (nm, dre, dim) in (("lin_f",fre,fim), ("lin_p",pre,pim))
    for (ti,t) in enumerate(times)
        for i in isub, j in jsub
            vre = dre[ti][2][i,j]; vim = dim[ti][2][i,j]
            pr(rpad(nm,6), "  ", ti, "  ", @sprintf("%.17e", t), "  ",
               rpad(i,2), "  ", rpad(j,2), "  ",
               @sprintf("% .17e", vre), "  ", @sprintf("% .17e", vim))
        end
    end
end
close(io)
println("WROTE ../../te_evolution_fixture.txt")
