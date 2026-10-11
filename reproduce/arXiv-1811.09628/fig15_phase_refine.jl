# MR19 Fig. 15: one phase-gradient secant refinement of the fiducial m=10 mode.
# All physics, conventions and diagnostic reductions live in AccretionDisc.Modes.
# Preregistered comparison: emri-disk evidence/mr19-fig15/TOLERANCES.json.
# julia --project=reproduce -t1 reproduce/arXiv-1811.09628/fig15_phase_refine.jl --check
include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
opts=ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts,:heavy)
using GrTools
using GrTools.AccretionDisc: mr19_phase_refinement
using DelimitedFiles: writedlm

println("[repro] MR19 m=10, 60000 points, one phase-gradient step"); flush(stdout)
t0=time()
figure=mr19_phase_refinement()
println("[repro] solve wall=$(time()-t0) s"); flush(stdout)
s=figure.result
idx=1:30:length(figure.r)
@repro_result "r" figure.r[idx]
@repro_result "before" figure.before[idx]
@repro_result "after" figure.after[idx]
@repro_result "C_before" collect(s.C_before)
@repro_result "C_after" collect(s.C_after)
@repro_result "X_before" collect(s.X_before)
@repro_result "X_after" collect(s.X_after)
@repro_result "X_ratio" collect(s.X_after ./ s.X_before)
@repro_result "relative_C_change" collect(s.relative_C_change)
@repro_result "relative_state_change" collect(s.relative_state_change)
@repro_result "boundary_indices" collect(s.indices)

# Optional full-resolution verification output; never a physics input or a golden substitute.
if haskey(ENV,"MR19_FIG15_EVIDENCE")
    dir=ENV["MR19_FIG15_EVIDENCE"]
    writedlm(joinpath(dir,"computed.csv"),hcat(figure.r,figure.before,figure.after),',')
    metrics=Dict("relative_C_change"=>collect(s.relative_C_change),
                 "relative_state_change"=>collect(s.relative_state_change),
                 "X_before"=>collect(s.X_before),"X_after"=>collect(s.X_after),
                 "C_before"=>[[real(c),imag(c)] for c in s.C_before],
                 "C_after"=>[[real(c),imag(c)] for c in s.C_after])
    write(joinpath(dir,"metrics.json"),ReproHarness.json_write(metrics))
end

# Presentation only: same touching panels, shared axes and labels as the source vector figure.
# Python is used for the Matplotlib Computer Modern mathtext and precise shared-panel layout.
if opts.mode === :regen
    mktemp() do table,io
        writedlm(io,hcat(figure.r,figure.before,figure.after),','); close(io)
        png=joinpath(@__DIR__,"fig15_phase_refine.png")
        presentation=raw"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter, MultipleLocator, FormatStrFormatter
plt.rcParams.update({'font.family':'serif','font.serif':['cmr10'],'mathtext.fontset':'cm','axes.unicode_minus':False,'axes.formatter.use_mathtext':True,'font.size':14.3462,'axes.linewidth':1.5,'xtick.direction':'in','ytick.direction':'in'})
a=np.loadtxt(sys.argv[1],delimiter=',')
fig,axs=plt.subplots(2,1,sharex=True,figsize=(3.875,4.4583),gridspec_kw={'hspace':0})
fig.subplots_adjust(left=47.4/279,right=277.75/279,top=1-1.45/321,bottom=1-289.4/321,hspace=0)
for i,ax in enumerate(axs):
    ax.plot(a[:,0],a[:,i+1],color='black',lw=1.0)
    ax.set(xscale='log',xlim=(.05,5),ylim=(-.045,.015))
    ax.set_yticks(np.arange(-.04,.011,.01))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_minor_locator(MultipleLocator(.005))
    ax.tick_params(which='both',top=True,right=True,length=3.15,width=1.5)
    ax.tick_params(which='minor',length=1.6)
    ax.xaxis.set_minor_locator(LogLocator(base=10,subs=np.arange(2,10)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.text(.02,.95,['(a) Before Refinement','(b) After Refinement'][i],transform=ax.transAxes,va='top',fontsize=11.9552)
axs[0].tick_params(labelbottom=False)
axs[1].set_xticks([.1,.3,1.,3.],['0.1','0.3','1.0','3.0'])
axs[1].set_xlabel(r'$r/r_\mathrm{p}$',labelpad=1)
fig.text(.028,.55,r'$\epsilon_{10}/k_{10}$',rotation=90,va='center',ha='center')
fig.savefig(sys.argv[2],dpi=144,bbox_inches='tight',pad_inches=.02)
"""
        run(`python3 -c $presentation $table $png`)
        filesize(png)<=200_000 || error("figure exceeds 200 KB")
        println("[repro] wrote $png ($(filesize(png)) bytes)")
    end
end
exit(ReproHarness.conclude!(opts;paper="1811.09628",id="fig15_phase_refine",script=@__FILE__,rtol=1e-7,atol=1e-10))
