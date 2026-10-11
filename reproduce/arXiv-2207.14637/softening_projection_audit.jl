# A small independent forcing audit for the candidate policies, not another bank.
include("softening_bank.jl")
using Test, DelimitedFiles
rows=[]
@testset "Softening: independent Gauss-E forcing audit" begin
    for (paper,e) in (("fr22",.5),("fr24",.12)), rule in ("fixed","moving","field")
        paper=="fr24" && rule=="field" && continue
        s=setup(paper,rule,e,80)
        modes=[F.NewtonianMode(10,8,1.),F.NewtonianMode(40,40,1.),F.NewtonianMode(80,84,1.)]
        radii=[1-e,1.,1+e]
        coarse=M.potential_bank(s.d,[s.o],radii,modes;softening=s.soft,Nt=1024,Nphi=1024)
        fine=M.potential_bank(s.d,[s.o],radii,modes;softening=s.soft,Nt=2048,Nphi=2048)
        for (j,mode) in enumerate(modes),(i,r) in enumerate(radii)
            oracle=collect(M.fr22_potential_quadrature(s.d,s.o,r,mode;softening=s.soft,panels=32,Nphi=2048))
            a=[getfield(coarse.caches[1,j],key)[i] for key in (:P,:dP,:d2P)]
            b=[getfield(fine.caches[1,j],key)[i] for key in (:P,:dP,:d2P)]
            error=maximum(abs.(b-oracle)./max.(abs.(oracle),1.))
            refinement=maximum(abs.(b-a)./max.(abs.(b),1.))
            println("AUDIT $paper $rule r=$r m=$(mode.m) l=$(mode.m+mode.k) oracle=$error refinement=$refinement")
            @test error<2e-8
            @test refinement<1e-4
            push!(rows,(paper,rule,r,mode.m,mode.m+mode.k,a...,
                (abs.(b-oracle)./max.(abs.(oracle),1.))...,
                (abs.(b-a)./max.(abs.(b),1.))...,error,refinement))
        end
    end
    # Preserve measurements even if a test above fails at testset completion.
    open(ARGS[1],"w") do io
        println(io,"paper,rule,r,m,l,Phi,dPhi,d2Phi,oracle_Phi_error,oracle_dPhi_error,oracle_d2Phi_error,refine_Phi,refine_dPhi,refine_d2Phi,oracle_scaled_error,angular_refinement_scaled")
        for row in rows;println(io,join(row,','));end
    end
    println("PROJECTION rows=$(length(rows)) max_oracle=$(maximum(row[end-1] for row in rows)) max_refinement=$(maximum(row[end] for row in rows))")
end
