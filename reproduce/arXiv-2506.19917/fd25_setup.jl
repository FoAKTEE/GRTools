using AccretionDisc

"FD25 2_linear_theory.tex:80-101,189: fixed softening, code units GM=a=n=Sigma_p=1."
function fd25_disc_setup(::Type{T}=Float64) where {T<:AbstractFloat}
    [AccretionDisc.NewtonianDisc{T}(h0=T(3)/50,p=T(p)/2,q=T(q),beta=zero(T),
        soft=T(9)/500,indirect=0) for (p,q) in ((3,1),(3,0),(1,0))]
end
