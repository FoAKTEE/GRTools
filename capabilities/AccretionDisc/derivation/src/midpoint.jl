"""
Independent C-64 certificate for arXiv:2608.17003 (47), (64)-(65).

Uses the existing R3 symbolic verifier and the actual SpiralWaves assembly and
solver. This module belongs to the standalone derivation environment only.
"""
module MidpointCertificate

using Symbolics
using LinearAlgebra
using AccretionDisc.SpiralWaves: midpoint_system, solve_midpoint

# The R3 verifier itself, from the R3Derivation instance the including module has
# already loaded (runtests.jl, certify_midpoint.jl, test_midpoint.jl): one alphabet,
# one Certificate type.
using ..R3Derivation: Certificate, certify_symbolic

export midpoint_certificates, midpoint_order_check, MIDPOINT_STATEMENT

const MIDPOINT_STATEMENT = raw"""
Sources: mirrored arXiv-2608.17003/src/main.tex:609-611 (47),
742-767 (64)-(65), label eq:discrete-master-equation; (65) defines h=r[i+1]-r[i].
Runtime: AccretionDisc/src/SpiralWaves/solver.jl:428-439 assembles the cells;
561-566 evaluates C,S at (r[i]+r[i+1])/2; 457-460 calls the actual linear solver.

Scheme: A_i y[i+1] = B_i y[i] + h_i S(m_i),
A_i=I-h_i C(m_i)/2, B_i=I+h_i C(m_i)/2, m_i=(r[i]+r[i+1])/2.
Thus G_i=A_i^(-1)B_i and g_i=A_i^(-1)h_i S(m_i).

For any finite-dimensional real or complex system y'=C(r)y+S(r), define the
integrated exact-solution defect rho_i=A_i y(r[i+1])-B_i y(r[i])-h_i S(m_i).
Taylor integration about m gives integral(y')=h y'(m)+h^3 y'''(m)/24+O(h^5).
The endpoint average is y(m)+h^2 y''(m)/8+O(h^4), giving
rho_i=h_i^3 D(m_i)+O(h_i^5), D=(y'''-3 C y'')/24.
The componentwise jet certificate is dimension-independent.
Using y'=Cy+S and the product rule without assuming C C'=C' C:
24 D=(C''+2 C' C-2 C C'-2 C^3)y+(2 C'-2 C^2)S-2 C S'+S''.
The one-step error EXACT MINUS NUMERICAL from the exact left value is
A_i^(-1)rho_i=h_i^3 D(m_i)+O(h_i^4). The opposite error convention negates D.
The derivative-form residual rho_i/h_i is O(h_i^2).

Assumptions: a fixed compact interval [a,b], C,S in C^4 with bounded derivatives
(hence y in C^5) for the stated O(h^5) centered remainder; C,S in C^2 suffice
for local O(h^3) and global O(H^2). Actual midpoint values are used, and roundoff
is excluded from the analytic statement. Constants need not be uniform as a
regularization parameter vanishes or a corotation singularity is approached.

Global IVP proof on arbitrary increasing meshes: let H=max h_i, M=sup ||C||,
L=b-a, h_i M/2 <= theta < 1, and ||rho_i|| <= K h_i^3. The Neumann bound gives
||A_i^(-1)|| <= 1/(1-theta), and
||G_i|| <= (1+h_i M/2)/(1-h_i M/2) <= exp(h_i M/(1-theta)).
The exact error recurrence e[i+1]=G_i e[i]+A_i^(-1)rho_i, e[0]=0, implies
max ||e[i]|| <= exp(M L/(1-theta))*K*L*H^2/(1-theta), since sum h_i^3<=L H^2.
No adjacent-step-ratio condition is needed. Initial error O(H^2) preserves order.
For the globally assembled two-point BVP the same order additionally requires
a mesh-uniform stability bound from cell defects in the sum norm and boundary
errors to the solution max norm, with exact or O(H^2) boundary data. This is a
conditional discretization theorem, not a certificate of physical WKB boundaries
or of uniform stability near singular BVPs. Manufactured tests use exact split BCs.

Code/paper identity: symbolic inputs are passed THROUGH midpoint_system; its
left block, right block, and source are compared with the paper before comparing
G and g. For the 2x2 update the only excluded locus is det(I-h C/2)=0.
Degree accounting over Q[h,C11,C12,C21,C22,S1,S2]: A,B entries degree<=2;
adj(A) degree<=2; det(A) degree<=4; numerators of G,g degree<=4.
A cross-multiplied update difference has total degree<=8 (two denominators).
Every residual is reduced symbolically to zero; no random identity test or
probabilistic fallback is used, so no sampling bound or seed is needed.
"""

"Symbolically derive the defect and compare the paper with the executed assembly."
function midpoint_certificates(; verbose::Bool = true)
    certs = Certificate[]
    @variables h s
    yjet = Symbolics.variables(:yjet, 0:4)
    cyjet = Symbolics.variables(:cyjet, 0:4)
    # One arbitrary component: cyjet[k+1]=(C(m)y^(k)(m))_j.
    yp = sum(yjet[k+1] * (h/2)^k / factorial(k) for k in 0:4)
    ym = sum(yjet[k+1] * (-h/2)^k / factorial(k) for k in 0:4)
    cyp = sum(cyjet[k+1] * (h/2)^k / factorial(k) for k in 0:4)
    cym = sum(cyjet[k+1] * (-h/2)^k / factorial(k) for k in 0:4)
    integrated = yp - ym
    push!(certs, certify_symbolic("midpoint quadrature through degree four", "(47),(64)",
        integrated - h*yjet[2] - h^3*yjet[4]/24; verbose))
    rho = Symbolics.substitute(integrated - h*(cyp+cym)/2 - h*s,
                               Dict(yjet[2] => cyjet[1]+s))
    leading = yjet[4]/24 - cyjet[3]/8
    # The h^5 term belongs to the remainder; it is known for these degree-4 jets.
    push!(certs, certify_symbolic("integrated defect h^3 D through degree four", "(64)",
        rho - h^3*leading + h^5*cyjet[5]/384; verbose))
    third = Symbolics.expand_derivatives((Differential(h)^3)(rho))
    defect = Symbolics.substitute(third, Dict(h => 0)) / 6

    C = Symbolics.variables(:C, 1:2, 1:2)
    C1 = Symbolics.variables(:C1, 1:2, 1:2)
    C2 = Symbolics.variables(:C2, 1:2, 1:2)
    S = Symbolics.variables(:S, 1:2)
    S1 = Symbolics.variables(:S1, 1:2)
    S2 = Symbolics.variables(:S2, 1:2)
    y = Symbolics.variables(:y, 1:2)
    dy = C*y+S
    d2y = C1*y+C*dy+S1
    d3y = C2*y+2C1*dy+C*d2y+S2
    expanded = ((C2+2C1*C-2C*C1-2C^3)*y+(2C1-2C^2)*S-2C*S1+S2)/24
    push!(certs, certify_symbolic("noncommuting coefficient leading term", "(47),(64)",
        (d3y-3C*d2y)/24-expanded; verbose))

    ident = Matrix{Num}(I, 2, 2)
    paper_A, paper_B, paper_b = ident-h*C/2, ident+h*C/2, h*S
    assembled, rhs = midpoint_system(Num[0,h], [C], [S], Num[1,0], Num(0),
                                     Num[0,1], Num(0))
    code_A = Matrix(assembled[2:3,3:4])
    code_B = -Matrix(assembled[2:3,1:2])
    code_b = rhs[2:3]
    push!(certs, certify_symbolic("actual assembly blocks equal paper", "(64)-(65)",
        vcat(vec(code_A-paper_A), vec(code_B-paper_B), code_b-paper_b); verbose))
    paper_G, paper_g = inv(paper_A)*paper_B, inv(paper_A)*paper_b
    code_G, code_g = inv(code_A)*code_B, inv(code_A)*code_b
    push!(certs, certify_symbolic("actual homogeneous update equals paper", "(64)",
        vec(code_G-paper_G); verbose))
    push!(certs, certify_symbolic("actual source update equals paper", "(64)",
        code_g-paper_g; verbose))
    push!(certs, certify_symbolic("propagation recurrence from assembled row", "(64)",
        vcat(vec(code_A*code_G-paper_B), code_A*code_g-paper_b); verbose))
    (; certificates = certs, defect, leading, quadrature_only = yjet[4]/24,
       source_update = code_g)
end

"""
    midpoint_order_check(::Type{T}=Float64; variable_steps=false)

Five nested meshes (16,32,64,128,256 cells) on [0,1], using the actual
`solve_midpoint` entry point and exact split boundary values. The variable mesh
starts with alternating widths 1/24,1/12, then bisects EVERY cell four times,
so H halves exactly and nonuniformity persists. BigFloat precision is inherited.
"""
function midpoint_order_check(::Type{T} = Float64; variable_steps::Bool = false) where {T<:AbstractFloat}
    Cf(r) = Complex{T}[r/5+im*T(1)/10 (1+r)/7+im*r/8;
                       -r/6+im*T(1)/9 -r/4-im*T(1)/11]
    yex(r) = Complex{T}[exp(r)+im*sin(r), cos(r)+im*exp(-r)]
    dyex(r) = Complex{T}[exp(r)+im*cos(r), -sin(r)-im*exp(-r)]
    base = variable_steps ? T[(3*(i÷2)+(isodd(i) ? 1 : 0))//24 for i in 0:16] :
                            T[i//16 for i in 0:16]
    rows = NamedTuple[]
    for level in 0:4
        subdivisions = 2^level
        r = T[base[i]+(base[i+1]-base[i])*T(j)/subdivisions
              for i in 1:16 for j in 0:subdivisions-1]
        push!(r, one(T))
        rm = [(r[i]+r[i+1])/2 for i in 1:length(r)-1]
        Cmid = Cf.(rm)
        Smid = [dyex(x)-C*yex(x) for (x,C) in zip(rm,Cmid)]
        row_in, row_out = Complex{T}[1,0], Complex{T}[0,1]
        yin, yout = yex(first(r))[1], yex(last(r))[2]
        sol = solve_midpoint(r, Cmid, Smid, row_in, yin, row_out, yout)
        err = maximum(abs, sol-hcat(yex.(r)...))
        A,b = midpoint_system(r, Cmid, Smid, row_in, yin, row_out, yout)
        algebraic_residual = norm(A*vec(sol)-b, Inf)/max(one(T), norm(b, Inf))
        order = isempty(rows) ? T(NaN) : log2(rows[end].error/err)
        push!(rows, (; cells = length(r)-1, hmax = maximum(diff(r)), error = err, order,
                      algebraic_residual, solution_type = eltype(sol), precision_bits = precision(T)))
    end
    rows
end

end # module MidpointCertificate
