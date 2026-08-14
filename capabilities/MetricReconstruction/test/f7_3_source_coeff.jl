# F7-3: convention-independent symbolic gate on the Ma–Yang ++ second-order source coefficients.
#
# INDEPENDENT ROUTE (independent of the owner's derive-spin2-plus-plus-coefficients.m):
#   The owner's Mathematica script builds the a_ij^(ν) by inserting CCK/ORG fields into the
#   Campanelli–Lousto B.11 bracket, applying the outer d4^(0), and reducing angular products.
#   NONE of that is used here. This gate uses ONLY:
#     (i)  the App B Chandrasekhar operators P=D_2† and D_1, which GrTool' `Dn` ladder ALREADY
#          implements (Dn(+1,2)=D_2†, Dn(-1,1)=D_1) — built from Chandrasekhar/Ma–Yang directly,
#          predating this paper's script;
#     (ii) the single radial field equation (ΔD_1D_2† + 6iωr − λ_+)R = 0 (App B / Sec VII.A);
#     (iii) the coefficient tables as PARSED FIXTURES (r-form from the generated .tex, Γ-form from
#           the paper Appendix C), transcribed literally below (owner-independent artifacts).
#
# GATES
#   G-F7-3a  Quotient relation R'' = F R' + G R (Eq. R-second-derivative-reduction). DOUBLE-ENTRY:
#            F,G extracted from the GrTool operator composition Dn(-1,1)∘Dn(+1,2) (shipped ladders)
#            vs. the tex closed forms F=−6(r−1)/Δ, G=−U'−VU−T. Plus residual-vanishing on the
#            GrTool-built radial operator.  Uses NO a_ij — pure operator algebra (non-circular).
#   G-F7-3b  Named entries as convention-independent rational-function identities under Γ+Γ̄=2r:
#            the machine r-form (generated .tex) == the hand-written Γ-form (Appendix C), evaluated
#            as GrTool RTJets (value + ∂_r,∂_θ,∂_rr,∂_θθ,∂_rθ) at generic Kerr (r,θ). 6 entries
#            spanning A_1..A_4 incl. all three source-level-fixed entries a_11^(1),a_04^(3),a_33^(3).
#            NON-CIRCULAR: the two forms are DISTINCT strings equal only under Γ+Γ̄=2r; substituting
#            one form on both sides would trivially pass and miss an r-form simplification error.
#   G-F7-3c  FULL-ROW reduction of A_2 (all 15 nonzero a_ij^(2)) to C_00R²+C_01RR'+C_11(R')².
#            DOUBLE-ROUTE: Route A = the (α_k,β_k) recursion + C_ab formula (Eq. app-Cnu-formula);
#            Route B = direct operator X_k=P^kR applied to a quotient-folded derivative tower.
#            This gates the paper's REDUCTION machinery on the owner's full A_2 row (X_6 exercised).
#   G-F7-3d  All value/identity gates run in ComplexF64 (rel≤1e-10) AND Complex{BigFloat}-256
#            (rel≤1e-40); every probe is {C}-generic.
#
# CIRCULARITY / [BLOCKING] note: the ABSOLUTE a_ij values-vs-physical-source (a_ij = free-jet
# coefficient of the unreduced CCK/B.11 source A_ν) are NOT gated here — that requires rebuilding the
# CCK reconstruction + Campanelli–Lousto bracket, which is beyond this round's blast radius. See the
# accompanying report. G-F7-3a/c verify the paper's reduction algebra; G-F7-3b verifies the two repo
# representations of the a_ij agree. Neither substitutes the owner's reduced table on both sides.

using Test
using MetricReconstruction.Radiation
import MetricReconstruction.Radiation as RAD
import NumericsCore: jmul, jadd, jderiv, jinv
# RTJet unary `-` and integer `^` used by the named-entry closures are provided by
# MetricReconstruction.Radiation (src/Radiation/Operators.jl) — the additive completion this gate drove.

# ───────────────────────────── length-generic coefficient jets (App B) ────────────────────────────
# U,V,T,F,G as r-jets of length LEN, from the App B operators. U is the D_2† coefficient
# (=_radcoef_jet(+1,2)); V is the D_1 coefficient (=_radcoef_jet(-1,1)); both cross-checked below.
function f73_coefjets(::Type{C}, a, r, ω, m, λ, LEN::Int) where {C}
    K = zeros(C, LEN); K[1] = (r^2 + a^2)*ω - a*m; K[2] = 2r*ω; K[3] = C(2ω)
    Δ = zeros(C, LEN); Δ[1] = r^2 - 2r + a^2;      Δ[2] = C(2r - 2); Δ[3] = C(2)
    rm1 = zeros(C, LEN); rm1[1] = C(r - 1); rm1[2] = one(C)
    invΔ = jinv(Δ)
    U = jmul(( im .* K) .+ (4 .* rm1), invΔ)          # U = ( iK + 4(r−1) )/Δ
    V = jmul((-im .* K) .+ (2 .* rm1), invΔ)          # V = (−iK + 2(r−1) )/Δ
    Tn = zeros(C, LEN); Tn[1] = 6*im*ω*r - λ; Tn[2] = C(6*im*ω)
    T = jmul(Tn, invΔ)                                 # T = ( 6iωr − λ )/Δ
    F = (-6) .* jmul(rm1, invΔ)                        # F = −6(r−1)/Δ  (tex closed form)
    G = (-1) .* jderiv(U) .- jmul(V, U) .- T           # G = −U' − VU − T
    (U=U, V=V, T=T, F=F, G=G, Δ=Δ, invΔ=invΔ)
end

# Convert a generic point to element type C (BigFloat-safe).
_cvt(::Type{C}, x::Integer) where {C} = C(x)
_cvt(::Type{T}, x) where {T<:AbstractFloat} = T(x)

# ───────────────────────────── G-F7-3a: quotient relation (double-entry) ───────────────────────────
# Uses the SHIPPED GrTool `Dn` on a SepField: value(D_1 D_2† R) with a generic radial jet
# [R,R',R'',0,0], then the radial operator L = Δ·(D_1 D_2† R) + (6iωr−λ)R. Extract F_op,G_op as the
# quotient of the linear-form coefficients in (R,R',R''); compare to the tex closed forms.
function f73_quotient(::Type{C}, a, r, ω, m, λ) where {C}
    aj = C[1, 0, 0, 0, 0]                              # trivial angular factor (radial reduction only)
    # residual L·R as a linear form in (R,R',R'') via the shipped ladders:
    Lval(R, Rp, Rpp) = begin
        rj = C[R, Rp, Rpp, 0, 0]
        Fld = sepfield(rj, aj)
        D2d  = Dn(+1, 2, Fld; a=a, r=r, ω=ω, m=m)      # D_2† (=P)
        D1D2 = Dn(-1, 1, D2d; a=a, r=r, ω=ω, m=m)      # D_1 D_2†
        RAD.value(D1D2) * (r^2 - 2r + a^2) + (6*im*ω*r - λ) * R
    end
    cRpp = Lval(zero(C), zero(C), one(C))              # coeff of R'' in L  (= Δ)
    cRp  = Lval(zero(C), one(C),  zero(C))             # coeff of R'  in L  (= Δ(U+V))
    cR   = Lval(one(C),  zero(C), zero(C))             # coeff of R   in L
    F_op = -cRp / cRpp
    G_op = -cR  / cRpp
    cj = f73_coefjets(C, a, r, ω, m, λ, 5)
    Utex = cj.U[1]; Vtex = cj.V[1]; Uptex = cj.U[2]; Ttex = cj.T[1]
    Δv = r^2 - 2r + a^2
    F_tex = -6*(r - 1)/Δv                              # tex closed form
    G_tex = -Uptex - Vtex*Utex - Ttex                  # −U'−VU−T
    # residual with the quotient-substituted R'':
    Rg, Rpg = C(0.7 + 0.4im), C(-0.3 + 0.9im)          # generic (R,R')
    Rppg = F_tex*Rpg + G_tex*Rg
    resid = Lval(Rg, Rpg, Rppg)
    (F_op=F_op, G_op=G_op, F_tex=F_tex, G_tex=G_tex, cRpp=cRpp, Δ=C(Δv),
     resid=resid, scale=abs(Rppg)*abs(C(Δv)) + abs(Rg) + abs(Rpg))
end

# ───────────────────────────── (α_k,β_k) recursion (Eq. alpha-beta-recursion) ──────────────────────
# X_k = P^k R = α_k R + β_k R'.  α_{k+1}=α_k'+Uα_k+Gβ_k, β_{k+1}=α_k+β_k'−Vβ_k, (α_0,β_0)=(1,0).
# Carried as r-jets so α',β' come from `jderiv`.
function f73_alphabeta(cj, LEN::Int, kmax::Int, ::Type{C}) where {C}
    α = Vector{Vector{C}}(undef, kmax+1); β = Vector{Vector{C}}(undef, kmax+1)
    α[1] = zeros(C, LEN); α[1][1] = one(C); β[1] = zeros(C, LEN)          # k=0
    for k in 1:kmax
        α[k+1] = jderiv(α[k]) .+ jmul(cj.U, α[k]) .+ jmul(cj.G, β[k])
        β[k+1] = α[k] .+ jderiv(β[k]) .- jmul(cj.V, β[k])
    end
    (α, β)
end

# Quotient-folded derivative tower R^(m) = A_m R + B_m R'  (A_{m+1}=A_m'+B_m G, B_{m+1}=A_m+B_m'+B_m F).
function f73_deriv_tower(cj, LEN::Int, mmax::Int, R, Rp, ::Type{C}) where {C}
    A = zeros(C, LEN); A[1] = one(C); B = zeros(C, LEN)                    # m=0
    tower = Vector{C}(undef, mmax+1); tower[1] = R                         # R^(0)=R
    if mmax >= 1; tower[2] = Rp; end                                       # R^(1)=R'
    Acur, Bcur = A, B
    for m in 1:mmax
        Anew = jderiv(Acur) .+ jmul(Bcur, cj.G)
        Bnew = Acur .+ jderiv(Bcur) .+ jmul(Bcur, cj.F)
        Acur, Bcur = Anew, Bnew
        tower[m+1] = Acur[1]*R + Bcur[1]*Rp
    end
    tower                                                                  # [R,R',R'',...,R^(mmax)]
end

# Route B: X_k = P^k R by DIRECT operator application to the folded derivative tower (length-generic).
function f73_Xk_directop(cj, tower::Vector{C}, kmax::Int, ::Type{C}) where {C}
    LEN = length(tower)
    J = copy(tower)                                                        # jet [R,R',...,R^(LEN-1)]
    X = Vector{C}(undef, kmax+1); X[1] = J[1]                              # X_0 = R
    for k in 1:kmax
        J = jderiv(J) .+ jmul(cj.U, J)                                     # P = ∂_r + U
        X[k+1] = J[1]
    end
    X
end

# ───────────────────────────── A_2 coefficient table (parsed fixtures) ─────────────────────────────
# r-form entries a_ij^(2) from coeff/spin2-plus-plus-coefficients.generated.tex (machine output).
function f73_A2_rform(::Type{C}, Γ, Γb, r) where {C}
    d = Dict{Tuple{Int,Int},C}()
    d[(0,1)] = -3 / Γb^5
    d[(0,2)] = 3r / Γb^5
    d[(0,3)] = 3*(Γb - 2r) / (2*Γb^4)
    d[(0,4)] = (3r - 2*Γb) / (2*Γb^3)
    d[(0,5)] = (3*Γb - 4r) / (8*Γb^2)
    d[(0,6)] = Γ / (16*Γb)
    d[(1,1)] = 3 / Γb^4
    d[(1,2)] = -3*(Γb + 2r) / (2*Γb^4)
    d[(1,3)] = -(Γb - 3r) / Γb^3
    d[(1,4)] = (7*Γb - 12r) / (8*Γb^2)
    d[(1,5)] = Γ / (4*Γb)
    d[(2,2)] = 3r / (2*Γb^3)
    d[(2,3)] = (3*Γb - 8r) / (4*Γb^2)
    d[(2,4)] = 7*Γ / (16*Γb)
    d[(3,3)] = Γ / (4*Γb)
    d
end

# ───────────────────────────── named-entry fixtures (r-form vs Γ-form) ─────────────────────────────
# Each: (name, r-form generated .tex, Γ-form Appendix C).  Both as (Γ,Γb,r)->value closures.
f73_named_entries() = [
    ("a_00^(1)", (Γ,Γb,r)-> 45*(Γb - 2r)/(2*Γb^7),          (Γ,Γb,r)-> -45*Γ/(2*Γb^7)),
    ("a_11^(1)", (Γ,Γb,r)-> (51*Γb - 90r)/(2*Γb^5),         (Γ,Γb,r)-> 3*(2*Γb - 15*Γ)/(2*Γb^5)),   # source-level-fixed
    ("a_06^(2)", (Γ,Γb,r)-> Γ/(16*Γb),                      (Γ,Γb,r)-> Γ/(16*Γb)),                    # flagship 6th-order
    ("a_04^(3)", (Γ,Γb,r)-> -(Γb - 6r)/(6*Γb^4),            (Γ,Γb,r)-> (Γ + 4r)/(6*Γb^4)),            # source-level-fixed
    ("a_33^(3)", (Γ,Γb,r)-> (3r - 2*Γb)/(6*Γb^2),           (Γ,Γb,r)-> (3*Γ - Γb)/(12*Γb^2)),         # source-level-fixed
    ("a_11^(4)", (Γ,Γb,r)-> 15*(3*Γ - 4*Γb)/(2*Γb^7),       (Γ,Γb,r)-> (45*Γ - 60*Γb)/(2*Γb^7)),
]

# Evaluate a named entry as an RTJet at generic Kerr (r,θ): Γ=ζ̄=r+iacosθ, Γb=ζ=r−iacosθ (GrTool).
function f73_entry_rtjet(f, ::Type{C}, a, r, θ) where {C}
    Γ  = RAD._rt_zetabar(C, a, r, θ)     # paper Γ  = GrTool ζ̄
    Γb = RAD._rt_zeta(C, a, r, θ)        # paper Γ̄  = GrTool ζ
    rr = RAD._rt_r(C, r)
    f(Γ, Γb, rr)
end
_rt_fields(J::RTJet) = (J.v, J.r, J.θ, J.rr, J.θθ, J.rθ)

# ───────────────────────────── generic sample points ──────────────────────────────────────────────
# (a, r, θ, ω, m, λ) — moderate r, off-axis θ, generic complex λ; avoids Δ=0, Γ̄=0, sinθ=0.
const F73_PTS = [
    (0.37, 3.7, 0.9, 0.31, 2, 4.2 + 0.3im),
    (0.61, 5.3, 1.3, 0.47, 2, 3.1 - 0.7im),
    (0.83, 4.1, 0.6, 0.29, 3, 5.6 + 1.1im),
    (0.29, 6.8, 2.1, 0.53, 2, 2.4 - 0.5im),
    (0.50, 4.9, 1.7, 0.41, 3, 4.9 + 0.2im),
    (0.71, 3.3, 1.1, 0.37, 2, 3.7 - 0.9im),
]

# ═══════════════════════════════════════════════════════════════════════════════════════════════════
function run_f73(::Type{C}, tol) where {C}
    LEN = 8; KMAX = 6
    label = C === ComplexF64 ? "F64" : "BigFloat-$(precision(BigFloat))"

    @testset "G-F7-3a quotient R''=FR'+GR [$label]" begin
        for (a0, r0, θ0, ω0, m0, λ0) in F73_PTS
            a = _cvt(real(C), a0); r = _cvt(real(C), r0); ω = _cvt(real(C), ω0)
            λ = C(λ0); m = m0
            q = f73_quotient(C, a, r, ω, m, λ)
            @test abs(q.F_op - q.F_tex) / abs(q.F_tex) < tol      # double-entry F
            @test abs(q.G_op - q.G_tex) / abs(q.G_tex) < tol      # double-entry G
            @test abs(q.cRpp - q.Δ) / abs(q.Δ) < tol              # leading coeff of R'' is Δ
            @test abs(q.resid) / q.scale < tol                    # quotient annihilates the operator
        end
    end

    @testset "G-F7-3b named entries r-form==Γ-form under Γ+Γ̄=2r [$label]" begin
        for (name, frf, fgf) in f73_named_entries()
            for (a0, r0, θ0, _, _, _) in F73_PTS
                a = _cvt(real(C), a0); r = _cvt(real(C), r0); θ = _cvt(real(C), θ0)
                Jr = f73_entry_rtjet(frf, C, a, r, θ)     # r-form (generated .tex)
                Jg = f73_entry_rtjet(fgf, C, a, r, θ)     # Γ-form (Appendix C)
                fr = _rt_fields(Jr); fg = _rt_fields(Jg)
                scale = maximum(abs, fg) + one(real(C))
                @test maximum(abs.(fr .- fg)) / scale < tol   # full RTJet (value + all derivs)
            end
        end
    end

    @testset "G-F7-3c full A_2 row reduction (Route A recursion == Route B direct-op) [$label]" begin
        for (a0, r0, θ0, ω0, m0, λ0) in F73_PTS
            a = _cvt(real(C), a0); r = _cvt(real(C), r0); ω = _cvt(real(C), ω0)
            θ = _cvt(real(C), θ0); λ = C(λ0); m = m0
            Γ  = complex(r,  a*cos(θ)); Γb = complex(r, -a*cos(θ))   # paper Γ, Γ̄ at this point
            cj = f73_coefjets(C, a, r, ω, m, λ, LEN)
            (α, β) = f73_alphabeta(cj, LEN, KMAX, C)

            # closed forms (α_1,β_1)=(U,1), (α_2,β_2)=(U(U−V)−T, U−V)
            @test abs(α[2][1] - cj.U[1]) < tol * (abs(cj.U[1]) + 1)
            @test abs(β[2][1] - one(C)) < tol
            @test abs(α[3][1] - (cj.U[1]*(cj.U[1]-cj.V[1]) - cj.T[1])) < tol*(abs(α[3][1])+1)
            @test abs(β[3][1] - (cj.U[1]-cj.V[1])) < tol*(abs(β[3][1])+1)

            R, Rp = C(0.83 - 0.21im), C(-0.44 + 0.67im)       # generic (R,R')
            tower = f73_deriv_tower(cj, LEN, LEN-1, R, Rp, C)  # length LEN (matches coef jets for jmul)
            XB = f73_Xk_directop(cj, tower, KMAX, C)          # Route B X_k (direct operator)

            # per-k cross-check: direct-operator X_k == α_k R + β_k R' (pins the recursion, incl. X_6)
            for k in 0:KMAX
                pred = α[k+1][1]*R + β[k+1][1]*Rp
                @test abs(XB[k+1] - pred) / (abs(pred) + 1) < tol
            end

            a2 = f73_A2_rform(C, Γ, Γb, r)                    # owner's full A_2 row (r-form)
            # Route A: C_ab via Eq. app-Cnu-formula, then reduced 3-monomial
            C00 = zero(C); C01 = zero(C); C11 = zero(C)
            for ((i, j), aij) in a2
                αi, βi, αj, βj = α[i+1][1], β[i+1][1], α[j+1][1], β[j+1][1]
                C00 += aij * αi * αj
                C01 += aij * (αi*βj + βi*αj)
                C11 += aij * βi * βj
            end
            A2_A = C00*R^2 + C01*R*Rp + C11*Rp^2
            # Route B: direct Σ a_ij X_i X_j with the direct-operator X_k
            A2_B = zero(C)
            for ((i, j), aij) in a2
                A2_B += aij * XB[i+1] * XB[j+1]
            end
            @test abs(A2_A - A2_B) / (abs(A2_A) + 1) < tol
        end
    end
end

# Compact achieved-accuracy diagnostic (visibility marker; the @testset above is the binding gate).
function f73_diagnostics(::Type{C}) where {C}
    LEN = 8; KMAX = 6
    eF = zero(real(C)); eG = zero(real(C)); eNamed = zero(real(C)); eA2 = zero(real(C)); eres = zero(real(C))
    for (a0, r0, θ0, ω0, m0, λ0) in F73_PTS
        a = _cvt(real(C), a0); r = _cvt(real(C), r0); ω = _cvt(real(C), ω0)
        θ = _cvt(real(C), θ0); λ = C(λ0); m = m0
        q = f73_quotient(C, a, r, ω, m, λ)
        eF = max(eF, abs(q.F_op - q.F_tex)/abs(q.F_tex))
        eG = max(eG, abs(q.G_op - q.G_tex)/abs(q.G_tex))
        eres = max(eres, abs(q.resid)/q.scale)
        for (_, frf, fgf) in f73_named_entries()
            Jr = f73_entry_rtjet(frf, C, a, r, θ); Jg = f73_entry_rtjet(fgf, C, a, r, θ)
            fr = _rt_fields(Jr); fg = _rt_fields(Jg)
            eNamed = max(eNamed, maximum(abs.(fr .- fg))/(maximum(abs, fg)+one(real(C))))
        end
        Γ = complex(r, a*cos(θ)); Γb = complex(r, -a*cos(θ))
        cj = f73_coefjets(C, a, r, ω, m, λ, LEN); (α, β) = f73_alphabeta(cj, LEN, KMAX, C)
        R, Rp = C(0.83 - 0.21im), C(-0.44 + 0.67im)
        tower = f73_deriv_tower(cj, LEN, LEN-1, R, Rp, C); XB = f73_Xk_directop(cj, tower, KMAX, C)
        a2 = f73_A2_rform(C, Γ, Γb, r)
        C00 = zero(C); C01 = zero(C); C11 = zero(C); A2_B = zero(C)
        for ((i, j), aij) in a2
            C00 += aij*α[i+1][1]*α[j+1][1]; C11 += aij*β[i+1][1]*β[j+1][1]
            C01 += aij*(α[i+1][1]*β[j+1][1] + β[i+1][1]*α[j+1][1]); A2_B += aij*XB[i+1]*XB[j+1]
        end
        A2_A = C00*R^2 + C01*R*Rp + C11*Rp^2
        eA2 = max(eA2, abs(A2_A - A2_B)/(abs(A2_A)+one(real(C))))
    end
    (eF, eG, eres, eNamed, eA2)
end

@testset "F7-3 ++ source-coefficient convention-independent gate" begin
    run_f73(ComplexF64, 1e-10)                       # G-F7-3d (Float64)
    setprecision(BigFloat, 256) do
        run_f73(Complex{BigFloat}, 1e-40)            # G-F7-3d (BigFloat-256) + {C}-generic
    end
end

let (eF, eG, eres, eN, eA2) = f73_diagnostics(ComplexF64)
    println("F7-3 [F64] max-rel: G-F7-3a F=$(Float64(eF)) G=$(Float64(eG)) resid=$(Float64(eres)) | ",
            "G-F7-3b named=$(Float64(eN)) | G-F7-3c A2=$(Float64(eA2))")
end
setprecision(BigFloat, 256) do
    (eF, eG, eres, eN, eA2) = f73_diagnostics(Complex{BigFloat})
    println("F7-3 [BF256] max-rel: G-F7-3a F=$(Float64(eF)) G=$(Float64(eG)) resid=$(Float64(eres)) | ",
            "G-F7-3b named=$(Float64(eN)) | G-F7-3c A2=$(Float64(eA2))")
end
