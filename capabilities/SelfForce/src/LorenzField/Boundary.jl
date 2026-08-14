# Homogeneous in/out boundary-condition series for the Lorenz-gauge radial ODE (port of
# ref-code/h1Lorenz/src/boundary_conditions.c, legacy H1Lorenz/boundary_conditions.jl). For each
# coupled parity sector, build the asymptotic power/recurrence series of the homogeneous fields at the
# outer boundary (outgoing e^{+iωr*}) and near-horizon (ingoing e^{-iωr*}). Each BC function takes
# a0s (the unit basis selector for the field) and returns
#   (R::Vector{Complex{T}}, R_rs_deriv::Vector{Complex{T}}, rs_boundary::T).
# These series SEED the integrator of the consumed MetricReconstruction.Lorenz.LorenzRadial RHS
# (even_ode!/odd_ode!/monopole_ode!/odd_dipole_ode!) — Boundary never re-owns those RHS. The
# Schwarzschild lapse `_f` and tortoise inverse `r_of_rs` are CONSUMED from Spacetime via Orbit.jl.
# Type-generic {T}; tolerances inherit precision so BigFloat reaches full mantissa.

using LinearAlgebra

# series truncation tolerance: precision-inherited (a hard 1e-14 would cap BigFloat at ~14 digits).
function _bc_tol(::Type{T}) where {T}
    if real(T) <: BigFloat
        bits = precision(real(T))
        return max(big(2.0)^(10 - bits), 8*eps(real(T)))
    end
    return T(1e-14)
end

# Hard cap on the asymptotic-series term count (the outgoing 1/r series only converges beyond a
# mode-dependent radius; in BigFloat there is no overflow to break a misseeded loop).
const _BC_JMAX_CAP = 600

# type-generic FULL null-space basis of a singular square matrix A: Gauss-Jordan to RREF (partial
# pivoting), then one basis vector per free column. Returns a Vector of vectors (handles degenerate
# eigenspaces of dimension > 1).
function _nullbasis(A::Matrix{T}) where {T}
    n = size(A, 1); M = copy(A)
    thresh = sqrt(eps(real(T))) * (1 + maximum(abs, M))
    pivcol = fill(0, n); colpiv = fill(0, n); r = 1
    for c in 1:n
        r > n && break
        pr = r; best = abs(M[r, c])
        for rr in (r+1):n
            v = abs(M[rr, c]); if v > best; best = v; pr = rr; end
        end
        best <= thresh && continue
        if pr != r; for cc in 1:n; M[r,cc], M[pr,cc] = M[pr,cc], M[r,cc]; end; end
        piv = M[r, c]; for cc in 1:n; M[r, cc] /= piv; end
        for rr in 1:n
            rr == r && continue
            f = M[rr, c]; f == 0 && continue
            for cc in 1:n; M[rr, cc] -= f*M[r, cc]; end
        end
        pivcol[r] = c; colpiv[c] = r; r += 1
    end
    freecols = [c for c in 1:n if colpiv[c] == 0]
    vecs = Vector{Vector{T}}()
    for fc in freecols
        x = zeros(T, n); x[fc] = one(T)
        for rr in 1:n
            pc = pivcol[rr]; pc == 0 && continue
            x[pc] = -M[rr, fc]
        end
        nrm = sqrt(sum(abs2, x)); nrm > 0 && (x ./= nrm)
        push!(vecs, x)
    end
    vecs
end

# eigensolve of a small real matrix M, sorted by |Re(λ)| desc; returns (λs, Vt, Vinv), rows of Vt the
# eigenvectors (matching the C gsl transpose), Vinv = inv(Vt). Float64: LAPACK eigen (validated path).
# BigFloat: LAPACK gives the order + the INTEGER eigenvalues (leading-power matrices have integer
# spectra), then eigenvectors are recomputed at full precision via nullspace(M-λI).
function _eig_sorted(M::Matrix{Tr}) where {Tr}
    F = LinearAlgebra.eigen(complex(Float64.(real.(M))))
    perm = sortperm(F.values; by = z -> abs(real(z)), rev = true)
    Vt64 = permutedims(F.vectors[:, perm])
    Tr <: Union{Float64,Complex{Float64}} && return (F.values[perm], Vt64, inv(Vt64))
    T = real(Tr); n = size(M, 1)
    λi = [round(Int, real(F.values[p])) for p in perm]
    Vt = zeros(Complex{T}, n, n); Mc = Complex{T}.(M)
    row = 1
    for λ in unique(λi)
        mult = count(==(λ), λi)
        basis = _nullbasis(Matrix(Mc - Complex{T}(λ)*I))
        for b in 1:mult
            v = b <= length(basis) ? basis[b] : zeros(Complex{T}, n)
            Vt[row, :] = v; row += 1
        end
    end
    (Complex{T}.(λi), Vt, inv(Vt))
end

# common pack helpers: multiply sums by exp(±iωr*) and return (R, R_rs_deriv, rs_boundary)
function _pack_out(sums, dsums, ω, rs_out)
    T = real(eltype(sums)); ex = exp(im*ω*rs_out)
    (Complex{T}[s*ex for s in sums], Complex{T}[d*ex for d in dsums], T(rs_out))
end
function _pack_in(sums, dsums, ω, rs_in)
    T = real(eltype(sums)); ex = exp(-im*ω*rs_in)
    (Complex{T}[s*ex for s in sums], Complex{T}[d*ex for d in dsums], T(rs_in))
end

# ───────── EVEN (l≥2): 5 coupled fields R1,R3,R5,R6,R7 ─────────
function even_out_bcs(a0s::Vector{T}, l, ω, p; rs_out=T(110000)) where {T}
    L = T(l*(l + 1)); lambda = T((l + 2)*(l - 1)); w = ω
    tol = _bc_tol(T)
    a_data = T[ -(2+L)  2      2          2      0;
                 2     -(2+L) -2         -2      0;
                 2L    -2L    -(4+L)     -2L     2;
                 2     -2     -2         -(2+L)  0;
                 0      0      2lambda    0     -(L-2) ]
    M = permutedims(a_data)
    λs, Vt, Vinv = _eig_sorted(M)
    highest_power = round(Int, real(λs[1])); middle_power = round(Int, real(λs[2])); lowest_power = round(Int, real(λs[5]))
    highest_l = round(Int, (-1 + sqrt(1 + 4*abs(highest_power)))/2)
    middle_l  = round(Int, (-1 + sqrt(1 + 4*abs(middle_power)))/2)
    lowest_l  = round(Int, (-1 + sqrt(1 + 4*abs(lowest_power)))/2)
    a13 = a33 = a53 = a63 = a73 = zero(Complex{T})
    if a0s[1] == 1
        sc = (ω*p)^(highest_l - middle_l)
        a13 = real(Vinv[1,1])*sc; a33 = real(Vinv[2,1])*sc; a53 = real(Vinv[3,1])*sc; a63 = real(Vinv[4,1])*sc; a73 = real(Vinv[5,1])*sc
    elseif a0s[2] == 1
        a13 = real(Vinv[1,2]); a33 = real(Vinv[2,2]); a53 = real(Vinv[3,2]); a63 = real(Vinv[4,2]); a73 = real(Vinv[5,2])
    elseif a0s[3] == 1
        a13 = real(Vinv[1,3]); a33 = real(Vinv[2,3]); a53 = real(Vinv[3,3]); a63 = real(Vinv[4,3]); a73 = real(Vinv[5,3])
    elseif a0s[4] == 1
        a13 = real(Vinv[1,4]); a33 = real(Vinv[2,4]); a53 = real(Vinv[3,4]); a63 = real(Vinv[4,4]); a73 = real(Vinv[5,4])
    elseif a0s[5] == 1
        sc = (ω*p)^(lowest_l - middle_l)
        a13 = real(Vinv[1,5])*sc; a33 = real(Vinv[2,5])*sc; a53 = real(Vinv[3,5])*sc; a63 = real(Vinv[4,5])*sc; a73 = real(Vinv[5,5])*sc
    end
    a0s[1] = real(a13); a0s[2] = real(a33); a0s[3] = real(a53); a0s[4] = real(a63); a0s[5] = real(a73)
    _even_series_out(a13, a33, a53, a63, a73, l, ω, w, L, lambda, tol; rs_start=T(rs_out))
end

# shared even-parity outer power series (used by even & even-dipole after the leading amplitudes)
function _even_series_out(a13, a33, a53, a63, a73, l, ω, w, L, lambda, tol; ndip::Bool=false, rs_start=nothing)
    Tc = typeof(a13); T = real(Tc)
    rs_out = rs_start === nothing ? T(110000) : T(rs_start); increase_level = 0
    while true
        r_out = r_of_rs(rs_out)
        j_max = 5
        a1 = Complex{T}[0,0,0,a13,0]; a3 = Complex{T}[0,0,0,a33,0]
        a5 = Complex{T}[0,0,0,a53,0]; a6 = Complex{T}[0,0,0,a63,0]; a7 = Complex{T}[0,0,0,a73,0]
        sum1=a1[4]; sum3=a3[4]; sum5=a5[4]; sum6=a6[4]; sum7=a7[4]
        ds1=im*w*a1[4]; ds3=im*w*a3[4]; ds5=im*w*a5[4]; ds6=im*w*a6[4]; ds7=im*w*a7[4]
        j = 4
        while j < j_max
            C1=(j-4)*(j-3)+4im*w*(j-4)-2-L
            D1=2*(5+L-2*(j-5)^2-3*(j-5)); D3=2*(2*(j-5)-8+4im*w)
            E1=4*((j-6)^2+2*(j-6)-3); E3=8*(5-2*(j-6)); F3=16*(j-7-2)
            C5=(j-4)*(j-3)+4im*w*(j-4)-4-L; D5=2*(12+L-2*(j-5)^2-3*(j-5)); E5=4*((j-6)^2+2*(j-6)-8)
            C7=(j-4)*(j-3)+4im*w*(j-4)+2-L; D7=2*(-3+L-2*(j-5)^2-3*(j-5)); E7=4*((j-6)^2+2*(j-6)+1)
            a14=(C1*a1[4]+(2-4im*w)*a3[4]+2*a5[4]+2*a6[4]+D1*a1[3]+D3*a3[3]-12*a5[3]-20*a6[3]+E1*a1[2]+E3*a3[2]+16*a5[2]+56*a6[2]+F3*a3[1]-48*a6[1])/(2im*w*(j-3))
            a34=(C1*a3[4]+2*(a1[4]-a5[4]-a6[4])+D1*a3[3]+4*(-a1[3]+a5[3]+3*a6[3])+E1*a3[2]-16*a6[2])/(2im*w*(j-3))
            a64=(C1*a6[4]+2*(a1[4]-a5[4]-a3[4])+D1*a6[3]+4*(-a1[3]+a5[3]+3*a3[3])+E1*a6[2]-16*a3[2])/(2im*w*(j-3))
            a54=(C5*a5[4]+2*L*(a1[4]-a3[4]-a6[4])+2*a7[4]+D5*a5[3]+2*L*(-2*a1[3]+4*a3[3]+5*a6[3])-10*a7[3]+E5*a5[2]+4*(-2*L*a3[2]-3*L*a6[2]+3*a7[2]))/(2im*w*(j-3))
            a74=(C7*a7[4]+2*lambda*a5[4]+D7*a7[3]-4*lambda*a5[3]+E7*a7[2])/(2im*w*(j-3))
            rp = r_out^(T(j)-3)
            sum1+=a14/rp; sum3+=a34/rp; sum5+=a54/rp; sum6+=a64/rp; sum7+=a74/rp
            ds1+=(im*w*a14/rp-(1-2/r_out)*(j-3)*a14/r_out^(T(j)-2))
            ds3+=(im*w*a34/rp-(1-2/r_out)*(j-3)*a34/r_out^(T(j)-2))
            ds5+=(im*w*a54/rp-(1-2/r_out)*(j-3)*a54/r_out^(T(j)-2))
            ds6+=(im*w*a64/rp-(1-2/r_out)*(j-3)*a64/r_out^(T(j)-2))
            ds7+=(im*w*a74/rp-(1-2/r_out)*(j-3)*a74/r_out^(T(j)-2))
            d1=abs(a14/rp); d3=abs(a34/rp); d5=abs(a54/rp); d6=abs(a64/rp); d7=abs(a74/rp)
            (d1>=tol||d3>=tol||d5>=tol||d6>=tol||d7>=tol) && (j_max < _BC_JMAX_CAP) && (j_max+=1)
            a1[1:3]=a1[2:4]; a1[4]=a14; a3[1:3]=a3[2:4]; a3[4]=a34
            a6[1:3]=a6[2:4]; a6[4]=a64; a5[1:3]=a5[2:4]; a5[4]=a54; a7[1:3]=a7[2:4]; a7[4]=a74
            j += 1
        end
        ω == 0 && return _pack_out([sum1,sum3,sum5,sum6,sum7],[ds1,ds3,ds5,ds6,ds7],ω,rs_out)
        if all(isfinite, (sum1,sum3,sum5,sum6,sum7,ds1,ds3,ds5,ds6,ds7))
            return _pack_out([sum1,sum3,sum5,sum6,sum7],[ds1,ds3,ds5,ds6,ds7],ω,rs_out)
        else
            rs_out += increase_level==0 ? T(500) : increase_level==1 ? T(1000) : T(5000)
            increase_level = min(increase_level+1, 2)
        end
    end
end

function even_in_bcs(a0s::Vector{T}, l, ω) where {T}
    rs_in = T(-50); w = ω; r_in = r_of_rs(rs_in); epsilon = r_in - 2
    L = T(l*(l + 1)); lambda = T((l + 2)*(l - 1)); tol = _bc_tol(T)
    b1=Complex{T}[0,0,0,a0s[1],0]; b3=Complex{T}[0,0,0,a0s[2],0]
    b5=Complex{T}[0,0,0,a0s[3],0]; b6=Complex{T}[0,0,0,a0s[4],0]; b7=Complex{T}[0,0,0,a0s[5],0]
    sum1=b1[4]; sum3=b3[4]; sum5=b5[4]; sum6=b6[4]; sum7=b7[4]
    ds1=-im*w*b1[4]; ds3=-im*w*b3[4]; ds5=-im*w*b5[4]; ds6=-im*w*b6[4]; ds7=-im*w*b7[4]
    j=4; j_max=5
    while j < j_max
        C1=4*(1-L+3*(j-4)^2-16im*w*(j-4)-(j-4)); C3=8*(2im*w-(j-4))
        D1=2*((j-5)*(3*(j-5)-2)-2*L-1-24im*w*(j-5)); D3=4*(4im*w-(j-5)-1)
        E1=(j-6)*(j-7)-16im*w*(j-6)-2-L; F1=2im*w*(j-7); F3=2im*w*(j-6)
        G3=2*((j-4)*(2*(j-4)-1)-L+1-12im*w*(j-4)); H3=(j-5)*(j-6)-L-2-12im*w*(j-5)
        I5=2*((j-4)*(2*(j-4)-1)-L+4-12im*w*(j-4)); J5=(j-5)*(j-6)-L-4-12im*w*(j-5); K7=(j-5)*(j-6)-L+2-12im*w*(j-5)
        den4=(j-3)*4*(4im*w-(j-3)); den8=(j-3)*8*(4im*w-(j-3))
        b14=(C1*b1[4]+C3*b3[4]-8*b5[4]+D1*b1[3]+D3*b3[3]-8*b6[3]+E1*b1[2]+(2+4im*w)*b3[2]+2*b5[2]+2*b6[2]-F1*b1[1])/den8
        b34=(G3*b3[4]+4*(b1[4]-b5[4]+b6[4])+H3*b3[3]+2*(b1[3]-b5[3]-b6[3])-F3*b3[2])/den4
        b64=(G3*b6[4]+4*(b1[4]-b5[4]+b3[4])+H3*b6[3]+2*(b1[3]-b5[3]-b3[3])-F3*b6[2])/den4
        b54=(I5*b5[4]+2*L*(2*b1[4]+b6[4])-2*b7[4]+J5*b5[3]+2*L*(b1[3]-b3[3]-b6[3])+2*b7[3]-F3*b5[2])/den4
        b74=(G3*b7[4]+4*lambda*b5[4]+K7*b7[3]+2*lambda*b5[3]-F3*b7[2])/den4
        ep=epsilon^(T(j)-3)
        sum1+=b14*ep; sum3+=b34*ep; sum5+=b54*ep; sum6+=b64*ep; sum7+=b74*ep
        ds1+=-im*w*b14*ep+(j-3)*b14*ep/r_in; ds3+=-im*w*b34*ep+(j-3)*b34*ep/r_in
        ds5+=-im*w*b54*ep+(j-3)*b54*ep/r_in; ds6+=-im*w*b64*ep+(j-3)*b64*ep/r_in; ds7+=-im*w*b74*ep+(j-3)*b74*ep/r_in
        d1=abs(b14*ep); d3=abs(b34*ep); d5=abs(b54*ep); d6=abs(b64*ep); d7=abs(b74*ep)
        (d1>=tol||d3>=tol||d5>=tol||d6>=tol||d7>=tol) && (j_max < _BC_JMAX_CAP) && (j_max+=1)
        b1[1:3]=b1[2:4]; b1[4]=b14; b3[1:3]=b3[2:4]; b3[4]=b34
        b6[1:3]=b6[2:4]; b6[4]=b64; b5[1:3]=b5[2:4]; b5[4]=b54; b7[1:3]=b7[2:4]; b7[4]=b74
        j += 1
    end
    _pack_in([sum1,sum3,sum5,sum6,sum7],[ds1,ds3,ds5,ds6,ds7],ω,rs_in)
end

# ───────── EVEN DIPOLE (l=1): 4 coupled fields R1,R3,R5,R6 ─────────
function even_dipole_out_bcs(a0s::Vector{T}, l, ω, p; rs_out=T(110000)) where {T}
    L = T(l*(l + 1)); lambda = T((l + 2)*(l - 1)); w = ω; tol = _bc_tol(T)
    if ω < T(1)/100
        if a0s[1]==1;     a13=Complex{T}(-1); a33=Complex{T}(1);  a53=Complex{T}(2);  a63=Complex{T}(-4)
        elseif a0s[2]==1; a13=Complex{T}(1);  a33=Complex{T}(-1); a53=Complex{T}(3);  a63=Complex{T}(-1)
        elseif a0s[3]==1; a13=Complex{T}(3);  a33=Complex{T}(2);  a53=Complex{T}(-1); a63=Complex{T}(2)
        else;             s=(ω*p)^2; a13=Complex{T}(s*1); a33=Complex{T}(s*-1); a53=Complex{T}(s*-2); a63=Complex{T}(s*-1)
        end
    else
        a13=Complex{T}(a0s[1]); a33=Complex{T}(a0s[2]); a53=Complex{T}(a0s[3]); a63=Complex{T}(a0s[4])
    end
    a73 = zero(Complex{T})
    R, dR, rsb = _even_series_out(a13, a33, a53, a63, a73, l, ω, w, L, lambda, tol; rs_start=T(rs_out))
    (R[1:4], dR[1:4], rsb)
end

function even_dipole_in_bcs(a0s::Vector{T}, l, ω) where {T}
    R, dR, rsb = even_in_bcs([a0s[1],a0s[2],a0s[3],a0s[4],zero(T)], l, ω)
    (R[1:4], dR[1:4], rsb)
end

# ───────── MONOPOLE (l=0): 3 coupled fields R1,R3,R6 ─────────
function monopole_out_bcs(a0s::Vector{T}, l, ω; rs_out=T(110000)) where {T}
    L = T(l*(l + 1)); lambda = T((l + 2)*(l - 1)); w = ω; tol = _bc_tol(T)
    rs_out = T(rs_out)
    while true
        r_out = r_of_rs(rs_out); j_max=5
        a1=Complex{T}[0,0,0,a0s[1],0]; a3=Complex{T}[0,0,0,a0s[2],0]
        a5=zeros(Complex{T},5); a6=Complex{T}[0,0,0,a0s[3],0]; a7=zeros(Complex{T},5)
        sum1=a1[4]; sum3=a3[4]; sum5=a5[4]; sum6=a6[4]; sum7=a7[4]
        ds1=im*w*a1[4]; ds3=im*w*a3[4]; ds5=im*w*a5[4]; ds6=im*w*a6[4]; ds7=im*w*a7[4]
        j=4
        while j < j_max
            C1=(j-4)*(j-3)+4im*w*(j-4)-2-L
            D1=2*(5+L-2*(j-5)^2-3*(j-5)); D3=2*(2*(j-5)-8+4im*w)
            E1=4*((j-6)^2+2*(j-6)-3); E3=8*(5-2*(j-6)); F3=16*(j-7-2)
            C5=(j-4)*(j-3)+4im*w*(j-4)-4-L; D5=2*(12+L-2*(j-5)^2-3*(j-5)); E5=4*((j-6)^2+2*(j-6)-8)
            C7=(j-4)*(j-3)+4im*w*(j-4)+2-L; D7=2*(-3+L-2*(j-5)^2-3*(j-5)); E7=4*((j-6)^2+2*(j-6)+1)
            a14=(C1*a1[4]+(2-4im*w)*a3[4]+2*a5[4]+2*a6[4]+D1*a1[3]+D3*a3[3]-12*a5[3]-20*a6[3]+E1*a1[2]+E3*a3[2]+16*a5[2]+56*a6[2]+F3*a3[1]-48*a6[1])/(2im*w*(j-3))
            a34=(C1*a3[4]+2*(a1[4]-a5[4]-a6[4])+D1*a3[3]+4*(-a1[3]+a5[3]+3*a6[3])+E1*a3[2]-16*a6[2])/(2im*w*(j-3))
            a64=(C1*a6[4]+2*(a1[4]-a5[4]-a3[4])+D1*a6[3]+4*(-a1[3]+a5[3]+3*a3[3])+E1*a6[2]-16*a3[2])/(2im*w*(j-3))
            a54=(C5*a5[4]+2*L*(a1[4]-a3[4]-a6[4])+2*a7[4]+D5*a5[3]+2*L*(-2*a1[3]+4*a3[3]+5*a6[3])-10*a7[3]+E5*a5[2]+4*(-2*L*a3[2]-3*L*a6[2]+3*a7[2]))/(2im*w*(j-3))
            a74=(C7*a7[4]+2*lambda*a5[4]+D7*a7[3]-4*lambda*a5[3]+E7*a7[2])/(2im*w*(j-3))
            rp=r_out^(T(j)-3)
            sum1+=a14/rp; sum3+=a34/rp; sum6+=a64/rp; sum5+=a54/rp; sum7+=a74/rp
            ds1+=(im*w*a14/rp-(1-2/r_out)*(j-3)*a14/r_out^(T(j)-2))
            ds3+=(im*w*a34/rp-(1-2/r_out)*(j-3)*a34/r_out^(T(j)-2))
            ds6+=(im*w*a64/rp-(1-2/r_out)*(j-3)*a64/r_out^(T(j)-2))
            ds5+=(im*w*a54/rp-(1-2/r_out)*(j-3)*a54/r_out^(T(j)-2))
            ds7+=(im*w*a74/rp-(1-2/r_out)*(j-3)*a74/r_out^(T(j)-2))
            d1=abs(a14/rp); d3=abs(a34/rp); d6=abs(a64/rp)
            (d1>=tol||d3>=tol||d6>=tol) && (j_max < _BC_JMAX_CAP) && (j_max+=1)
            a1[1:3]=a1[2:4]; a1[4]=a14; a3[1:3]=a3[2:4]; a3[4]=a34
            a5[1:3]=a5[2:4]; a5[4]=a54; a6[1:3]=a6[2:4]; a6[4]=a64; a7[1:3]=a7[2:4]; a7[4]=a74
            j += 1
        end
        ω == 0 && return _pack_out([sum1,sum3,sum6],[ds1,ds3,ds6],ω,rs_out)
        if all(isfinite,(sum1,sum3,sum5,sum6,sum7,ds1,ds3,ds5,ds6,ds7))
            return _pack_out([sum1,sum3,sum6],[ds1,ds3,ds6],ω,rs_out)
        else
            rs_out += T(5000)
        end
    end
end

function monopole_in_bcs(a0s::Vector{T}, l, ω) where {T}
    rs_in=T(-50); w=ω; r_in=r_of_rs(rs_in); epsilon=r_in-2
    L=T(l*(l+1)); lambda=T((l+2)*(l-1)); tol=_bc_tol(T)
    b1=Complex{T}[0,0,0,a0s[1],0]; b3=Complex{T}[0,0,0,a0s[2],0]
    b5=zeros(Complex{T},5); b6=Complex{T}[0,0,0,a0s[3],0]; b7=zeros(Complex{T},5)
    sum1=b1[4]; sum3=b3[4]; sum5=b5[4]; sum6=b6[4]; sum7=b7[4]
    ds1=-im*w*b1[4]; ds3=-im*w*b3[4]; ds5=-im*w*b5[4]; ds6=-im*w*b6[4]; ds7=-im*w*b7[4]
    j=4; j_max=5
    while j < j_max
        C1=4*(1-L+3*(j-4)^2-16im*w*(j-4)-(j-4)); C3=8*(2im*w-(j-4))
        D1=2*((j-5)*(3*(j-5)-2)-2*L-1-24im*w*(j-5)); D3=4*(4im*w-(j-5)-1)
        E1=(j-6)*(j-7)-16im*w*(j-6)-2-L; F1=2im*w*(j-7); F3=2im*w*(j-6)
        G3=2*((j-4)*(2*(j-4)-1)-L+1-12im*w*(j-4)); H3=(j-5)*(j-6)-L-2-12im*w*(j-5)
        I5=2*((j-4)*(2*(j-4)-1)-L+4-12im*w*(j-4)); J5=(j-5)*(j-6)-L-4-12im*w*(j-5); K7=(j-5)*(j-6)-L+2-12im*w*(j-5)
        den4=(j-3)*4*(4im*w-(j-3)); den8=(j-3)*8*(4im*w-(j-3))
        b14=(C1*b1[4]+C3*b3[4]-8*b5[4]+D1*b1[3]+D3*b3[3]-8*b6[3]+E1*b1[2]+(2+4im*w)*b3[2]+2*b5[2]+2*b6[2]-F1*b1[1])/den8
        b34=(G3*b3[4]+4*(b1[4]-b5[4]+b6[4])+H3*b3[3]+2*(b1[3]-b5[3]-b6[3])-F3*b3[2])/den4
        b64=(G3*b6[4]+4*(b1[4]-b5[4]+b3[4])+H3*b6[3]+2*(b1[3]-b5[3]-b3[3])-F3*b6[2])/den4
        b54=(I5*b5[4]+2*L*(2*b1[4]+b6[4])-2*b7[4]+J5*b5[3]+2*L*(b1[3]-b3[3]-b6[3])+2*b7[3]-F3*b5[2])/den4
        b74=(G3*b7[4]+4*lambda*b5[4]+K7*b7[3]+2*lambda*b5[3]-F3*b7[2])/den4
        ep=epsilon^(T(j)-3)
        sum1+=b14*ep; sum3+=b34*ep; sum6+=b64*ep; sum5+=b54*ep; sum7+=b74*ep
        ds1+=-im*w*b14*ep+(j-3)*b14*ep/r_in; ds3+=-im*w*b34*ep+(j-3)*b34*ep/r_in
        ds6+=-im*w*b64*ep+(j-3)*b64*ep/r_in; ds5+=-im*w*b54*ep+(j-3)*b54*ep/r_in; ds7+=-im*w*b74*ep+(j-3)*b74*ep/r_in
        d1=abs(b14*ep); d3=abs(b34*ep); d6=abs(b64*ep)
        (d1>=tol||d3>=tol||d6>=tol) && (j_max < _BC_JMAX_CAP) && (j_max+=1)
        b1[1:3]=b1[2:4]; b1[4]=b14; b3[1:3]=b3[2:4]; b3[4]=b34
        b5[1:3]=b5[2:4]; b5[4]=b54; b6[1:3]=b6[2:4]; b6[4]=b64; b7[1:3]=b7[2:4]; b7[4]=b74
        j += 1
    end
    _pack_in([sum1,sum3,sum6],[ds1,ds3,ds6],ω,rs_in)
end

# ───────── ODD (l≥2): 2 coupled fields R9,R10 ─────────
function odd_out_bcs(a0s::Vector{T}, l, ω, p; rs_out=T(110000)) where {T}
    L=T(l*(l+1)); lambda=T((l+2)*(l-1)); w=ω; tol=_bc_tol(T)
    rs_out=T(rs_out)
    while true
        r_out=r_of_rs(rs_out); j_max=5
        F=T(2); G=T(-10); H=T(12); N=2lambda; P=-4lambda
        a9=Complex{T}[0,0,a0s[1],0]; a10=Complex{T}[0,0,a0s[2],0]
        if ω != 0
            a_data = T[ -(L+4)  2;  2lambda  -(L-2) ]
            M = permutedims(a_data)
            λs, Vt, Vinv = _eig_sorted(M)
            hp=round(Int,real(λs[1])); lp=round(Int,real(λs[2]))
            hl=round(Int,(-1+sqrt(1+4*abs(hp)))/2); ll=round(Int,(-1+sqrt(1+4*abs(lp)))/2)
            if a0s[1]==1
                sc=(ω*p)^(hl-ll); a9[3]=real(Vinv[1,1])*sc; a10[3]=real(Vinv[2,1])*sc
            elseif a0s[2]==1
                a9[3]=real(Vinv[1,2]); a10[3]=real(Vinv[2,2])
            end
        end
        a0s[1]=real(a9[3]); a0s[2]=real(a10[3])
        sum9=a9[3]; sum10=a10[3]; ds9=im*w*a9[3]; ds10=im*w*a10[3]
        j = ω == 0 ? j_max + 1 : 3
        while j <= j_max
            k=j-2
            C=(4im*w*(k-1)+(k-1)*k-L-4); D=-6*(k-2)-4*(k-2)^2+24+2L; E=4*((k-3)^2+2*(k-3)-8)
            Ii=4im*w*(k-1)+(k-1)*k-L+2; J=-6*(k-2)-4*(k-2)^2-6+2L; K=4*((k-3)^2+2*(k-3)+1)
            a93=(C*a9[3]+D*a9[2]+E*a9[1]+F*a10[3]+G*a10[2]+H*a10[1])/(2im*w*k)
            a103=(Ii*a10[3]+J*a10[2]+K*a10[1]+N*a9[3]+P*a9[2])/(2im*w*k)
            rp=r_out^(T(j)-2)
            sum9+=a93/rp; ds9+=(im*w*a93/rp-(1-2/r_out)*(j-2)*a93/r_out^(T(j)-1))
            sum10+=a103/rp; ds10+=(im*w*a103/rp-(1-2/r_out)*(j-2)*a103/r_out^(T(j)-1))
            d9=abs(a93/r_out^T(k)); d10=abs(a93/r_out^T(k))
            (d9>tol||d10>tol) && (j_max < _BC_JMAX_CAP) && (j_max+=1)
            a9[1:2]=a9[2:3]; a9[3]=a93; a10[1:2]=a10[2:3]; a10[3]=a103
            j += 1
        end
        w == 0 && return _pack_out([sum9,sum10],[ds9,ds10],ω,rs_out)
        if all(isfinite,(sum9,sum10,ds9,ds10))
            return _pack_out([sum9,sum10],[ds9,ds10],ω,rs_out)
        else
            rs_out += T(5000)
        end
    end
end

function odd_in_bcs(a0s::Vector{T}, l, ω) where {T}
    rs_in=T(-50); w=ω; r_in=r_of_rs(rs_in); epsilon=r_in-2
    L=T(l*(l+1)); lambda=T((l+2)*(l-1)); tol=_bc_tol(T)
    b9=Complex{T}[0,0,a0s[1],0]; b10=Complex{T}[0,0,a0s[2],0]
    sum9=b9[3]; sum10=b10[3]; ds9=-im*w*b9[3]; ds10=-im*w*b10[3]
    j=3; j_max=5
    while j < j_max
        k=j-2
        C=2*((k-1)+12im*w*(k-1)-2*(k-1)^2+L-4); D=(4+12im*w*(k-2)+L-(k-2)*(k-3)); E=2im*w*(k-3)
        F=T(2); G=T(-2); H=2*(k-1+12im*w*(k-1)-2*(k-1)^2+L-1); J=(-2+12im*w*(k-2)+L-(k-2)*(k-3))
        N=-4lambda; P=-2lambda
        b93=(C*b9[3]+D*b9[2]+E*b9[1]+F*b10[3]+G*b10[2])/(k*4*(k-4im*w))
        b103=(H*b10[3]+J*b10[2]+E*b10[1]+N*b9[3]+P*b9[2])/(4*k*(k-4im*w))
        ep=epsilon^(T(j)-2)
        sum9+=b93*ep; sum10+=b103*ep
        ds9+=-im*w*b93*ep+(j-2)*b93*ep/r_in; ds10+=-im*w*b103*ep+(j-2)*b103*ep/r_in
        d9=abs(b93*ep); d10=abs(b103*ep)
        (d9>tol||d10>tol) && (j_max < _BC_JMAX_CAP) && (j_max+=1)
        b9[1:2]=b9[2:3]; b9[3]=b93; b10[1:2]=b10[2:3]; b10[3]=b103
        j += 1
    end
    _pack_in([sum9,sum10],[ds9,ds10],ω,rs_in)
end

# ───────── ODD DIPOLE (l=1): 1 coupled field R9 ─────────
function odd_dipole_out_bcs(a0s::Vector{T}, l, ω, p; rs_out=T(110000)) where {T}
    L=T(l*(l+1)); lambda=T((l+2)*(l-1)); w=ω; tol=_bc_tol(T)
    rs_out=T(rs_out)
    while true
        r_out=r_of_rs(rs_out); j_max=5
        F=T(2); G=T(-10); H=T(12); N=2lambda; P=-4lambda
        a9=Complex{T}[0,0,a0s[1],0]; a10=zeros(Complex{T},4)
        sum9=a9[3]; sum10=a10[3]; ds9=im*w*a9[3]; ds10=im*w*a10[3]
        j = ω == 0 ? j_max + 1 : 3
        while j <= j_max
            k=j-2
            C=(4im*w*(k-1)+(k-1)*k-L-4); D=-6*(k-2)-4*(k-2)^2+24+2L; E=4*((k-3)^2+2*(k-3)-8)
            Ii=4im*w*(k-1)+(k-1)*k-L+2; J=-6*(k-2)-4*(k-2)^2-6+2L; K=4*((k-3)^2+2*(k-3)+1)
            a93=(C*a9[3]+D*a9[2]+E*a9[1]+F*a10[3]+G*a10[2]+H*a10[1])/(2im*w*k)
            a103=(Ii*a10[3]+J*a10[2]+K*a10[1]+N*a9[3]+P*a9[2])/(2im*w*k)
            rp=r_out^(T(j)-2)
            sum9+=a93/rp; ds9+=(im*w*a93/rp-(1-2/r_out)*(j-2)*a93/r_out^(T(j)-1))
            sum10+=a103/rp; ds10+=(im*w*a103/rp-(1-2/r_out)*(j-2)*a103/r_out^(T(j)-1))
            d9=abs(a93/r_out^T(k)); d10=abs(a93/r_out^T(k))
            (d9>tol||d10>tol) && (j_max < _BC_JMAX_CAP) && (j_max+=1)
            a9[1:2]=a9[2:3]; a9[3]=a93; a10[1:2]=a10[2:3]; a10[3]=a103
            j += 1
        end
        w == 0 && return _pack_out([sum9],[ds9],ω,rs_out)
        if all(isfinite,(sum9,ds9))
            return _pack_out([sum9],[ds9],ω,rs_out)
        else
            rs_out += T(5000)
        end
    end
end

function odd_dipole_in_bcs(a0s::Vector{T}, l, ω) where {T}
    rs_in=T(-50); w=ω; r_in=r_of_rs(rs_in); epsilon=r_in-2
    L=T(l*(l+1)); lambda=T((l+2)*(l-1)); tol=_bc_tol(T)
    b9=Complex{T}[0,0,a0s[1],0]; b10=zeros(Complex{T},4)
    sum9=b9[3]; sum10=b10[3]; ds9=-im*w*b9[3]; ds10=-im*w*b10[3]
    j=3; j_max=5
    while j < j_max
        k=j-2
        C=2*((k-1)+12im*w*(k-1)-2*(k-1)^2+L-4); D=(4+12im*w*(k-2)+L-(k-2)*(k-3)); E=2im*w*(k-3)
        F=T(2); G=T(-2); H=2*(k-1+12im*w*(k-1)-2*(k-1)^2+L-1); J=(-2+12im*w*(k-2)+L-(k-2)*(k-3))
        N=-4lambda; P=-2lambda
        b93=(C*b9[3]+D*b9[2]+E*b9[1]+F*b10[3]+G*b10[2])/(k*4*(k-4im*w))
        b103=(H*b10[3]+J*b10[2]+E*b10[1]+N*b9[3]+P*b9[2])/(4*k*(k-4im*w))
        ep=epsilon^(T(j)-2)
        sum9+=b93*ep; sum10+=b103*ep
        ds9+=-im*w*b93*ep+(j-2)*b93*ep/r_in; ds10+=-im*w*b103*ep+(j-2)*b103*ep/r_in
        d9=abs(b93*ep); d10=abs(b103*ep)
        (d9>tol||d10>tol) && (j_max < _BC_JMAX_CAP) && (j_max+=1)
        b9[1:2]=b9[2:3]; b9[3]=b93; b10[1:2]=b10[2:3]; b10[3]=b103
        j += 1
    end
    _pack_in([sum9],[ds9],ω,rs_in)
end

# ───────── EVEN STATIC (l+m even, m=0, n=0 ⇒ ω=0): 3 coupled fields R1,R3,R5 (log-augmented) ─────────
function even_static_out_bcs(a0s::Vector{T}, l, ω; rs_out=T(110000)) where {T}
    rs_out=T(rs_out); r_out=r_of_rs(rs_out)
    sum1=zero(T); sum3=zero(T); sum5=zero(T); ds1=zero(T); ds3=zero(T); ds5=zero(T)
    M=one(T)
    a1=zeros(T,4); a3=zeros(T,4); a5=zeros(T,4); a1l=zeros(T,4); a3l=zeros(T,4); a5l=zeros(T,4)
    c1=a0s[1]; c2=a0s[2]; c3=a0s[3]; lf=T(l)
    a3[1]=c1; a5[1]=c2; a1[1]=a3[1]+a5[1]/(lf+1)
    a1[2]=((-1+lf)*(c2+c1*lf))/lf
    a3[2]=c1+c2/lf+c1*lf
    a5[2]=(2*c1*lf*(1+lf)+c2*(-2-lf^2+lf^3))/(lf*(1+lf))
    a1[3]=-c3/(2lf)+(2*c1*lf*(-6-2lf+12lf^2+9lf^3+4lf^4+5lf^5+2lf^6)+c2*(12+8lf+6lf^2-15lf^3-8lf^4+15lf^5+6lf^6))/(2*lf^2*(1+lf)*(1+2lf)*(3+2lf))
    a3[3]=c3/(2lf)+(c2*(-12-8lf+6lf^2+51lf^3+42lf^4+7lf^5-2lf^6)+2*c1*lf*(6+6lf+6lf^2+23lf^3+28lf^4+13lf^5+2lf^6))/(2*lf^2*(1+lf)*(1+2lf)*(3+2lf))
    a5[3]=c3
    a1l[3]=(-2*(2+lf)*(-(c2*(-3+lf))+2*c1*(1+lf)))/((1+lf)*(1+2lf)*(3+2lf))
    a3l[3]=-a1l[3]; a5l[3]=2lf*a3l[3]
    for k in 0:2
        kk=k+1
        sum1+=(a1[kk]+a1l[kk]*log(r_out))*r_out^(-(lf+k))
        sum3+=(a3[kk]+a3l[kk]*log(r_out))*r_out^(-(lf+k))
        sum5+=(a5[kk]+a5l[kk]*log(r_out))*r_out^(-(lf+k))
        ds1+=_f(r_out)*(-(lf+k)*a1[kk]+(1-(lf+k)*log(r_out))*a1l[kk])/r_out^(lf+k+1)
        ds3+=_f(r_out)*(-(lf+k)*a3[kk]+(1-(lf+k)*log(r_out))*a3l[kk])/r_out^(lf+k+1)
        ds5+=_f(r_out)*(-(lf+k)*a5[kk]+(1-(lf+k)*log(r_out))*a5l[kk])/r_out^(lf+k+1)
    end
    for k in (l+3):(l+10)
        C1=(lf*(lf+1)+1-T(k)^2); C5=(lf*(lf+1)+k*(1-k)); D5=2k-one(T)
        j=k-1; D1=T(j)*(j-1); D3=2*(j+one(T)); E1=1-2*T(j); G3=(lf*(lf+1)+2-2*T(j)^2); G1=-2*T(j); H3=4*T(j); E5=(T(j)*(1-j)+2); D52=2*T(j)-1
        n=k-2; F3=n+one(T); I3=T(n)^2-1; J3=-2*T(n); L=lf*(lf+1)
        Mat=T[ C1 2k -(k+1) 1 -1 0;
               0  C1 0 -(k+1) 0 -1;
               -(k+1) 1 C1 2k 1 0;
               0 -(k+1) 0 C1 0 1;
               -2L 0 2L 0 C5 D5;
               0 -2L 0 2L 0 C5 ]
        A=-2M*(D1*a1[3]+D3*a3[3]+a5[3]+E1*a1l[3]-2*a3l[3])+4M^2*(F3*a3[2]-a3l[2])
        B=-2M*(D1*a1l[3]+D3*a3l[3]+a5l[3])+4M^2*F3*a3l[2]
        C=2M*(G3*a3[3]+G1*a1[3]+H3*a3l[3]+2*a1l[3])+4M^2*(I3*a3[2]+J3*a3l[2])
        D=2M*(G3*a3l[3]+G1*a1l[3])+4M^2*I3*a3l[2]
        E=2M*(E5*a5[3]+2L*a3[3]+D52*a5l[3])
        Ff=2M*(E5*a5l[3]+2L*a3l[3])
        Vec=T[A,B,C,D,E,Ff]
        row=Mat\Vec
        a1l[4]=row[2]; a3l[4]=row[4]; a5l[4]=row[6]; a1[4]=row[1]; a3[4]=row[3]; a5[4]=row[5]
        sum1+=(a1[4]+a1l[4]*log(r_out))/r_out^k; sum3+=(a3[4]+a3l[4]*log(r_out))/r_out^k; sum5+=(a5[4]+a5l[4]*log(r_out))/r_out^k
        ds1+=_f(r_out)*(-k*a1[4]+(1-k*log(r_out))*a1l[4])/r_out^(k+1)
        ds3+=_f(r_out)*(-k*a3[4]+(1-k*log(r_out))*a3l[4])/r_out^(k+1)
        ds5+=_f(r_out)*(-k*a5[4]+(1-k*log(r_out))*a5l[4])/r_out^(k+1)
        a1[1:3]=a1[2:4]; a3[1:3]=a3[2:4]; a5[1:3]=a5[2:4]; a1l[1:3]=a1l[2:4]; a3l[1:3]=a3l[2:4]; a5l[1:3]=a5l[2:4]
    end
    if l >= 18
        sf=r_out^l; sum1*=sf; sum3*=sf; sum5*=sf; ds1*=sf; ds3*=sf; ds5*=sf
    end
    (Complex{T}[sum1,sum3,sum5], Complex{T}[ds1,ds3,ds5], rs_out)
end

function even_static_in_bcs(a0s::Vector{T}, l, ω) where {T}
    rs_in=T(-30); r_in=r_of_rs(rs_in); eps_=r_in-2
    L=T(l*(l+1)); lf=T(l)
    c1=a0s[1]; c2=a0s[2]; c3=a0s[3]
    b1=zeros(T,4); b3=zeros(T,4); b5=zeros(T,4)
    b3[1]=c1; b3[2]=c2; b5[2]=c3
    b5[3]=T(1)/4*(-4*b5[2]+lf*(b3[1]+b5[2])+lf^2*(b3[1]+b5[2]))
    b1[3]=(lf*(lf+1)*b3[1]+2*b5[2])/8
    sum3=c1+c2*eps_; ds3=b3[2]*eps_/r_in
    sum1=zero(T); ds1=zero(T); sum5=zero(T); ds5=zero(T)
    for k in 0:2
        kk=k+1
        sum1+=b1[kk]*eps_^k; sum5+=b5[kk]*eps_^k
        ds1+=k*b1[kk]*eps_^k/r_in; ds5+=k*b5[kk]*eps_^k/r_in
    end
    k=3; k_max=6; tol=_bc_tol(T)
    while k <= k_max
        C5=2L-4*((k-1)^2+1); D5=L-(k-2)*(k-1); C3=4*(k-1)*(2-k); C1=-4k
        D3=2*(L+(k-2)*(1-2*(k-2))); D1=2*(k-4); E3=L+1-(k-3)^2; E1=T(k-4)
        F1=4*(L+1+(k-1)*(4-3*(k-1))); G1=2*(2L+2+(k-2)*(2-3*(k-2))); G3=2*(k-2)
        b53=(C5*b5[3]-4L*b1[3]+D5*b5[2]+2L*(b3[2]-b1[2]))/(4*k*(k-1))
        b13=(F1*b1[3]+G1*b1[2]+G3*b3[2]-2*b5[2]+E3*b1[1]+E1*b3[1]-b5[1])/(8*k*(k-2))
        b32=(8*k*b13-C1*b1[3]-4*b5[3]-D3*b3[2]-D1*b1[2]-4*b5[2]-E3*b3[1]-E1*b1[1]-b5[1])/C3
        d1=abs(b13*eps_^k); d3=abs(b32*eps_^(k-1)); d5=abs(b53*eps_^k)
        (d1>=tol||d3>=tol||d5>=tol) && (k_max < _BC_JMAX_CAP) && (k_max+=1)
        sum1+=b13*eps_^k; sum3+=b32*eps_^(k-1); sum5+=b53*eps_^k
        ds1+=k*b13*eps_^k/r_in; ds3+=(k-1)*b32*eps_^(k-1)/r_in; ds5+=k*b53*eps_^k/r_in
        b1[1:3]=b1[2:4]; b1[3]=b13; b3[1]=b3[2]; b3[2]=b32; b5[1:3]=b5[2:4]; b5[3]=b53
        k += 1
    end
    (Complex{T}[sum1,sum3,sum5], Complex{T}[ds1,ds3,ds5], rs_in)
end
