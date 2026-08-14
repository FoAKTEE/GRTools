# Variation-of-parameters / two-point-shoot integrator that turns the green inputs (CONSUMED
# LorenzRadial RHS + Orbit R1..R10 sources + Boundary in/out BC series) into the radial
# homogeneous+particular solution, assembles the algebraic gauge field, builds NModeData, and WRITES
# the h1 corpus. CLOSES LorenzField. Port of H1Lorenz integrator.jl/scaling.jl/gauge_fields.jl/
# static.jl(R8)/assemble.jl/driver.jl/io_hdf5.jl, REDESIGNED to consume NumericsCore's canonical DP5/
# GBS integrators (`integrate_to_targets!`/`gbs_integrate_to_targets!`) + `blas_single_thread` guard +
# `Data.H1Schema`/`write_mode_h5` (the storage contract), instead of the legacy vendored copies.
# Type-generic {T}: BigFloat far-field series-fill preserved.

using LinearAlgebra: BLAS
using NumericsCore: integrate_to_targets!, gbs_integrate_to_targets!, gbs_kmax, blas_single_thread,
                    legendre_P_thetaderiv
import Data
using Data: H1Schema

# ── integrator tolerance + far-field seeding (precision-inherited; H1Lorenz integrator.jl/boundary) ──

"RK step/solve tolerance, precision-inherited (the C path used a hard 1e-12)."
function _rk_tol(::Type{T}) where {T}
    T <: BigFloat && return big(2.0)^(20 - precision(BigFloat))
    return T(1e-12)
end

# Optional high-order GBS path for BigFloat (off by default — DP5 is the validated Float64 oracle path
# and the measured-faster BigFloat default; see H1Lorenz integrator.jl notes).
const USE_GBS = Ref{Bool}(false)
const GBS_KMAX = Ref{Any}(nothing)
_gbs_kmax(::Type{T}) where {T} = GBS_KMAX[] !== nothing ? Int(GBS_KMAX[]) : gbs_kmax(T)

function _integrate_targets(f!, y0::Vector{T}, x0::T, targets, dir::Int; tol::T) where {T}
    if T <: BigFloat && USE_GBS[]
        return gbs_integrate_to_targets!(f!, y0, x0, targets, dir; tol=tol, kmax=_gbs_kmax(T))
    end
    return integrate_to_targets!(f!, y0, x0, targets, dir; tol=tol)
end

# outer-boundary r* at which to seed the homogeneous outgoing series (Float64: exact C value 110000 ⇒
# validated oracle path byte-identical; BigFloat: just beyond the outermost grid r* so the empty
# 1e5-M leg collapses without changing the convergent series).
const BC_RS_OUT = Ref{Any}(nothing)
function _far_field_rs_out(::Type{T}, grid_rs_max::Real) where {T}
    BC_RS_OUT[] !== nothing && return T(BC_RS_OUT[])
    T <: BigFloat && return T(grid_rs_max) + T(20)
    return T(110000)
end

# match radius for the BigFloat analytic far-field FILL (∝1/ω, clamped); Float64/ω=0 ⇒ Inf so the
# validated single-sweep integration path is byte-identical.
const BC_MATCH_R = Ref{Any}(nothing)
const BC_MATCH_COEF = Ref{Any}(nothing)
function _far_field_match_r(::Type{T}, ω) where {T}
    BC_MATCH_R[] !== nothing && return T(BC_MATCH_R[])
    (T <: BigFloat && ω != 0) || return T(Inf)
    coef = BC_MATCH_COEF[] === nothing ? T(25) : T(BC_MATCH_COEF[])
    clamp(coef/abs(T(ω)), T(60), T(400))
end

# ── homogeneous-field integration (two-point shoot, EHS) ────────────────────────────────────────

"""
    integrate_field_equations(orbit, l, ω, cset) -> NamedTuple

Integrate the homogeneous fields inward from the outer boundary and outward from the horizon, one
basis vector `d` per coupled field, landing exactly on each grid r*. Returns
`hom_out`/`hom_in` :: `Array{T,4}` `[k,d,part(1..4),grid]`, the matching `R_in`/`R_out`(+rs-derivs)
at r0, and the asymptotic amplitudes `C_in_hom`/`C_out_hom`. `part`: 1=Re(R),2=Im(R),3=Re(dR/dr),
4=Im(dR/dr). Consumes `cset.ode!` (the LorenzRadial RHS) + `cset.out_bc!`/`in_bc!` (Boundary series).
"""
function integrate_field_equations(orbit::OrbitalParams{T}, l, ω, cset) where {T}
    nf = cset.num_coupled
    grid = orbit.grid; gs = length(grid)
    ra = orbit.ra_index; r0i = orbit.r0_index; rb = orbit.rb_index
    tol = _rk_tol(T); p = orbit.p
    grid_rs_max = rs_of_r(grid[gs])
    rs_out_seed = _far_field_rs_out(T, grid_rs_max)

    hom_out = zeros(T, nf, nf, 4, gs - ra)
    hom_in  = zeros(T, nf, nf, 4, rb + 1)
    R_out = zeros(Complex{T}, nf, nf); R_out_d = zeros(Complex{T}, nf, nf)
    R_in  = zeros(Complex{T}, nf, nf); R_in_d  = zeros(Complex{T}, nf, nf)
    C_out_hom = zeros(T, nf, nf); C_in_hom = zeros(T, nf, nf)

    r_match = ω == 0 ? T(Inf) : _far_field_match_r(T, ω)
    all_idxs = collect((gs-1):-1:ra)
    far_idxs  = [i for i in all_idxs if grid[i+1] > r_match]
    near_idxs = [i for i in all_idxs if grid[i+1] <= r_match]

    # ----- OUTER fields (integrate inward, decreasing rs) -----
    for d in 1:nf
        a0s = zeros(T, nf); a0s[d] = 1
        for i in far_idxs
            af = zeros(T, nf); af[d] = 1
            Rf, dRf, _ = cset.out_bc!(af, l, ω, p; rs_out=rs_of_r(grid[i+1]))
            ef = _f(grid[i+1]); store = gs - 1 - i
            for k in 1:nf
                hom_out[k, d, 1, store+1] = real(Rf[k]); hom_out[k, d, 2, store+1] = imag(Rf[k])
                hom_out[k, d, 3, store+1] = real(dRf[k])/ef; hom_out[k, d, 4, store+1] = imag(dRf[k])/ef
            end
            i == r0i && (for k in 1:nf; R_out[k, d] = Rf[k]; R_out_d[k, d] = dRf[k]; end)
        end
        seed_rs = isempty(far_idxs) ? rs_out_seed : rs_of_r(r_match)
        R, dR, rs_b = cset.out_bc!(a0s, l, ω, p; rs_out=seed_rs)
        for k in 1:nf; C_out_hom[k, d] = a0s[k]; end
        isempty(near_idxs) && continue
        # power scaling for static (ω=0) modes: scale up the tiny BCs (Runge-Kutta.c L57-72)
        psc = zero(T)
        if ω == 0
            for k in 1:nf
                v = -log10(abs(real(R[k]))); v > psc && (psc = v)
            end
        end
        scl = T(10)^psc
        y = zeros(T, 4nf)
        for k in 1:nf
            y[1+4(k-1)] = real(R[k])*scl;  y[2+4(k-1)] = real(dR[k])*scl
            y[3+4(k-1)] = imag(R[k])*scl;  y[4+4(k-1)] = imag(dR[k])*scl
        end
        f!(g, rs, yy) = cset.ode!(g, rs, yy, l, ω)
        targets = [rs_of_r(grid[i+1]) for i in near_idxs]
        Y = _integrate_targets(f!, y, T(rs_b), targets, -1; tol=tol)
        for (row, i) in enumerate(near_idxs)
            ef = _f(grid[i+1]); store = gs - 1 - i
            for k in 1:nf
                hom_out[k, d, 1, store+1] = Y[row, 1+4(k-1)]
                hom_out[k, d, 2, store+1] = Y[row, 3+4(k-1)]
                hom_out[k, d, 3, store+1] = Y[row, 2+4(k-1)]/ef
                hom_out[k, d, 4, store+1] = Y[row, 4+4(k-1)]/ef
            end
            if i == r0i
                for k in 1:nf
                    R_out[k, d]   = Y[row, 1+4(k-1)] + im*Y[row, 3+4(k-1)]
                    R_out_d[k, d] = Y[row, 2+4(k-1)] + im*Y[row, 4+4(k-1)]
                end
            end
        end
    end

    # ----- INNER fields (integrate outward from rs_in, increasing rs) -----
    for d in 1:nf
        a0s = zeros(T, nf); a0s[d] = 1
        R, dR, rs_b = cset.in_bc!(a0s, l, ω)
        for k in 1:nf; C_in_hom[k, d] = a0s[k]; end
        y = zeros(T, 4nf)
        for k in 1:nf
            y[1+4(k-1)] = real(R[k]);  y[2+4(k-1)] = real(dR[k])
            y[3+4(k-1)] = imag(R[k]);  y[4+4(k-1)] = imag(dR[k])
        end
        f!(g, rs, yy) = cset.ode!(g, rs, yy, l, ω)
        idxs = 0:rb
        targets = [rs_of_r(grid[i+1]) for i in idxs]
        Y = _integrate_targets(f!, y, T(rs_b), targets, +1; tol=tol)
        for (row, i) in enumerate(idxs)
            ef = _f(grid[i+1])
            for k in 1:nf
                hom_in[k, d, 1, i+1] = Y[row, 1+4(k-1)]
                hom_in[k, d, 2, i+1] = Y[row, 3+4(k-1)]
                hom_in[k, d, 3, i+1] = Y[row, 2+4(k-1)]/ef
                hom_in[k, d, 4, i+1] = Y[row, 4+4(k-1)]/ef
            end
            if i == r0i
                for k in 1:nf
                    R_in[k, d]   = Y[row, 1+4(k-1)] + im*Y[row, 3+4(k-1)]
                    R_in_d[k, d] = Y[row, 2+4(k-1)] + im*Y[row, 4+4(k-1)]
                end
            end
        end
    end

    (hom_out=hom_out, hom_in=hom_in, R_out=R_out, R_out_rs_deriv=R_out_d,
     R_in=R_in, R_in_rs_deriv=R_in_d, C_out_hom=C_out_hom, C_in_hom=C_in_hom)
end

# ── scaling coefficients (EHS Wronskian solve, χ=π/2 collapse) ──────────────────────────────────

"""
    calculate_scaling_coefficients(orbit, l, m, integ, cset) -> (C_in, C_out)

Method of Extended Homogeneous Solutions: build the `2nf×2nf` Wronskian matrix from the in/out
fields + r*-derivatives at r0, RHS = source jump vector (χ=π/2 collapse), solve `Φ\\J` (BLAS-guarded
to a single thread for the deadlock-safe, bit-deterministic tiny-matrix LU). `C = (1/(2 f_p))·sol`.
"""
function calculate_scaling_coefficients(orbit::OrbitalParams{T}, l, m, integ, cset) where {T}
    nf = cset.num_coupled
    R_in  = integ.R_in;  R_in_d  = integ.R_in_rs_deriv
    R_out = integ.R_out; R_out_d = integ.R_out_rs_deriv
    Φ = zeros(Complex{T}, 2nf, 2nf)
    for i in 1:nf, j in 1:nf
        Φ[i, j]       = -R_in[i, j];   Φ[i, j+nf]    =  R_out[i, j]
        Φ[i+nf, j]    = -R_in_d[i, j]; Φ[i+nf, j+nf] =  R_out_d[i, j]
    end
    J = zeros(Complex{T}, 2nf)
    for i in 1:nf; J[i+nf] = cset.sources[i](l, m, orbit); end
    sol = blas_single_thread(() -> Φ \ J)
    pref = 1/(2*_f(orbit.p))
    (Complex{T}[pref*sol[i] for i in 1:nf], Complex{T}[pref*sol[nf+i] for i in 1:nf])
end

# ── algebraic gauge-field reconstruction on the grid (H1Lorenz gauge_fields.jl) ─────────────────

# EVEN / EVEN-DIPOLE / MONOPOLE: construct R2 (and R4 for l≥1) from the integrated coupled fields.
function construct_R2_R4_gauge_fields_on_grid!(inhom, orbit::OrbitalParams{T}, l, ω, num_coupled) where {T}
    grid = orbit.grid; gs = length(grid); r0i = orbit.r0_index
    ω == 0 && return inhom
    for i in 0:gs
        r_p = i == gs ? grid[r0i+1] : grid[i+1]
        f_p = _f(r_p); df = _df(r_p); idx = i+1
        cplx(k) = inhom[k,1,idx] + im*inhom[k,2,idx]
        cplxd(k) = inhom[k,3,idx] + im*inhom[k,4,idx]
        R1 = cplx(1); R1d = cplxd(1); R3 = cplx(2); R3d = cplxd(2)
        R5 = zero(Complex{T}); R5d = zero(Complex{T}); R6 = zero(Complex{T}); R6d = zero(Complex{T})
        R7 = zero(Complex{T}); R7d = zero(Complex{T})
        if l == 0
            R6 = cplx(3); R6d = cplxd(3)
        else
            R5 = cplx(3); R5d = cplxd(3); R6 = cplx(4); R6d = cplxd(4)
            l >= 2 && (R7 = cplx(5); R7d = cplxd(5))
        end
        V_l = (one(T)/4)*f_p*(df/r_p + l*(l+1)/r_p^2) - (one(T)/4)*ω^2
        R1dd = 1/f_p^2*(4*V_l*R1 + 4/r_p^2*f_p^2*R3d + 2*f_p/r_p^2*(1-4/r_p)*(R1-R5-f_p*R3) - 2*f_p^2/r_p^2*(1-6/r_p)*R6 - f_p*df*R1d)
        R3dd = 1/f_p^2*(4*V_l*R3 - 2*f_p/r_p^2*(R1-R5-(1-4/r_p)*(R3+R6)) - f_p*df*R3d)
        R5dd = 1/f_p^2*(4*V_l*R5 + 4*f_p/r_p^2*((1-4.5/r_p)*R5 - (one(T)/2)*l*(l+1)*(R1-f_p*R3) + (one(T)/2)*(1-3/r_p)*(l*(l+1)*R6-R7)) - f_p*df*R5d)
        R2  = im*f_p/ω*(R1d - f_p*R3d + 1/r_p*(R1-R5-f_p*R3-2*f_p*R6))
        R2d = im*df/ω*(R1d - f_p*R3d + 1/r_p*(R1-R5-f_p*R3-2*f_p*R6)) + im*f_p/ω*(R1dd - f_p*R3dd - df*R3d + R1d/r_p - R1/r_p^2 - R5d/r_p + R5/r_p^2 - f_p/r_p*R3d - df/r_p*R3 + f_p/r_p^2*R3 - 2*f_p/r_p*R6d - 2*df*R6/r_p + 2*f_p*R6/r_p^2)
        R4  = im*f_p/ω*(R5d + 2*R5/r_p + l*(l+1)*R6/r_p - R7/r_p)
        R4d = im*df/ω*(R5d + 2*R5/r_p + l*(l+1)*R6/r_p - R7/r_p) + im*f_p/ω*(R5dd + 2*R5d/r_p - 2*R5/r_p^2 + l*(l+1)*(R6d/r_p - R6/r_p^2) - R7d/r_p + R7/r_p^2)
        inhom[num_coupled+1,1,idx]=real(R2); inhom[num_coupled+1,2,idx]=imag(R2)
        inhom[num_coupled+1,3,idx]=real(R2d); inhom[num_coupled+1,4,idx]=imag(R2d)
        if l >= 1
            inhom[num_coupled+2,1,idx]=real(R4); inhom[num_coupled+2,2,idx]=imag(R4)
            inhom[num_coupled+2,3,idx]=real(R4d); inhom[num_coupled+2,4,idx]=imag(R4d)
        end
    end
    inhom
end

# ODD / ODD-DIPOLE: construct R8 from R9,R10 via the gauge condition (ω≠0); static m=0 ⇒ series match.
function construct_R8_gauge_field_on_grid!(inhom, orbit::OrbitalParams{T}, l, m, ω, num_coupled) where {T}
    grid = orbit.grid; gs = length(grid); r0i = orbit.r0_index
    (m == 0 && ω == 0) && return _construct_R8_static_on_grid!(inhom, orbit, l, num_coupled)
    for i in 0:gs
        r_p = i == gs ? grid[r0i+1] : grid[i+1]
        f_p = _f(r_p); df = _df(r_p); idx = i+1
        R9  = inhom[1,1,idx] + im*inhom[1,2,idx]; R9d = inhom[1,3,idx] + im*inhom[1,4,idx]
        R10 = zero(Complex{T}); R10d = zero(Complex{T})
        if num_coupled != 1
            R10  = inhom[2,1,idx] + im*inhom[2,2,idx]; R10d = inhom[2,3,idx] + im*inhom[2,4,idx]
        end
        V_l = (one(T)/4)*f_p*(df/r_p + l*(l+1)/r_p^2)
        R9dd = 1/f_p^2*((4*V_l + 4*f_p/r_p^2*(1-4.5/r_p) - ω^2)*R9 - 2*f_p/r_p^2*(1-3/r_p)*R10 - f_p*df*R9d)
        R8  = im*f_p/ω*(R9d + 2*R9/r_p - R10/r_p)
        R8d = im/ω*(f_p*R9dd + df*R9d + 2*f_p*(R9d/r_p - R9/r_p^2) + 2*df*R9/r_p - f_p*(R10d/r_p - R10/r_p^2) - df*R10/r_p)
        inhom[num_coupled+1,1,idx]=real(R8); inhom[num_coupled+1,2,idx]=imag(R8)
        inhom[num_coupled+1,3,idx]=real(R8d); inhom[num_coupled+1,4,idx]=imag(R8d)
    end
    inhom
end

# static R8 (m=0,n=0): series solutions + the Wronskian matching integral collapsed to χ=π/2.
function _construct_R8_static_on_grid!(inhom, orbit::OrbitalParams{T}, l, num_coupled) where {T}
    grid = orbit.grid; gs = length(grid); r0i = orbit.r0_index
    R8coeff = -4*16*T(π)*orbit.L*legendre_P_thetaderiv(l, 0, T(π)/2; normalized=true)
    function matchint(piece::Symbol)
        r_p = orbit.p; f_p = _f(r_p)
        Ro = R8_static_out_solution(l, r_p); Rod = R8_static_out_rs_deriv_solution(l, r_p)
        Ri = R8_static_in_solution(l, r_p);  Rid = R8_static_in_rs_deriv_solution(l, r_p)
        wronsk = Ri*Rod - Ro*Rid
        T(π)*((piece === :in ? Ro : Ri)*f_p/(wronsk*r_p^2)*(1/(2*T(π))))
    end
    c_out = R8coeff*matchint(:out); c_in = R8coeff*matchint(:in)
    for i in 0:gs
        r_p = i == gs ? grid[r0i+1] : grid[i+1]
        f_p = _f(r_p); idx = i+1
        if i <= r0i
            inhom[num_coupled+1,1,idx] = c_in*R8_static_in_solution(l, r_p)
            inhom[num_coupled+1,3,idx] = c_in*R8_static_in_rs_deriv_solution(l, r_p)/f_p
        else
            inhom[num_coupled+1,1,idx] = c_out*R8_static_out_solution(l, r_p)
            inhom[num_coupled+1,3,idx] = c_out*R8_static_out_rs_deriv_solution(l, r_p)/f_p
        end
        inhom[num_coupled+1,2,idx]=0; inhom[num_coupled+1,4,idx]=0
        inhom[1,1,idx]=0; inhom[1,2,idx]=0; inhom[1,3,idx]=0; inhom[1,4,idx]=0
        if l > 1
            inhom[2,1,idx]=0; inhom[2,2,idx]=0; inhom[2,3,idx]=0; inhom[2,4,idx]=0
        end
    end
    inhom
end

# EVEN STATIC: construct R6,R7 (from R1,R3,R5) + zero R2,R4 gauge slots.
function construct_R6_R7_R2_R4_gauge_fields_on_grid!(inhom, orbit::OrbitalParams{T}, l, num_coupled) where {T}
    grid = orbit.grid; gs = length(grid); r0i = orbit.r0_index
    for i in 0:gs
        r_p = i == gs ? grid[r0i+1] : grid[i+1]
        f_p = _f(r_p); df = _df(r_p); idx = i+1
        R1  = inhom[1,1,idx] + im*inhom[1,2,idx]; R1d = inhom[1,3,idx] + im*inhom[1,4,idx]
        R3  = inhom[2,1,idx] + im*inhom[2,2,idx]; R3d = inhom[2,3,idx] + im*inhom[2,4,idx]
        R5  = inhom[3,1,idx] + im*inhom[3,2,idx]; R5d = inhom[3,3,idx] + im*inhom[3,4,idx]
        V_l = (one(T)/4)*f_p*(df/r_p + l*(l+1)/r_p^2)
        R6 = 1/(2*f_p)*(r_p*R1d - r_p*f_p*R3d + R1 - R5 - f_p*R3)
        R1dd = 1/f_p^2*(4*V_l*R1 + 4/r_p^2*f_p^2*R3d + 2*f_p/r_p^2*(1-4/r_p)*(R1-R5-f_p*R3) - 2*f_p^2/r_p^2*(1-6/r_p)*R6 - f_p*df*R1d)
        R3dd = 1/f_p^2*(4*V_l*R3 - 2*f_p/r_p^2*(R1-R5-(1-4/r_p)*(R3+R6)) - f_p*df*R3d)
        R6d = 1/f_p*(-1/(2*f_p)*(R1*df - R5*df - 2*f_p*R1d + r_p*df*R1d + 2*f_p^2*R3d + f_p*R5d - r_p*f_p*R1dd + r_p*f_p^2*R3dd))
        R7 = r_p*R5d + 2*R5 + l*(l+1)*R6
        R5dd = 1/f_p^2*(4*V_l*R5 + 4*f_p/r_p^2*((1-4.5/r_p)*R5 - (one(T)/2)*l*(l+1)*(R1-f_p*R3) + (one(T)/2)*(1-3/r_p)*(l*(l+1)*R6-R7)) - f_p*df*R5d)
        R7d = r_p*R5dd + R5d + 2*R5d + l*(l+1)*R6d
        inhom[num_coupled+1,1,idx]=real(R6); inhom[num_coupled+1,2,idx]=imag(R6)
        inhom[num_coupled+1,3,idx]=real(R6d); inhom[num_coupled+1,4,idx]=imag(R6d)
        inhom[num_coupled+2,1,idx]=real(R7); inhom[num_coupled+2,2,idx]=imag(R7)
        inhom[num_coupled+2,3,idx]=real(R7d); inhom[num_coupled+2,4,idx]=imag(R7d)
        inhom[num_coupled+3,1,idx]=0; inhom[num_coupled+3,2,idx]=0; inhom[num_coupled+3,3,idx]=0; inhom[num_coupled+3,4,idx]=0
        inhom[num_coupled+4,1,idx]=0; inhom[num_coupled+4,2,idx]=0; inhom[num_coupled+4,3,idx]=0; inhom[num_coupled+4,4,idx]=0
    end
    inhom
end

# ── R8 static series solutions (H1Lorenz static.jl libSch.c R8_static_*) ─────────────────────────
_factT(::Type{T}, x::Integer) where {T} = x <= 1 ? one(T) : prod(T(2):T(x))
_R8_alk(::Type{T}, l, k) where {T} = l*(l+1)*_factT(T, l+k-1)/(_factT(T, l-k+1)*_factT(T, k+1)*_factT(T, k))

function R8_static_in_solution(l, r_p::T) where {T}
    l == 1 && return r_p^2
    x = r_p/2 - 1; s = zero(T)
    for k in 0:(l+1); s += x/(1+x)*_R8_alk(T, l, k)*x^k; end
    s
end
function R8_static_in_rs_deriv_solution(l, r_p::T) where {T}
    fp = _f(r_p)
    l == 1 && return 2*r_p*fp
    x = r_p/2 - 1; s = zero(T)
    for k in 0:(l+1); s += 1/(1+x)*_R8_alk(T, l, k)*(k + 1/(1+x))*x^k; end
    fp/2*s
end
function _R8_out_series(r_p::T, l) where {T}
    c0 = one(T); sum8 = c0/r_p^l; Dsum8 = -l*c0/r_p^(l+1)
    k = 1; k_max = 1
    while k <= k_max
        c1 = 2*c0*(((k-1)+l)*((k-1)+l+1)-2)/(k*(k+2l+1))
        sum8 += c1/r_p^(l+k); Dsum8 += -(l+k)*c1/r_p^(l+k+1)
        rtol = T <: BigFloat ? max(big(2.0)^(10 - precision(T)), 8*eps(T)) : T(1e-16)
        (abs(c1/r_p^(l+k)/sum8) > rtol) && (k_max < 600) && (k_max += 1)
        c0 = c1; k += 1
    end
    (sum8, Dsum8)
end
R8_static_out_solution(l, r_p::T) where {T} = l == 1 ? 1/r_p : _R8_out_series(r_p, l)[1]
R8_static_out_rs_deriv_solution(l, r_p::T) where {T} =
    l == 1 ? -1/r_p^2*_f(r_p) : _f(r_p)*_R8_out_series(r_p, l)[2]

# ── per-(l,m) mode assembly + corpus write (H1Lorenz assemble.jl/driver.jl/io_hdf5.jl) ──────────

"""
    NModeData{T}

One assembled `(l,m)` tensor-harmonic mode: `inhom[field,part(1..4),grididx(1..gs+1)]` (the .h5
inhom_left/right source; `nf = num_coupled+num_gauge`, last grid index = right-limit at r0),
`C_inf`/`C_horiz` asymptotic amplitudes (2·num_coupled, re/im interleaved), and the cset `kind`.
"""
struct NModeData{T}
    l::Int
    m::Int
    omega::T
    inhom::Array{T,3}
    nf::Int
    num_coupled::Int
    C_inf::Vector{T}
    C_horiz::Vector{T}
    cset_kind::Symbol
end

"""
    compute_mode(orbit, l, m) -> NModeData

Compute one `(l,m)` mode for the circular Schwarzschild orbit `orbit` (radius p), ω = m·Ω_φ.
Handles every parity/special case: the l=m=0 analytic static monopole, the (l+m)even & m=0
even-static swap, EHS rescale to `inhom[field,part,grididx]`, the algebraic gauge fields, and the
l=1,m=0 odd-dipole pure-gauge h̄9 admixture.
"""
function compute_mode(orbit::OrbitalParams{T}, l::Int, m::Int) where {T}
    ω = m*orbit.Omega_phi
    gs = length(orbit.grid); r0i = orbit.r0_index

    if l == 0 && m == 0
        inhom = construct_static_monopole(orbit)
        nf = size(inhom, 1)
        return NModeData{T}(l, m, ω, inhom, nf, 3, zeros(T, 6), zeros(T, 6), :monopole)
    end

    cset = select_coupled_set(l, m; static=true)
    nc = cset.num_coupled; ng = cset.num_gauge; nf = nc + ng

    integ = integrate_field_equations(orbit, l, ω, cset)
    C_in, C_out = calculate_scaling_coefficients(orbit, l, m, integ, cset)

    inhom = zeros(T, nf, 4, gs + 1)
    hom_in = integ.hom_in; hom_out = integ.hom_out
    for i in 0:gs, k in 1:nc
        R = zero(Complex{T}); dR = zero(Complex{T})
        for d in 1:nc
            if i <= r0i
                j = i; C = C_in[d]
                Rv  = hom_in[k,d,1,j+1] + im*hom_in[k,d,2,j+1]
                dRv = hom_in[k,d,3,j+1] + im*hom_in[k,d,4,j+1]
            elseif i == gs
                j = gs - 1 - r0i; C = C_out[d]
                Rv  = hom_out[k,d,1,j+1] + im*hom_out[k,d,2,j+1]
                dRv = hom_out[k,d,3,j+1] + im*hom_out[k,d,4,j+1]
            else
                j = gs - 1 - i; C = C_out[d]
                Rv  = hom_out[k,d,1,j+1] + im*hom_out[k,d,2,j+1]
                dRv = hom_out[k,d,3,j+1] + im*hom_out[k,d,4,j+1]
            end
            R += C*Rv; dR += C*dRv
        end
        inhom[k,1,i+1]=real(R); inhom[k,2,i+1]=imag(R)
        inhom[k,3,i+1]=real(dR); inhom[k,4,i+1]=imag(dR)
    end

    l > 0 && cset.gauge!(inhom, orbit, l, m, ω, nc)

    # odd-dipole l=1,m=0 pure-gauge h̄9 fix (GSF_ecc.c L358-367)
    if l == 1 && m == 0
        A = -T(256)/3*sqrt(3*T(π))*orbit.E*orbit.Omega_phi/(orbit.p - 2)
        for i in 0:gs
            r = i == gs ? orbit.p : orbit.grid[i+1]
            inhom[1,1,i+1] += A/r^2
            inhom[1,3,i+1] += -2*A/r^3
        end
    end

    C_inf = zeros(T, 2nc); C_horiz = zeros(T, 2nc)
    for i in 1:nc, j in 1:nc
        C_inf[2i-1]   += real(C_out[j]*integ.C_out_hom[i,j])
        C_inf[2i]     += imag(C_out[j]*integ.C_out_hom[i,j])
        C_horiz[2i-1] += real(C_in[j]*integ.C_in_hom[i,j])
        C_horiz[2i]   += imag(C_in[j]*integ.C_in_hom[i,j])
    end
    NModeData{T}(l, m, ω, inhom, nf, nc, C_inf, C_horiz, cset.kind)
end

# build a single mode at the requested precision (per-task setprecision for BigFloat thread-safety)
function _compute_one(p, grid, ra, r0, rb, l, m, bits)
    if bits === nothing
        orbit = OrbitalParams(Float64(p), Float64.(grid), ra, r0, rb)
        return compute_mode(orbit, l, m), orbit
    end
    setprecision(BigFloat, bits) do
        orbit = OrbitalParams(BigFloat(p), BigFloat.(grid), ra, r0, rb)
        compute_mode(orbit, l, m), orbit
    end
end

"""
    compute_modes(p, lmax; grid, indices, outdir=nothing, bits=nothing, parallel=:threads) -> Dict

Compute all `(l,m)` modes (l=0..lmax, m=0..l) for a circular orbit at radius `p` on `grid` with
`indices=(ra,r0,rb)` (0-based). Threads.@threads sweeps the flat γ list into lock-free per-index
slots; OpenBLAS is pinned to one thread for the whole sweep (deadlock-safe + serial==threaded
bit-identical). `outdir` writes one `h1-l{l}m{m}.h5` per mode; `bits` selects BigFloat precision.
"""
function compute_modes(p, lmax::Int; grid, indices, outdir=nothing, bits=nothing,
                       parallel::Symbol=:threads)
    ra, r0, rb = indices
    outdir !== nothing && mkpath(outdir)
    modes = mode_list(lmax); n = length(modes)
    slots = Vector{Any}(undef, n)
    work = function (idx)
        l, m = modes[idx]
        nmode, orbit = _compute_one(p, grid, ra, r0, rb, l, m, bits)
        slots[idx] = nmode
        outdir !== nothing && write_h1_mode(joinpath(outdir, "h1-l$(l)m$(m).h5"), nmode, orbit)
        nothing
    end
    blas_single_thread() do
        if parallel === :threads
            Threads.@threads for idx in 1:n; work(idx); end
        else
            for idx in 1:n; work(idx); end
        end
    end
    Dict{Tuple{Int,Int},Any}(modes[i] => slots[i] for i in 1:n)
end

"""
    compute_one_mode(p, l, m; grid, indices, bits=nothing) -> (NModeData, OrbitalParams)

Single-mode entry (no file output). Same grid/precision options as `compute_modes`.
"""
function compute_one_mode(p, l, m; grid, indices, bits=nothing)
    ra, r0, rb = indices
    _compute_one(p, grid, ra, r0, rb, l, m, bits)
end

# ── h1 corpus write (build a Data.H1Schema in the loader layout + Data.write_mode_h5) ───────────

"""
    write_h1_mode(path, nmode, orbit) -> path

Write one `(l,m)` mode to `path` by building a `Data.H1Schema` in the loader `(4nf, grid)` layout
(inhom_left columns 0..r0index, inhom_right = right-limit-at-r0 + columns r0index+1..gs-1) and
calling `Data.write_mode_h5` (the versioned storage contract; HDF5 weakdep extension).
"""
function write_h1_mode(path, nmode::NModeData{T}, orbit::OrbitalParams{T}) where {T}
    schema = h1_schema(nmode, orbit)
    Data.write_mode_h5(path, schema)
end

"Assemble the `Data.H1Schema{Float64}` (loader layout) for one mode (shared by write + round-trip)."
function h1_schema(nmode::NModeData{T}, orbit::OrbitalParams{T}) where {T}
    grid = orbit.grid; gs = length(grid); r0i = orbit.r0_index
    nf = nmode.nf; inhom = nmode.inhom
    nL = r0i + 1
    L = Array{Float64}(undef, 4nf, nL)
    for i in 0:r0i, k in 1:nf, part in 1:4
        L[(k-1)*4 + part, i+1] = Float64(inhom[k, part, i+1])
    end
    nR = gs - r0i
    R = Array{Float64}(undef, 4nf, nR)
    for k in 1:nf, part in 1:4
        R[(k-1)*4 + part, 1] = Float64(inhom[k, part, gs+1])
    end
    for r in 1:(nR-1), k in 1:nf, part in 1:4
        R[(k-1)*4 + part, r+1] = Float64(inhom[k, part, r0i + r + 1])
    end
    H1Schema{Float64}(Float64.(grid), L, R, Float64.(nmode.C_inf), Float64.(nmode.C_horiz),
                      r0i, orbit.ra_index, orbit.rb_index)
end
