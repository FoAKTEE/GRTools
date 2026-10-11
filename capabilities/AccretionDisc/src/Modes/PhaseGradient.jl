# MR19 Appendix A / mirrored MR cooling code: phase-gradient-error secant correction.
# Reuses the shooting basis; changing a boundary coefficient needs no ODE integration.

_phase_scale_height(d, r) = csiso(d, r) / Omega(d, r)
_mr19_scale_height(d, r) = csadi(d, r) / Omega(d, r)

"""
    phase_gradient_error(d, f, x, m)

Real phase-gradient error `epsilon = Im(conj(δh)*δh′)/abs2(δh) - Re(k)`.
No phase unwrapping or finite differencing of the sampled phase is needed. Zeros of δh
produce NaN; callers must keep boundary windows inside the nonzero solved domain.
"""
function phase_gradient_error(d, f, x, m)
    [imag(conj(f[i,1]) * f[i,2]) / abs2(f[i,1]) - real(kwave(d,x[i],m)) for i in eachindex(x)]
end

"""
    phase_error_diagnostic(d, f, x, m, i_in, i_out; scale_height)

Peak-to-peak epsilon in one scale height of each boundary. The default uses `csiso/Omega`
and nearest-grid cutoffs with half-open windows, exactly as mirrored `main.py:489–523`.
MR19 uses the adiabatic scale height; pass its callable as `scale_height` for that convention.
"""
function phase_error_diagnostic(d, f, x, m, i_in, i_out; scale_height = _phase_scale_height)
    cut_in = argmin(abs.(x .- (x[i_in] + scale_height(d,x[i_in]))))
    cut_out = argmin(abs.(x .- (x[i_out] - scale_height(d,x[i_out]))))
    cut_in > i_in && cut_out < i_out || throw(ArgumentError("phase windows need at least one grid interval"))
    errors = map((i_in:cut_in-1, cut_out:i_out-1)) do window
        lo, hi = extrema(imag(conj(f[i,1])*f[i,2])/abs2(f[i,1]) - real(kwave(d,x[i],m)) for i in window)
        hi - lo
    end
    errors
end

"""
    phase_gradient_refine(d, x, m, fh1, fh2, fih, Ci, Co, i_in, i_out;
                          baseline, delta=1e-6, scale_height,
                          safeguard=false, max_backtracks=6)

ONE secant step in `(Ci,Co)` using the two existing homogeneous solutions and particular
solution. The Jacobian uses `C*(1±delta)` and complex coefficient differences, as in mirrored
`main.py:641–697`. Only Co is updated for m=1. `baseline` is the pre-boundary-update mode
whose error drives the step (the ordering in main.py:639), not a silently rebuilt combination.
The returned X is that baseline error. Floating-point precision is inherited from the disc.
The default returns `(mode, C, X, jacobian)` and preserves the unguarded MR step exactly.

For the FR22 consumer, `safeguard=true` requires finite, strictly smaller X at both edges
for m≥2. For m=1 only X_out gates descent: C_in is fixed, and the relative L2 change
of enthalpy on [2r_in, r_corotation] must not exceed 0.1. An empty or nonfinite
resolved inner interval rejects the step. X_in remains diagnostic only because its
inner edge layer is not mesh-converged. This 10% bound is a conservative acceptance
policy, not a claimed physical error estimate. Try the full step, then halve its common
fraction at most `max_backtracks` times. Only Co varies at m=1. On failure, retain the
pre-step coefficients and recombine their mode. The supplied pre-update `baseline` still
drives the secant and acceptance comparison; it need not equal that recombination when
adaptive boundaries just moved. The extra `safeguard` result records `status` (:accepted,
:backtracked or :rejected), `fraction`, `backtracks`, `X_before` and the actual `X_after`.
`inner_response_change` records the m=1 bound (zero for other m). Integer and
harmonic modes use the same azimuthal-number branch. No ODE integrations are added.
"""
function phase_gradient_refine(d::NewtonianDisc{T}, x, m, fh1, fh2, fih, Ci, Co, i_in, i_out;
                               baseline, delta = T(1)/10^6, scale_height = _phase_scale_height,
                               safeguard::Bool = false, max_backtracks::Integer = 6) where {T}
    delta > 0 || throw(ArgumentError("delta must be positive"))
    max_backtracks >= 0 || throw(ArgumentError("max_backtracks must be nonnegative"))
    m1 = _azimuthal_number(m) == 1
    combine(ci,co) = begin
        a1,a2 = homogeneous_coefficients(fh1,fh2,fih,ci,co,i_in,i_out)
        a1 .* fh1 .+ a2 .* fh2 .+ fih
    end
    diagnostic(f) = phase_error_diagnostic(d,f,x,m,i_in,i_out; scale_height)
    X = diagnostic(baseline)
    jac = Matrix{Complex{T}}(undef,2,2)
    for j in 1:2
        cp,cm = (j == 1 ? Ci : Co) .* (1+delta,1-delta)
        fp = j == 1 ? combine(cp,Co) : combine(Ci,cp)
        fm = j == 1 ? combine(cm,Co) : combine(Ci,cm)
        ep,em = diagnostic(fp),diagnostic(fm)
        for i in 1:2
            jac[i,j] = (ep[i]-em[i])/(cp-cm)
        end
    end
    C_before = (Ci,Co)
    # The unresolved m=1 inner edge layer must not gate a phase correction.
    # Bound the change of enthalpy on the resolved inner disc instead (F11 of
    # the solver-generalisation certificate: stable for r >= 2 r_in).
    inner = searchsortedfirst(x, 2x[i_in]):corotation_index(d,x,m)
    inner_change(f) = begin
        isempty(inner) && return T(Inf)
        scale = sum(abs2, @view baseline[inner,1])
        change = sum(abs2(f[i,1]-baseline[i,1]) for i in inner)
        iszero(scale) ? (iszero(change) ? zero(T) : T(Inf)) : sqrt(change/scale)
    end
    if m1
        Co -= X[2]/jac[2,2]
    else
        det = jac[1,1]*jac[2,2]-jac[1,2]*jac[2,1]
        Ci -= (jac[2,2]*X[1]-jac[1,2]*X[2])/det
        Co -= (-jac[2,1]*X[1]+jac[1,1]*X[2])/det
    end
    mode = combine(Ci,Co)
    if safeguard
        C_full = (Ci,Co)
        fraction = one(T)
        for backtracks in 0:max_backtracks
            X_after = diagnostic(mode)
            response_change = m1 ? inner_change(mode) : zero(T)
            descent = m1 ? (isfinite(X[2]) && isfinite(X_after[2]) && X_after[2] < X[2] &&
                            isfinite(response_change) && response_change <= T(1)/10) :
                           (all(isfinite, X) && all(isfinite, X_after) && all(X_after .< X))
            if descent
                status = backtracks == 0 ? :accepted : :backtracked
                audit = (; status, fraction, backtracks, X_before=X, X_after,
                           inner_response_change=response_change)
                return (; mode, C=(Ci,Co), X, jacobian=jac, safeguard=audit)
            end
            backtracks == max_backtracks && break
            fraction /= 2
            Ci = m1 ? C_before[1] : C_before[1] + fraction*(C_full[1]-C_before[1])
            Co = C_before[2] + fraction*(C_full[2]-C_before[2])
            mode = combine(Ci,Co)
        end
        mode = combine(C_before...)
        audit = (; status=:rejected, fraction=zero(T), backtracks=max_backtracks,
                   X_before=X, X_after=diagnostic(mode),
                   inner_response_change=m1 ? inner_change(mode) : zero(T))
        return (; mode, C=C_before, X, jacobian=jac, safeguard=audit)
    end
    (; mode, C=(Ci,Co), X, jacobian=jac)
end

"MR19 App. A outgoing coefficient with D_S (the cooling code uses D instead)."
function mr19_boundary_coeffs(d, x, i_in, i_out, m)
    map((i_in,i_out)) do i
        r=x[i]; DS=Dfreq(d,r,m)+Nr2(d,r); dDS=dDdr(d,r,m)+dNr2dr(d,r)
        C=im*kwave(d,r,m)+dDS/(4DS)-dSigmadr(d,r)/(2Sigma(d,r))+dcsadidr(d,r)/(2csadi(d,r))-1/(2r)
        _azimuthal_number(m) == 1 && i == i_in ? C-2im*kwave(d,r,m) : C
    end
end

"""
    mr19_phase_refinement([T=Float64]; npts=60000, beta=1e8, lap_nphi=800,
                          rtol=1e-10, atol=1e-12)

MR19 Fig. 15 experiment: m=10, hp=0.1 (ADIABATIC sound-speed convention), p=q=1,
gamma=7/5, softening=0.06, pole displacement=1e-6, r=0.05…5, two re-seeds and one
phase step. Returns radius, before/after fractional errors, and shooting diagnostics.
Finite beta approximates no cooling; compare beta=1e7 for the approximation error.
The full published grid is the default. Uses the shared direct solver, never a second solver.
"""
function mr19_phase_refinement(::Type{T}=Float64; npts=60000, beta=T(10)^8, lap_nphi=800,
                               rtol=T(1)/10^10, atol=T(1)/10^12) where {T<:AbstractFloat}
    gamma=T(7)/5; hp=T(1)/10
    d=NewtonianDisc{T}(h0=hp/sqrt(gamma),q=one(T),p=one(T),gamma=gamma,beta=T(beta),
                       soft=T(6)/100,indirect=0,pole_disp=T(1)/10^6,lap_nphi=lap_nphi)
    x=loggrid(d,T(1)/20,T(5),npts)
    result=solve_m_mode(d,10,x; reint_steps=2,phgrad_steps=1,diagnostics=true,
                        boundary_coeffs=mr19_boundary_coeffs,phase_scale_height=_mr19_scale_height,
                        rtol,atol)
    k=real.([kwave(d,r,10) for r in x])
    fractional(f) = [iszero(k[i]) ? T(NaN) : e/k[i] for (i,e) in enumerate(phase_gradient_error(d,f,x,10))]
    (; r=x, before=fractional(result.before), after=fractional(result.mode), result)
end
