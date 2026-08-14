# =============================================================================================================
# SchwarzschildRadial — a=0 (Schwarzschild) FULL RADIAL RECONSTRUCTION: the physical radial functions
# P_{-2}(r),P_{-1}(r),h(r),κ(r) from the validated jumps (legacy KerrLorenzCirc/radial.jl).
#
# P_{-2},P_{-1},h: homogeneous physical solutions with a jump at r0, via variation-of-parameters on the
# Perturbation R^up/R^in (TeukolskyRadial). κ: INHOMOGENEOUS (𝒟Δ𝒟†−2iωr−λ0)κ=½r²h, solved by Green's-function
# VoP using the SAME s=0 R^up/R^in, with Abel-summed conditionally-convergent tails (Iup_inf/Iin_hor) cached.
# The JUMPS come from the consolidated closure (closure_a0) + jumps.jl (jumps_pm, hj1).
#
# This is the a=0 standalone surface AND the a→0 regression oracle the spin-general KerrComponents module
# needs. BigFloat inherits from BigFloat a,ω (Monodromy ν, MST series); the VoP quadratures are precision-aware
# via atol/rtol; the κ chain converges 256→512.
# =============================================================================================================

import Perturbation: TeukolskyRadial, teukolsky_radial_deriv
import QuadGK: quadgk
import Printf: @printf

export reconstruct, HomField, KappaField, val, der, der2, κval, κder, jump_data_radial, radial_gates

const _AA = 0.0   # this module reconstructs the a=0 metric

# ---- jumps per (r0,l): {J0,J1} of P_{-2},P_{-1},h,κ (a=0). Reuse the consolidated closure. ----
# Type-generic {T} (§9): element type inherits from r0 via orbit_background (BigFloat in ⇒ BigFloat out);
# the upstream chain (orbit_background/jumps_pm/hj1/closure_a0) is already {T}. Float64 in ⇒ values
# identical to the pre-{T} Float64 path (the conversions below are exact for T=Float64).
"""
    jump_data_radial(r0, l, m) -> NamedTuple

a=0 jump data at `(r0, l, m)`: the `{J0, J1}` jumps at the orbit radius ``r_0`` of the four radial
fields ``P_{-2}, P_{-1}, h, κ`` — `Pm2`/`Pm1`/`h`/`κ` as `(J0, J1)` tuples — from `jumps_pm`/`hj1`
plus the consolidated banded closure `closure_a0`, together with `ω, Δ0` and the eigenvalues
`Λ=(l-1)(l+2)`, `λ=λ0=l(l+1)`. Element type inherits from `r0` (BigFloat in ⇒ BigFloat out).
Legacy: `KerrLorenzCirc/radial.jl`.
"""
function jump_data_radial(r0, l, m)
    bg = orbit_background(0, r0, m)
    T = typeof(bg.r0); CT = Complex{T}
    ω = bg.ω; Δ0 = bg.Δ0
    J0m2,J1m2 = jumps_pm(bg, l, -1)
    J1h = hj1(bg, l)
    sol = closure_a0(r0, m; ls=l:l)[l]      # [J0(P_{-1}),J1(P_{-1}),J0(κ),J1(κ)]
    Λ = T((l-1)*(l+2)); λ = T(l*(l+1)); λ0 = T(l*(l+1))
    (ω=ω, Δ0=Δ0,
     Pm2=(CT(J0m2), CT(J1m2)), Pm1=(CT(sol[1]), CT(sol[2])),
     h=(zero(CT), CT(J1h)), κ=(CT(sol[3]), CT(sol[4])),
     Λ=Λ, λ=λ, λ0=λ0)
end

# generic c=iK/Δ, c'(r) (a=0)
cfun(r,ω,m)  = (Δ=r^2-2r+_AA^2; K=ω*(r^2+_AA^2)-_AA*m; im*K/Δ)
cpfun(r,ω,m) = (Δ=r^2-2r+_AA^2; K=ω*(r^2+_AA^2)-_AA*m; Kp=2ω*r; im*(Kp*Δ-K*(2r-2))/Δ^2)

# R'' from each homogeneous ODE, given R,R' (a=0).
function Rpp_sm2(r,ω,m,Λ,R,Rp)
    Δ=r^2-2r+_AA^2; Δp=2r-2; c=cfun(r,ω,m); cp=cpfun(r,ω,m)
    (-Δ*(-cp*R - c^2*R - (Δp/Δ)*(Rp - c*R)) + 6im*ω*r*R + Λ*R)/Δ
end
function Rpp_sm1(r,ω,m,λ,R,Rp)
    Δ=r^2-2r+_AA^2; c=cfun(r,ω,m); cp=cpfun(r,ω,m)
    (-Δ*(-cp*R - c^2*R) + 2im*ω*r*R + λ*R)/Δ
end
function Rpp_s0(r,ω,m,λ0,R,Rp)
    Δ=r^2-2r+_AA^2; Δp=2r-2; K=ω*(r^2+_AA^2)-_AA*m; Kp=2ω*r
    (-(Δp*Rp + (im*Kp + K^2/Δ)*R - 2im*ω*r*R - λ0*R))/Δ
end

# Build a homogeneous solution branch (value/deriv/2nd-deriv) for spin s. a=0 is passed at the
# precision of ω so a BigFloat ω builds a BigFloat TeukolskyRadial (precision inherits from a,ω).
function branch(s, l, m, ω, eig)
    R = TeukolskyRadial(s, l, m, oftype(float(real(ω)), _AA), ω; lambda=eig)
    Rpp = s==-2 ? Rpp_sm2 : s==-1 ? Rpp_sm1 : Rpp_s0
    up(r)  = R.Up(r);  upd(r)=teukolsky_radial_deriv(R.Up,r); updd(r)=Rpp(r,ω,m,eig,R.Up(r),teukolsky_radial_deriv(R.Up,r))
    inn(r) = R.In(r);  ind(r)=teukolsky_radial_deriv(R.In,r); indd(r)=Rpp(r,ω,m,eig,R.In(r),teukolsky_radial_deriv(R.In,r))
    (R=R, up=up, upd=upd, updd=updd, inn=inn, ind=ind, indd=indd, eig=eig)
end

# 2×2 jump solve ⇒ (A_up, A_in).
function jump_amps(br, r0, J0, J1)
    Rup=br.up(r0); Rin=br.inn(r0); Rupd=br.upd(r0); Rind=br.ind(r0)
    det = Rup*(-Rind) - (-Rin)*Rupd
    Aup = (J0*(-Rind) - (-Rin)*J1)/det
    Ain = (Rup*J1 - Rupd*J0)/det
    (Aup, Ain)
end

# {T}-parametric (S3-b): mirrors KerrRadial's HomFieldK{T}/KappaFieldK{T}. Float64 in ⇒ HomField{Float64}
# (identical to the pre-{T} behavior); BigFloat in ⇒ Complex{BigFloat} amplitudes, no silent demotion.
"""
    HomField{T}

a=0 jumped homogeneous radial field: ``A^{\\rm in} R^{\\rm in}(r)`` for ``r < r_0`` and
``A^{\\rm up} R^{\\rm up}(r)`` for ``r ≥ r_0``, on the `Perturbation.TeukolskyRadial` legs, with the
amplitudes fixed by the 2×2 jump solve at ``r_0``. Evaluate with [`val`](@ref) / [`der`](@ref) /
[`der2`](@ref). `Complex{T}` amplitudes (BigFloat-inherited). Legacy: `KerrLorenzCirc/radial.jl`.
"""
struct HomField{T<:Real}
    r0::T; Aup::Complex{T}; Ain::Complex{T}; br
end
"""
    val(f::HomField, r)

Value of the jumped homogeneous field: the In leg ``A^{\\rm in} R^{\\rm in}(r)`` for ``r < r_0``,
else the Up leg ``A^{\\rm up} R^{\\rm up}(r)``.
"""
val(f::HomField,r)  = r<f.r0 ? f.Ain*f.br.inn(r)  : f.Aup*f.br.up(r)
"""
    der(f::HomField, r)

First radial derivative of [`val`](@ref) (analytic, via `teukolsky_radial_deriv` on the active leg).
"""
der(f::HomField,r)  = r<f.r0 ? f.Ain*f.br.ind(r)  : f.Aup*f.br.upd(r)
"""
    der2(f::HomField, r)

Second radial derivative of [`val`](@ref), closed by the field's own homogeneous ODE
(``R''`` expressed through ``R, R'``).
"""
der2(f::HomField,r) = r<f.r0 ? f.Ain*f.br.indd(r) : f.Aup*f.br.updd(r)

"""
    KappaField{T}

a=0 inhomogeneous κ field: the solution of ``(𝒟Δ𝒟^† − 2iωr − λ_0)κ = \\tfrac12 r^2 h`` as the
jumped homogeneous piece plus a Green's-function variation-of-parameters particular integral on the
same s=0 ``R^{\\rm up}/R^{\\rm in}`` legs; the conditionally-convergent ``r→∞``/horizon tails are
Abel-regularized and cached (`Iup_inf`/`Iin_hor`). Evaluate with [`κval`](@ref) / [`κder`](@ref).
Legacy: `KerrLorenzCirc/radial.jl`.
"""
struct KappaField{T<:Real}
    r0::T; hom::HomField{T}; brs0; hfield::HomField{T}; W̃::Complex{T}; Rin_lo::T; Rout_hi::T
    Tup_tail::Complex{T}; Tin_head::Complex{T}; Iup_inf::Complex{T}; Iin_hor::Complex{T}; rtol::Float64
end
src(kf::KappaField, s) = 0.5*s^2*val(kf.hfield, s)
function kappa_I(kf::KappaField, r)
    brs0=kf.brs0; W̃=kf.W̃; r0=kf.r0; rt=kf.rtol
    integ_up(s) = brs0.up(s) * (0.5*s^2*val(kf.hfield,s)) / W̃
    integ_in(s) = brs0.inn(s) * (0.5*s^2*val(kf.hfield,s)) / W̃
    if r ≤ r0
        I1 = first(quadgk(integ_up, r, r0; rtol=rt, atol=1e-12)) + kf.Tup_tail
        if r ≥ kf.Rin_lo
            I2 = kf.Iin_hor + first(quadgk(integ_in, kf.Rin_lo, r; rtol=rt, atol=1e-12))
        else
            # Direct horizon→r Abel-regularized forward integral (avoids the cancellation of the
            # Iin_hor + ∫_{Rin_lo}^{r} difference for r<Rin_lo; same ε-regulator Richardson as abel0_in).
            ε1 = 4e-3; ε2 = 2e-3
            rh = one(kf.r0) + sqrt(one(kf.r0) - _AA^2)
            lo = rh + oftype(float(real(kf.r0)), 1e-6)
            fI(s, ε) = integ_in(s) * exp(-ε*(s - rh))
            t1 = first(quadgk(s -> fI(s, ε1), lo, r; rtol=1e-8, atol=0.0))
            t2 = first(quadgk(s -> fI(s, ε2), lo, r; rtol=1e-8, atol=0.0))
            I2 = 2t2 - t1
        end
    else
        I1 = first(quadgk(integ_up, r, kf.Rout_hi; rtol=rt, atol=1e-12)) + kf.Iup_inf
        I2 = kf.Tin_head + first(quadgk(integ_in, r0, r; rtol=rt, atol=1e-12))
    end
    (I1, I2)
end
kappa_p(kf::KappaField, r)       = let (I1,I2)=kappa_I(kf,r); kf.brs0.inn(r)*I1 + kf.brs0.up(r)*I2 end
kappa_p_deriv(kf::KappaField, r) = let (I1,I2)=kappa_I(kf,r); kf.brs0.ind(r)*I1 + kf.brs0.upd(r)*I2 end
"""
    κval(kf::KappaField, r)

``κ(r)``: the jumped homogeneous piece plus the Green's-function particular integral
(VoP on the s=0 legs, Abel-regularized tails).
"""
κval(kf::KappaField, r) = val(kf.hom, r) + kappa_p(kf, r)
"""
    κder(kf::KappaField, r)

First radial derivative ``κ'(r)`` (analytic: the VoP integrals differentiate onto the leg
derivatives only).
"""
κder(kf::KappaField, r) = der(kf.hom, r) + kappa_p_deriv(kf, r)

"""
    reconstruct(r0, l; m=2, Rin_lo=nothing, Rout_hi=nothing) -> NamedTuple

a=0 reconstructed radial fields at (r0,l): `rec.Pm2,rec.Pm1,rec.h` are `HomField` (use `val/der/der2`), `rec.κ`
is a `KappaField` (use `κval/κder`). `rec.jd` carries ω,Δ0,Λ,λ,λ0 and the jumps.
"""
function reconstruct(r0, l; m=2, Rin_lo=nothing, Rout_hi=nothing)
    jd = jump_data_radial(r0, l, m); ω=jd.ω
    T = typeof(float(real(ω)))          # element type inherited from r0 (via jump_data_radial's bg)
    brm2 = branch(-2, l, m, ω, jd.Λ); brm1 = branch(-1, l, m, ω, jd.λ); br0 = branch(0, l, m, ω, jd.λ0)
    Pm2 = HomField{T}(r0, jump_amps(brm2,r0,jd.Pm2...)..., brm2)
    Pm1 = HomField{T}(r0, jump_amps(brm1,r0,jd.Pm1...)..., brm1)
    h   = HomField{T}(r0, jump_amps(br0, r0,jd.h...)...,  br0)
    κhom = HomField{T}(r0, jump_amps(br0, r0, jd.κ...)..., br0)
    rW=T(r0); ΔW=rW^2-2rW+_AA^2
    W̃ = ΔW*(br0.inn(rW)*br0.upd(rW) - br0.up(rW)*br0.ind(rW))
    Rin_lo = Rin_lo===nothing ? rW - T(2.6) : T(Rin_lo)     # T(2.6)==2.6 exactly for T=Float64
    Rout_hi = Rout_hi===nothing ? rW + T(11) : T(Rout_hi)
    rtol=1e-9
    iup(s)=br0.up(s) *(0.5*s^2*val(h,s))/W̃
    iin(s)=br0.inn(s)*(0.5*s^2*val(h,s))/W̃
    abel_tail(integ, a, ε) = first(quadgk(s->integ(s)*exp(-ε*(s-a)), a, a+30.0/ε; rtol=1e-8, atol=1e-14))
    function abel0(integ, a)
        ε1=4e-3; ε2=2e-3
        t1=abel_tail(integ,a,ε1); t2=abel_tail(integ,a,ε2)
        2t2 - t1
    end
    Iup_inf = abel0(iup, Rout_hi)
    rh = one(T) + sqrt(one(T) - _AA^2)
    function abel0_in(integ, b)
        ε1=4e-3; ε2=2e-3
        f(s,ε)=integ(s)*exp(-ε*(s-rh))
        lo = rh + T(1e-6)
        t1=first(quadgk(s->f(s,ε1), lo, b; rtol=1e-8, atol=1e-14))
        t2=first(quadgk(s->f(s,ε2), lo, b; rtol=1e-8, atol=1e-14))
        2t2 - t1
    end
    Iin_hor = abel0_in(iin, Rin_lo)
    Tup_tail=first(quadgk(iup, rW, Rout_hi; rtol=rtol, atol=1e-12)) + Iup_inf
    Tin_head=first(quadgk(iin, Rin_lo, rW; rtol=rtol, atol=1e-12)) + Iin_hor
    κf = KappaField{T}(rW, κhom, br0, h, Complex{T}(W̃), Rin_lo, Rout_hi,
                       Complex{T}(Tup_tail), Complex{T}(Tin_head), Complex{T}(Iup_inf), Complex{T}(Iin_hor), rtol)
    (jd=jd, Pm2=Pm2, Pm1=Pm1, h=h, κ=κf)
end

# ---- residuals ----
res_sm2(r,ω,m,Λ,R,Rp,Rpp)=(Δ=r^2-2r+_AA^2;Δp=2r-2;c=cfun(r,ω,m);cp=cpfun(r,ω,m); Δ*(Rpp-cp*R-c^2*R-(Δp/Δ)*(Rp-c*R))-6im*ω*r*R-Λ*R)
res_sm1(r,ω,m,λ,R,Rp,Rpp)=(Δ=r^2-2r+_AA^2;c=cfun(r,ω,m);cp=cpfun(r,ω,m); Δ*(Rpp-cp*R-c^2*R)-2im*ω*r*R-λ*R)
res_s0(r,ω,m,λ0,R,Rp,Rpp)=(Δ=r^2-2r+_AA^2;Δp=2r-2;K=ω*(r^2+_AA^2)-_AA*m;Kp=2ω*r; Δ*Rpp+Δp*Rp+(im*Kp+K^2/Δ)*R-2im*ω*r*R-λ0*R)
function res_kappa(kf, r, m)
    jd_ω = kf.brs0.R.Up.omega; λ0 = kf.brs0.eig; brs0=kf.brs0
    local I1,I2; I1,I2 = kappa_I(kf, r)
    S = 0.5*r^2*val(kf.hfield,r)
    κ  = val(kf.hom,r)  + brs0.inn(r)*I1  + brs0.up(r)*I2
    κp = der(kf.hom,r)  + brs0.ind(r)*I1  + brs0.upd(r)*I2
    κpp = der2(kf.hom,r) + brs0.indd(r)*I1 + brs0.updd(r)*I2 + (-brs0.ind(r)*brs0.up(r) + brs0.upd(r)*brs0.inn(r))*S/kf.W̃
    Δ=r^2-2r+_AA^2; Δp=2r-2; K=jd_ω*(r^2+_AA^2)-_AA*m; Kp=2jd_ω*r
    Δ*κpp + Δp*κp + (im*Kp + K^2/Δ)*κ - 2im*jd_ω*r*κ - λ0*κ - S
end

"""
    radial_gates(r0; ls=2:8, m=2) -> (maxjump, maxres_h, maxjump_k, maxres_k)

Self-consistency gates for the a=0 radial reconstruction at orbit r0 (Printf output): P_{-2},P_{-1},h jumps +
continuity, the h s=0-ODE residual, and the κ jumps + inhomogeneous-ODE residual.
"""
function radial_gates(r0; ls=2:8, m=2)
    println("\n############################################################################")
    @printf("################  RADIAL RECONSTRUCTION GATES — a=0, m=2, r0=%.0f  ################\n", r0)
    println("############################################################################")
    ε=1e-7
    grid = [r0+1.5, r0+3.0, r0+6.0, r0+10.0, r0-1.0, r0-2.0]
    grid = filter(r->r>2.2, grid)
    maxjump=0.0; maxres_h=0.0; maxres_k=0.0; maxjump_k=0.0
    for l in ls
        rec = reconstruct(r0, l; m=m); jd=rec.jd; ω=jd.ω
        flush(stdout)
        @printf("\n=== l=%d  (ω=%.5f, Λ=%.1f, λ=%.1f, λ0=%.1f) ===\n", l, ω, jd.Λ, jd.λ, jd.λ0)
        for (nm,f,J) in (("P_{-2}",rec.Pm2,jd.Pm2), ("P_{-1}",rec.Pm1,jd.Pm1), ("h     ",rec.h,jd.h))
            jv = val(f,r0+ε)-val(f,r0-ε); jd1=der(f,r0+ε)-der(f,r0-ε)
            e0 = abs(J[1])>1e-8 ? abs(jv-J[1])/abs(J[1]) : abs(jv)
            e1 = abs(J[2])>1e-8 ? abs(jd1-J[2])/abs(J[2]) : abs(jd1)
            maxjump=max(maxjump,e0,e1)
            @printf("  GATE1 %s: [f]=%+.4e%+.4ei (J0 rel %.1e)  [f']=%+.4e%+.4ei (J1 rel %.1e)\n",
                    nm, real(jv),imag(jv),e0, real(jd1),imag(jd1),e1)
        end
        relh=maximum(let R=val(rec.h,r),Rp=der(rec.h,r),Rpp=der2(rec.h,r); abs(res_s0(r,ω,m,jd.λ0,R,Rp,Rpp))/(abs(jd.λ0*R)+1e-300) end for r in grid)
        maxres_h=max(maxres_h,relh)
        @printf("  GATE2 h ODE residual (max over grid) = %.2e\n", relh)
        kf=rec.κ
        jv = κval(kf,r0+ε)-κval(kf,r0-ε); jd1=κder(kf,r0+ε)-κder(kf,r0-ε)
        e0 = abs(jd.κ[1])>1e-8 ? abs(jv-jd.κ[1])/abs(jd.κ[1]) : abs(jv)
        e1 = abs(jd.κ[2])>1e-8 ? abs(jd1-jd.κ[2])/abs(jd.κ[2]) : abs(jd1)
        maxjump_k=max(maxjump_k,e0,e1)
        @printf("  GATE3 κ jumps: [κ]=%+.4e%+.4ei (J0=%+.4f rel %.1e)  [κ']=%+.4e%+.4ei (J1=%+.4f rel %.1e)\n",
                real(jv),imag(jv), real(jd.κ[1]),e0, real(jd1),imag(jd1), real(jd.κ[2]),e1)
        relk=0.0; rworst=0.0
        for r in grid
            res=res_kappa(kf,r,m); S=0.5*r^2*val(rec.h,r); scale=abs(S)+abs(jd.λ0*κval(kf,r))+1e-300
            rel=abs(res)/scale; if rel>relk; relk=rel; rworst=r; end
            maxres_k=max(maxres_k,rel)
        end
        @printf("  GATE3 κ inhomog ODE residual |L[κ]−½r²h|/scale (max over grid) = %.2e (at r=%.1f)\n", relk, rworst)
    end
    @printf("\n---- r0=%.0f SUMMARY ----\n", r0)
    @printf("  GATE1 max jump rel.err (P_{-2},P_{-1},h)      = %.2e  ⟹ %s\n", maxjump, maxjump<1e-4 ? "PASS" : "CHECK")
    @printf("  GATE2 max h s=0-ODE residual                  = %.2e  ⟹ %s\n", maxres_h, maxres_h<1e-8 ? "PASS" : "CHECK")
    @printf("  GATE3 max κ jump rel.err                      = %.2e  ⟹ %s\n", maxjump_k, maxjump_k<1e-4 ? "PASS" : "CHECK")
    @printf("  GATE3 max κ inhomog-ODE residual              = %.2e  ⟹ %s\n", maxres_k, maxres_k<1e-6 ? "PASS" : "CHECK")
    (maxjump,maxres_h,maxjump_k,maxres_k)
end
