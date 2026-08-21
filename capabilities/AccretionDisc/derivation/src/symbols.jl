# symbols.jl — the symbol alphabet of the route-R3 derivation and the radial derivation `ddr`.
#
# DESIGN (binding, see derivation/README.md):
#   • the background is GENERIC: ν, ψ, μ, ω̃ of the Page–Thorne line element (4) are FREE functions of
#     r, represented by the *exponentiated* symbols enu = e^ν, epsi = e^ψ, emu = e^μ plus ω̃ = `omt`
#     and by INDEPENDENT first-derivative symbols nu_r, psi_r, mu_r, omt_r.  The Kerr closed forms
#     (A1a–e) are NEVER used here — they live in `Spacetime.kerr_equatorial_potentials` and enter only
#     at evaluation time.  Exponentiating keeps every expression a RATIONAL FUNCTION of the alphabet,
#     which is what makes the exact-rational certificates of certify.jl exact (no transcendentals, no
#     rounding) and what matches the emission convention `enu^a*epsi^b*emu^c` of the R1/R2 generators.
#   • the Lorentz factor W = 1/√(1−v²) (8) is an INDEPENDENT symbol subject to the single algebraic
#     relation W²(1−v²) = 1, so that √ never appears.  Certificate points are drawn from the rational
#     parametrisation v = 2t/(1+t²), W = (1+t²)/(1−t²) which satisfies the relation EXACTLY.
#   • Ω is NOT a symbol: Ω = ω̃ + v e^{ν−ψ} (10) is an expression, so d/dr and the Newtonian limit are
#     automatic.  ω = mΩ• (46) IS a constant symbol (`omega`); the Doppler frequency ω_d = ω − mΩ (51)
#     is introduced only at the very END (`to_omega_d`), exactly as route R2 did.
#   • d/dr is the chain-rule derivation `ddr` built on `Symbolics.derivative` (never a hand-written
#     tree walker): ddr(f) = Σ_s (∂f/∂s) · s′ over the symbols actually occurring in f.
#
# Second-derivative symbols (`*_rr`) exist so that `ddr` is total; the reduction must not produce any
# (route R2 certified "no second derivatives"), and `assert_no_second_derivatives` is the gate.

using Symbolics

# ---------------------------------------------------------------------------------------------
# alphabet
# ---------------------------------------------------------------------------------------------

# background geometry (4): e^ν, e^ψ, e^μ, ω̃ and their r-derivatives (ν′, ψ′, μ′, ω̃′)
@variables enu epsi emu omt
@variables nu_r psi_r mu_r omt_r
@variables nu_rr psi_rr mu_rr omt_rr

# background fluid: ZAMO speed v (10), Lorentz factor W = 1/√(1−v²) (8), equation of state
@variables v v_r v_rr W
@variables eps eps_r eps_rr p p_r p_rr Gamma1 Gamma1_r

# mode labels: azimuthal number m, forcing frequency ω = mΩ• (46), Doppler frequency ω_d (51)
@variables m omega omega_d

# Eulerian perturbation amplitudes (45a): δp̂, δε̂, δû^t, δû^r, δû^φ  (hats dropped in the names)
@variables dp dp_r dp_rr
@variables deps deps_r
@variables dut dut_r dur dur_r dur_rr duphi duphi_r

# master variables (49): y0 = δÛ = δp̂/(ε+p), y1 = −i δû^r
@variables y0 y0_r y1 y1_r

# COVARIANT Fourier metric-perturbation modes ĥ_μν = δĝ_μν (14), (60) — index positions per ERRATUM
# E-M2.  ĥ_tz = ĥ_rz = ĥ_φz = 0 by the equatorial reflection symmetry of the background.
@variables htt htr htphi hrr hrphi hphiphi hzz
@variables htt_r htr_r htphi_r hrr_r hrphi_r hphiphi_r hzz_r
@variables htt_rr htr_rr htphi_rr hrr_rr hrphi_rr hphiphi_rr hzz_rr

"Coordinate order of every 4-index object in this project: (t, r, φ, z)."
const COORDS = (:t, :r, :φ, :z)
const IT, IR, IPHI, IZ = 1, 2, 3, 4

# ---------------------------------------------------------------------------------------------
# the radial derivation d/dr
# ---------------------------------------------------------------------------------------------

"""
    Ω_expr, Ω_r_expr

Ω = ω̃ + v e^{ν−ψ} (10) and its r-derivative, as EXPRESSIONS in the alphabet (never symbols).
"""
const Ω_expr = omt + v * enu / epsi
const Ω_r_expr = omt_r + enu * (v_r + v * (nu_r - psi_r)) / epsi

# d/dr rule table.  Keys are the underlying SymbolicUtils symbols (`Symbolics.value`).
const _DDR = Dict{Any,Num}()
_rule!(s, rhs) = (_DDR[Symbolics.value(s)] = Num(rhs); nothing)

for (s, sr) in ((enu, enu * nu_r), (epsi, epsi * psi_r), (emu, emu * mu_r), (omt, omt_r),
                (nu_r, nu_rr), (psi_r, psi_rr), (mu_r, mu_rr), (omt_r, omt_rr),
                (v, v_r), (v_r, v_rr), (W, W^3 * v * v_r),
                (eps, eps_r), (eps_r, eps_rr), (p, p_r), (p_r, p_rr), (Gamma1, Gamma1_r),
                (m, 0), (omega, 0), (omega_d, -m * Ω_r_expr),
                (dp, dp_r), (dp_r, dp_rr), (deps, deps_r),
                (dur, dur_r), (dur_r, dur_rr), (dut, dut_r), (duphi, duphi_r),
                (y0, y0_r), (y1, y1_r),
                (htt, htt_r), (htr, htr_r), (htphi, htphi_r), (hrr, hrr_r),
                (hrphi, hrphi_r), (hphiphi, hphiphi_r), (hzz, hzz_r),
                (htt_r, htt_rr), (htr_r, htr_rr), (htphi_r, htphi_rr), (hrr_r, hrr_rr),
                (hrphi_r, hrphi_rr), (hphiphi_r, hphiphi_rr), (hzz_r, hzz_rr))
    _rule!(s, sr)
end

# δû^t and δû^φ are ELIMINATED algebraically, so their r-derivatives must not survive: they DO get
# formed inside `four_acceleration_up` (the term u^ν∂_νδû^μ) but are multiplied by the background
# u^r = 0 and drop.  `dut_r`, `duphi_r` therefore exist, and `ELIMINATED_DERIVS` is the gate that they
# are absent from the reduced equations — an assertion, not an assumption.
const ELIMINATED_DERIVS = [dut_r, duphi_r, deps_r, dur_rr, dp_rr]

"""
    ddr(x)

Radial derivative d/dr of a `Num`/`Complex{Num}`/number, by the chain rule over the alphabet.
Errors on a symbol with no rule — a new symbol must declare its derivative.
"""
function ddr(x::Num)
    vars = Symbolics.get_variables(x)
    isempty(vars) && return Num(0)
    acc = Num(0)
    for s in vars
        haskey(_DDR, s) || error("ddr: no derivation rule for symbol $(s)")
        rule = _DDR[s]
        isequal(rule, Num(0)) && continue
        acc += Symbolics.derivative(x, Num(s)) * rule
    end
    acc
end
ddr(x::Complex{Num}) = ddr(real(x)) + im * ddr(imag(x))
ddr(x::Number) = zero(x) * Num(1)

"""
    SECOND_DERIVS

The `*_rr` symbols.  No expression of the reduction may contain one.
"""
const SECOND_DERIVS = [nu_rr, psi_rr, mu_rr, omt_rr, v_rr, eps_rr, p_rr, Gamma1_r, dp_rr, dur_rr,
                       htt_rr, htr_rr, htphi_rr, hrr_rr, hrphi_rr, hphiphi_rr, hzz_rr]

"`true` iff `x` contains none of [`SECOND_DERIVS`](@ref)."
function no_second_derivatives(x)
    vs = Set(Symbolics.get_variables(x isa Complex ? real(x) + imag(x) : x))
    all(s -> !(Symbolics.value(s) in vs), SECOND_DERIVS)
end

"""
    to_omega_d(x)

Replace the constant forcing frequency ω by ω_d + mΩ (51), i.e. rewrite in terms of the Doppler
frequency — the final variable of the master system, exactly as route R2 did.
"""
to_omega_d(x::Num) = Symbolics.substitute(x, Dict(omega => omega_d + m * Ω_expr))
to_omega_d(x::Complex{Num}) = to_omega_d(real(x)) + im * to_omega_d(imag(x))
to_omega_d(x::AbstractArray) = map(to_omega_d, x)
