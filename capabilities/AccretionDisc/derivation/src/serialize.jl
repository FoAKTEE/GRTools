# serialize.jl — the hand-off artefacts of stages 2 and 3.
#
# Route R3 stops at SYMBOLIC expressions; turning them into `master_derived.jl` / `r2_derived.jl` with
# the `_master_entries(::MasterInputs)` signature and the R1/R2 emission conventions (n//d rationals,
# hoisted `enu^a*epsi^b*emu^c`, BigFloat-inheriting) is the next node, `EGD::jd-codegen`.  What is
# written here is the machine-readable INPUT of that node: one expression per line, in Symbolics'
# printed form, which `Meta.parse` reads back verbatim.
#
# Complex expressions are split into `.re` and `.im` because the alphabet is real; recombine as
# `re + im*im_part`.

"""
    expr_lines(name, e) -> Vector{String}

`name = <expr>` for a real expression; `name.re = …`, `name.im = …` for a complex one.
"""
expr_lines(name, e::Num) = [string(name, " = ", string(e))]
expr_lines(name, e::Complex{Num}) = [string(name, ".re = ", string(real(e))),
                                     string(name, ".im = ", string(imag(e)))]

"The named expressions of stage 2 (the master system and the eliminated quantities)."
function master_exports()
    ex = Pair{String,Any}[]
    for (n, e) in zip(MASTER_NAMES, MASTER_ONSHELL)
        push!(ex, n => e)
    end
    for (n, e) in zip(("C00", "C01", "C10", "C11", "S0", "S1"), MASTER_RAW)
        push!(ex, string(n, "_raw") => e)
    end
    push!(ex, "dut" => ELIMINATED_ONSHELL.dut)
    push!(ex, "duphi_up" => ELIMINATED_ONSHELL.duphi_up)
    push!(ex, "deps" => ELIMINATED_ONSHELL.deps)
    push!(ex, "duphi_low" => ELIMINATED_ONSHELL.duphi_low)
    ex
end

"The named expressions of stage 3 (the diagnostics integrands)."
function diagnostics_exports()
    ex = Pair{String,Any}[]
    lbl = ("t", "r", "phi", "z")
    for α in 1:4, β in α:4
        iszero_entry(DXI_T[α, β]) && continue
        push!(ex, string("dxiT_", lbl[α], lbl[β]) => to_omega_d(onshell_nu(DXI_T[α, β])))
    end
    push!(ex, "amf_mode" => to_omega_d(onshell_nu(AMF_MODE)))
    push!(ex, "torque_mode" => to_omega_d(onshell_nu(TORQUE_MODE)))
    push!(ex, "xi_r" => to_omega_d(onshell_nu(XI_UP[IR])))
    ex
end

const _MANIFEST = """
# R3 symbol manifest — arXiv:2608.17003, route R3 (Julia-native derivation)
#
# Every exported expression is a RATIONAL FUNCTION of the symbols below.  The alphabet carries exactly
# one algebraic relation, W^2 (1 - v^2) = 1, and no transcendental function: e^nu, e^psi, e^mu appear
# only through enu, epsi, emu, which is what makes the emitted code inherit any float precision.
#
# background geometry, line element (4)   ds^2 = -e^{2nu} dt^2 + e^{2psi}(dphi - omt dt)^2 + e^{2mu} dr^2 + dz^2
#   enu, epsi, emu   e^nu, e^psi, e^mu                          (SQRT(-g) = enu*epsi*emu)
#   omt              omega-tilde(r)                              (frame dragging)
#   nu_r, psi_r, mu_r, omt_r    d/dr of nu, psi, mu, omega-tilde
# background fluid, (8)-(12)
#   v                ZAMO speed (10),  Omega = omt + v*enu/epsi
#   v_r              dv/dr                                       (Omega' follows; see below)
#   W                Lorentz factor 1/sqrt(1-v^2) (8),  u^t = W/enu
#   eps, eps_r       energy density and its r-derivative
#   p, p_r           height-integrated pressure and its r-derivative
#   Gamma1           adiabatic index (22)
# mode labels, (45)-(51)
#   m                azimuthal number;  omega_d = omega - m*Omega  Doppler frequency (51)
# master variables, (48),(49)
#   y0 = delta-U-hat = delta-p-hat/(eps+p);   y1 = -i delta-u-hat^r
#   y0c, y1c         complex conjugates (diagnostics bilinears only)
# metric perturbation modes, (14),(60) -- COVARIANT (ERRATUM E-M2)
#   htt, htr, htphi, hrr, hrphi, hphiphi, hzz and their _r derivatives
#   *c, *c_r         complex conjugates (diagnostics bilinears only)
#   htr = hrphi = 0 for a circular equatorial orbit; h_tz = h_rz = h_phiz = 0 by equatorial symmetry
# gauge
#   Fg, Fg_r, Fgc, Fgc_r    arbitrary trivial displacement xi -> xi + F u (diagnostics only; drops out)
#
# CONVERSIONS the consumer may need
#   Omega   = omt + v*enu/epsi                                     (10)
#   Omega_r = omt_r + enu*(v_r + v*(nu_r - psi_r))/epsi            (d/dr of the above)
#   v_r     = (Omega_r - omt_r)*epsi/enu - v*(nu_r - psi_r)        (the inverse; `MasterInputs` carries
#                                                                   Omega_r, `dv_dr` carries v_r)
#   ON-SHELL form: nu_r eliminated by the hydrostatic balance (11),
#     nu_r = v^2*psi_r - epsi*v*omt_r/enu - p_r/((eps+p)*W^2)
#   the *_raw entries keep nu_r AND p_r independent (off-shell); on a background satisfying (11) the
#   two forms agree.
#
# AVERAGING conventions of the diagnostics (73)
#   F_J   = sqrt(-g) * 4*pi * sum_{m>=1} Re[amf_mode]
#   <T>   = 4*pi * sum_{m>=1} Re[torque_mode]
"""

"""
    write_exports(path, exports; header)

Write `name = expr` lines plus the symbol manifest.  Deterministic: the same expressions produce a
byte-identical file, which is what makes the regeneration check of `EGD::jd-codegen` meaningful.
"""
function write_exports(path::AbstractString, exports; header::AbstractString = "")
    open(path, "w") do io
        println(io, header)
        print(io, _MANIFEST)
        println(io)
        for (n, e) in exports
            for line in expr_lines(n, e)
                println(io, line)
            end
        end
    end
    path
end
