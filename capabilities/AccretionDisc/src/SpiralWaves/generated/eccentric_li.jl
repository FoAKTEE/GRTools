# GENERATED: closure = li; do not edit by hand.
# generator: derivation/src/eccentric_emission.jl
# input: ECC_li.jls.txt sha256 e3edf9a3cac8af44bb4dcc1c883f5e503576e0d5f09c928a7a0e2585512e219c
# derivation: sha256 f6728339670d32802255cd9f5676bd7453c99b2c9507479851040c80c02539c0
# Inputs use the R3 alphabet (enu, epsi, emu are exponentials).
# htr/hrphi and barred amplitudes are independent complex inputs.
# Physical conjugacy is supplied by the caller; no reality projection is applied.
# The background is on shell; omega_d != 0 and all algebraic pivots must be nonzero.
# Scalar kernels accept a NamedTuple or a property-compatible input struct.

function _eccentric_li_C00(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    eps_r = inp.eps_r
    epsi = inp.epsi
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    return (((-eps_r - p_r) / (emu^2) - (eps + p)*((eps*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2)))) / p + ((-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*((W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)) / ((emu^2)*enu) + (-2W*(epsi^2)*(omt + (enu*v) / epsi)*psi_r) / ((emu^2)*enu))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + ((m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) - (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*((W*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)) / ((emu^2)*enu) + (W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*(omt + (enu*v) / epsi)) / ((emu^2)*enu))) / (((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + ((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((eps + p)^3)) / ((emu^2)*((((eps + p)^2) / (emu^2))^2))
end

function _eccentric_li_C01(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    v_r = inp.v_r
    return (((((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*((W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)) / ((emu^2)*enu) + (-2W*(epsi^2)*(omt + (enu*v) / epsi)*psi_r) / ((emu^2)*enu))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + ((-((W*(omt + (enu*v) / epsi)*p_r) / enu) - ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*((W*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)) / ((emu^2)*enu) + (W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*(omt + (enu*v) / epsi)) / ((emu^2)*enu))) / (((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))*((eps + p)^3)) / ((emu^2)*((((eps + p)^2) / (emu^2))^2))
end

function _eccentric_li_C10(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    p = inp.p
    v = inp.v
    return ((-eps - p)*((eps + p)^3)*((eps*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) / p + (m*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + ((omega_d + m*(omt + (enu*v) / epsi))*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))) / ((emu^4)*((((eps + p)^2) / (emu^2))^2))
end

function _eccentric_li_C11(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    eps_r = inp.eps_r
    epsi = inp.epsi
    m = inp.m
    mu_r = inp.mu_r
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    v_r = inp.v_r
    return ((-eps - p)*(eps_r + ((-((W*(omt + (enu*v) / epsi)*p_r) / enu) - ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*m) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + ((-((W*(omt + (enu*v) / epsi)*p_r) / enu) - ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(omega_d + m*(omt + (enu*v) / epsi))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (eps + p)*(mu_r + (1//2)*((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r) / (enu^2) + (((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2)) + (1//2)*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))))*((eps + p)^2)) / ((emu^4)*((((eps + p)^2) / (emu^2))^2))
end

function _eccentric_li_S0(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hphiphi = inp.hphiphi
    hphiphi_r = inp.hphiphi_r
    hrphi = inp.hrphi
    hrr = inp.hrr
    htphi = inp.htphi
    htphi_r = inp.htphi_r
    htr = inp.htr
    htt = inp.htt
    htt_r = inp.htt_r
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    return ((((hrr*p_r) / (emu^4) + ((eps + p)*((W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)) / ((emu^2)*enu) + (-2W*(epsi^2)*(omt + (enu*v) / epsi)*psi_r) / ((emu^2)*enu))*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + ((-eps - p)*((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2))*((W*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)) / ((emu^2)*enu) + (W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*(omt + (enu*v) / epsi)) / ((emu^2)*enu))) / ((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu) - (eps + p)*(((1//2)*(W^2)*((-htt_r) / (emu^2) + ((2(epsi^2)*omt*omt_r - 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) + 2(epsi^2)*(omt^2)*psi_r)*hrr) / (emu^4))) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)*((-htphi_r) / (emu^2) + ((-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*hrr) / (emu^4))) / (enu^2) + ((1//2)*(W^2)*((omt + (enu*v) / epsi)^2)*((-hphiphi_r) / (emu^2) + (2(epsi^2)*hrr*psi_r) / (emu^4))) / (enu^2)))*((eps + p)^3)) / ((emu^2)*((((eps + p)^2) / (emu^2))^2))) + im * (((((-((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2)) - (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)) / ((emu^2)*enu) + (-2W*(epsi^2)*(omt + (enu*v) / epsi)*psi_r) / ((emu^2)*enu))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*((W*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)) / ((emu^2)*enu) + (W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*(omt + (enu*v) / epsi)) / ((emu^2)*enu))) / (((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) - (eps + p)*((-(W^2)*htr*(omega_d + m*(omt + (enu*v) / epsi))) / ((emu^2)*(enu^2)) + ((W^2)*hrphi*m*((omt + (enu*v) / epsi)^2)) / ((emu^2)*(enu^2)) + ((W^2)*(htr*m - hrphi*(omega_d + m*(omt + (enu*v) / epsi)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((eps + p)^3)) / ((emu^2)*((((eps + p)^2) / (emu^2))^2)))
end

function _eccentric_li_S1(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hphiphi = inp.hphiphi
    hrphi = inp.hrphi
    hrr = inp.hrr
    htphi = inp.htphi
    htr = inp.htr
    htt = inp.htt
    hzz = inp.hzz
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    return (((-eps - p)*(((-eps - p)*m*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + ((-eps - p)*(omega_d + m*(omt + (enu*v) / epsi))*((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2))) / ((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu) + (eps + p)*(((1//2)*W*((htt*(omega_d + m*(omt + (enu*v) / epsi))) / (enu^2) + ((htt*m + 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2))) / enu + ((1//2)*W*(omt + (enu*v) / epsi)*((-htt*m) / (enu^2) + (hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2))) / enu + (-W*hrr*(omega_d + m*(omt + (enu*v) / epsi))) / (2(emu^2)*enu) + (W*hrr*m*(omt + (enu*v) / epsi)) / (2(emu^2)*enu) + ((-1//2)*W*hzz*(omega_d + m*(omt + (enu*v) / epsi))) / enu + ((1//2)*W*hzz*m*(omt + (enu*v) / epsi)) / enu + ((1//2)*W*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + ((1//2)*W*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu))*((eps + p)^2)) / ((emu^4)*((((eps + p)^2) / (emu^2))^2))) + im * (((((-((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2)) - (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*m) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + ((-((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2)) - (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(omega_d + m*(omt + (enu*v) / epsi))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (eps + p)*(((1//2)*W*((-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-hrphi*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) / (emu^2) + (htr*omt) / ((emu^2)*(enu^2))) + (2(epsi^2)*omt*omt_r - 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) + 2(epsi^2)*(omt^2)*psi_r)*(htr / ((emu^2)*(enu^2)) + (hrphi*omt) / ((emu^2)*(enu^2))))) / enu + (W*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*(htr / (enu^2) + (hrphi*omt) / (enu^2))) / (2(emu^2)*enu) + (W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*enu) + (-W*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)*psi_r) / ((emu^2)*enu) + ((1//2)*W*(2(epsi^2)*psi_r*((-hrphi*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) / (emu^2) + (htr*omt) / ((emu^2)*(enu^2))) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*(htr / ((emu^2)*(enu^2)) + (hrphi*omt) / ((emu^2)*(enu^2))))*(omt + (enu*v) / epsi)) / enu + (W*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*(omt + (enu*v) / epsi)*(htr / (enu^2) + (hrphi*omt) / (enu^2))) / (2(emu^2)*enu)))*((eps + p)^3)) / ((emu^4)*((((eps + p)^2) / (emu^2))^2)))
end

function _eccentric_li_dut(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hphiphi = inp.hphiphi
    hrphi = inp.hrphi
    htphi = inp.htphi
    htr = inp.htr
    htt = inp.htt
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    v_r = inp.v_r
    y0 = inp.y0
    y1 = inp.y1
    return (((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2)) / ((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu) + ((m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) - (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*y0*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*y1*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + im * (((-((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2)) - (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_duphi_up(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hphiphi = inp.hphiphi
    hrphi = inp.hrphi
    htphi = inp.htphi
    htr = inp.htr
    htt = inp.htt
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    v_r = inp.v_r
    y0 = inp.y0
    y1 = inp.y1
    return ((((-1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((-1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + (-(W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2)) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + ((-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*y0) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + ((-((W*(omt + (enu*v) / epsi)*p_r) / enu) - ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*y1) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + im * (((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2)))) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_deps(inp)
    eps = inp.eps
    p = inp.p
    y0 = inp.y0
    return ((eps + p)*eps*y0) / p
end

function _eccentric_li_duphi_low(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hphiphi = inp.hphiphi
    hrphi = inp.hrphi
    htphi = inp.htphi
    htr = inp.htr
    htt = inp.htt
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    v_r = inp.v_r
    y0 = inp.y0
    y1 = inp.y1
    return ((W*htphi) / enu + (W*hphiphi*(omt + (enu*v) / epsi)) / enu + (-(epsi^2)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(epsi^2)*omt*((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2))) / ((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu) - y1*((((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(epsi^2)) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(epsi^2)*omt*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + (eps + p)*y0*(((epsi^2)*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + ((epsi^2)*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*omt*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))) + im * ((((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(epsi^2)) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(epsi^2)*omt*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end
