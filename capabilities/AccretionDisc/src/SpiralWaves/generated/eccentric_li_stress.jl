# GENERATED: closure = li_stress; do not edit by hand.
# generator: derivation/src/eccentric_emission.jl
# input: ECC_li_stress.jls.txt sha256 54c2b8c446b83413744ac5f85bc5cb00391a520758e56c81b2a750d17b8f6f5c
# derivation: sha256 f6728339670d32802255cd9f5676bd7453c99b2c9507479851040c80c02539c0
# Inputs use the R3 alphabet (enu, epsi, emu are exponentials).
# htr/hrphi and barred amplitudes are independent complex inputs.
# Physical conjugacy is supplied by the caller; no reality projection is applied.
# The background is on shell; omega_d != 0 and all algebraic pivots must be nonzero.
# Scalar kernels accept a NamedTuple or a property-compatible input struct.

function _eccentric_li_stress_T0_tt(inp)
    W = inp.W
    enu = inp.enu
    eps = inp.eps
    p = inp.p
    return (-p) / (enu^2) + ((W^2)*(eps + p)) / (enu^2)
end

function _eccentric_li_stress_dT_tt(inp)
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
    return ((2W*(eps + p)*((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2))) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) + (2W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*y1*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + p*(((-htt) / (enu^2) + (-htphi*omt) / (enu^2)) / (enu^2) + (-omt*(htphi / (enu^2) + (hphiphi*omt) / (enu^2))) / (enu^2)) + (eps + p)*y0*((W^2) / (enu^2) + ((W^2)*eps) / ((enu^2)*p) + (-2W*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + -1 / (enu^2))) + im * ((-2W*((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_stress_dTbar_tt(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hphiphic = inp.hphiphic
    hrphic = inp.hrphic
    htphic = inp.htphic
    htrc = inp.htrc
    httc = inp.httc
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    v_r = inp.v_r
    y0c = inp.y0c
    y1c = inp.y1c
    return ((2W*(eps + p)*((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((httc*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-httc*m - 2htphic*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphic*m - hphiphic*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphic*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-httc*m*omt) / (enu^2) - hphiphic*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*httc) / (enu^2) + (-2(W^2)*htphic*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphic*((omt + (enu*v) / epsi)^2)) / (enu^2))) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) + (2W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*y1c*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + p*(((-httc) / (enu^2) + (-htphic*omt) / (enu^2)) / (enu^2) + (-omt*(htphic / (enu^2) + (hphiphic*omt) / (enu^2))) / (enu^2)) + (eps + p)*y0c*((W^2) / (enu^2) + ((W^2)*eps) / ((enu^2)*p) + (-2W*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + -1 / (enu^2))) + im * ((2W*((((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_stress_dT_tr(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    hrphi = inp.hrphi
    htr = inp.htr
    omt = inp.omt
    p = inp.p
    y1 = inp.y1
    return ((p*(htr / (enu^2) + (hrphi*omt) / (enu^2))) / (emu^2)) + im * ((W*(eps + p)*y1) / enu)
end

function _eccentric_li_stress_dTbar_tr(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    hrphic = inp.hrphic
    htrc = inp.htrc
    omt = inp.omt
    p = inp.p
    y1c = inp.y1c
    return ((p*(htrc / (enu^2) + (hrphic*omt) / (enu^2))) / (emu^2)) + im * ((-W*(eps + p)*y1c) / enu)
end

function _eccentric_li_stress_T0_tphi(inp)
    W = inp.W
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    omt = inp.omt
    p = inp.p
    v = inp.v
    return (-omt*p) / (enu^2) + ((W^2)*(eps + p)*(omt + (enu*v) / epsi)) / (enu^2)
end

function _eccentric_li_stress_dT_tphi(inp)
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
    return ((-W*(eps + p)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (W*(eps + p)*(omt + (enu*v) / epsi)*((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2))) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) - y1*((W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + (eps + p)*y0*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2) + ((W^2)*eps*(omt + (enu*v) / epsi)) / ((enu^2)*p) + (W*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + p*((-omt*(htt / (enu^2) + (htphi*omt) / (enu^2))) / (enu^2) + (htphi / (enu^2) + (hphiphi*omt) / (enu^2))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) + im * ((W*((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_stress_dTbar_tphi(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hphiphic = inp.hphiphic
    hrphic = inp.hrphic
    htphic = inp.htphic
    htrc = inp.htrc
    httc = inp.httc
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    v_r = inp.v_r
    y0c = inp.y0c
    y1c = inp.y1c
    return ((-W*(eps + p)*(((1//2)*(W^2)*((httc*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-httc*m - 2htphic*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphic*m - hphiphic*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphic*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-httc*m*omt) / (enu^2) - hphiphic*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (W*(eps + p)*(omt + (enu*v) / epsi)*((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((httc*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-httc*m - 2htphic*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphic*m - hphiphic*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphic*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-httc*m*omt) / (enu^2) - hphiphic*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*httc) / (enu^2) + (-2(W^2)*htphic*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphic*((omt + (enu*v) / epsi)^2)) / (enu^2))) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) - y1c*((W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + (eps + p)*y0c*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2) + ((W^2)*eps*(omt + (enu*v) / epsi)) / ((enu^2)*p) + (W*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + p*((-omt*(httc / (enu^2) + (htphic*omt) / (enu^2))) / (enu^2) + (htphic / (enu^2) + (hphiphic*omt) / (enu^2))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) + im * ((-W*((((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (W*((((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_stress_T0_rr(inp)
    emu = inp.emu
    p = inp.p
    return p / (emu^2)
end

function _eccentric_li_stress_dT_rr(inp)
    emu = inp.emu
    eps = inp.eps
    hrr = inp.hrr
    p = inp.p
    y0 = inp.y0
    return ((eps + p)*y0) / (emu^2) + (-hrr*p) / (emu^4)
end

function _eccentric_li_stress_dTbar_rr(inp)
    emu = inp.emu
    eps = inp.eps
    hrrc = inp.hrrc
    p = inp.p
    y0c = inp.y0c
    return ((eps + p)*y0c) / (emu^2) + (-hrrc*p) / (emu^4)
end

function _eccentric_li_stress_dT_rphi(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hrphi = inp.hrphi
    htr = inp.htr
    omt = inp.omt
    p = inp.p
    v = inp.v
    y1 = inp.y1
    return (p*((-hrphi*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) / (emu^2) + (htr*omt) / ((emu^2)*(enu^2)))) + im * ((W*(eps + p)*(omt + (enu*v) / epsi)*y1) / enu)
end

function _eccentric_li_stress_dTbar_rphi(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hrphic = inp.hrphic
    htrc = inp.htrc
    omt = inp.omt
    p = inp.p
    v = inp.v
    y1c = inp.y1c
    return (p*((-hrphic*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) / (emu^2) + (htrc*omt) / ((emu^2)*(enu^2)))) + im * ((-W*(eps + p)*(omt + (enu*v) / epsi)*y1c) / enu)
end

function _eccentric_li_stress_T0_phiphi(inp)
    W = inp.W
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    omt = inp.omt
    p = inp.p
    v = inp.v
    return ((W^2)*(eps + p)*((omt + (enu*v) / epsi)^2)) / (enu^2) + p*((-(omt^2)) / (enu^2) + 1 / (epsi^2))
end

function _eccentric_li_stress_dT_phiphi(inp)
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
    return ((-2W*(eps + p)*(omt + (enu*v) / epsi)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-2W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(omt + (enu*v) / epsi)*y1) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (eps + p)*y0*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*eps*((omt + (enu*v) / epsi)^2)) / ((enu^2)*p) + (2W*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*(omt + (enu*v) / epsi)) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + 1 / (epsi^2)) + (((-((htt*omt) / (enu^2)) - htphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*omt) / (enu^2) + ((htphi*omt) / (enu^2) + hphiphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*p) + im * ((2W*((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(omt + (enu*v) / epsi)) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_stress_dTbar_phiphi(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hphiphic = inp.hphiphic
    hrphic = inp.hrphic
    htphic = inp.htphic
    htrc = inp.htrc
    httc = inp.httc
    m = inp.m
    omega_d = inp.omega_d
    omt = inp.omt
    omt_r = inp.omt_r
    p = inp.p
    p_r = inp.p_r
    psi_r = inp.psi_r
    v = inp.v
    v_r = inp.v_r
    y0c = inp.y0c
    y1c = inp.y1c
    return ((-2W*(eps + p)*(omt + (enu*v) / epsi)*(((1//2)*(W^2)*((httc*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-httc*m - 2htphic*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphic*m - hphiphic*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphic*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-httc*m*omt) / (enu^2) - hphiphic*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-2W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(omt + (enu*v) / epsi)*y1c) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (eps + p)*y0c*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*eps*((omt + (enu*v) / epsi)^2)) / ((enu^2)*p) + (2W*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*(omt + (enu*v) / epsi)) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + 1 / (epsi^2)) + (((-((httc*omt) / (enu^2)) - htphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*omt) / (enu^2) + ((htphic*omt) / (enu^2) + hphiphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*p) + im * ((-2W*((((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htrc*omt) / (enu^2) + hrphic*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(omt + (enu*v) / epsi)) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_stress_T0_zz(inp)
    p = inp.p
    return p
end

function _eccentric_li_stress_dT_zz(inp)
    eps = inp.eps
    hzz = inp.hzz
    p = inp.p
    y0 = inp.y0
    return (eps + p)*y0 - hzz*p
end

function _eccentric_li_stress_dTbar_zz(inp)
    eps = inp.eps
    hzzc = inp.hzzc
    p = inp.p
    y0c = inp.y0c
    return (eps + p)*y0c - hzzc*p
end

function _eccentric_li_stress_dJt_t(inp)
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
    return ((((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2))*((2W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) / ((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu) + (W*(eps + p)*(epsi^2)*omt*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + htphi*((-omt*p) / (enu^2) + ((W^2)*(eps + p)*(omt + (enu*v) / epsi)) / (enu^2)) + htt*((-p) / (enu^2) + ((W^2)*(eps + p)) / (enu^2)) - y1*(((-((W*(omt + (enu*v) / epsi)*p_r) / enu) - ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*((2W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(epsi^2)*omt) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + (-(enu^2) + (epsi^2)*(omt^2))*p*(((-htt) / (enu^2) + (-htphi*omt) / (enu^2)) / (enu^2) + (-omt*(htphi / (enu^2) + (hphiphi*omt) / (enu^2))) / (enu^2)) - (epsi^2)*omt*p*((-omt*(htt / (enu^2) + (htphi*omt) / (enu^2))) / (enu^2) + (htphi / (enu^2) + (hphiphi*omt) / (enu^2))*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) + (eps + p)*((eps*(((W^2)*(-(enu^2) + (epsi^2)*(omt^2))) / (enu^2) + (-(W^2)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / (enu^2))) / p + ((m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) - (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*((2W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*(epsi^2)*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*omt) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) - (epsi^2)*omt*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)) + (-(enu^2) + (epsi^2)*(omt^2))*((W^2) / (enu^2) + -1 / (enu^2)))*y0) + im * (((-((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2)) - (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*((2W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(epsi^2)*omt) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_stress_dJt_r(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hrphi = inp.hrphi
    htr = inp.htr
    omt = inp.omt
    p = inp.p
    v = inp.v
    y1 = inp.y1
    return ((htr*p) / (emu^2) - (epsi^2)*omt*p*((-hrphi*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) / (emu^2) + (htr*omt) / ((emu^2)*(enu^2))) + (-(enu^2) + (epsi^2)*(omt^2))*p*(htr / ((emu^2)*(enu^2)) + (hrphi*omt) / ((emu^2)*(enu^2)))) + im * (y1*((W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu))
end

function _eccentric_li_stress_dJt_phi(inp)
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
    return (((-(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2)) - (((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2)) - (((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2)))*((W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))*(omt + (enu*v) / epsi)*((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2))) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) + htt*((-omt*p) / (enu^2) + ((W^2)*(eps + p)*(omt + (enu*v) / epsi)) / (enu^2)) - y1*((((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*((W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(-(enu^2) + (epsi^2)*(omt^2))*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + htphi*(((W^2)*(eps + p)*((omt + (enu*v) / epsi)^2)) / (enu^2) + p*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) + (-(enu^2) + (epsi^2)*(omt^2))*p*(((-htt*omt) / (enu^2) - htphi*((omt^2) / (enu^2) + -1 / (epsi^2))) / (enu^2) + ((-((htphi*omt) / (enu^2)) - hphiphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*omt) / (enu^2)) + (eps + p)*((eps*(((W^2)*(-(enu^2) + (epsi^2)*(omt^2))*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*(epsi^2)*omt*((omt + (enu*v) / epsi)^2)) / (enu^2))) / p + ((-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*((W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*(-(enu^2) + (epsi^2)*(omt^2))*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) - (epsi^2)*omt*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (-(enu^2) + (epsi^2)*(omt^2))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*y0 - (epsi^2)*(((-((htt*omt) / (enu^2)) - htphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*omt) / (enu^2) + ((htphi*omt) / (enu^2) + hphiphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*omt*p) + im * ((((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((W*(eps + p)*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (-W*((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(-(enu^2) + (epsi^2)*(omt^2))*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_stress_dJphi_t(inp)
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
    return ((((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2))*((W*(eps + p)*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(eps + p)*(epsi^2)*omt) / enu)) / ((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu) + (-W*(eps + p)*(epsi^2)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + hphiphi*((-omt*p) / (enu^2) + ((W^2)*(eps + p)*(omt + (enu*v) / epsi)) / (enu^2)) + htphi*((-p) / (enu^2) + ((W^2)*(eps + p)) / (enu^2)) - y1*(((-((W*(omt + (enu*v) / epsi)*p_r) / enu) - ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*((W*(eps + p)*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(eps + p)*(epsi^2)*omt) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(epsi^2)) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) - (epsi^2)*omt*p*(((-htt) / (enu^2) + (-htphi*omt) / (enu^2)) / (enu^2) + (-omt*(htphi / (enu^2) + (hphiphi*omt) / (enu^2))) / (enu^2)) + (epsi^2)*p*((-omt*(htt / (enu^2) + (htphi*omt) / (enu^2))) / (enu^2) + (htphi / (enu^2) + (hphiphi*omt) / (enu^2))*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) + (eps + p)*((eps*(((W^2)*(epsi^2)*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*(epsi^2)*omt) / (enu^2))) / p + ((m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) - (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*((W*(eps + p)*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(eps + p)*(epsi^2)*omt) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (W*(epsi^2)*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (epsi^2)*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)) - (epsi^2)*omt*((W^2) / (enu^2) + -1 / (enu^2)))*y0) + im * (((-((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2)) - (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*((W*(eps + p)*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(eps + p)*(epsi^2)*omt) / enu)) / ((eps + p)*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (W*((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(epsi^2)) / (enu*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end

function _eccentric_li_stress_dJphi_r(inp)
    W = inp.W
    emu = inp.emu
    enu = inp.enu
    eps = inp.eps
    epsi = inp.epsi
    hrphi = inp.hrphi
    htr = inp.htr
    omt = inp.omt
    p = inp.p
    v = inp.v
    y1 = inp.y1
    return ((hrphi*p) / (emu^2) + (epsi^2)*p*((-hrphi*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) / (emu^2) + (htr*omt) / ((emu^2)*(enu^2))) - (epsi^2)*omt*p*(htr / ((emu^2)*(enu^2)) + (hrphi*omt) / ((emu^2)*(enu^2)))) + im * (y1*((W*(eps + p)*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-W*(eps + p)*(epsi^2)*omt) / enu))
end

function _eccentric_li_stress_dJphi_phi(inp)
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
    return (((-(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2)) - (((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2)) - (((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2)))*((-W*(eps + p)*(epsi^2)*omt) / enu + (2W*(eps + p)*(epsi^2)*(omt + (enu*v) / epsi)) / enu)) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-W*(eps + p)*(epsi^2)*omt*(omt + (enu*v) / epsi)*((((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)*(((1//2)*(W^2)*((htt*(omega_d + m*(omt + (enu*v) / epsi))*omt) / (enu^2) + (-htt*m - 2htphi*(omega_d + m*(omt + (enu*v) / epsi)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (enu^2) + ((1//2)*(W^2)*(((-2htphi*m - hphiphi*(omega_d + m*(omt + (enu*v) / epsi)))*omt) / (enu^2) + hphiphi*m*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)) / (enu^2) + ((W^2)*((-htt*m*omt) / (enu^2) - hphiphi*(omega_d + m*(omt + (enu*v) / epsi))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (enu^2))) / ((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu) + (-(W^2)*htt) / (enu^2) + (-2(W^2)*htphi*(omt + (enu*v) / epsi)) / (enu^2) + (-(W^2)*hphiphi*((omt + (enu*v) / epsi)^2)) / (enu^2))) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)) + htphi*((-omt*p) / (enu^2) + ((W^2)*(eps + p)*(omt + (enu*v) / epsi)) / (enu^2)) - y1*((((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*((-W*(eps + p)*(epsi^2)*omt) / enu + (2W*(eps + p)*(epsi^2)*(omt + (enu*v) / epsi)) / enu)) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (W*((W*(omt + (enu*v) / epsi)*p_r) / enu + ((W*omt_r) / enu + (W*v_r) / epsi + (-W*psi_r*v) / epsi + ((W^3)*(omt + (enu*v) / epsi)*v*v_r) / enu + (W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / (2enu) + ((1//2)*W*(((-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*omt) / (enu^2) + (-(epsi^2)*omt_r - 2(epsi^2)*omt*psi_r)*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))) / enu + (W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / (2enu) + ((1//2)*W*((((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*omt) / (enu^2) + 2(epsi^2)*psi_r*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*(omt + (enu*v) / epsi)) / enu + enu*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2))*((-W*(omt + (enu*v) / epsi)) / (enu^2) + (W*v) / (enu*epsi)))*(eps + p))*(epsi^2)*omt*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu))) + hphiphi*(((W^2)*(eps + p)*((omt + (enu*v) / epsi)^2)) / (enu^2) + p*((-(omt^2)) / (enu^2) + 1 / (epsi^2))) - (epsi^2)*omt*p*(((-htt*omt) / (enu^2) - htphi*((omt^2) / (enu^2) + -1 / (epsi^2))) / (enu^2) + ((-((htphi*omt) / (enu^2)) - hphiphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*omt) / (enu^2)) + (eps + p)*((eps*(((W^2)*(epsi^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + (-(W^2)*(epsi^2)*omt*(omt + (enu*v) / epsi)) / (enu^2))) / p + ((-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*((-W*(eps + p)*(epsi^2)*omt) / enu + (2W*(eps + p)*(epsi^2)*(omt + (enu*v) / epsi)) / enu)) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (W*(epsi^2)*(-m*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) + (omega_d + m*(omt + (enu*v) / epsi))*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*omt*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (epsi^2)*((-(omt^2)) / (enu^2) + ((W^2)*((omt + (enu*v) / epsi)^2)) / (enu^2) + 1 / (epsi^2)) - (epsi^2)*omt*((-omt) / (enu^2) + ((W^2)*(omt + (enu*v) / epsi)) / (enu^2)))*y0 + (epsi^2)*(((-((htt*omt) / (enu^2)) - htphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*omt) / (enu^2) + ((htphi*omt) / (enu^2) + hphiphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((-(omt^2)) / (enu^2) + 1 / (epsi^2)))*p) + im * ((((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*((-W*(eps + p)*(epsi^2)*omt) / enu + (2W*(eps + p)*(epsi^2)*(omt + (enu*v) / epsi)) / enu)) / ((eps + p)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)) + (W*((((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*p_r) / (emu^2) + (eps + p)*(((W^2)*(-2(epsi^2)*omt*omt_r + 2(enu^2)*((-p_r) / ((W^2)*(eps + p)) + (-epsi*omt_r*v) / enu + psi_r*(v^2)) - 2(epsi^2)*(omt^2)*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))) / (2(emu^2)*(enu^2)) + (-(W^2)*(epsi^2)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*((omt + (enu*v) / epsi)^2)*psi_r) / ((emu^2)*(enu^2)) + ((W^2)*((epsi^2)*omt_r + 2(epsi^2)*omt*psi_r)*((htr*omt) / (enu^2) + hrphi*((omt^2) / (enu^2) + -1 / (epsi^2)))*(omt + (enu*v) / epsi)) / ((emu^2)*(enu^2))))*(epsi^2)*omt*(omt + (enu*v) / epsi)*((2W*(epsi^2)*(omt + (enu*v) / epsi)) / enu + (-2W*(epsi^2)*omt) / enu)) / (enu*((2W*(-(enu^2) + (epsi^2)*(omt^2))) / enu + (-2W*(epsi^2)*omt*(omt + (enu*v) / epsi)) / enu)*((-W*(omega_d + m*(omt + (enu*v) / epsi))) / enu + (W*m*(omt + (enu*v) / epsi)) / enu)))
end
