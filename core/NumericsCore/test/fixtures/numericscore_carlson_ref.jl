# AUTO-GENERATED one-shot fixture — DO NOT EDIT.
# Baked KerrGeodesics complex-Carlson / complete+incomplete-Legendre / Jacobi oracle
# outputs on the EXACT argument grid of test_carlson.jl. Captured once (legacy still
# resolvable) by scratchpad/gen_fixtures.jl:
#   KG._RC/_RF/_RD/_RJ, KG._ellipK/_ellipE/_ellipPi, KG._ellipFinc/_ellipEinc/
#   _ellipPiinc, KG._phi_reduce, KG._jacobiAM/_jacobiSN/_jacobiCN/_jacobiDN
# (packages/KerrGeodesics). Preserves the Phase-7-B characterization at the SAME rtol=1e-12.

const REF_RC_cxy   = 0.7070125013945383 - 0.020205824601440135im
const REF_RF_cxyz  = 0.6552033345708749 - 0.08715220718497983im
const REF_RD_cxyz  = 0.22519915301402793 - 0.07891058770501491im
const REF_RJ_cxyzp = 0.45460250827019827 - 0.2290146623111023im

const REF_RC_23   = 0.6154797086703871
const REF_RF_123  = 0.7269459354689083
const REF_RD_123  = 0.29046028102899063
const REF_RJ_1234 = 0.23984809974956778

const REF_ellipK = Dict(0.001 => 1.5711892469233446 - 0.0im, 0.3 => 1.7138894481787907 - 0.0im, 0.7 => 2.075363135292469 + 0.0im, 0.999 => 4.841132560550298 + 0.0im, 2.4972 => 1.125028950544552 - 1.2333995228464811im)
const REF_ellipE = Dict(0.001 => 1.5704035540514238 - 0.0im, 0.3 => 1.445363064412665 - 0.0im, 0.7 => 1.2416705679458229 + 0.0im, 0.999 => 1.0021707908344464 + 0.0im, 2.4972 => 0.5266654431524656 + 0.8188299688332643im)
const REF_ellipPi = Dict((0.2, 0.3) => 1.925728340271729 + 0.0im, (0.6, 0.7) => 3.5111846140642324 + 0.0im)
const REF_ellipPi_nC = 2.420252487792902 + 0.4962240697902826im

const REF_ellipK_mod = Dict(0.3 => 1.608048619930513, 0.6 => 1.7507538029157528, 0.9 => 2.2805491384227703)
const REF_ellipE_mod = Dict(0.3 => 1.5348334649232493, 0.6 => 1.4180833944487246, 0.9 => 1.171697052781614)

const REF_ellipFinc = Dict(0.7 => 0.7165177159853932, 4.0 => 4.315943388399585, -1.3 => -1.3918738208303512)
const REF_ellipEinc = Dict(0.7 => 0.6841406078067003, 4.0 => 3.721065268268602, -1.3 => -1.2176345070733994)
const REF_ellipPiinc = Dict(0.7 => 0.7392911132848218, 4.0 => 4.781240367044864, -1.3 => -1.5255973416148914)
const REF_phi_reduce_4 = (1.0, 0.8584073464102069)

const REF_jacobiAM = Dict((0.8, 0.5) => 0.7627814197609588, (1.7, 0.2) => 1.6069120400784802, (2.3, 0.9) => 1.4818337945150484)
const REF_jacobiSN = Dict((0.8, 0.5) => 0.6909348508664388, (1.7, 0.2) => 0.999347898512061, (2.3, 0.9) => 0.9960454430973229)
const REF_jacobiCN = Dict((0.8, 0.5) => 0.7229170297192977, (1.7, 0.2) => -0.03610786257212465, (2.3, 0.9) => 0.08884523220217044)
const REF_jacobiDN = Dict((0.8, 0.5) => 0.8725276591198047, (1.7, 0.2) => 0.8945729459065401, (2.3, 0.9) => 0.32726766989201955)
