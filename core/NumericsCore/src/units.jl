# Split-named physical constants.
#
# The two legacy forks DISAGREE on more than naming: KerrBackwards carries the canonical LAL/global.h
# #define values, while FastEMRIWaveforms re-derives them from lisaconstants 2.0.2 (a slightly
# different GM_SUN), so MTSUN/MRSUN/YRSID differ at ~1e-7. We take the KerrBackwards (LAL) values as
# canonical. Separately, the symbol `Gpc` collides between packages (meters in FEW, but seconds in the
# C `GPCINSEC`); we resolve that by SPLITTING into distinct named axes Gpc_in_m / Gpc_in_s.

"Mass-to-time conversion G M_sun / c^3 in seconds (LAL/global.h canonical value)."
const MTSUN_SI = 4.925491025873693e-06

"Sidereal year in seconds (LAL/global.h)."
const YRSID_SI = 31558149.763545603

"Solar-mass Schwarzschild radius G M_sun / c^2 in metres = MTSUN_SI*c (LAL/global.h)."
const MRSUN_SI = 1476.6250615036158

"Astronomical unit as a light-travel time in seconds."
const AUsec = 499.004783836156412

"One gigaparsec in metres (the `Gpc` meters axis)."
const Gpc_in_m = 3.0856775814913674e25

"One gigaparsec as a light-travel time in seconds (the `Gpc` seconds axis; C `GPCINSEC`)."
const Gpc_in_s = 1.029271251e17

"Legacy alias for the seconds-axis gigaparsec (KerrBackwards `GPCINSEC`)."
const GPCINSEC = Gpc_in_s

"π to upstream literal precision (Float64 round-trip)."
const PI = 3.141592653589793238462643383279502884

"π/2."
const PI_2 = 1.570796326794896619231321691639751442

"π/4."
const PI_4 = 0.785398163397448309615660845819875721
