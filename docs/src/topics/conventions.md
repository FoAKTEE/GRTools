# [Conventions](@id topic-conventions)

The cross-cutting convention surface of GrTools: the unified negative-``m`` parity flag,
the dimensionless MST/EMRI convention helpers, and the split-named physical constants.
This page is the **single home** of the topic — package pages
([Harmonics](../capabilities/harmonics.md), [NumericsCore](../core/numericscore.md),
[EMRIWaveforms](../capabilities/emriwaveforms.md)) link here rather than restating it.

Everything below lives in exactly two places:

| Home | Symbols |
|---|---|
| `NumericsCore` (`conventions.jl`, `units.jl`, `legendre.jl`) | `p_to_y`, `eps_to_omega`/`omega_to_eps`, `q_to_a`/`a_to_q`, `imagunit`, `MTSUN_SI`, `YRSID_SI`, `MRSUN_SI`, `AUsec`, `Gpc_in_m`/`Gpc_in_s`/`GPCINSEC`, `PI`/`PI_2`/`PI_4`, and the bare seed `spin_weighted_Y` |
| `Harmonics` (`SphericalHarmonics.jl`) | `spin_weighted_spherical_Y`, `angular_parity`, `SpinWeightedSphericalHarmonicY`, and the legacy shims `ylm`/`ylm_few`/`ylms_few` |

## [The `m<0` parity surface](@id conventions-parity)

The legacy packages disagreed on how to *package* the negative-``m`` spin-weighted
spherical harmonic: FEW's `GetYlms` and KerrBackwards' `ylm_few` apply a ``(-1)^\ell``
operator, the scalar/real-``Y`` users (DiscModes, Punctures) build ``m<0`` from the
reality relation, and BHPT-style code uses the bare seed. GrTools collapses these forks
into **one parity-flagged front door**:

```julia
spin_weighted_spherical_Y(s, l, m, θ, φ; convention=:none)   # Harmonics
angular_parity(l, m; convention=:none)                        # the sign helper
```

For ``m \ge 0`` every convention is a no-op — they differ **only** in the ``m<0`` packaging:

| `convention` | ``m<0`` rule | Legacy fork it replaces |
|---|---|---|
| `:none` (default) | bare seed ``{}_sY_{\ell m}`` — the **BHPT** convention, bit-identical to the legacy `SpinWeightedSphericalHarmonicY` | BHPT `SpinWeightedSpheroidalHarmonics`, FEW `ylm` |
| `:few` | ``(-1)^\ell\, {}_sY_{\ell m}`` | FEW `GetYlms`, KerrBackwards `ylm_few`/`ylms_few` |
| `:reality` | built via the **s-flip rule** below (value-equivalent, indeed bit-identical, to the bare seed) | DiscModes real-``Y``, Punctures ``(-1)^m\,\mathrm{conj}`` |

`angular_parity(l, m; convention)` returns only the multiplicative *sign*: `+1` for
``m\ge 0`` or `:none`, ``(-1)^\ell`` for `:few`, and the ``m``-part ``(-1)^m`` for
`:reality`. It takes `(l, m)` because ``(-1)^\ell`` genuinely needs the degree. An unknown
convention symbol throws an `ArgumentError`.

### [The s-flip rule](@id conventions-s-flip)

`:reality` constructs ``m<0`` from the general spin-weighted reality relation

```math
{}_sY_{\ell,-|m|} \;=\; (-1)^{s+|m|}\,\overline{{}_{-s}Y_{\ell,|m|}}\,,
```

i.e. flip the spin weight ``s \to -s``, evaluate the *positive*-``m`` seed, conjugate, and
apply the full phase ``(-1)^{s+m}``. In code, `angular_parity(...; convention=:reality)`
supplies the ``(-1)^m`` factor and the front door adds the ``(-1)^s`` factor and the
``s\to-s`` conjugated seed. For ``s=0`` this reduces to the scalar
``{}_0Y_{\ell,-|m|} = (-1)^{|m|}\,\overline{{}_0Y_{\ell,|m|}}``.

### Named entry points and shims

- `SpinWeightedSphericalHarmonicY(s, l, m, θ, φ; parity=:none)` — the BHPT-named alias;
  `parity` routes straight to `convention`.
- `ylm(s, l, m, θ, φ)` and the 2-argument `ylm(l, m, θ, φ)` (fixes ``s=-2``, the
  gravitational FEW default) — the bare-seed shim.
- `ylm_few(l, m, θ, φ)` / vectorised `ylms_few(l_arr, m_arr, θ, φ)` — the
  ``s=-2``, `convention=:few` shims (FEW `few.utils.ylm.GetYlms`).

All of these evaluate the single canonical seed `NumericsCore.spin_weighted_Y(s,l,m,θ,φ)`
(half-angle binomial sum; exact-integer BigFloat path at high ``\ell`` or BigFloat
colatitude, so precision is inherited — see [Precision](precision.md)). The real
scalar-``Y`` path `Ylm_real(l, m, θ)` / `Ylm_eq(l, m, θ)` also seeds in `NumericsCore` and
is re-exported by `Harmonics`.

## [ε↔ω and q↔a convention helpers](@id conventions-eps-omega)

`NumericsCore.conventions` is the single source of truth for the dimensionless
MST/EMRI variable conventions (the MST notebook's ``\epsilon = 2M\omega`` and
``q = a/M``). All helpers are ``M``-parametrized with default ``M=1``, so callers in
either convention agree:

| Helper | Map |
|---|---|
| `omega_to_eps(omega, M=1)` | ``\omega \mapsto \epsilon = 2M\omega`` |
| `eps_to_omega(eps, M=1)` | ``\epsilon \mapsto \omega = \epsilon/(2M)`` |
| `q_to_a(q, M=1)` | dimensionless spin ``q = a/M \mapsto`` Boyer–Lindquist ``a = qM`` |
| `a_to_q(a, M=1)` | ``a \mapsto q = a/M`` |
| `imagunit()` / `imagunit(T)` | the numeric imaginary unit; `imagunit(Float64) == ComplexF64(0,1)` etc. for precision control |

!!! note "The numeric and the symbolic imaginary unit"
    `NumericsCore.imagunit` is *numeric*. `PostNewtonian.MST` exports the symbolic
    `imagunit_sym()` — the Symbolics atom `𝕚` carried inside PN coefficient polynomials —
    whose substitution `imagunit_sym() => im` evaluates to the numeric one. They are
    deliberately distinct symbols in distinct packages, and (since the export-collision
    cleanup) distinct NAMES, so `using NumericsCore, PostNewtonian.MST` is unambiguous.

## [The `p_to_y` map](@id conventions-p-to-y)

```julia
p_to_y(p, e) = log(p - 2e - 21//10)       # NumericsCore
```

maps the semilatus rectum ``p`` to the amplitude-grid coordinate
``y = \log(p - 2e - 21/10)``; the shift keeps ``y`` finite and clear of the Schwarzschild
separatrix ``p_s = 6 + 2e``. It is the canonical home of FEW's
`schwarzecc_p_to_y` and KerrBackwards' `p_to_y` (deduplicated to one definition): the
`EMRIWaveforms` amplitude interpolants (`AmplitudeModels.jl`) evaluate on the ``(y, e)``
grid through this one definition, and `EMRIWaveforms.Mappings` re-exports it with no new
body.

## [Units: the split-named constants](@id conventions-units)

The legacy forks disagreed on more than naming: KerrBackwards carries the canonical
**LAL** (`global.h`) values while FastEMRIWaveforms re-derives from `lisaconstants 2.0.2`
(a slightly different ``GM_\odot``), so `MTSUN`/`MRSUN`/`YRSID` differed at ``\sim 10^{-7}``
relative. GrTools takes the **LAL values as canonical**. Separately, the bare symbol
`Gpc` *collided* between packages — metres in FEW but light-seconds in the C `GPCINSEC` —
so it is **split into two distinct named axes**:

| Constant (`NumericsCore`) | Value | Meaning |
|---|---|---|
| `MTSUN_SI` | `4.925491025873693e-6` | ``GM_\odot/c^3`` in seconds (LAL) |
| `MRSUN_SI` | `1476.6250615036158` | ``GM_\odot/c^2`` in metres (`= MTSUN_SI·c`) |
| `YRSID_SI` | `31558149.763545603` | sidereal year in seconds (LAL) |
| `AUsec` | `499.004783836156412` | astronomical unit as light-travel time, s |
| `Gpc_in_m` | `3.0856775814913674e25` | one Gpc in **metres** (the FEW axis) |
| `Gpc_in_s` | `1.029271251e17` | one Gpc as **light-travel time in seconds** (the C `GPCINSEC` axis) |
| `GPCINSEC` | `= Gpc_in_s` | legacy KerrBackwards alias |
| `PI`, `PI_2`, `PI_4` | — | ``\pi, \pi/2, \pi/4`` to upstream literal precision |

Pick the axis by what the formula divides by:

- **metres axis** — `EMRIWaveforms.Generate` forms the dimensionless luminosity distance
  as `(dist * Gpc_in_m) / (mu * MRSUN_SI)` (both metres);
- **seconds axis** — the AAK kernel (`Kernels/src/aak.jl`) forms the amplitude scale as
  `zeta = (mu * MTSUN_SI) / dist / GPCINSEC` (both seconds);
- **time conversions** — `EMRIWaveforms.Trajectory` converts observation years to
  dimensionless time via `T * YRSID_SI / (Mtot * MTSUN_SI)`.

## [Hyperboloidal charts in GrTools](@id conventions-hyperboloidal)

The time-domain solvers use **three distinct hyperboloidal charts**. "Hyperboloidal" is an
*umbrella property*, not one chart: the constant-time slices are spacelike yet reach future null
infinity ``\mathcal{I}^+``, so the boost ``H \equiv (dh/dr)/(dr_*/dr) \to 1`` as ``r\to\infty`` (the
outgoing coordinate speed stays finite). All three below share that property but fix **different**
(foliation, radial coordinate, field rescaling) triples — so they are genuinely different charts and,
for RT vs HPHC, genuinely different **foliations**, not relabelled radii. Each writes the new time as
``T = t - h(r)`` (Boyer–Lindquist ``t``, ``M=1``):

| chart | radial map | ``\mathcal{I}^+`` | horizon | offset ``h(r)`` | rescaling | source | home |
|---|---|---|---|---|---|---|---|
| **RT** | ``r=2R/(1-R^2)`` | ``R=1`` | penetrating | ``G_{RT}(R)+r-r_*`` | ``\Psi=\Delta^{-s}r^{-1}\psi`` | arXiv:1301.1591 | `TimeDomain/rtchart.jl` |
| **HPHC** | ``R=c_\ell^2/r`` | ``R=0`` | penetrating | ``2r+4M\ln r-r_*`` | peel ``\times R^{-\text{falloff}}`` | arXiv:2010.00162 §V.C | `NPGHP/hphc.jl` |
| **Layer** | ``r_*=\rho/\Omega`` | ``\rho=s`` | **no** (Cauchy interior) | ``r_*-\rho`` | ``1/(1+H)`` combos | arXiv:2307.01349 §II.D | `TimeDomain/hyperboloidal.jl` |

with ``G_{RT}(R)=-4M\log|1-R^2|+(1+R^2)/(1-R^2)``, tortoise ``dr_*/dr=(r^2+a^2)/\Delta``, and the layer
function ``H=1-\Omega^2/(\Omega-\rho\Omega')`` (= ``1-d\rho/dr_*``, so ``h=r_*-\rho``).

### [RT vs HPHC are different foliations (F0)](@id conventions-hyperboloidal-F0)

RT (its Eq. (3)) and HPHC (its ``dv=dt+dr_*-dr``) use the **same** intermediate ingoing-Kerr time
``\tau\equiv v=t-r+r_*``, then bend to ``\mathcal{I}^+`` with **different** heights. Two charts share a
foliation iff ``h_1-h_2=\text{const}``, i.e. ``dh_1/dr\equiv dh_2/dr``. Numerically,

```math
\frac{dh_{\rm HPHC}}{dr}-\frac{dh_{\rm RT}}{dr}=1+\frac{4M}{r}-\frac{dG_{RT}}{dr}\neq 0 ,
```

running from ``+0.92`` (``a=0``) / ``+1.19`` (``a=0.7``) at the horizon to ``+1.1\times10^{-4}`` at
``r=200``, vanishing only as ``r\to\infty``. Since the common ``r_*`` cancels, the difference is
**``a``-independent** (identical for ``a=0`` and ``a=0.7`` to ``9\times10^{-16}``). **Verdict: RT and
HPHC are different foliations** — the [HYPOTHESIS] is confirmed. The layer chart differs from both: it
is Cauchy (``H=0``) in the strong field and so is not even horizon-penetrating. The verdict rests on
tabulated ``dh/dr`` differences, with the analytic derivatives checked against finite differences.

Coefficients and maps live in the module docstrings: `rt_r_of_R`/`rt_F_of_R` (RT), `hphc_R_of_r`/
`hphc_horizon_R` (HPHC), and the `Layer` struct with `Hfun`/`rstar_of_ρ` (layer).
