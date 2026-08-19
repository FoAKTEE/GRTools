# Provenance — `Wigner/test` data assets

| field | value |
|---|---|
| immediate source | predecessor package `Punctures`, `test/wignerD_reference.csv` in the `emri-disk` monorepo (byte-identical, `cmp`-verified) |
| capture commit | `2fb2c909d821bb717a6087aafe43061d12d6e00f` (emri-disk, `packages/Punctures/test`, 2026-06-28) |
| generator | `packages/Punctures/oracle/gen_wignerD.wl` — `wolframscript -file`, 30 digits, closed form |
| upstream physics | BlackHolePerturbationToolkit / **Punctures**, <https://github.com/BlackHolePerturbationToolkit/Punctures.git> @ `fd4ee50003cd571e0629c832e204cefbb565a706` (2025-09-16) |
| recorded | 2026-08-19 |

Content: 287 rows of `l, m, mp, Re|Im` at the FIXED Euler angles (pi, pi/2, pi/2) — the oracle for
`Wigner.wignerD_fixed` (`test/test_wignerD.jl`). The identical file is also vendored at
`capabilities/SelfForce/test/fixtures/puncture/wignerD_reference.csv` (same sha256): the two
packages are firewalled, neither may reach into the other's test tree.

| file | bytes | sha256 |
|---|---:|---|
| `wignerD_reference.csv` | 9320 | `3548c69e0e788c2c5cae6bb2c3430df7f0bcd6e1c1be7b61235dc4eb99e5ff20` |
