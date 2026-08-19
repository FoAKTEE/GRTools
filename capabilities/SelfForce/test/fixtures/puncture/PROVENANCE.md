# Provenance — `SelfForce/test/fixtures/puncture`

Captured `wolframscript` oracle output for the puncture / regularized-sum / worldline-self-field /
mode-mixing gates. Byte-identical to the predecessor capture (verified with `cmp`).

| field | value |
|---|---|
| immediate source | predecessor package `Punctures`, directory `test/` in the `emri-disk` monorepo |
| capture commit | `2fb2c909d821bb717a6087aafe43061d12d6e00f` (emri-disk, `packages/Punctures/test`, 2026-06-28) |
| generators | `packages/Punctures/oracle/gen_*.wl` — Mathematica, `wolframscript -file`, 30 digits, deterministic |
| upstream physics | BlackHolePerturbationToolkit / **Punctures**, <https://github.com/BlackHolePerturbationToolkit/Punctures.git> @ `fd4ee50003cd571e0629c832e204cefbb565a706` (2025-09-16) |
| recorded | 2026-08-19 |

The generators are deterministic closed-form Mathematica evaluations (e.g. the fixed-angle
`WignerD` table at 30 digits for `wignerD_reference.csv`), so re-running one reproduces its CSV
byte-for-byte under the same Mathematica version.

| file | bytes | sha256 | generator (`oracle/`) |
|---|---:|---|---|
| `D2hR1_reference.csv` | 10220 | `390ae91ccfbe5580b2531fc1615710fc17f088c7d75e9f04cca23028a66e5746` | `gen_worldline.wl` |
| `DhR1_reference.csv` | 2328 | `847b6a3d84555390a175fffee7b47d30f3d8f0a9d2e6af03e3bd03cdf96ec89a` | `gen_worldline.wl` |
| `alphaint_reference.csv` | 17095 | `c45db2a917491be43c7dcb2676e4465f391c6c5cb4bb6f35f88ba1613916ac63` | `gen_alphaint.wl` |
| `coord_reference.csv` | 26311 | `1c13bb373e1fc122cf4dd5a059761e252c700196951a3ffae2edd3617ce347e0` | `gen_tensor_basis.wl` |
| `ddhret_reference.csv` | 4743 | `ca5640aa2ccfe1e17223b78dd2cbd51d208c6007b52bca5f9470aae876b99cde` | `gen_ddhret.wl` |
| `dhP_reference.csv` | 767 | `be1f583ef9a43d7d88c9943652356bc808f71e2a2de20c211ad9284cad6e69ee` | `gen_dhP.wl` |
| `dhR1Data_reference.csv` | 774 | `24a169a51856ee855b02b1bfaff0d7ebf4283133b97bf3ba2bf0f8662db5a8ff` | `gen_worldline_assembly.wl` |
| `hP_reference.csv` | 2532 | `2d161c6c801280e050c370eebce82a462fa7722640f6fcc338ae1782c8c589d8` | `gen_hP_reference.wl / gen_hP_tail.wl` |
| `hR1Data_reference.csv` | 196 | `f3ea9ba8c3d48da762c09869ebb1a124ad3d99fed44ae6ed8306907fd05d2ddb` | `gen_worldline_assembly.wl` |
| `h_reference.csv` | 1057 | `a6b4b494da5fc1e79ab62e18fa777b5f40ee31ecec167ac28d886c9002a30a30` | `gen_h_reference.wl` |
| `regsum_reference.csv` | 1099 | `0f49e6f495ec32ec3b0c84d88bee20b8279c4c2f0f7b44105faac743294fef86` | `gen_regsum.wl` |
| `tensor_basis_reference.csv` | 59741 | `a1f9ce4939acc72280d1c09e9f53082e75f6fecc6c58235981ccacab6d35be7f` | `gen_tensor_basis.wl` |
| `wignerD_reference.csv` | 9320 | `3548c69e0e788c2c5cae6bb2c3430df7f0bcd6e1c1be7b61235dc4eb99e5ff20` | `gen_wignerD.wl` |
