# Provenance — `ModifiedGravity/test/data`

Frozen source-mode tables (12-column tab-separated `Sourcefiletotable`/`mywrite` output) used as the
ORACLE by `test/test_source_tables.jl`. Vendored VERBATIM (all 3 files byte-identical to upstream;
verified with `cmp` against `emri-disk/ref-code/EMRI-modified-Teukolsky-equation`).

| field | value |
|---|---|
| source | Lukas Polcar / **EMRI-modified-Teukolsky-equation** (Polcar & Witzany 2025, arXiv:2507.15720) |
| URL | <https://github.com/LukasPolcar/EMRI-modified-Teukolsky-equation.git> |
| commit | `652d6c8c7b9c6c2f75bcc3f2c247c714bdc19628` |
| upstream date | 2026-05-16 |
| vendored path | repository root (`*.dat`), unmodified |
| local mirror | `emri-disk/ref-code/EMRI-modified-Teukolsky-equation` |
| recorded | 2026-08-19 |

File naming: `i<inclination>_p<semi-latus-rectum>_e<eccentricity>_r<ring-radius>.dat`. Columns:
integers `m, n, k`, the mode frequency, then the four projected source amplitudes as interleaved
Re/Im pairs (`Data.read_source_table`); the ring radius and ring mass are caller-supplied and are
NOT stored in the file.

| file | bytes | sha256 |
|---|---:|---|
| `i0.5_p12.0_e0_r20.0.dat` | 92863 | `1f98ae2c06d1adb69b21081d23120782048600641d9d28499bbd0623195f38e9` |
| `i0.5_p12.0_e0_r400.0.dat` | 32948 | `5d48c4583a8b4675331db34a4bb3c24dc2977645a942f0bf6ca7b2313db5b34b` |
| `i1.0_p12.0_e0_r20.0.dat` | 1712 | `57a2fef8e7d1527079a70bfb125fef74b801c7b96409d329d69209b3a91e8156` |
