# runtime_config fixture provenance

`ScalarSelfForce.yaml` — BYTE-EXACT copy (sha256
`480aba9f8d6a7f14541b136cd2d7ba07d0a3d7cb022eba866978972f6a61b45d`, 6562 bytes) of the
READ-ONLY mirror `ref-paper/arxiv-2606.04998/src/anc/ScalarSelfForce.yaml` — the CRP
production spectre input file behind the figures of arXiv:2606.04998 (Vu et al.). It is a
Jinja2 template; the production example parameters are in its header (lines 9-12): spin=0,
orbital_radius=10, m_mode=2, inner_radius=2, outer_radius=1e4, worldtube_radius=4.8,
worldtube_dz2=0.25, P=4, L=0.

DO NOT EDIT: `test/test_runtime_config.jl` pins this file's sha256 (and cross-checks the
ref-paper mirror when the repo checkout is present). Any change here or in the mirror MUST
force a review of the `crp_paper_config()` transcription in `src/Domain/config.jl` and a
re-record of the hash.
