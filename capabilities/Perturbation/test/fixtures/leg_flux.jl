# BAKED test-oracle fixture, captured ONCE from the LIVE legacy Teukolsky.energy_flux / angular_momentum_flux.
# Frozen so Perturbation tests need no packages/* at test time (Phase 7-B).
# DO NOT hand-edit; regenerate via scratchpad gen.jl against the same legacy commit.
const REF_FLUX_LEG = Dict(
  (-2, 2, 2, 0.0, 0.1) => (eI=12.175353146529993, eH=1.3619605978551405e-5, jI=243.50706293059986, jH=0.0002723921195710281),
  (-2, 2, 1, 0.5, 0.15) => (eI=6.012520072360492, eH=9.496837249519675e-8, jI=40.083467149069946, jH=6.331224833013117e-7),
  (-1, 2, 2, 0.6, 0.2) => (eI=0.1543802947991385, eH=-3.953191063405665e-6, jI=1.543802947991385, jH=-3.953191063405665e-5),
  (0, 2, 2, 0.5, 0.2) => (eI=0.01273239544735163, eH=-4.0360076312958e-5, jI=0.1273239544735163, jH=-0.00040360076312958),
  (1, 2, 1, 0.9, 0.2) => (eI=1.0078324643245411e-5, eH=-1.7838973446488136e-6, jI=5.039162321622705e-5, jH=-8.919486723244068e-6),
  (2, 3, 3, 0.7, 0.25) => (eI=2.1587866875331188e-9, eH=-4.8683223823463415e-9, jI=2.5905440250397425e-8, jH=-5.8419868588156095e-8),
)
