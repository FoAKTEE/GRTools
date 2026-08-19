# Provenance — `PostNewtonian/SeriesData`

Upstream PN self-force series catalog, vendored VERBATIM (all 94 files byte-identical to upstream;
verified with `cmp` against the mirror at `emri-disk/ref-code/PostNewtonianSelfForce/SeriesData`).

| field | value |
|---|---|
| source | BlackHolePerturbationToolkit / **PostNewtonianSelfForce** |
| URL | <https://github.com/BlackHolePerturbationToolkit/PostNewtonianSelfForce.git> |
| commit | `a9d83f1b8b775ac1d6d90741ecabd9903d9195b5` ("Fix glitch in documentation guide page") |
| upstream date | 2025-07-16 |
| vendored path | `SeriesData/` (whole subtree, unmodified) |
| local mirror | `emri-disk/ref-code/PostNewtonianSelfForce` |
| recorded | 2026-08-19 |
| license | upstream `LICENSE` (MIT); upstream `CITATION.cff` names the authors to cite |

Consumers: `PostNewtonian.Catalog` (`.m` discovery + `read_mma_inputform`) and
`PostNewtonian.SelfForce` (`postnewtonian_coefficient`, `ResummedSeriesData`).

| file | bytes | sha256 |
|---|---:|---|
| `Kerr/Circular/Flux/Energy/Horizon.m` | 383244 | `c4034c6d97352b9200c88acd109d348b33ba781d40adf59b36d5522113ccfc9b` |
| `Kerr/Circular/Flux/Energy/Infinity.m` | 127851 | `ef9197c396f46b3bafd4cc16db32fdfdc0463981915f842d11a3570adcd6b61e` |
| `Kerr/Circular/Local/RedShift.m` | 5213 | `0fb154af529c591d50c54722b506cd291883c1ea5771a502c407da0f814781a3` |
| `Kerr/Circular/Local/Spin.m` | 3173 | `c68359584519d897dd8c877f8f271c2752d09eb1a380ac1183b33cee01d486a5` |
| `Kerr/Circular/Memory/h-l10m0.m` | 900 | `1502c9f6260f50209a62a4d9c9a2f56f645f220993099adf02feff0143cf7cca` |
| `Kerr/Circular/Memory/h-l12m0.m` | 769 | `b9c5db0672fdc760d7d5107c34e04f0e42cf43c87d2028b625dc1fa3214d32d3` |
| `Kerr/Circular/Memory/h-l14m0.m` | 724 | `d512fa0419a337b388b2874463842308fea72d7bd24fc351c52ecc41c7ddcf0c` |
| `Kerr/Circular/Memory/h-l2m0.m` | 3609 | `56ab30281e0595d8da6febfaf8303621f643345799d29e049a70641a87abecd8` |
| `Kerr/Circular/Memory/h-l4m0.m` | 3390 | `01cea07748c82c1d63ff2f0af1f719cfe247e4ce7b0702f4b5d3baba687a9bf0` |
| `Kerr/Circular/Memory/h-l6m0.m` | 4718 | `eff1740acacd5790479fdb14bb8ad2b2e78df83eaf6fe564d223c959cf0fbb36` |
| `Kerr/Circular/Memory/h-l8m0.m` | 6826 | `8a7eee18724f788c8c2e030698eb42829bc373cb50046e0e5823404e2a0a1ea5` |
| `Kerr/Circular/Memory/hdot-l10m0.m` | 932 | `3c8e659137515d686fb9a3d25f2e4193b4ad94a6b90597737644e61b32d8d565` |
| `Kerr/Circular/Memory/hdot-l12m0.m` | 781 | `c8dc40c4223e94abf773b13abdbc7793899ab2739e4d28b74505654fd280407e` |
| `Kerr/Circular/Memory/hdot-l14m0.m` | 740 | `fc19b28da16cf070892c28f0de70bd9f87794c0aac2c1df6ab75ccca9e2cbb68` |
| `Kerr/Circular/Memory/hdot-l2m0.m` | 2662 | `2cf3e58f43002cb93b964dfefc7a4c7c23ecbeaddb8c505b73bc40e11a0b0c62` |
| `Kerr/Circular/Memory/hdot-l4m0.m` | 2686 | `49d37eb0547494b3abc0f442e969e812c75b83005894b7eafc46dd40856b2b80` |
| `Kerr/Circular/Memory/hdot-l6m0.m` | 3417 | `ab23950a788227a30c73bcd21ae7204ebade40892c932a23402c819d8431c540` |
| `Kerr/Circular/Memory/hdot-l8m0.m` | 4523 | `16cc287beeea200878556cfee71c7e302d7eecc6f48d857b9b2af35e72c410bc` |
| `Kerr/Circular/README.md` | 1455 | `fdac5712f9f26df6709b3475500a09db6079eebee6259d5132f4f3cf36b99df6` |
| `Kerr/Eccentric/Local/Redshift.m` | 35900 | `18d43ec76929226b264cf3146463c89385be1988dcf9e973599cbc4e6ed27339` |
| `Kerr/Generic/Flux/EInfinity-v.m` | 4172 | `076b3cd850225f55f973964f90aa78456a68e226709a45b0d89c79ab25e53bdb` |
| `Kerr/Spherical/Flux/AngularMomentum/Horizon.m` | 3189620 | `497cb86954a401a9d686364def063ddeb95c6204367eec5de897ef07f31c6cfc` |
| `Kerr/Spherical/Flux/AngularMomentum/Infinity.m` | 273262 | `f9373dc3b4e74b3eda1060809e8ff926c4e915c55c1a652e8f51a98c933b152d` |
| `Kerr/Spherical/Flux/Energy/Horizon.m` | 3395691 | `e75cff3699c4cbf894104d3d4ab4d81afe0cf1382cbc749a4c794b40812e2ae9` |
| `Kerr/Spherical/Flux/Energy/Infinity.m` | 250315 | `f40c52999c81ecd562427cd054116cee85791557ff6f3eeb910cc61101cd9057` |
| `Kerr/Spherical/Orbit/AngularMomentum.m` | 3754 | `fa4c28b6124faa2291c38834a014d298f197dc2aac8820a56515576d8259dba0` |
| `Kerr/Spherical/Orbit/AzimuthalFrequency.m` | 3305 | `1965afeabe9d7a314152d0daf50cfaf5957a6c525ac3f77bcff52942d0382132` |
| `Kerr/Spherical/Orbit/CarterConstant.m` | 3476 | `7e8a4f181db02c2dd5f5e6f45986fb59efb377d45cb7dfc63268c003d825d60d` |
| `Kerr/Spherical/Orbit/Energy.m` | 2965 | `b5e81ddbb28b724e0424ad80129b57db010b09a445be338904e1b18ab1503d6b` |
| `Kerr/Spherical/Orbit/PolarFrequency.m` | 2339 | `bb6f4169e97a050d7a473cfb539885a43c0302d6b44f783f7511e4e933003f73` |
| `Kerr/Spherical/ScalarFlux/AngularMomentum/Horizon.m` | 3593992 | `320b4b60023fbc64996926c177052b3d419cc0fd67fede570b86989842f0e3f5` |
| `Kerr/Spherical/ScalarFlux/AngularMomentum/Infinity.m` | 837499 | `560e768bc7f661bbe3601349daec5f509315349a441f2faa72d9a094b81164fb` |
| `Kerr/Spherical/ScalarFlux/Energy/Horizon.m` | 3735737 | `8e4fb77d318e6f4af5e05e6c437dd6fbee32a727a7a3c4612387faf9fe37aeea` |
| `Kerr/Spherical/ScalarFlux/Energy/Infinity.m` | 1018015 | `bd8933920c40074025176f06eaaf79ee31fd2bed762d52504b57c721b416ea47` |
| `Schwarzschild/Circular/Flux/Horizon.m` | 121099 | `1e86c103f4b160690eef90e8f2b6dab3801669d7f3357b4f8a5bcfc142be4cff` |
| `Schwarzschild/Circular/Flux/Infinity-l2m1.m` | 363827 | `198483693d358ca0db1b86ee2e1e83a0a0f4edbf14e7b1da691c6cdbe92443c1` |
| `Schwarzschild/Circular/Flux/Infinity-l2m2.m` | 338471 | `f070f929e94b6f55ae7388213fa76deb56858a3bf284704750808f1f44e687f1` |
| `Schwarzschild/Circular/Flux/Infinity-l3m1.m` | 370678 | `fb0784583eb767e34254f8fbc92ecfcc8d241ccf05526037ab2622b739c69c53` |
| `Schwarzschild/Circular/Flux/Infinity-l3m2.m` | 270260 | `5c801ad521ba6e4b792321f65e42bcae2d5faeb8c81633a1ac2997c0a296f835` |
| `Schwarzschild/Circular/Flux/Infinity-l3m3.m` | 709697 | `b43b89af5c2c297c397d512f1858bc6c6f0c773dba9ab60123892c2a8b2956b1` |
| `Schwarzschild/Circular/Flux/Infinity-l4m1.m` | 231130 | `3250bd970ba3b947ef86e19a3683f7eb353aca015285c4c243adb2520c340249` |
| `Schwarzschild/Circular/Flux/Infinity-l4m2.m` | 310320 | `daadee50900f039e677bd3f652a7fbf185d114b6923bbeeab7abf5c42204f025` |
| `Schwarzschild/Circular/Flux/Infinity-l4m3.m` | 425823 | `4396744f8e047b89c0196f1868633b72f6f224e0a5fdc01017f491823ef02f78` |
| `Schwarzschild/Circular/Flux/Infinity-l4m4.m` | 322404 | `bf973c63d208a95e6b276fc665f7a89af2d63be153684cc873c0d902668ffa03` |
| `Schwarzschild/Circular/Flux/Infinity-l5m1.m` | 242556 | `b16114a66d86657537ea07b0dc5c554e2dcb58b5eee294f170f28ee322c7b6fc` |
| `Schwarzschild/Circular/Flux/Infinity-l5m2.m` | 184130 | `380dfb2f79b66e82813937bb44df7e73f92b46123c1b90485453838a95fe6fd3` |
| `Schwarzschild/Circular/Flux/Infinity-l5m3.m` | 445793 | `5be3509de304de9559f628501bc3978c3221620d7dcc2fca48406bad5aa1aa89` |
| `Schwarzschild/Circular/Flux/Infinity-l5m4.m` | 191436 | `b3e660d3a41908d0ce2258da12e750ee21717e2335da00635654dff99e94060d` |
| `Schwarzschild/Circular/Flux/Infinity-l5m5.m` | 443313 | `bb23e248bc7508a0596b256b8e1f164715b333bd0068c328ed5c7e4fee418aa3` |
| `Schwarzschild/Circular/Flux/Infinity.m` | 1220602 | `c1ec5cfb431fc413c46e95f0529cd60fbdbbab630429966d101170004ce199e2` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l2m1.m` | 1588 | `8301d962bf50e6f893b4e07803f847d4b144d8da80835f4c73a390f63eead1fa` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l2m2.m` | 1179 | `ee1c6b0291c5a94e0fc230b72a9bfdc602f1ae23eed949640b7c2beab17999ff` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l3m1.m` | 1255 | `d86261e86a4b0a94469ea647bcad7206f9f8952a0cea8e1e162bcda2662b8ed3` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l3m2.m` | 1638 | `f8df8499e9944cf305fb79a9dfcf374853bf33aa9963004a2dd24cd07cac6b9c` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l3m3.m` | 1336 | `73c127ac04b3c75fca093758744b2301eeb104c6c4c9e03ef469b887e1240468` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l4m1.m` | 1748 | `1038d97135e85e19ad74b73900339682710469ae8337d6445081a3dfd1ee6ddf` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l4m2.m` | 1361 | `395d8b7072a653cf88edeed3cb592562086b3361bfb076ed367147342d22c56a` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l4m3.m` | 1908 | `b33fc56073f2da446c18a56fac821cdf4e33cae4f9e403f4c5a820b6bff6261d` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l4m4.m` | 1378 | `47dc2648ea74225b51a2385a63182c16211f5430f2ba88a169a75319e5e2bc45` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l5m1.m` | 1354 | `9156dd70ba5f641e7d5e317488f0b422c443bb7863e9bdde4dd09787826ffed5` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l5m2.m` | 1870 | `0b1f32a74eba029793bb67bccd258ffadc1b198bab20e5981d236543750b92cd` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l5m3.m` | 1431 | `71b9968bf5e303fd6cc9212717c0ea62134f95c9d12cddbcd0b87d189facc52a` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l5m4.m` | 1957 | `c2fb57f529aa38846264c8edc29bb8dcbc4ff66ae1ecfdb76b97dac2fdee469f` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l5m5.m` | 1461 | `f0545d3f1ee33272679f41efb1730bcc5f24ecef278c72e3f4fedf14e3b8ecd2` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l6m1.m` | 1875 | `08b89cff93fc12ca6d857b10c92b290f305a89245d54336fb1bfe389225cc114` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l6m2.m` | 1413 | `e256808e72da6f536ca7919733458410becab84e618adec457794c188bb01f84` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l6m3.m` | 2105 | `4e2beb7106d0084eaf8a262b94dcf71d7bcd755f216d01f63c2794715d42f933` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l6m4.m` | 1450 | `f2e271e43c71e2b2cd68d0c0854b8f034e753e24d17009078bbaccdcb62780c7` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l6m5.m` | 2255 | `ce66e716ade872a5677f76e5d8f9e37e2d7d3eb1af45e9249c67dd578cffe7f1` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l6m6.m` | 1494 | `8382bc7d858c2e494012ad50947b6d5625624d33016c9052bd756b9a27a65da8` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l7m1.m` | 1457 | `025850655b9ee08824653f722fc87ccf3d1a4d846dbf3910a28fc87b22397e02` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l7m2.m` | 2026 | `45d93691df41e94460c07f1a72dd0e65c498cc01dbacd3ecdb828787b859a383` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l7m3.m` | 1584 | `0ba89bfd9de38a8529307ec00cf278e677f7d947b9bfea5b7eecf944ad794da1` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l7m4.m` | 2196 | `bbcc9531ad5f995ba008300ae41b1d8eb8f17465810e743e7b409c1662b550de` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l7m5.m` | 1615 | `211a24e00a5ba00e9f3b77afe47b000ba2619c5c4c7691e3116054c0270be135` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l7m6.m` | 2288 | `9b9b2c5b7b70b6212ed0b46cd8bbb28468e5c58cf53c552745cba0b3ce87f20d` |
| `Schwarzschild/Circular/Flux/Spin/Infinity-l7m7.m` | 1629 | `a027760d42aeebda2437d65a458b11583b08dcb4b5250a8ba7704dc7a1820998` |
| `Schwarzschild/Circular/Local/Redshift.m` | 695928 | `e1e13a6eb6eb41c257655b63c91e02d042baf44bbfe807093c3f2244d1be3b05` |
| `Schwarzschild/Circular/Local/Spin.m` | 929861 | `aa293853a1372ec8d4308c19fe9cff9cba76e6c0a5a284edad52e7d02e8bbe16` |
| `Schwarzschild/Circular/Local/Spin/dEdt.m` | 1806 | `3f87888cc6921ba87f278e2a5996a5c3cbb8dc90d05dc3e32b564e46378284ba` |
| `Schwarzschild/Circular/Local/Tidal-lambdaB.m` | 913142 | `ef5205257fb73488875de255904c50c61f5844dd08fc96c65ad37f023593a6a5` |
| `Schwarzschild/Circular/Local/Tidal-lambdaE1.m` | 979648 | `d5081b76c36fc6a993b650f46dc6d9aa048af2929cb8215817787ee84f571a7f` |
| `Schwarzschild/Circular/Local/Tidal-lambdaE2.m` | 980686 | `3c5907e48e79f304c05ecfba3d12cc09f56608d5365c4b80f6486a1264c2d299` |
| `Schwarzschild/Circular/Local/Tidal-lambdaE3.m` | 986306 | `9bec17f584c0e5112f5d986037493032337487867f18f3637bdcafa15556d1bc` |
| `Schwarzschild/Circular/README.md` | 6896 | `6228671e68f99c296e0b5c3051559fd6f108fcd94ad4a11b579a605e60caf44b` |
| `Schwarzschild/Eccentric/Flux/EHorizon.m` | 1336017 | `cb6eca430a8e0be06f56e94ba4a28653454fa333cc80a5136f04fa363069eada` |
| `Schwarzschild/Eccentric/Flux/EInfinity.m` | 3728867 | `581d5299d9b04fc6c1a076cc80b627625b69d5c01f259a6c796cdfd6e582e4e4` |
| `Schwarzschild/Eccentric/Flux/LHorizon.m` | 1106699 | `0fa4f8154dc5815cc12d744c5c75a34c6cf5bf5fa2b07b7e0df4df2e455ec8fb` |
| `Schwarzschild/Eccentric/Flux/LInfinity.m` | 3713026 | `fdcebbcd0ca0173fb425eace9d9347a17adf2d9034aac3707ee5d7b89178dd9d` |
| `Schwarzschild/Eccentric/Local/Redshift-p.m` | 119741 | `4c37f096402720f68264ad26c5b50cea4ec71561eee5c50c5b8e3f63e5f7a6b1` |
| `Schwarzschild/Eccentric/Local/Redshift-y.m` | 120569 | `47f2d44818b51efa102128c53ed9b5fcfcab37ebea3cee319b3a0ed6b43b712f` |
| `Schwarzschild/Eccentric/Local/Spin-p.m` | 57894 | `ee7ae3bd485164bd93f4c311f5ab37a952680497ab9d706f62bba469d5bae95d` |
| `Schwarzschild/Eccentric/Local/Spin-y.m` | 57855 | `ed28f5d584a9d985a0c5e5429ebf69807eb8ff33855f90fe2108e31753df8fd3` |
| `Schwarzschild/Eccentric/README.md` | 3807 | `f3adf39a1126879c8ec6e838221a2d16bdc263c2a228c0a52b5dad5d8ab1472d` |
