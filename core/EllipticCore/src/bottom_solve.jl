# bottom_solve.jl — coarsest-level matrix by column probes -> SparseArrays.sparse ->
# UMFPACK lu (ComplexF64), cached across V-cycles and AMR iterations; dense-LU generic
# fallback for exotic T.
# DAG node: multigrid-bottom-solver. Built TDD — no code lands before its failing tests.
