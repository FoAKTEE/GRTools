# amr.jl — h-split (halve element in both/either dims) + p-increment mesh
# transformations, 2:1 enforcement by splitting offenders, Lagrange interpolation of
# the current solution onto the new mesh; criteria are caller-supplied predicates.
# DAG nodes: amr-h-refinement, amr-p-refinement. Built TDD — no code lands before
# their failing tests.
