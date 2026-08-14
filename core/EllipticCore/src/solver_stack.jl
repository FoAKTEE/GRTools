# solver_stack.jl — composition API solve!(u, A, b; pc = MultigridV(...;
# smoother=SchwarzSmoother(...; passes=2), bottom=CachedSparseLU())); right-
# preconditioned GMRES; residual monitor.
# DAG node: solver-stack-core. Built TDD — no code lands before its failing tests.
