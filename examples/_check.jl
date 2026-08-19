# Shared assertion helper for the example scripts.
#
# Every example PRINTS labeled cross-checks against exact limits or identities; each
# printed check is also ASSERTED through `check`/`check_close`, so a regression aborts the
# script with a nonzero exit code instead of quietly printing a wrong number (CI runs all
# of examples/ as a job). Tolerances are the ones each script's own print format implies:
# a quantity printed to k digits alongside its exact value is gated at that scale, and a
# residual printed as "%.2e" is gated a few orders above the value the check reports.

"""
    check(name, ok) -> nothing

Assert an example cross-check. Throws an `ErrorException` naming the failed check when
`ok` is `false`.
"""
check(name::AbstractString, ok::Bool) = ok ? nothing :
    error("example cross-check FAILED: ", name)

"""
    check_close(name, got, want; atol = 0, rtol = 0) -> nothing

`check` on `|got - want| ≤ atol + rtol*|want|` (works for real and complex `got`/`want`).
At least one of `atol`, `rtol` must be given, since both default to zero (exact equality).
"""
function check_close(name::AbstractString, got, want; atol = 0.0, rtol = 0.0)
    d = abs(got - want)
    tol = atol + rtol * abs(want)
    check("$name — got $got, want $want, |Δ| = $d > tol = $tol", d ≤ tol)
end
