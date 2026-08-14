# Convention-reconciliation helpers shared across the EMRI/MST/PN entry points.
#
# Single source of truth for the dimensionless EMRI conventions: the MST notebook variables
# `eps = 2*M*omega` and `q = a/M` (KerrBackwards `p_to_y`; ScatteringPhaseMST `_nu_lambda_omega`;
# TeukolskyPN `imagunit`). All are M-parametrized (default M=1) so callers in either convention agree.

"""`p_to_y(p, e)` — map semilatus rectum `p` to the amplitude-grid coordinate `y = log(p - 2e - 21/10)`
(the shift keeps `y` finite and clear of the Schwarzschild separatrix). (KerrBackwards `p_to_y`.)"""
p_to_y(p, e) = log(p - 2*e - 21//10)

"""`eps_to_omega(eps, M=1)` — MST frequency parameter `eps = 2*M*omega` -> `omega = eps/(2M)`."""
eps_to_omega(eps, M=1) = eps/(2*M)

"""`omega_to_eps(omega, M=1)` — `omega` -> the MST frequency parameter `eps = 2*M*omega`."""
omega_to_eps(omega, M=1) = 2*M*omega

"""`q_to_a(q, M=1)` — dimensionless spin `q = a/M` -> Boyer-Lindquist `a = q*M`."""
q_to_a(q, M=1) = q*M

"""`a_to_q(a, M=1)` — Boyer-Lindquist spin `a` -> dimensionless `q = a/M`."""
a_to_q(a, M=1) = a/M

"""`imagunit()` / `imagunit(::Type{T})` — the numeric imaginary unit (the value TeukolskyPN's symbolic
`imagunit()` substitutes to); type-parametrized form returns `Complex{T}(0,1)` for precision control."""
imagunit() = im
imagunit(::Type{T}) where {T<:AbstractFloat} = Complex{T}(0, 1)
