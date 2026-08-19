# EllipticCore building blocks: spectral reference elements and complex-generic GMRES.
#
# Run from the repository root:
#     julia --project=. examples/06_elliptic_core.jl
#
# EllipticCore is the physics-free hp-DG elliptic solver stack (reference elements,
# multi-block meshes, matrix-free DG operators, GMRES/Schwarz/multigrid). This
# example exercises its two most reusable primitives directly.

using EllipticCore
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "_check.jl"))   # check / check_close: every printed cross-check is asserted

# --- Legendre-Gauss reference element: exact quadrature and spectral differentiation ---
N = 8
el = refelement(Float64, N)                       # N nodes on (-1, 1)
@printf "N=%d Legendre-Gauss element:  sum(w) = %.15f (exact 2)\n" N sum(el.weights)
check_close("Gauss weights sum to the measure of (-1,1) == 2 (exact)", sum(el.weights), 2.0; atol = 1e-13)

# Gauss quadrature integrates polynomials up to degree 2N-1 exactly.
integ = dot(el.weights, el.nodes .^ (2N - 2))     # int x^(2N-2) dx = 2/(2N-1)
@printf "quadrature of x^%d: %.15f (exact %.15f)\n" (2N - 2) integ 2 / (2N - 1)
check_close("Gauss quadrature is exact at degree 2N-2", integ, 2 / (2N - 1); rtol = 1e-13)

# The differentiation matrix D is exact on the polynomial space.
u = el.nodes .^ 5
du = el.D * u
maxerr = maximum(abs, du .- 5 .* el.nodes .^ 4)
@printf "spectral derivative of x^5: max error = %.2e\n" maxerr
check("differentiation matrix is exact on the polynomial space (max err = $maxerr)", maxerr < 1e-12)

# --- Complex-generic GMRES on a Helmholtz-like shifted system ---------------------------
# gmres(apply!, b) takes a matrix-free applier (out, v) -> out = A*v.
n = 200
T = ComplexF64
A = Matrix{T}(5.0I, n, n)
for i in 1:n-1
    A[i, i+1] = -1.0 + 0.3im                      # non-hermitian complex off-diagonals
    A[i+1, i] = -1.0 - 0.2im
end
b = normalize!(T[sin(0.7k) + im * cos(0.3k) for k in 1:n])

res = gmres((out, v) -> mul!(out, A, v), b)
@printf "\nGMRES on %dx%d complex system: converged = %s in %d iterations\n" n n res.converged res.iterations
truer = norm(b - A * res.x)
@printf "true residual |b - A x| = %.2e (reported %.2e)\n" truer res.residual
check("GMRES reports convergence", res.converged)
check("GMRES true residual is at solver tolerance (|b - A x| = $truer)", truer < 1e-10)
check_close("GMRES reported residual matches the recomputed one", res.residual, truer; rtol = 0.05, atol = 1e-14)
