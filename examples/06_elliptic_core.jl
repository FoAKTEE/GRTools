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

# --- Legendre-Gauss reference element: exact quadrature and spectral differentiation ---
N = 8
el = refelement(Float64, N)                       # N nodes on (-1, 1)
@printf "N=%d Legendre-Gauss element:  sum(w) = %.15f (exact 2)\n" N sum(el.weights)

# Gauss quadrature integrates polynomials up to degree 2N-1 exactly.
integ = dot(el.weights, el.nodes .^ (2N - 2))     # int x^(2N-2) dx = 2/(2N-1)
@printf "quadrature of x^%d: %.15f (exact %.15f)\n" (2N - 2) integ 2 / (2N - 1)

# The differentiation matrix D is exact on the polynomial space.
u = el.nodes .^ 5
du = el.D * u
@printf "spectral derivative of x^5: max error = %.2e\n" maximum(abs, du .- 5 .* el.nodes .^ 4)

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
@printf "true residual |b - A x| = %.2e (reported %.2e)\n" norm(b - A * res.x) res.residual
