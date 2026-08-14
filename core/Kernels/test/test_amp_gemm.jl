# Characterization of the ROMAN dense-NN forward GEMM. NO live heavy dep (FEW/KerrBackwards forbidden).
# Grounding:
#   - leaky_relu: closed form x>=0 ? x : 0.2x, at x ∈ {-2, -1e-30, 0, 3}.
#   - roman_forward: HAND-COMPUTED 2-layer dense pass with tiny literal W/b/reduced_basis.
#       x=[1,2]; layer1 W1=I, b1=[0,-3] ⇒ W1x+b1=[1,-1] ⇒ leaky_relu=[1,-0.2];
#       layer2 (last, no relu) W2=[1 0;0 1;1 1;1 -1], b2=0 ⇒ out=[1,-0.2,0.8,1.2]; bi=2;
#       tfi=2 ⇒ nn=[(1+0.8i)2, (-0.2+1.2i)2]=[2+1.6i, -0.4+2.4i];
#       reduced_basis=[1 0; i 1] ⇒ modes = [2+1.6i,  i(2+1.6i)+(-0.4+2.4i)=-2+4.4i].
#   - PARALLEL == SERIAL bit-identity (===) over output rows.

using Test
using Kernels

@testset "amp_gemm" begin
    @testset "leaky_relu closed form" begin
        @test leaky_relu(-2.0) == -0.4
        @test leaky_relu(-1e-30) == 0.2 * (-1e-30)
        @test leaky_relu(0.0) === 0.0
        @test leaky_relu(3.0) === 3.0
    end

    @testset "roman_forward: hand-computed 2-layer pass" begin
        W1 = [1.0 0.0; 0.0 1.0]
        W2 = [1.0 0.0; 0.0 1.0; 1.0 1.0; 1.0 -1.0]
        weights = [W1, W2]
        biases = [[0.0, -3.0], [0.0, 0.0, 0.0, 0.0]]
        reduced_basis = ComplexF64[1 0; im 1]
        tfi = 2.0
        x = [1.0, 2.0]
        modes = roman_forward(weights, biases, reduced_basis, tfi, x)
        @test modes ≈ ComplexF64[2 + 1.6im, -2 + 4.4im] atol = 1e-12
    end

    @testset "roman_forward: parallel === serial (bit-identity)" begin
        W1 = [0.7 -0.3 0.1; 0.2 0.5 -0.4]          # 2x3 (n_out=2, n_in=3)
        W2 = [reshape([0.3 * (i + j) - 0.5 for i in 1:6, j in 1:2], 6, 2)]   # placeholder slot
        weights = [W1, W2[1]]                       # last layer 6 outputs ⇒ bi=3
        biases = [[0.05, -0.1], zeros(6)]
        x = [0.4, -1.2, 0.9]
        nrows = 40
        reduced_basis = ComplexF64[cis(0.13 * (r + c)) + 0.5im for r in 1:nrows, c in 1:3]
        tfi = 1.7
        ms = roman_forward(weights, biases, reduced_basis, tfi, x; threaded = false)
        mt = roman_forward(weights, biases, reduced_basis, tfi, x; threaded = true)
        @test all(mt[r] === ms[r] for r in 1:nrows)
    end
end
