# Puncture.TensorBasis — Martel-Poisson tensor-harmonic angular basis (11 functions) + coordinate
# reconstruction of hR (htt..hφφ) and its deriv[nt,nr,nθ,nφ] from per-(l,m) Lorenz amplitudes on the
# worldline. The base seed θ-derivative `_dthY` is supplied by ForwardDiff over the consumed
# Harmonics.spin_weighted_spherical_Y (s=0), replacing the unexposed legacy SWSH `_y_dtheta`.
#
# Grounding:
#   • CAPTURED upstream oracle: packages/Punctures/test/tensor_basis_reference.csv (11 basis fns on
#     an (l,m,θ,φ) grid, rtol 1e-9/atol 1e-11) + coord_reference.csv (htt..hφφ + deriv on synthetic
#     modes, rtol 1e-8/atol 1e-10) from test_tensor_basis.jl, copied into test/fixtures/puncture/.
#   • Julia-first FD cross-check: the analytic θ-derivative machinery (ForwardDiff _dthY + the AngTerm
#     _dtheta term lists) equals central finite differences of the consumed Harmonics seed and of the
#     reconstructed coordinate component.
#   • reality/parity rule hR[i][l,-m] = (-1)^m conj + odd-in-m vanishing at m=0, enumerated.

using Test
using SelfForce
import Harmonics

const PT = SelfForce.Puncture
const TBREF = joinpath(@__DIR__, "fixtures", "puncture", "tensor_basis_reference.csv")
const COREF = joinpath(@__DIR__, "fixtures", "puncture", "coord_reference.csv")

# synthetic base amplitudes (m >= 0); MUST match gen_tensor_basis.wl.
_pt_synH(i, l, m) = (i + 1 + l - m) / 7 + im * (2 + i - l + 2m) / 5
_pt_synDH(i, l, m) = (3 + i + 2l + m) / 11 + im * (1 + i + l - m) / 4
_pt_synDDH(i, l, m) = (5 + 2i - l + m) / 13 + im * (4 + i - 2l + m) / 3
_pt_base(i, l, m) = (ComplexF64(_pt_synH(i, l, m)),
                     ComplexF64(_pt_synDH(i, l, m)),
                     ComplexF64(_pt_synDDH(i, l, m)))

function _pt_load_basis(path)
    rows = Tuple{Symbol,Int,Int,Float64,Float64,Float64,Float64}[]
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        isempty(strip(line)) && continue
        f = split(line, ',')
        push!(rows, (Symbol(f[1]), parse(Int, f[2]), parse(Int, f[3]),
                     parse(Float64, f[4]), parse(Float64, f[5]),
                     parse(Float64, f[6]), parse(Float64, f[7])))
    end
    return rows
end

function _pt_load_coord(path)
    rows = Tuple{Symbol,Int,Int,Int,Int,Int,Float64}[]
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        isempty(strip(line)) && continue
        f = split(line, ',')
        push!(rows, (Symbol(f[1]), parse(Int, f[2]), parse(Int, f[3]),
                     parse(Int, f[4]), parse(Int, f[5]), parse(Int, f[6]),
                     parse(Float64, f[7])))
    end
    return rows
end

@testset "Puncture.TensorBasis" begin

    @testset "basis functions vs captured oracle (l,m,θ,φ grid)" begin
        brows = _pt_load_basis(TBREF)
        @test length(brows) > 400
        for (key, l, m, θ, φ, re, im_) in brows
            v = PT.tensor_harmonic(key, l, m, θ, φ)
            @test isapprox(real(v), re; rtol = 1e-9, atol = 1e-11)
            @test isapprox(imag(v), im_; rtol = 1e-9, atol = 1e-11)
        end
    end

    @testset "coordinate components + deriv vs captured oracle (synthetic modes)" begin
        crows = _pt_load_coord(COREF)
        @test length(crows) > 500
        for (comp, l, nt, nr, nθ, nφ, val) in crows
            got = PT.coord_deriv(comp, l, _pt_base, 10.0, 1.0;
                                 nt = nt, nr = nr, nθ = nθ, nφ = nφ)
            @test isapprox(got, val; rtol = 1e-8, atol = 1e-10)
        end
    end

    @testset "ForwardDiff _dthY === central-FD of the consumed Harmonics seed" begin
        seed(l, m, θ) = real(Harmonics.spin_weighted_spherical_Y(0, l, m, θ, 0.0))
        h = 1e-5
        for l in (1, 2, 3), m in -l:l, θ in (0.6, 1.1, 2.0)
            # 1st derivative
            fd1 = (seed(l, m, θ + h) - seed(l, m, θ - h)) / (2h)
            @test isapprox(PT._dthY(l, m, θ, 1), fd1; rtol = 1e-5, atol = 1e-7)
            # 2nd derivative
            fd2 = (seed(l, m, θ + h) - 2seed(l, m, θ) + seed(l, m, θ - h)) / h^2
            @test isapprox(PT._dthY(l, m, θ, 2), fd2; rtol = 1e-4, atol = 1e-5)
            # 0th = seed itself
            @test PT._dthY(l, m, θ, 0) ≈ seed(l, m, θ)
        end
    end

    @testset "analytic nθ machinery === central-FD of coord component (Julia-first)" begin
        h = 1e-5
        for comp in (:htt, :htth, :hthth), l in (2, 3)
            for θ0 in (0.9, 1.3)
                ana = PT.coord_deriv(comp, l, _pt_base, 10.0, 1.0; nθ = 1, θ = θ0)
                fdp = PT.coord_deriv(comp, l, _pt_base, 10.0, 1.0; nθ = 0, θ = θ0 + h)
                fdm = PT.coord_deriv(comp, l, _pt_base, 10.0, 1.0; nθ = 0, θ = θ0 - h)
                @test isapprox(ana, (fdp - fdm) / (2h); rtol = 1e-5, atol = 1e-6)
            end
        end
    end

    @testset "reality identity Y_{l,-m} = (-1)^m conj(Y_{lm}) (oracle-independent)" begin
        keys = (:Ythth, :Zth, :Zph, :Xth, :Xph, :Zthth, :Zthph, :Zphph,
                :Xthth, :Xthph, :Xphph)
        for key in keys, l in (2, 3, 5), m in 1:min(l, 3)
            θ = 0.7; φ = 0.4
            lhs = PT.tensor_harmonic(key, l, -m, θ, φ)
            rhs = (-1)^m * conj(PT.tensor_harmonic(key, l, m, θ, φ))
            @test isapprox(lhs, rhs; rtol = 1e-10, atol = 1e-12)
        end
    end

    @testset "m=0 vanishing of odd-in-m basis (oracle-independent)" begin
        for key in (:Zph, :Xth, :Zthph, :Xthth, :Xphph), l in (2, 3, 5)
            @test abs(PT.tensor_harmonic(key, l, 0, 0.9, 0.3)) < 1e-12
        end
    end
end
