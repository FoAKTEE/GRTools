using Test
using NumericsCore
using LinearAlgebra

# Independent O(N^2) definition: reduced phases avoid argument-reduction error
# growing with N. No FFT butterfly, permutation, or twiddle table is reused.
function fft_test_dft(x::AbstractVector{Complex{T}}; inverse=false) where {T}
    n = length(x)
    direction = inverse ? 1 : -1
    y = similar(x)
    for k in 0:n-1
        total = zero(Complex{T})
        for j in 0:n-1
            phase = T(direction * 2 * mod(j * k, n)) / T(n)
            total += x[j+1] * complex(cospi(phase), sinpi(phase))
        end
        y[k+1] = inverse ? total / T(n) : total
    end
    return y
end

fft_test_input(T, n) = [complex(T(mod(17j + 3, 101) - 50) / T(51),
                               T(mod(29j + 11, 103) - 51) / T(53)) for j in 0:n-1]

function fft_test_accuracy(::Type{T}) where {T}
    max_error = zero(T)
    max_fraction = zero(T)
    for n in (1, 2, 4, 8, 16, 32, 64, 128, 256)
        x = fft_test_input(T, n)
        original = copy(x)
        bound = T(10) * eps(T) * T(log2(n)) * norm(x)
        forward = fft_radix2(x)
        error = norm(forward - fft_test_dft(x))
        @test error <= bound
        @test norm(fft_radix2(x; inverse=true) - fft_test_dft(x; inverse=true)) <= bound
        @test norm(fft_radix2(forward; inverse=true) - x) <= bound
        @test abs(sum(abs2, forward) / T(n) - sum(abs2, x)) <= bound * norm(x)
        @test x == original
        @test eltype(forward) === Complex{T}
        if T === BigFloat
            @test all(z -> precision(real(z)) == precision(T) &&
                           precision(imag(z)) == precision(T), forward)
        end
        max_error = max(max_error, error)
        n > 1 && (max_fraction = max(max_fraction, error / bound))
    end
    println("FFT_ACCURACY type=", T, " bits=", precision(T),
            " N=1:2^8 max_L2_error=", max_error,
            " max_fraction_of_bound=", max_fraction)
end

@testset "radix-2 FFT" begin
    @testset "DFT, inverse, Parseval, precision" begin
        fft_test_accuracy(Float64)
        for bits in (128, 256)
            setprecision(BigFloat, bits) do
                fft_test_accuracy(BigFloat)
            end
        end
    end

    @testset "phase, single bin, real symmetry, edge signals" begin
        for T in (Float64, BigFloat)
            n = 32
            bound = T(10) * eps(T) * T(log2(n))
            for bin in (0, 1, 7, n ÷ 2, n - 1)
                x = [complex(cospi(T(2 * mod(bin*j, n)) / T(n)),
                             sinpi(T(2 * mod(bin*j, n)) / T(n))) for j in 0:n-1]
                expected = zeros(Complex{T}, n)
                expected[bin+1] = T(n)
                @test norm(fft_radix2(x) - expected) <= bound * norm(x)
            end
            impulse = zeros(Complex{T}, n)
            impulse[4] = one(T)
            @test norm(fft_radix2(impulse) - fft_test_dft(impulse)) <= bound
            @test fft_radix2(zeros(Complex{T}, n)) == zeros(Complex{T}, n)
            @test fft_radix2(ones(Complex{T}, n)) == vcat(Complex{T}(n), zeros(Complex{T}, n-1))
            real_input = real.(fft_test_input(T, n))
            y = fft_radix2(real_input)
            @test eltype(y) === Complex{T}
            @test norm(y - fft_test_dft(complex.(real_input))) <= bound * norm(real_input)
            @test abs(imag(y[1])) <= bound * norm(real_input)
            @test abs(imag(y[n÷2+1])) <= bound * norm(real_input)
            @test all(k -> abs(y[k+1] - conj(y[n-k+1])) <= bound * norm(real_input), 1:n-1)
        end
        @test fft_radix2([1, 2, 3, 4]) == ComplexF64[10, -2+2im, -2, -2-2im]
        @test eltype(fft_radix2(Float32[1, 2])) === ComplexF32
    end

    @testset "selected axes and strided 2-D composition" begin
        for T in (Float64, BigFloat)
            a = reshape(fft_test_input(T, 32), 4, 8)
            by_rows = reduce(vcat, [permutedims(fft_test_dft(a[i, :])) for i in axes(a, 1)])
            reference = reduce(hcat, [fft_test_dft(by_rows[:, j]) for j in axes(a, 2)])
            bound = T(10) * eps(T) * T(log2(length(a))) * norm(a)
            @test norm(fft_radix2(a; dims=2) - by_rows) <= bound
            @test norm(fft_radix2(a; dims=(1, 2)) - reference) <= bound
            @test norm(fft_radix2(a; dims=(2, 1)) - reference) <= bound
            @test norm(fft_radix2(fft_radix2(a); inverse=true) - a) <= bound
            @test fft_radix2(a; dims=()) == a
            parent = fill(Complex{T}(99), 8, 16)
            strided = @view parent[1:2:8, 2:2:16]
            strided .= a
            @test fft_radix2!(strided; dims=(1, 2)) === strided
            @test norm(strided - reference) <= bound
            @test all(==(Complex{T}(99)), parent[2:2:8, :])
            @test all(==(Complex{T}(99)), parent[:, 1:2:16])
            row = @view parent[1, 2:2:16]
            row .= a[1, :]
            @test fft_radix2!(row) === row
            @test norm(row - by_rows[1, :]) <= bound
            # An untransformed batch dimension need not be a power of two.
            batch = reshape(fft_test_input(T, 24), 4, 3, 2)
            transformed = fft_radix2(batch; dims=[1, 3])
            @test norm(fft_radix2(transformed; dims=(1, 3), inverse=true) - batch) <=
                  T(10) * eps(T) * T(3) * norm(batch)
        end
    end

    @testset "validation before mutation" begin
        for n in (0, 3, 6, 12)
            x = fft_test_input(Float64, n)
            original = copy(x)
            @test_throws ArgumentError fft_radix2!(x)
            @test x == original
        end
        x = reshape(fft_test_input(Float64, 12), 4, 3)
        original = copy(x)
        for dims in ((1, 2), (1, 1), (0,), (3,), (1, 1.5))
            @test_throws ArgumentError fft_radix2!(x; dims)
            @test x == original
        end
        @test fft_radix2!(ComplexF64[2+3im]) == [2+3im]
    end

    @testset "concurrent calls, deterministic and precision-preserving" begin
        # Interleaving needs >= 2 threads (CI sets JULIA_NUM_THREADS=4). A serial run
        # still checks determinism and precision but records the race gate as skipped.
        Threads.nthreads() >= 2 || @test_skip "concurrent FFT race gate needs -t2"
        original_precision = precision(BigFloat)
        # Change precision only before spawning, and join every task before
        # leaving the block. No task mutates the global precision setting.
        for bits in (128, 256, 128)
            setprecision(BigFloat, bits) do
                float_input = reshape(fft_test_input(Float64, 128), 8, 16)
                big_input = reshape(fft_test_input(BigFloat, 128), 8, 16)
                cases = [(a, dims, inverse) for a in (float_input, big_input),
                         dims in ((1,), (2,), (1, 2)), inverse in (false, true)]
                serial = map(c -> fft_radix2(c[1]; dims=c[2], inverse=c[3]), cases)
                matched = fill(false, length(cases))
                Threads.@threads for i in eachindex(cases)
                    a, dims, inverse = cases[i]
                    matched[i] = all(1:8) do _
                        result = fft_radix2(a; dims, inverse)
                        result == serial[i] && precision(BigFloat) == bits &&
                            (eltype(result) !== Complex{BigFloat} ||
                             all(z -> precision(real(z)) == bits && precision(imag(z)) == bits, result))
                    end
                end
                @test all(matched)
                @test precision(BigFloat) == bits
            end
        end
        @test precision(BigFloat) == original_precision
        println("FFT concurrency: threads=", Threads.nthreads(), " bits=128,256,128 PASS")
    end

    @testset "large transforms and bounded scratch allocation" begin
        for n in (1024, 4096)
            x = fft_test_input(Float64, n)
            y = fft_radix2(x)
            bound = 10eps(Float64) * log2(n) * norm(x)
            @test norm(fft_radix2(y; inverse=true) - x) <= bound
            @test abs(sum(abs2, y) / n - sum(abs2, x)) <= bound * norm(x)
            work = copy(x)
            fft_radix2!(work) # compile before measuring
            copyto!(work, x)
            allocated = @allocated fft_radix2!(work)
            # One half-length ComplexF64 twiddle table, with generous room for
            # fixed bookkeeping. Reject full-size scratch copies/per-stage tables.
            @test allocated <= 8n + 4096
            println("FFT allocation: N=", n, " bytes=", allocated)
        end
    end
end
