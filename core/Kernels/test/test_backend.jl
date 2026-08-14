# Characterization of the backend seam. Grounding:
#   - select_backend / backend_array / backend_collect: verbatim vs NumericsCore.backend_iface (LIGHT live
#     dep). gpu=false and the gpu=true fallback (no GPU ext loaded) both resolve to CPUBackend(); the host
#     array/collect helpers are identity on CPUBackend.
#   - ka_device(CPUBackend()) is the concrete KernelAbstractions CPU device (identified by name, KA being a
#     transitive dep not a test dep).
#   - with_kernel_threads: closed-form coverage — applies f to every index 1:n exactly once, and the
#     bounded-threaded path is BIT-IDENTICAL (===) to the serial path because each work-item writes a
#     disjoint output (the §8 parallelism contract / bit-equality gate).

using Test
using Kernels

@testset "backend" begin
    @testset "select_backend + host helpers (vs NumericsCore)" begin
        @test select_backend() === CPUBackend()
        @test select_backend(; gpu = false) === CPUBackend()
        # No GPU extension is loaded, so the gpu=true request falls back to the host backend.
        @test select_backend(; gpu = true) === CPUBackend()
        @test select_backend(; gpu = true, kind = :cuda) === CPUBackend()

        b = CPUBackend()
        a = backend_array(b, Float64, 3)
        @test a isa Array{Float64,1} && length(a) == 3
        src = [1.0, 2.0, 3.0]
        @test backend_collect(b, src) == src && backend_collect(b, src) isa Array
    end

    @testset "ka_device" begin
        @test nameof(typeof(ka_device(CPUBackend()))) === :CPU
    end

    @testset "with_kernel_threads: covers 1:n once, parallel===serial" begin
        for n in (0, 1, 5, 23, 100)
            # full coverage, each index visited exactly once
            seen = zeros(Int, n)
            with_kernel_threads(n) do i
                seen[i] += 1
            end
            @test all(==(1), seen)

            # serial vs bounded-threaded must be BIT-IDENTICAL (disjoint per-index writes)
            xs = Float64[sin(i) * 1e3 + cos(3i) for i in 1:n]
            out_s = Vector{Float64}(undef, n)
            out_t = Vector{Float64}(undef, n)
            # a nontrivial serial inner reduction inside each work-item
            kern = (dst, i) -> begin
                acc = 0.0
                for k in 1:7
                    acc += xs[i] / (k + 0.5) - xs[i] * k
                end
                dst[i] = acc
            end
            with_kernel_threads(n; threaded = false) do i
                kern(out_s, i)
            end
            with_kernel_threads(n; threaded = true) do i
                kern(out_t, i)
            end
            @test out_s === out_s
            @test out_t == out_s
            @test all(out_t[i] === out_s[i] for i in 1:n)
        end
    end
end
