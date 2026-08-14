# Characterization of the scoped BLAS-thread guard (genuinely-new primitive).
#
# No legacy oracle: ground by (1) the restore invariant — BLAS.get_num_threads() returns to the prior
# count after the block, including on exception (try/finally); and (2) parallel==serial bit-identity —
# a threaded LAPACK region (eigvals) run inside blas_single_thread is bit-for-bit identical to the
# serial result, because pinning OpenBLAS to one thread removes the reduction-order nondeterminism.

using NumericsCore
using LinearAlgebra
using LinearAlgebra: BLAS

@testset "blas_guard" begin
    @testset "restore invariant" begin
        n0 = BLAS.get_num_threads()
        BLAS.set_num_threads(4)
        got = NumericsCore.with_blas_threads(2) do
            BLAS.get_num_threads()
        end
        @test got == 2
        @test BLAS.get_num_threads() == 4              # restored
        # restored even when the block throws
        @test_throws ErrorException NumericsCore.blas_single_thread() do
            error("boom")
        end
        @test BLAS.get_num_threads() == 4
        gots = NumericsCore.blas_single_thread() do
            BLAS.get_num_threads()
        end
        @test gots == 1
        @test BLAS.get_num_threads() == 4
        BLAS.set_num_threads(n0)
    end

    @testset "overlapping leases, conflicting serialization, and nesting" begin
        n0 = BLAS.get_num_threads()
        BLAS.set_num_threads(4)
        try
            # Two same-count leases overlap. Releasing the first must not restore the process-global
            # count while the second task still holds its lease.
            entered = Channel{Int}(2)
            release1 = Channel{Nothing}(1)
            release2 = Channel{Nothing}(1)
            t1 = @async NumericsCore.with_blas_threads(2) do
                put!(entered, BLAS.get_num_threads())
                take!(release1)
            end
            @test take!(entered) == 2
            t2 = @async NumericsCore.with_blas_threads(2) do
                put!(entered, BLAS.get_num_threads())
                take!(release2)
            end
            @test take!(entered) == 2
            put!(release1, nothing)
            wait(t1)
            @test BLAS.get_num_threads() == 2
            put!(release2, nothing)
            wait(t2)
            @test BLAS.get_num_threads() == 4

            # A conflicting request cannot enter until the current lease drains.
            held = Channel{Nothing}(1)
            release_held = Channel{Nothing}(1)
            t3 = @async NumericsCore.with_blas_threads(2) do
                put!(held, nothing)
                take!(release_held)
            end
            take!(held)
            started = Channel{Nothing}(1)
            entered3 = Channel{Int}(1)
            t4 = @async begin
                put!(started, nothing)
                NumericsCore.with_blas_threads(3) do
                    put!(entered3, BLAS.get_num_threads())
                end
            end
            take!(started)
            for _ in 1:4
                yield()
            end
            @test !isready(entered3)
            @test BLAS.get_num_threads() == 2
            put!(release_held, nothing)
            wait(t3)
            @test take!(entered3) == 3
            wait(t4)
            @test BLAS.get_num_threads() == 4

            # Same-task, same-count nesting is reentrant. A conflicting nested count is rejected
            # immediately (rather than waiting forever for the task's own outer lease).
            nested = NumericsCore.with_blas_threads(2) do
                inner = NumericsCore.with_blas_threads(2) do
                    BLAS.get_num_threads()
                end
                @test BLAS.get_num_threads() == 2
                @test_throws ArgumentError NumericsCore.with_blas_threads(3) do
                    nothing
                end
                @test BLAS.get_num_threads() == 2
                inner
            end
            @test nested == 2
            @test BLAS.get_num_threads() == 4

            @test_throws ErrorException NumericsCore.with_blas_threads(2) do
                NumericsCore.with_blas_threads(2) do
                    error("nested boom")
                end
            end
            @test BLAS.get_num_threads() == 4
        finally
            BLAS.set_num_threads(n0)
        end
    end

    @testset "parallel == serial bit-identity (single-thread BLAS)" begin
        # deterministic symmetric positive-definite matrix
        A = [1.0/(i + j - 1) + (i == j ? 5.0 : 0.0) for i in 1:40, j in 1:40]
        A = Symmetric((A + A')/2)
        ref = NumericsCore.blas_single_thread() do
            eigvals(A)
        end
        results = Vector{Vector{Float64}}(undef, 8)
        NumericsCore.blas_single_thread() do
            Threads.@threads for i in 1:8
                results[i] = eigvals(A)
            end
        end
        # every per-task LAPACK result is bit-for-bit identical to the serial reference
        @test all(r -> r == ref, results)
    end
end
