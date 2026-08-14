# Scoped BLAS-thread guard: OpenBLAS inside a Threads.@threads region deadlocks on its internal
# locks, and multi-thread BLAS reductions are nondeterministic in the last bits — pin for the
# scope, restore on exit/exception.

using LinearAlgebra: BLAS

# BLAS' thread count is process-global.  A plain get/set/restore sequence is therefore racy: one
# task can restore its old value while another task still relies on the scoped value.  Leases with
# the same requested count may overlap, while a different count waits for the current lease group
# to drain.  Per-task depth makes same-count nesting reentrant; a conflicting nested request is an
# error because waiting for one's own outer lease can never make progress.
const _blas_lease_lock = ReentrantLock()
const _blas_lease_changed = Threads.Condition(_blas_lease_lock)
const _blas_active_threads = Ref{Union{Nothing,Int}}(nothing)
const _blas_previous_threads = Ref(0)
const _blas_lease_count = Ref(0)
const _blas_task_depth = IdDict{Task,Int}()

function _acquire_blas_threads(n::Int)
    n > 0 || throw(ArgumentError("BLAS thread count must be positive, got $n"))
    task = current_task()
    lock(_blas_lease_lock)
    try
        depth = get(_blas_task_depth, task, 0)
        if depth > 0
            active = _blas_active_threads[]
            active == n || throw(ArgumentError(
                "cannot request $n BLAS threads inside this task's active $active-thread lease"))
            _blas_task_depth[task] = depth + 1
            _blas_lease_count[] += 1
            return nothing
        end

        while _blas_active_threads[] !== nothing && _blas_active_threads[] != n
            wait(_blas_lease_changed)
        end
        if _blas_active_threads[] === nothing
            previous = Int(BLAS.get_num_threads())
            BLAS.set_num_threads(n)
            _blas_previous_threads[] = previous
            _blas_active_threads[] = n
        end
        _blas_task_depth[task] = 1
        _blas_lease_count[] += 1
        return nothing
    finally
        unlock(_blas_lease_lock)
    end
end

function _release_blas_threads()
    task = current_task()
    lock(_blas_lease_lock)
    try
        depth = get(_blas_task_depth, task, 0)
        depth > 0 || error("BLAS lease released by a task that does not hold one")
        if depth == 1
            delete!(_blas_task_depth, task)
        else
            _blas_task_depth[task] = depth - 1
        end
        _blas_lease_count[] -= 1
        if _blas_lease_count[] == 0
            try
                BLAS.set_num_threads(_blas_previous_threads[])
            finally
                _blas_active_threads[] = nothing
                _blas_previous_threads[] = 0
                notify(_blas_lease_changed; all = true)
            end
        end
        return nothing
    finally
        unlock(_blas_lease_lock)
    end
end

"""`with_blas_threads(f, n)` — run `f()` with BLAS pinned to `n` threads under a process-global
lease. Same-`n` scopes may overlap and are reference-counted; conflicting scopes serialize. A
same-`n` nested scope is reentrant, while a conflicting nested request throws instead of deadlocking.
The thread count preceding the first lease is restored after the final scope, including on exception."""
function with_blas_threads(f, n::Integer)
    ni = Int(n)
    _acquire_blas_threads(ni)
    try
        return f()
    finally
        _release_blas_threads()
    end
end

"""`blas_single_thread(f)` — run `f()` with OpenBLAS pinned to a single thread (the deadlock-safe,
bit-deterministic setting for tiny-matrix LAPACK inside a threaded sweep)."""
blas_single_thread(f) = with_blas_threads(f, 1)
