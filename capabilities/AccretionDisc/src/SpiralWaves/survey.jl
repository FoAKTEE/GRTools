# SpiralWaves/survey.jl — the per-m and per-configuration torque survey driver of arXiv:2608.17003
# Sec. III: the azimuthal modes of the master system are INDEPENDENT (the forcing (60) is diagonal in
# m and the disc response is linear), so the mode sum T = Σ_m T_m, T_L = Σ_m T_{L,m} of (73)-(75) is an
# embarrassingly parallel loop. Each task builds its own per-mode adaptive grid (App. C, the domain
# (C1)-(C3) depends on m), solves (64) and reduces to the four scalars of one mode; the SUM is then
# taken SERIALLY in a fixed m order, so the survey is BITWISE reproducible and independent of the
# thread count.
#
# COMPUTE POLICY (binding, mission `compute-usage` skill): at most 24 CPUs per job and never a GPU;
# per-process thread counts are capped by [`compute_budget`](@ref), which the driver reports (and
# which `torque_survey` refuses to exceed silently). Process-level fan-out over configurations is NOT
# used here — `survey_configs` loops serially and lets the m-loop use the threads.
# Type-generic `{T}`; BigFloat inherits. Additive (EGD::sw-survey-parallel).

"""
    compute_budget(; nthreads = Threads.nthreads(), cap = 24) -> (threads, cap, used, ok, note)

Compute-usage self-check for a survey run: `used` = min(`nthreads`, `cap`) is the number of worker
threads the survey will actually keep busy, `ok` is false when the process was started with more
threads than the mission cap (24 CPUs per agent job, no GPU). Reported by [`torque_survey`](@ref) when
`verbose = true` and returned so callers can log it.
"""
function compute_budget(; nthreads::Integer = Threads.nthreads(), cap::Integer = 24)
    used = min(Int(nthreads), Int(cap))
    (threads = Int(nthreads), cap = Int(cap), used = used, ok = nthreads ≤ cap,
     note = "SpiralWaves survey: $(used) worker thread(s) of $(nthreads) available (mission cap $(cap) CPUs, no GPU)")
end

"""
    ModeTorque{T}

Reduction of one solved azimuthal mode: the total torque `T` (74), the Lindblad torque `T_L` (75), the
corotation torque `T_c` = Δ_c𝔗 − Δ_cF_J, the AMF at the two domain edges `F_in`, `F_out`, the grid size
`nodes`, the discrete residual of (64) and the wall time of that mode.
"""
struct ModeTorque{T<:Real}
    m::Int
    T::T
    T_L::T
    T_c::T
    F_in::T
    F_out::T
    nodes::Int
    residual::T
    seconds::Float64
end

"""
    TorqueSurvey{T}

Result of [`torque_survey`](@ref): the mode list `ms`, the per-mode reductions `modes`, the fixed-order
sums `T_total`, `T_L_total`, the optional per-mode radial `profiles` (`(r, F_J, density)` NamedTuples,
kept only with `keep_profiles = true` — they are large), the wall time and the compute budget used.
"""
struct TorqueSurvey{T<:Real}
    ms::Vector{Int}
    modes::Vector{ModeTorque{T}}
    T_total::T
    T_L_total::T
    profiles::Union{Nothing,Vector{NamedTuple{(:r, :F_J, :density),Tuple{Vector{T},Vector{T},Vector{T}}}}}
    seconds::Float64
    budget::NamedTuple{(:threads, :cap, :used, :ok, :note),Tuple{Int,Int,Int,Bool,String}}
end

"Per-mode torques `T_m` in the order of `s.ms` (the input of [`m_cut`](@ref))."
mode_torques(s::TorqueSurvey) = [mt.T for mt in s.modes]
"Per-mode Lindblad torques `T_{L,m}` in the order of `s.ms`."
mode_lindblad_torques(s::TorqueSurvey) = [mt.T_L for mt in s.modes]

# one mode, start to finish (pure: identical inputs ⇒ identical bits, on any thread)
function _survey_one(bg::DiscBackground{T}, m::Integer; grid_kwargs, solve_kwargs, halfwidth,
                     keep_profiles::Bool) where {T}
    t0 = time()
    g = build_grid(bg, m; grid_kwargs...)
    sol = solve_mode(bg, m; grid = g, threads = false, solve_kwargs...)
    F = amf(bg, sol)
    D = torque_density(bg, sol)
    sp = corotation_split(bg, sol; halfwidth = halfwidth, flux = F, density = D)
    mt = ModeTorque{T}(Int(m), sp.T, sp.T_L, sp.T_c, F[1], F[end], length(sol.r), sol.residual, time() - t0)
    prof = keep_profiles ? (r = sol.r, F_J = F, density = D) : nothing
    (mt, prof)
end

"""
    torque_survey(bg::DiscBackground, ms; threads = true, keep_profiles = false, verbose = false,
                  halfwidth = 3bg.σε, grid_kwargs = (;), solve_kwargs = (;)) -> TorqueSurvey

Solve every azimuthal mode in `ms` (each on its OWN adaptive grid, `grid_kwargs` forwarded to
[`build_grid`](@ref), `solve_kwargs` to [`solve_mode`](@ref)) and reduce each to a [`ModeTorque`](@ref);
sum the per-mode torques SERIALLY in the order of `ms`.

`threads = true` distributes the modes with `Threads.@threads :static`; each task owns everything it
touches (its grid, its matrices, its source modes — the inner source projection is forced serial to
avoid nested threading) and writes one disjoint slot, so the result is BITWISE identical to
`threads = false`. `verbose = true` prints the [`compute_budget`](@ref) note before the loop.
"""
function torque_survey(bg::DiscBackground{T}, ms::AbstractVector{<:Integer};
                       threads::Bool = true, keep_profiles::Bool = false, verbose::Bool = false,
                       halfwidth = 3 * bg.σε, grid_kwargs = NamedTuple(),
                       solve_kwargs = NamedTuple()) where {T}
    isempty(ms) && throw(ArgumentError("torque_survey: `ms` must be non-empty"))
    all(m -> m ≥ 1, ms) || throw(ArgumentError("torque_survey: every m must be ≥ 1"))
    budget = compute_budget()
    verbose && @info budget.note
    t0 = time()
    n = length(ms)
    mts = Vector{ModeTorque{T}}(undef, n)
    profs = keep_profiles ?
        Vector{NamedTuple{(:r, :F_J, :density),Tuple{Vector{T},Vector{T},Vector{T}}}}(undef, n) : nothing
    if threads
        Threads.@threads :static for i in 1:n
            mt, pr = _survey_one(bg, ms[i]; grid_kwargs = grid_kwargs, solve_kwargs = solve_kwargs,
                                 halfwidth = halfwidth, keep_profiles = keep_profiles)
            @inbounds mts[i] = mt
            keep_profiles && (@inbounds profs[i] = pr)
        end
    else
        for i in 1:n
            mt, pr = _survey_one(bg, ms[i]; grid_kwargs = grid_kwargs, solve_kwargs = solve_kwargs,
                                 halfwidth = halfwidth, keep_profiles = keep_profiles)
            @inbounds mts[i] = mt
            keep_profiles && (@inbounds profs[i] = pr)
        end
    end
    Ttot = zero(T); TLtot = zero(T)                       # FIXED-ORDER serial reduction
    for mt in mts
        Ttot += mt.T; TLtot += mt.T_L
    end
    TorqueSurvey{T}(collect(Int.(ms)), mts, Ttot, TLtot, profs, time() - t0, budget)
end

"""
    survey_configs(configs, ms; kwargs...) -> Vector{TorqueSurvey}

Run [`torque_survey`](@ref) for every disc configuration in `configs` (a vector of `DiscBackground`s,
or of NamedTuples of `DiscBackground` keyword arguments). The configurations are looped SERIALLY and
the threads are spent inside the m-loop, which keeps the process within the mission compute cap
(≤ 24 CPUs, no GPU, no process fan-out) and keeps every result bitwise reproducible.
"""
function survey_configs(configs::AbstractVector, ms::AbstractVector{<:Integer}; kwargs...)
    [torque_survey(cfg isa DiscBackground ? cfg : DiscBackground(; cfg...), ms; kwargs...) for cfg in configs]
end
