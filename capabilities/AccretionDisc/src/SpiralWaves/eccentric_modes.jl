# End-to-end eccentric response. Projection, emitted physics, radiation rows,
# midpoint assembly and banded LU remain with their existing owners.

"""
    EccentricModeSource

Callable table of complex covariant source jets at explicitly sampled radii.
No interpolation or circular projection is performed. `source(r)` rejects an
unsampled radius; construct a new table when changing the integration grid.
"""
struct EccentricModeSource{T,H}
    mode::HarmonicMode{T}
    r::Vector{T}
    jets::Vector{H}
end

function (s::EccentricModeSource)(r)
    i = searchsortedfirst(s.r,r)
    i <= length(s.r) && s.r[i] == r ||
        throw(ArgumentError("source radius was not sampled; rebuild for this grid"))
    s.jets[i]
end

"""
    mode_source(bank::SourceAzimuthBank, mode)

Project one mode from an existing all-azimuth bank. The bank radii must be
strictly ascending. The returned table owns its radii and projected jets and
does not retain the larger azimuth bank. Banks can be reused for other modes.
"""
function mode_source(bank::SourceAzimuthBank, mode::HarmonicMode)
    all(diff(bank.rs) .> 0) || throw(ArgumentError("source radii must be strictly ascending"))
    EccentricModeSource(mode,copy(bank.rs),vec(source_modes(bank,(mode,);derivatives=true)))
end

"""
    mode_source(context, perturber, grid; radial_batch_size=32, bank_options=(;))

Project the eccentric source at every node AND cell midpoint, using
`source_modes(source_azimuth_bank(...))`. Radial batches bound the retained
azimuth workspace; each bank also enforces its own memory cap. `bank_options`
forwards Nq, Nphi, blocksize and max_memory_bytes to the source owner. Source
quadrature convergence remains a separate obligation from BVP refinement.
"""
function mode_source(c::ModeContext, p::EccentricPerturber, grid::EccentricGrid;
                     radial_batch_size::Integer=32, bank_options::NamedTuple=(;))
    # Reuse the context owner's spin, precision and orbital-frequency checks.
    ModeContext(c.disc,p,c.mode;closure=c.closure,domain=c.domain)
    grid.context_key == mode_cache_key(c) || throw(ArgumentError("grid belongs to a different context"))
    radial_batch_size > 0 || throw(ArgumentError("radial_batch_size must be positive"))
    r = grid.r
    length(r)>=2 && first(r)==c.domain[1] && last(r)==c.domain[2] && all(diff(r).>0) ||
        throw(ArgumentError("grid must span the domain in strictly ascending order"))
    samples = sort!(unique!(vcat(r,[(r[i]+r[i+1])/2 for i in 1:length(r)-1])))
    options = merge((mmax=c.mode.m,),bank_options)
    chunks = map(collect(Iterators.partition(samples,radial_batch_size))) do radii
        bank = source_azimuth_bank(c.disc,p,radii;options...)
        mode_source(bank,c.mode).jets
    end
    EccentricModeSource(c.mode,samples,reduce(vcat,chunks))
end

function _eccentric_boundary_row(::Val{Adiabatic}, c::ModeContext, r; outward)
    x = master_inputs(c,r;hmodes=_zero_source(c),landau=false)
    fields = (;(n => getfield(x,n) for n in fieldnames(typeof(x)))...)
    h = sqrt(eps(typeof(x.enu)))*abs(real(x.omega_d))
    C = _context_operator(Val(c.closure),x)
    Cω = imag.(_context_operator(Val(c.closure),merge(fields,(omega_d=x.omega_d+im*h,))))/h
    last(outgoing_branch(C,Cω;outward))
end

_eccentric_boundary_row(::Val{LocallyIsothermalOptionA},c::ModeContext,r;outward) =
    li_boundary(c,r;outward).row

"""
    mode_boundary_rows(context)

Homogeneous radiation conditions at the physical inner/outer edges. Adiabatic
rows use the certified C-66 total-frequency group-velocity/decay selector on
the unregularized emitted operator; OPTION A uses `li_boundary`, including
its first thermal amplitude/impedance correction. These conditions retain their
exterior-truncation error; grid refinement alone cannot certify it away.
An edge at corotation or a degenerate turning point is rejected by the owners.
"""
function mode_boundary_rows(c::ModeContext{T}) where {T}
    map(((first(c.domain),false),(last(c.domain),true))) do (r,outward)
        (_eccentric_boundary_row(Val(c.closure),c,r;outward),zero(Complex{T}))
    end
end

"""
    solve_mode(context, perturber; grid=build_grid(context),
               radial_batch_size=32, bank_options=(;),
               boundary_rows=mode_boundary_rows(context), kwargs...)

End-to-end eccentric projection and response at omega=m*Omega_phi+k*Omega_r.
Returns the existing `EccentricModeSolution`, retaining a sampled source for
reconstruction. `kwargs` forwards landau, linsolve and cache to the existing
context solve; the default backend there is the certified generic banded LU.
The solution's `seconds` measures BVP work only, excluding source projection.
Supply explicit boundary_rows for an independently defined exterior model.
"""
function solve_mode(c::ModeContext,p::EccentricPerturber;
                    grid::EccentricGrid=build_grid(c), radial_batch_size::Integer=32,
                    bank_options::NamedTuple=(;), boundary_rows=mode_boundary_rows(c), kwargs...)
    source = mode_source(c,p,grid;radial_batch_size,bank_options)
    solve_mode(c;source,grid,boundary_rows,kwargs...)
end

"""
    solve_modes(context_builder, perturber, modes; bank=nothing,
                bank_options=(;), grid_builder=build_grid, transform=identity,
                boundary_builder=mode_boundary_rows, kwargs...)

Solve an ordered, nonempty mode collection with one caller-owned source bank.
`context_builder(mode)` supplies each context; `grid_builder(context)` supplies
its grid. By default, build a `SourceRadialBank` on the union of the context
domains, forwarding `bank_options` to `source_radial_bank`. Its Chebyshev source
interpolation uses the documented 1e-8 |ref| + 1e-12 U source budget, not a
bound on response error. Quadrature convergence remains the caller's duty.

An explicit `SourceAzimuthBank` instead uses the direct projector and must
contain every grid node and midpoint. Identical samples and quadrature give
bitwise-identical results to separate solves. An explicit bank must belong to
the supplied perturber instance; no global cache or circular dispatch is used.

`transform(solution)` runs immediately after each solve, so survey reductions
need not retain every full solution. Results preserve input order, including
duplicates. With an explicit bank, contexts are also constructed one at a time;
without one, their domain metadata is collected before building the bank.
`kwargs` forwards to the existing context solver. The single-mode
`solve_mode(context, perturber)` interface and defaults are unchanged.
"""
function solve_modes(context_builder, p::EccentricPerturber, modes;
                     bank=nothing, bank_options::NamedTuple=(;),
                     grid_builder=build_grid, transform=identity,
                     boundary_builder=mode_boundary_rows, kwargs...)
    selected = collect(modes)
    isempty(selected) && throw(ArgumentError("at least one mode is required"))
    checked_context = function(mode)
        c = context_builder(mode)
        c.mode == mode || throw(ArgumentError("context builder changed the requested mode"))
        ModeContext(c.disc,p,c.mode;closure=c.closure,domain=c.domain)
    end
    if bank === nothing
        contexts = map(checked_context,selected)
        domain = (minimum(c.domain[1] for c in contexts),
                  maximum(c.domain[2] for c in contexts))
        bank = source_radial_bank(first(contexts).disc,p,domain,
                                  unique(c.mode for c in contexts);bank_options...)
    else
        bank isa Union{SourceAzimuthBank,SourceRadialBank} ||
            throw(ArgumentError("expected an eccentric source bank"))
        bank.perturber === p || throw(ArgumentError("bank belongs to another perturber"))
        isempty(bank_options) || throw(ArgumentError("bank_options require bank=nothing"))
        contexts = Iterators.map(checked_context,selected)
    end
    map(contexts) do c
        source = mode_source(bank,c.mode)
        grid = grid_builder(c)
        rows = boundary_builder(c)
        transform(solve_mode(c;source,grid,boundary_rows=rows,kwargs...))
    end
end

"Reconstruct closure-tagged Eulerian fluid amplitudes at all solved grid nodes."
function reconstruct_fluid(sol::EccentricModeSolution)
    [reconstruct_fluid(sol.context,r,view(sol.y,:,i);hmodes=sol.source(r),
                       landau=sol.params.landau) for (i,r) in enumerate(sol.r)]
end
