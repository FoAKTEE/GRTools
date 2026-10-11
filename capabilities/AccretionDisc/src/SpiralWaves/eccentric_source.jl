# Leading equal-time G26 source, extended to a signed eccentric worldline.
# Algebra: derivation/src/eccentric_source.jl:29-53 (distance_jet, field_jet).
# Model limitations: leading singular field, fixed BL component transport,
# Plummer softening; no retarded, subleading or far-zone accuracy claim.
using NumericsCore: fft_radix2!, cheb_lobatto, bary_weights, bary_row
using Spacetime: worldline_rates, orbit_frequency

const ECCENTRIC_SOURCE_COMPONENTS = (:htt, :htphi, :hphiphi, :hrr, :hzz, :htr, :hrphi)
const ECCENTRIC_SOURCE_JETS = (:htt, :htt_r, :htphi, :htphi_r, :hphiphi, :hphiphi_r,
                              :hrr, :hrr_r, :hzz, :hzz_r, :htr, :htr_r, :hrphi, :hrphi_r)

"Worldline tensor and positive-square distance factors; independent of field radius."
struct EccentricSourceEvent{T<:Real}
    rp::T
    soft::T
    grr::T
    gtan::T
    vr::T
    vtan::T
    tensor::NTuple{7,T}
end

function _eccentric_source_event(rp, soft, mass, gtt, gtp, gpp, grr, ut, ur, uphi)
    vt, vr, vp = gtt*ut+gtp*uphi, grr*ur, gtp*ut+gpp*uphi
    tensor = (2mass*(gtt+2vt^2), 2mass*(gtp+2vt*vp), 2mass*(gpp+2vp^2),
              2mass*(grr+2vr^2), 2mass, 4mass*vt*vr, 4mass*vr*vp)
    EccentricSourceEvent(rp,soft,grr,gpp/rp^2,vr,vp/rp,tensor)
end

function _eccentric_source_event(p::EccentricPerturber{T}, s) where {T}
    p.softening isa FixedSoftening ||
        throw(ArgumentError("eccentric GR source currently requires FixedSoftening"))
    g = kerr_equatorial_potentials(p.orbit.a,s.r)
    gpp, grr = exp(2g.ψ), exp(2g.μ)
    _eccentric_source_event(s.r,T(p.softening.eps),p.mass,
        -exp(2g.ν)+gpp*g.ω̃^2,-gpp*g.ω̃,gpp,grr,s.ut,s.ur,s.uphi)
end

# Positive squares avoid cancellation near the worldline; the velocity square
# retains the odd cross term 2 u_r u_phi dR dPhi. No azimuth parity is assumed.
@inline function _eccentric_distance(K::EccentricSourceEvent, r, psi)
    s, c = sincos(psi)
    _eccentric_distance(K,r,s,c,sin(psi/2)^2)
end

@inline function _eccentric_distance(K::EccentricSourceEvent, r, s, c, halfsin2)
    x, y = (r-K.rp)-2r*halfsin2, r*s
    v = K.vr*x+K.vtan*y
    vr = K.vr*c+K.vtan*s
    s2 = K.soft^2+K.grr*x^2+K.gtan*y^2+v^2
    s2_r = 2*(K.grr*x*c+K.gtan*y*s+v*vr)
    (;s2,s2_r)
end

@inline function _eccentric_kernel(K, r, psi)
    _eccentric_kernel(_eccentric_distance(K,r,psi))
end

@inline function _eccentric_kernel(d)
    I = inv(sqrt(d.s2))
    (I,-d.s2_r*I^3/2)
end

@inline _eccentric_source_jet(K, I, Ir) = NamedTuple{ECCENTRIC_SOURCE_JETS}(
    ntuple(i -> K.tensor[(i+1)÷2]*(isodd(i) ? I : Ir), Val(14)))

"""
    singular_field(perturber::EccentricPerturber, r, psi, q)

Covariant leading singular field and analytic first radial derivatives, at
relative azimuth psi and radial Mino phase q. All seven components, including
htr/hrphi, follow the eccentric path even at e=0. FixedSoftening only.
"""
function singular_field(p::EccentricPerturber, r, psi, q)
    K = _eccentric_source_event(p,worldline_rates(p.orbit,q))
    _eccentric_source_jet(K,_eccentric_kernel(K,r,psi)...)
end

"""
    strip_halfwidth(event::EccentricSourceEvent, r)

Conservative zero-free half-strip for the eccentric distance. Let Q be the
positive 2x2 spatial quadratic in (r cos(psi)-rp,r sin(psi)), with eigenvalues
lo,hi. Set D=soft²+lo*(r-rp)²/2 and delta=asinh(sqrt(D/(4hi*r²))).
For z=x+iy, Young's inequality gives Re(s²(z)) >= D - hi*r² *
((cosh(y)-1)²+sinh(y)²) >= D/2 for |y|<=delta. Thus both the kernel and
its radial derivative are analytic there. This is a lower bound, not the
nearest zero; at e=0 it is checked against source.jl's exact circular strip.
"""
function strip_halfwidth(K::EccentricSourceEvent, r)
    lo,hi = _eccentric_spatial_eigenvalues(K)
    D = K.soft^2+lo*(r-K.rp)^2/2
    D > 0 || throw(DomainError(r,"unsoftened source on the real azimuth axis"))
    asinh(sqrt(D/hi)/ (2r))
end

function _eccentric_spatial_eigenvalues(K)
    a, b, c = K.grr+K.vr^2, K.vr*K.vtan, K.gtan+K.vtan^2
    hi = (a+c+hypot(a-c,2b))/2
    lo = (K.grr*K.gtan+K.grr*K.vtan^2+K.gtan*K.vr^2)/hi
    lo,hi
end

function _source_orbit_data(p::EccentricPerturber{T}, Nq) where {T}
    q = [2T(pi)*j/Nq for j in 0:Nq-1]
    states = worldline_rates.(Ref(p.orbit),q)
    origin = first(states)
    (;q,events=_eccentric_source_event.(Ref(p),states),
      tw=[s.dt_dlambda/(Nq*p.orbit.Γ) for s in states],
      dt=[s.Δtr-origin.Δtr for s in states],dp=[s.Δφr-origin.Δφr for s in states])
end

function _eccentric_nphi(delta, mmax, ::Type{T}) where {T}
    # Same precision/strip safety factors as the circular evaluator, with a
    # conservative strip and power-of-two rounding. Refinement is still a gate.
    need = 2mmax+ceil(Int,(-log(eps(T))+max(zero(T),-log(delta))+8)/delta)
    nextpow(2,max(256,8mmax+64,need))
end

"""
Caller-owned, factorised azimuth bank. `kernels[m+1,qindex,rindex,d]` holds
FFT/N of 1/s or its radial derivative (d=1,2); `events[qindex].tensor` expands
it to all seven components. FFTs are shared across components and all m.
`nphi`/`strips` are per radial block; `blocksize` identifies those blocks.
No bank/cache is global. Numeric workspace bound excludes Julia/JIT overhead.
"""
struct SourceAzimuthBank{T,P}
    perturber::P
    rs::Vector{T}
    q::Vector{T}
    time_weights::Vector{T}
    delta_t::Vector{T}
    delta_phi::Vector{T}
    events::Vector{EccentricSourceEvent{T}}
    kernels::Array{Complex{T},4}
    nphi::Vector{Int}
    strips::Vector{T}
    blocksize::Int
    workspace_bound_bytes::Int
end

"""
    source_azimuth_bank(model, perturber, rs; Nq=1024, Nphi=nothing,
                        mmax=256, blocksize=4, max_memory_bytes=512*1024^2)

Endpoint-exclusive uniform Mino phases; source orbit/rates evaluated once.
Each radial block uses the minimum certified strip over its radii and all
sampled events. Nphi is an optional power-of-two LOWER FLOOR (never a way to
bypass strip adaptation). Double Nphi beyond max(bank.nphi) to audit azimuth
convergence; rebuild with doubled Nq independently to audit time convergence.
Stores modes 0:mmax. Only FixedSoftening is supported; no implicit scale-height
policy is chosen. Spin and precision must match; rs may be convertible to T.

The Float64 numeric storage bound includes retained arrays, a single maximum
block workspace and the FFT twiddle table; a small object-overhead allowance
is included. Oversized requests fail before allocating the large arrays.
Use successive calls on subsets of rs for bounded radial streaming. BigFloat
inherits ambient precision for small unit tests; the byte cap is a Float64
production guarantee, not a BigFloat heap bound.
"""
function source_azimuth_bank(d::DiscModel{T}, p::EccentricPerturber, rs;
        Nq::Int=1024, Nphi::Union{Nothing,Int}=nothing, mmax::Int=256,
        blocksize::Int=4, max_memory_bytes::Int=512*1024^2,
        _orbit_data=nothing) where {T}
    typeof(p.orbit.a) === T && p.orbit.a == d.χ ||
        throw(ArgumentError("disc and orbit spin/type must match"))
    Nq >= 2 && mmax >= 0 && blocksize > 0 && max_memory_bytes > 0 ||
        throw(ArgumentError("require Nq>=2, mmax>=0, positive blocksize and memory cap"))
    Nphi === nothing || (ispow2(Nphi) && Nphi>=2) ||
        throw(ArgumentError("Nphi must be a power of two >=2"))
    rr = T.(collect(rs))
    !isempty(rr) && all(r -> isfinite(r) && r>0,rr) ||
        throw(ArgumentError("radii must be finite, positive and nonempty"))
    n = length(rr)
    # BigInt arithmetic here prevents a huge requested bank overflowing its cap.
    retained = big(2)*sizeof(T)*(mmax+big(1))*Nq*n*2 + big(17)*sizeof(T)*Nq + sizeof(T)*big(n)
    retained+4096 <= max_memory_bytes || throw(ArgumentError("source bank exceeds memory cap"))
    data = _orbit_data === nothing ? _source_orbit_data(p,Nq) : _orbit_data
    q,events,tw,dt,dp = data.q,data.events,data.tw,data.dt,data.dp
    blocks = [i:min(i+blocksize-1,n) for i in 1:blocksize:n]
    strips = [minimum(strip_halfwidth(K,r) for K in events for r in rr[b]) for b in blocks]
    nphi = [max(something(Nphi,2),_eccentric_nphi(delta,mmax,T)) for delta in strips]
    peak = maximum(nphi)
    # Batch orbit phases too: NumericsCore then shares the FFT twiddle table
    # across those phases. Reduce the batch before exceeding the byte cap.
    # Includes states temporarily retained during construction (11 reals/event).
    fixed = retained+big(sizeof(T))*(4big(peak)+11big(Nq))+4096+64big(length(blocks))
    perphase = big(sizeof(T))*4big(peak)*min(blocksize,n)
    qbatch = Int(min(Nq,8,max(big(1),(max_memory_bytes-fixed)÷perphase)))
    bound = fixed+qbatch*perphase
    bound <= max_memory_bytes || throw(ArgumentError("strip-adapted FFT workspace exceeds memory cap"))
    kernels = Array{Complex{T}}(undef,mmax+1,Nq,n,2)
    # One buffer survives all blocks, so different block lengths do not leave
    # multiple large dead buffers awaiting collection.
    scratch = Matrix{Complex{T}}(undef,peak,2min(blocksize,n)*qbatch)
    nodes = Matrix{T}(undef,peak,3)
    for (ib,b) in enumerate(blocks)
        N = nphi[ib]
        for l in 0:N-1
            psi = 2T(pi)*l/N
            nodes[l+1,1],nodes[l+1,2] = sincos(psi)
            nodes[l+1,3] = sin(psi/2)^2
        end
        for qb in Iterators.partition(eachindex(q),qbatch)
            work = view(scratch,1:N,1:2length(b)*length(qb))
            for (jq,iq) in enumerate(qb)
                K = events[iq]
                for (j,ir) in enumerate(b), l in 0:N-1
                    I, Ir = _eccentric_kernel(_eccentric_distance(K,rr[ir],
                        nodes[l+1,1],nodes[l+1,2],nodes[l+1,3]))
                    column = 2*((jq-1)*length(b)+j-1)
                    work[l+1,column+1], work[l+1,column+2] = I, Ir
                end
            end
            fft_radix2!(work;dims=1)
            for (jq,iq) in enumerate(qb), (j,ir) in enumerate(b), derivative in 1:2, m in 0:mmax
                column = 2*((jq-1)*length(b)+j-1)+derivative
                kernels[m+1,iq,ir,derivative] = work[m+1,column]/N
            end
        end
    end
    SourceAzimuthBank(p,rr,q,tw,dt,dp,events,kernels,nphi,strips,blocksize,Int(bound))
end

"A Chebyshev panel with an analytic interpolation truncation bound for each jet."
struct SourceRadialPanel{T}
    lo::T
    hi::T
    indices::UnitRange{Int}
    nodes::Vector{T}
    weights::Vector{T}
    error_bound::NTuple{14,T}
end

struct SourceRadialGrid{T}
    panels::Vector{SourceRadialPanel{T}}
    edges::Vector{T}
    rs::Vector{T}
end

# On a complex-r ellipse z=u+iv, Re(s²) >= soft² + lambda_lo*
# dist(u,[rperi,rapo])² - lambda_hi*v². The rho=2 Bernstein ellipse has
# real semi-axis 5h/4 and imaginary semi-axis 3h/4. Require its subtraction
# <= 3/4 of the positive part. Hence it contains no kernel branch point.
# |I| <= B^(-1/2). For v=(z cos(psi)-rp,z sin(psi)) and n=(cos,sin),
# Cauchy-Schwarz in Q gives |n'Qv|² <= hi*v'Qv, where the first prime is
# transpose and v' in the latter product is conjugate transpose. Moreover
# v'Qv = Re(s²)-soft²+2*A*Im(z)² <= |s²|+2*hi*Im(z)².
# Thus |Ir| <= sqrt(hi)*sqrt(1+2*hi*Im(z)²/B)/B, without a large |z| factor.
# Chebyshev-Lobatto degree n error <= 4*M*rho^(-n)/(rho-1).
# Tensor/phase/time factors are independent of r, so their L1 sum propagates
# these bounds to every m,k of a fixed discrete orbit/azimuth projection or
# exact azimuth integral. Changes of Nphi between nodes contribute the
# separate azimuth-quadrature error, not this interpolation truncation bound.
function _source_radial_grid(p::EccentricPerturber{T}, domain, data) where {T}
    left,right = T.(domain)
    isfinite(left) && isfinite(right) && 0<left<right ||
        throw(ArgumentError("source domain must have finite 0<lo<hi"))
    soft = minimum(K.soft for K in data.events)
    soft>0 || throw(ArgumentError("radial interpolation requires positive fixed softening"))
    eigs = _eccentric_spatial_eigenvalues.(data.events)
    low,high = minimum(first,eigs),maximum(last,eigs)
    peri,apo = extrema(K.rp for K in data.events)
    amplitude = ntuple(c -> sum(abs(data.tw[i]*data.events[i].tensor[c]) for i in eachindex(data.q)),7)
    units = _source_units(p,true)
    panels,edges,rs = SourceRadialPanel{T}[],T[left],T[]
    while left<right
        distance = max(peri-left,left-apo,zero(T))
        # Keep the ellipse's real part positive, including very small domains.
        width = min(right-left,4left,2sqrt((soft^2+low*distance^2)/high))
        local bound,major,imaginary
        while true
            half = width/2
            center = left+half
            major,imaginary = 5half/4,3half/4
            distance = max(peri-(center+major),(center-major)-apo,zero(T))
            positive = soft^2+low*distance^2
            bound = positive-high*imaginary^2
            bound>=positive/4 && break
            width /= 2
        end
        stop = min(right,left+width)
        stop>left || throw(ArgumentError("radial panels cannot be resolved in this precision"))
        I = inv(sqrt(bound))
        Ir = sqrt(high)*sqrt(1+2high*imaginary^2/bound)/bound
        M = ntuple(j -> amplitude[(j+1)÷2]*(iseven(j) ? Ir : I),14)
        ratio = maximum(4M[j]/(T(2.5e-13)*units[j]) for j in 1:14)
        degree = max(8,ceil(Int,log2(max(one(T),ratio))))
        nodes = cheb_lobatto(T,degree;lo=left,hi=stop)
        nodes[1],nodes[end] = stop,left
        indices = length(rs)+1:length(rs)+length(nodes)
        weights = bary_weights((nodes .- (left+stop)/2) ./ ((stop-left)/2))
        weights ./= maximum(abs,weights)
        push!(panels,SourceRadialPanel(left,stop,indices,nodes,weights,
            ntuple(j -> 4M[j]*T(2)^(-degree),14)))
        append!(rs,nodes)
        push!(edges,stop)
        left = stop
    end
    SourceRadialGrid(panels,edges,rs)
end

"Caller-owned projected radial bank; all modes and both disc closures share panels."
struct SourceRadialBank{T,P,H}
    perturber::P
    grid::SourceRadialGrid{T}
    modes::Vector{HarmonicMode{T}}
    jets::Matrix{H}
end

"""
    source_radial_bank(disc, perturber, domain, modes; Nq=1024, Nphi=nothing,
                      method=:fft, radial_batch_size=4, max_memory_bytes=512*1024^2)

Build a source-adapted Chebyshev bank, dense across the peri/apo interval on
the softening scale and geometric in the tails. Every radial sample builds
ONE azimuth bank for all requested modes. Orbit samples are shared across
radial batches; only projected jets survive. No BVP grid enters construction.
Use `mode_source(bank, mode)` as the `source` of the existing context solver.
The callable supplies both field and analytic radial derivative interpolants;
the derivative is interpolated independently, never obtained by differencing.

Each panel's rho=2 ellipse bounds radial truncation below 2.5e-13 U, reserving
the rest of the design budget for quadrature/roundoff. This analyticity bound
is for the discrete source; it does NOT certify Nq/Nphi convergence or a BVP
response error. Recorded Float64 parity remains the floating-point gate.

The cap covers retained jets plus the azimuth batch's numeric workspace,
with 64 KiB reserved for object overhead. FFT temporary/output storage is
also reserved. Very large all-mode tables require an explicit larger cap;
never silently split modes and rebuild a bank per mode. Storage and work are
independent of the number of nodes on each consumer's BVP grid.
"""
function source_radial_bank(d::DiscModel{T},p::EccentricPerturber,domain,modes;
        Nq::Int=1024,Nphi::Union{Nothing,Int}=nothing,method::Symbol=:fft,
        radial_batch_size::Int=4,max_memory_bytes::Int=512*1024^2) where {T}
    typeof(p.orbit.a) === T && p.orbit.a == d.χ || throw(ArgumentError("disc and orbit spin/type must match"))
    Nq>=2 && radial_batch_size>0 && max_memory_bytes>0 || throw(ArgumentError("invalid radial bank sizes"))
    selected = collect(modes)
    !isempty(selected) && all(mode -> mode isa HarmonicMode{T} && mode.m>=0,selected) ||
        throw(ArgumentError("nonempty matching-precision modes with m>=0 required"))
    data = _source_orbit_data(p,Nq)
    grid = _source_radial_grid(p,domain,data)
    proto = NamedTuple{ECCENTRIC_SOURCE_JETS}(ntuple(_->zero(Complex{T}),14))
    retained = big(sizeof(proto))*length(grid.rs)*length(selected)+Base.summarysize(grid)+65536
    kmax = maximum(abs(mode.k) for mode in selected)
    # Nq/order dependent FFT reserve, including moment corrections and weights.
    fftN = nextpow(2,max(Nq,2kmax+1,2))
    # x<=pi/2 at this FFT size: 32 moments is conservative for Float64.
    largest_group = maximum(count(mode -> mode.m==m,selected) for m in unique(mode.m for mode in selected))
    reserve = big(sizeof(T))*(4big(fftN)*max(32,precision(T))+(8+2big(largest_group))*Nq) +
              big(sizeof(proto))*radial_batch_size*length(selected)
    available = max_memory_bytes-retained-reserve
    available>65536 || throw(ArgumentError("projected radial source bank exceeds memory cap"))
    jets = Matrix{typeof(proto)}(undef,length(grid.rs),length(selected))
    mmax = maximum(mode.m for mode in selected)
    for indices in Iterators.partition(eachindex(grid.rs),radial_batch_size)
        bank = source_azimuth_bank(d,p,grid.rs[indices];Nq,Nphi,mmax,
            max_memory_bytes=Int(available),_orbit_data=data)
        jets[indices,:] = source_modes(bank,selected;method)
    end
    SourceRadialBank(p,grid,selected,jets)
end

struct RadialModeSource{B}
    bank::B
    column::Int
end

"Select a mode without copying its shared radial bank; exact mode identity is checked."
function mode_source(b::SourceRadialBank{T},mode::HarmonicMode) where {T}
    mode isa HarmonicMode{T} || throw(ArgumentError("mode and radial bank precision types must match"))
    column = findfirst(==(mode),b.modes)
    column === nothing && throw(ArgumentError("mode is absent or has a different frequency"))
    RadialModeSource(b,column)
end

function _source_radial_row(grid,r)
    isfinite(r) && first(grid.edges)<=r<=last(grid.edges) || throw(DomainError(r,"source radius outside sampled domain"))
    index = min(searchsortedlast(grid.edges,r),length(grid.panels))
    panel = grid.panels[index]
    panel.indices,bary_row(panel.nodes,panel.weights,r)
end

@inline function _source_radial_value(b,indices,row,column)
    values = ntuple(Val(14)) do j
        z = zero(eltype(row))+0im
        @inbounds for i in eachindex(row)
            z += row[i]*b.jets[indices[i],column][j]
        end
        z
    end
    NamedTuple{ECCENTRIC_SOURCE_JETS}(values)
end

function (s::RadialModeSource)(r)
    indices,row = _source_radial_row(s.bank.grid,r)
    _source_radial_value(s.bank,indices,row,s.column)
end

"Interpolate selected modes on arbitrary radii; one barycentric row per radius."
function source_modes(b::SourceRadialBank,modes,rs)
    sources = [mode_source(b,mode) for mode in modes]
    radii = collect(rs)
    out = Matrix{eltype(b.jets)}(undef,length(radii),length(sources))
    for (i,r) in enumerate(radii)
        indices,row = _source_radial_row(b.grid,r)
        for (j,s) in enumerate(sources)
            out[i,j] = _source_radial_value(b,indices,row,s.column)
        end
    end
    out
end

"Expand one bank entry to covariant component jets (indices refer to stored q and rs)."
function singular_field_mode(b::SourceAzimuthBank, m::Integer, ri::Integer, qi::Integer)
    _eccentric_source_jet(b.events[qi],b.kernels[m+1,qi,ri,1],b.kernels[m+1,qi,ri,2])
end

"""
    singular_field_mode(perturber, m, r, q; N=1024)

Direct full-period trapezoid oracle at one event. Uses exp(-im psi), evaluates
both azimuth halves, and applies compensated sums. N is explicit; refine it
for accuracy. Production uses source_azimuth_bank and shares its FFTs.
"""
function singular_field_mode(p::EccentricPerturber{T}, m::Integer, r, q; N::Int=1024) where {T}
    N>=2 && abs(m)<N÷2 || throw(ArgumentError("oracle requires N>=2 and |m|<N/2"))
    K = _eccentric_source_event(p,worldline_rates(p.orbit,q))
    z = zero(Complex{T})
    a,c,ar,cr = z,z,z,z
    for j in 0:N-1
        psi = 2T(pi)*j/N
        I, Ir = _eccentric_kernel(K,r,psi)
        phase = cis(-m*psi)
        a,c = _kbn_add(a,c,I*phase)
        ar,cr = _kbn_add(ar,cr,Ir*phase)
    end
    _eccentric_source_jet(K,(a+c)/N,(ar+cr)/N)
end

# A component of a coordinate-time projection, using the circular source's
# compensated adder separately on real/imaginary parts. Componentwise ordering
# matters when one part nearly cancels while the other stays large.
function _eccentric_harmonic_sum(b, m, weights, ir, component, derivative)
    T = eltype(b.q)
    ar,cr,ai,ci = zero(T),zero(T),zero(T),zero(T)
    @inbounds for iq in eachindex(b.q)
        z = weights[iq]*b.events[iq].tensor[component]*b.kernels[m+1,iq,ir,derivative]
        ar,cr = _kbn_add(ar,cr,real(z))
        ai,ci = _kbn_add(ai,ci,imag(z))
    end
    complex(ar+cr,ai+ci)
end

"""
    source_modes(bank::SourceAzimuthBank, modes; derivatives=true)

Return a matrix indexed `[radius_index, mode_index]` of covariant component
NamedTuples. With derivatives=true each tuple includes all seven first radial
jets (same names as singular_field); otherwise it has seven field components.
Every k shares the stored orbit, tensor and azimuth bank. Frequencies must be
m*Omega_phi+k*Omega_r in the bank's floating type, and 0<=m<=bank.mmax.
No sector is silently dropped, including axisymmetric or secular modes.

Implements Spacetime's C-HARM convention (Geodesics/harmonics.jl:73-102):
sum[(dt/dlambda)/(Nq*Gamma) H_m exp(i(k*q+omega*Delta_t-m*Delta_phi))].
Real and imaginary sums are separately compensated. Returned complex Fourier
coefficients are not real-amplitude cosine coefficients: FR22's Newtonian
potential amplitude is -real(htt), since Phi=-h_tt/2 and the reality-pair
amplitude is 2*Phi_hat (E-ECC-2). No enforced reality projection is applied.

Nq is caller-controlled: compare independent rebuilt banks at Nq,2Nq,4Nq.
A single bank carries no convergence claim. No circular-source dispatch occurs.
"""
function _source_modes_direct(b::SourceAzimuthBank{T}, modes; derivatives::Bool=true) where {T}
    selected = collect(modes)
    names = derivatives ? ECCENTRIC_SOURCE_JETS : ECCENTRIC_SOURCE_COMPONENTS
    proto = NamedTuple{names}(ntuple(_ -> zero(Complex{T}),length(names)))
    out = Matrix{typeof(proto)}(undef,length(b.rs),length(selected))
    weights = Vector{Complex{T}}(undef,length(b.q))
    for (imode,mode) in enumerate(selected)
        omega = _source_mode_frequency(b,mode)
        for iq in eachindex(b.q)
            phase = mode.k*b.q[iq]+omega*b.delta_t[iq]-mode.m*b.delta_phi[iq]
            weights[iq] = b.time_weights[iq]*cis(phase)
        end
        for ir in eachindex(b.rs)
            values = ntuple(length(names)) do j
                component = derivatives ? (j+1)÷2 : j
                derivative = derivatives && iseven(j) ? 2 : 1
                _eccentric_harmonic_sum(b,mode.m,weights,ir,component,derivative)
            end
            out[ir,imode] = NamedTuple{names}(values)
        end
    end
    out
end

function _source_mode_frequency(b::SourceAzimuthBank{T}, mode) where {T}
    mode isa HarmonicMode || throw(ArgumentError("source modes must be HarmonicMode values"))
    0<=mode.m<size(b.kernels,1) || throw(ArgumentError("requested m is outside the azimuth bank"))
    typeof(mode.omega) === T || throw(ArgumentError("mode and source bank precision types must match"))
    omega = orbit_frequency(b.perturber.orbit,mode.m,mode.k)
    scale = max(abs(omega),abs(mode.m*b.perturber.orbit.Ωφ),abs(mode.k*b.perturber.orbit.Ωr))
    abs(mode.omega-omega)<=64eps(T)*scale ||
        throw(ArgumentError("mode frequency does not match m*Omega_phi+k*Omega_r"))
    omega
end

"Optional caller-owned operation counters, independent of wall-clock load."
Base.@kwdef mutable struct SourceProjectionWork
    spread_terms::Int = 0
    fft_butterflies::Int = 0
    direct_terms::Int = 0
end

# Component units used by the existing source parity gate (DESIGN section 4).
function _source_units(p, derivatives)
    r = p.orbit.p
    scales = (one(r),r,r^2,one(r),one(r),one(r),r)
    ntuple(derivatives ? 14 : 7) do j
        p.mass*scales[derivatives ? (j+1)÷2 : j]/r/(derivatives && iseven(j) ? r : one(r))
    end
end

# q is uniform in Mino time, NOT coordinate time. The exact k phase is
# k*theta, theta=q+Omega_r*Delta_t. Spread Taylor moments about the nearest
# uniform node: exp(i*k*delta)=sum_j (i*kmax*delta)^j/j! * (k/kmax)^j.
# The absolute omitted tail is <= exp(x)*x^J/J!, x=max|k*delta|.
function _source_k_plan(b::SourceAzimuthBank{T}, kmax) where {T}
    nq = length(b.q)
    if ispow2(nq) && 2kmax<nq && all(iszero,b.delta_t) &&
            all(b.q[i]==2T(pi)*(i-1)/nq for i in eachindex(b.q))
        return (;N=nq,indices=collect(1:nq),factors=ones(Complex{T},nq,1),order=1,tail=zero(T))
    end
    # Only the requested band must fit the FFT, not every input orbit node.
    # A modest oversampling keeps the Taylor expansion short at large Nq.
    N = nextpow(2,max(min(nq,8max(kmax,1)),2kmax+1,2))
    theta = b.q .+ b.perturber.orbit.Ωr .* b.delta_t
    nearest = round.(Int,theta .* (N/(2T(pi))))
    delta = theta .- nearest .* (2T(pi)/N)
    x = maximum(abs,delta)*kmax
    order,term = 1,x
    while exp(x)*term > eps(T)/16
        order += 1
        term *= x/order
    end
    factors = ones(Complex{T},nq,order)
    for power in 2:order, iq in 1:nq
        factors[iq,power] = factors[iq,power-1]*(im*kmax*delta[iq]/(power-1))
    end
    (;N,indices=mod.(nearest,N).+1,factors,order,tail=exp(x)*term)
end

"""
    source_modes(bank, modes; derivatives=true, method=:direct, work=nothing)

`:direct` preserves the componentwise compensated reference reduction.
`:fft` evaluates the SAME discrete coordinate-time sum with a nearest-lattice
Taylor FFT (one FFT per component for uniform coordinate-time nodes; a small
batch of moments otherwise). It does not replace the orbit quadrature.
Modes may be signed, repeated, unordered or sparse. No Nyquist aliasing is
introduced: the FFT length is at least `2maximum(abs(k))+1`.

Each output is screened using its input L1 norm, Taylor remainder and a
conservative floating-point allowance. If that estimate exceeds one quarter
of `1e-8*abs(h)+1e-12*U`, use the original compensated sum for that component.
This is a roundoff screen, not an interval-arithmetic certificate. In
particular an ordinary FFT cannot promise exact recovery of 1e16+1-1e16.
The direct fallback retains that guarantee. Optional `SourceProjectionWork`
counts actual spread terms, FFT butterflies and fallback terms.
"""
function source_modes(b::SourceAzimuthBank, modes; derivatives::Bool=true,
                      method::Symbol=:direct, work=nothing)
    method in (:direct,:fft) || throw(ArgumentError("method must be :direct or :fft"))
    selected = collect(modes)
    if method == :direct || isempty(selected)
        work === nothing || (work.direct_terms += length(b.q)*length(b.rs)*length(selected)*(derivatives ? 14 : 7))
        return _source_modes_direct(b,selected;derivatives)
    end
    _source_modes_fft(b,selected,derivatives,work)
end

function _source_modes_fft(b::SourceAzimuthBank{T}, selected, derivatives, counts) where {T}
    foreach(mode -> _source_mode_frequency(b,mode),selected)
    names = derivatives ? ECCENTRIC_SOURCE_JETS : ECCENTRIC_SOURCE_COMPONENTS
    nc = length(names)
    proto = NamedTuple{names}(ntuple(_ -> zero(Complex{T}),nc))
    out = Matrix{typeof(proto)}(undef,length(b.rs),length(selected))
    units = _source_units(b.perturber,derivatives)
    dtmax = maximum(abs,b.delta_t)
    mphase = abs(b.perturber.orbit.Ωφ)*dtmax+maximum(abs,b.delta_phi)
    kphase = maximum(abs,b.q)+abs(b.perturber.orbit.Ωr)*dtmax
    planned_kmax = maximum(abs(mode.k) for mode in selected)
    plan = _source_k_plan(b,planned_kmax)
    for m in unique(mode.m for mode in selected)
        columns = findall(mode -> mode.m==m,selected)
        kmax = maximum(abs(selected[j].k) for j in columns)
        if kmax != planned_kmax
            plan = _source_k_plan(b,kmax)
            planned_kmax = kmax
        end
        moments = zeros(Complex{T},plan.N,plan.order)
        corrections = similar(moments)
        values = Matrix{Complex{T}}(undef,nc,length(columns))
        weights = b.time_weights .* cis.(m .* (b.perturber.orbit.Ωφ .* b.delta_t .- b.delta_phi))
        fallback_weights = Dict{Int,Vector{Complex{T}}}()
        for ir in eachindex(b.rs)
            for j in 1:nc
                component,derivative = derivatives ? ((j+1)÷2,iseven(j) ? 2 : 1) : (j,1)
                fill!(moments,0); fill!(corrections,0)
                norm1 = zero(T)
                @inbounds for iq in eachindex(b.q)
                    z = weights[iq]*b.events[iq].tensor[component]*b.kernels[m+1,iq,ir,derivative]
                    norm1 += abs(z)
                    row = plan.indices[iq]
                    for power in 1:plan.order
                        a,c = moments[row,power],corrections[row,power]
                        term = z*plan.factors[iq,power]
                        ar,cr = _kbn_add(real(a),real(c),real(term))
                        ai,ci = _kbn_add(imag(a),imag(c),imag(term))
                        moments[row,power],corrections[row,power] = complex(ar,ai),complex(cr,ci)
                    end
                end
                moments .+= corrections
                fft_radix2!(moments;dims=1)
                # Include actual phase magnitudes: frequency ratios can grow
                # near the separatrix, so |m|+|k| alone is not a phase bound.
                allowance = norm1*(plan.tail+eps(T)*(32*(log2(plan.N)+plan.order)+
                    8*(abs(m)*mphase+kmax*kphase)))
                for (l,col) in enumerate(columns)
                    mode = selected[col]
                    row = mod(-mode.k,plan.N)+1
                    ratio = iszero(kmax) ? zero(T) : T(mode.k)/kmax
                    value = moments[row,end]
                    for power in plan.order-1:-1:1
                        value = value*ratio+moments[row,power]
                    end
                    if allowance > (T(1e-8)*abs(value)+T(1e-12)*units[j])/4
                        direct_weights = get!(fallback_weights,mode.k) do
                            omega = orbit_frequency(b.perturber.orbit,m,mode.k)
                            [b.time_weights[iq]*cis(mode.k*b.q[iq]+omega*b.delta_t[iq]-m*b.delta_phi[iq]) for iq in eachindex(b.q)]
                        end
                        value = _eccentric_harmonic_sum(b,m,direct_weights,ir,component,derivative)
                        counts === nothing || (counts.direct_terms += length(b.q))
                    end
                    values[j,l] = value
                end
            end
            for (l,col) in enumerate(columns)
                out[ir,col] = NamedTuple{names}(ntuple(j -> values[j,l],nc))
            end
        end
        if counts !== nothing
            counts.spread_terms += nc*length(b.rs)*length(b.q)*plan.order
            counts.fft_butterflies += nc*length(b.rs)*plan.order*(plan.N÷2)*trailing_zeros(plan.N)
        end
    end
    out
end
