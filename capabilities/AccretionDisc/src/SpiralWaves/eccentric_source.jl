# Leading equal-time G26 source, extended to a signed eccentric worldline.
# Algebra: derivation/src/eccentric_source.jl:29-53 (distance_jet, field_jet).
# Model limitations: leading singular field, fixed BL component transport,
# Plummer softening; no retarded, subleading or far-zone accuracy claim.
using NumericsCore: fft_radix2!
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
    x, y = (r-K.rp)-2r*sin(psi/2)^2, r*s
    v = K.vr*x+K.vtan*y
    vr = K.vr*c+K.vtan*s
    s2 = K.soft^2+K.grr*x^2+K.gtan*y^2+v^2
    s2_r = 2*(K.grr*x*c+K.gtan*y*s+v*vr)
    (;s2,s2_r)
end

@inline function _eccentric_kernel(K, r, psi)
    d = _eccentric_distance(K,r,psi)
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
    a, b, c = K.grr+K.vr^2, K.vr*K.vtan, K.gtan+K.vtan^2
    hi = (a+c+hypot(a-c,2b))/2
    lo = (K.grr*K.gtan+K.grr*K.vtan^2+K.gtan*K.vr^2)/hi
    D = K.soft^2+lo*(r-K.rp)^2/2
    D > 0 || throw(DomainError(r,"unsoftened source on the real azimuth axis"))
    asinh(sqrt(D/hi)/ (2r))
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
        blocksize::Int=4, max_memory_bytes::Int=512*1024^2) where {T}
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
    q = [2T(pi)*j/Nq for j in 0:Nq-1]
    states = worldline_rates.(Ref(p.orbit),q)
    events = _eccentric_source_event.(Ref(p),states)
    origin = first(states)
    tw = [s.dt_dlambda/(Nq*p.orbit.Γ) for s in states]
    dt = [s.Δtr-origin.Δtr for s in states]
    dp = [s.Δφr-origin.Δφr for s in states]
    blocks = [i:min(i+blocksize-1,n) for i in 1:blocksize:n]
    strips = [minimum(strip_halfwidth(K,r) for K in events for r in rr[b]) for b in blocks]
    nphi = [max(something(Nphi,2),_eccentric_nphi(delta,mmax,T)) for delta in strips]
    peak = maximum(nphi)
    # Includes states temporarily retained during construction (11 reals/event).
    bound = retained+big(sizeof(T))*(4big(peak)*min(blocksize,n)+peak+11big(Nq)) + 4096+64big(length(blocks))
    bound <= max_memory_bytes || throw(ArgumentError("strip-adapted FFT workspace exceeds memory cap"))
    kernels = Array{Complex{T}}(undef,mmax+1,Nq,n,2)
    # One buffer survives all blocks, so different block lengths do not leave
    # multiple large dead buffers awaiting collection.
    scratch = Array{Complex{T}}(undef,peak,2,min(blocksize,n))
    for (ib,b) in enumerate(blocks)
        N = nphi[ib]
        work = view(scratch,1:N,:,1:length(b))
        for iq in eachindex(q)
            K = events[iq]
            for (j,ir) in enumerate(b), l in 0:N-1
                I, Ir = _eccentric_kernel(K,rr[ir],2T(pi)*l/N)
                work[l+1,1,j], work[l+1,2,j] = I, Ir
            end
            fft_radix2!(work;dims=1)
            for (j,ir) in enumerate(b), derivative in 1:2, m in 0:mmax
                kernels[m+1,iq,ir,derivative] = work[m+1,derivative,j]/N
            end
        end
    end
    SourceAzimuthBank(p,rr,q,tw,dt,dp,events,kernels,nphi,strips,blocksize,Int(bound))
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
