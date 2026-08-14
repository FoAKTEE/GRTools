# Analytic first-order metric perturbations (mass/spin monopole+dipole) and the Lorenz-gauge
# second radial derivative. Pure closed-form, type-generic {T}, no green-dep. Transcribed from
# ref-code/SecondOrderRicci/src/{h1_M.cc, h1_S.cc, h1.cc}.

const _MassField{T} = @NamedTuple{h::Vector{Complex{T}}, dh::Vector{Complex{T}}, ddh::Vector{Complex{T}}}

"""
    h1_mass([T=Float64,] r; M=1) -> Dict{NTuple{3,Int}, @NamedTuple{h, dh, ddh}}

First-order mass metric perturbation on radial grid `r`. Transcription of `h1_M`
(`ref-code/SecondOrderRicci/src/h1_M.cc` L46-100): the only nonzero components are
`(i,l,m) = (1,0,0)`, `(3,0,0)`, `(6,0,0)`, given by log-bearing closed forms in `r` and
`f = 1 - 2M/r`, with `dh`, `ddh` the closed-form first/second radial derivatives. Returns a sparse
map keyed by `(i,l,m)`; each value is a NamedTuple of `Complex{T}` vectors `(h, dh, ddh)`.
"""
function h1_mass(::Type{T}, r::AbstractVector; M::Real=1) where {T<:AbstractFloat}
    rr = T.(r); n = length(rr); Mt = T(M)
    keys3 = (1, 3, 6)
    h = Dict(i => Vector{Complex{T}}(undef, n) for i in keys3)
    dh = Dict(i => Vector{Complex{T}}(undef, n) for i in keys3)
    ddh = Dict(i => Vector{Complex{T}}(undef, n) for i in keys3)
    for j in 1:n
        rj = rr[j]; fj = 1 - 2Mt / rj
        c = _h1_mass_point(Mt, rj, fj)
        for i in keys3
            h[i][j] = c[i][1]; dh[i][j] = c[i][2]; ddh[i][j] = c[i][3]
        end
    end
    return Dict{NTuple{3,Int},_MassField{T}}(
        (i, 0, 0) => (h = h[i], dh = dh[i], ddh = ddh[i]) for i in keys3
    )
end
h1_mass(r::AbstractVector; M::Real=1) = h1_mass(Float64, r; M=M)

# Scalar kernel: (h, dh, ddh) for i ∈ {1,3,6} at one (r, f). Faithful transcription of h1_M.cc
# (pow→^, M_PI→pi, M kept symbolic).
function _h1_mass_point(M::T, r::T, f::T) where {T}
    sq2pi = sqrt(2 * T(pi))
    L = log((r / M) / 2); Lr = log(r / M); Lf = log(f)
    L4 = log(T(4)); L256 = log(T(256))

    h1 = (2 * M^-1 * sq2pi * r^-4 * (16 * L * (M - r) * M^3 -
        2 * (4 * M^2 - r^2) * (-6 * M * r + 2 * f * (2 * M - r) * r + 5 * M^2 + 2 * r^2) +
        Lf * (-16 * r * M^3 + 24 * M^4 + M * f * r * (2 * M * r + 4 * M^2 + r^2) +
            3 * M * r^3 - r^4))) / 3
    dh1 = (-2 * sq2pi * f^-1 * r^-6 * (32 * r * M^3 - 48 * M^4 - 6 * M * r^3 +
        f * r * (8 * f * r * (-2 * M + r) * (6 * M + r) + 16 * Lr * (4 * M - 3 * r) * M^2 +
            4 * (31 + 6 * L4) * r * M^2 - 8 * (15 + L256) * M^3 - 30 * M * r^2 - 4 * r^3 +
            Lf * (-52 * r * M^2 + 88 * M^3 - 2 * M * r^2 +
                f * r * (4 * M * r + 12 * M^2 + r^2) + 3 * r^3)) + 2 * r^4)) / 3
    ddh1 = (4 * sq2pi * f^-2 * r^-8 * (2 * M * (16 * r * M^3 - 24 * M^4 - 3 * M * r^3 + r^4) +
        f * r * (f * r * (32 * L * (5 * M - 3 * r) * M^2 + 236 * r * M^2 -
            248 * M^3 + Lf * (2 * M - r) * (-2 * M * r + 104 * M^2 - 3 * r^2) - 54 * M * r^2 +
            f * r * (8 * (6 * M * r - 24 * M^2 + r^2) +
                Lf * (6 * M * r + 24 * M^2 + r^2)) - 4 * r^3) +
            2 * (66 * r * M^3 - 116 * M^4 + M^2 * r^2 - 6 * M * r^3 + r^4)))) / 3

    h3 = (2 * M^-1 * sq2pi * f^-1 * r^-4 *
        (16 * L * (-2 * M + r) * M^3 +
            f * r * (64 * M^3 + M * Lf * (2 * M * r + 4 * M^2 + r^2) - 2 * r^3) +
            Lf * (16 * r * M^3 - 24 * M^4 - 3 * M * r^3 + r^4) +
            2 * (3 * M * r + 10 * M^2 + r^2) * (-2 * M + r)^2)) / 3
    dh3 = (-2 * sq2pi * f^-2 * r^-6 * (f^2 * r^2 *
        (192 * M^2 + Lf * (4 * M * r + 12 * M^2 + r^2)) +
        2 * (2 * M - r) * (-10 * r * M^2 + (52 + 8 * L4) * M^3 - 16 * Lr * M^3 -
            3 * M * r^2 - r^3) + f * r *
        (-4 * (47 + 6 * L4) * r * M^2 + 16 * Lr * (-8 * M + 3 * r) * M^2 +
            8 * (43 + 8 * L4) * M^3 + 6 * M * r^2 - 2 * r^3 -
            3 * Lf * (-16 * r * M^2 + 32 * M^3 + r^3)) +
        2 * Lf * (16 * r * M^3 - 24 * M^4 - 3 * M * r^3 + r^4))) / 3
    ddh3 = (4 * sq2pi * f^-3 * r^-8 * (2 * M * (2 * M - r) *
        (-22 * r * M^2 + 4 * (29 + L256) * M^3 - 32 * Lr * M^3 - 7 * M * r^2 -
            r^3) + f^3 * (384 * M^2 + Lf * (6 * M * r + 24 * M^2 + r^2)) *
        r^3 + f^2 * r^2 * (-4 * (101 + 12 * L4) * r * M^2 +
            32 * Lr * (-10 * M + 3 * r) * M^2 + 16 * (57 + 10 * L4) * M^3 + 8 * M * r^2 -
            2 * r^3 - 3 * Lf * (-32 * r * M^2 + 80 * M^3 + r^3)) +
        4 * M * Lf * (16 * r * M^3 - 24 * M^4 - 3 * M * r^3 + r^4) +
        2 * f * r * (-2 * (153 + 16 * L4) * r * M^3 + 32 * Lr * (-5 * M + 2 * r) * M^3 +
            (548 + 80 * L4) * M^4 + 11 * M^2 * r^2 + 2 * M * r^3 + r^4 +
            Lf * (64 * r * M^3 - 120 * M^4 - 6 * M * r^3 + r^4)))) / 3

    h6 = (-4 * M^-1 * sq2pi * r^-3 * (-4 * r * M^2 + 4 * (5 + L4) * M^3 -
        8 * Lr * M^3 - M * r^2 + 4 * f * r * (2 * M * r + 4 * M^2 + r^2) -
        4 * r^3 + Lf * (-8 * M^3 + r^3))) / 3
    dh6 = (4 * sq2pi * f^-1 * r^-5 * (16 * M^3 +
        f * r * (8 * f * r * (4 * M + r) + 3 *
            (-8 * M * r + 4 * (3 + L4) * M^2 - 8 * (Lf + Lr) * M^2 - 3 * r^2)) -
        2 * r^3)) / 3
    ddh6 = (-8 * sq2pi * f^-2 * r^-7 * (16 * M^4 +
        f * r * (64 * M^3 + f * r * (-44 * M * r + 8 * f * r * (6 * M + r) + 4 * (13 + 6 * L4) * M^2 -
            48 * (Lf + Lr) * M^2 - 9 * r^2) - 2 * r^3) - 2 * M * r^3)) / 3

    return Dict(
        1 => (complex(h1), complex(dh1), complex(ddh1)),
        3 => (complex(h3), complex(dh3), complex(ddh3)),
        6 => (complex(h6), complex(dh6), complex(ddh6)),
    )
end

"""
    h1_spin([T=Float64,] r) -> Dict{NTuple{3,Int}, @NamedTuple{h, dh, ddh}}

First-order angular-momentum (spin) metric perturbation on radial grid `r` (`M = 1`). Transcription
of `h1_S` (`ref-code/SecondOrderRicci/src/h1_S.cc` L46-55): the only nonzero components are
`(8,1,0)` and `(9,1,0)`, `h[8] = -8√(π/3) r⁻²`, `h[9] = -16√(π/3) r⁻³`, with `dh = ∂_r h`,
`ddh = ∂_r dh` the exact radial derivatives.
"""
function h1_spin(::Type{T}, r::AbstractVector) where {T<:AbstractFloat}
    c = sqrt(T(pi) / 3)
    rr = T.(r)
    h8   = @. -8c * rr^-2
    dh8  = @. 16c * rr^-3
    ddh8 = @. -48c * rr^-4
    h9   = @. -16c * rr^-3
    dh9  = @. 48c * rr^-4
    ddh9 = @. -192c * rr^-5
    return Dict{NTuple{3,Int},_MassField{T}}(
        (8, 1, 0) => (h = complex.(h8), dh = complex.(dh8), ddh = complex.(ddh8)),
        (9, 1, 0) => (h = complex.(h9), dh = complex.(dh9), ddh = complex.(ddh9)),
    )
end
h1_spin(r::AbstractVector) = h1_spin(Float64, r)

"""
    d2hdr2([T=Float64,] i, l, m, r, f, fp, r0, hbar, dhbar) -> Complex{T}

Second radial derivative `∂²_r h̄ᵢ` of the trace-reversed first-order field from the Lorenz-gauge
field equations on a circular orbit. Transcription of `d2hdr2` (`ref-code/SecondOrderRicci/src/
h1.cc` L59-122): given the field `hbar` and its first derivative `dhbar` (1-based length-10 complex
vectors), the field equation for basis index `i ∈ 1:10` fixes `∂²_r h̄ᵢ` in terms of the other
components, `l`, `m`, the local `r,f=1-2M/r,fp=df/dr`, and the orbital frequency `Ω=√M r0^{-3/2}`
(via `dt=-i m Ω`). `M=1`.
"""
function d2hdr2(::Type{T}, i::Integer, l::Integer, m::Integer,
                r, f, fp, r0,
                hbar::AbstractVector, dhbar::AbstractVector) where {T<:AbstractFloat}
    M = one(T)
    rT = T(r); fT = T(f); fpT = T(fp); r0T = T(r0)
    Omega = sqrt(M) * r0T^(-T(3) / 2)
    dt = -im * T(m) * Omega
    ld = T(l); ll1 = ld * (ld + 1)
    m2Ω2 = -Omega^2 * T(m)^2
    hb = Complex{T}.(hbar); db = Complex{T}.(dhbar)
    r2 = rT^2; r3 = r2 * rT; f2 = fT^2

    if i == 1
        return -(fpT / fT) * db[1] + hb[1] * m2Ω2 / f2 +
               hb[1] * ((2M) / rT + ll1) / (fT * r2) +
               4 / f2 * (T(1)/2 * fT * fpT * db[1] - T(1)/2 * fpT * dt * hb[2] +
               f2 / (2 * r2) * (hb[1] - fT * hb[3] - hb[5] - fT * hb[6]))
    elseif i == 2
        return -((db[2] * fpT) / fT) +
               (4 * (-(dt * hb[1] * fpT) / 2 + (db[2] * fT * fpT) / 2 +
                (f2 * (hb[2] - hb[4])) / (2 * r2))) / f2 +
               (hb[2] * ((2M) / r3 + ll1 / r2)) / fT +
               hb[2] * m2Ω2 / f2
    elseif i == 3
        return -((db[3] * fpT) / fT) + (hb[3] * ((2M) / r3 + ll1 / r2)) / fT +
               hb[3] * m2Ω2 / f2 -
               (2 * (hb[1] - hb[5] - (hb[3] + hb[6]) * (1 - (4M) / rT))) / (fT * r2)
    elseif i == 4
        return -((db[4] * fpT) / fT) + (hb[4] * ((2M) / r3 + ll1 / r2)) / fT +
               (4 * (-(dt * hb[5] * fpT) / 4 + (db[4] * fT * fpT) / 4 -
                (fT * hb[2] * ll1) / (2 * r2) - (3 * fT * fpT * hb[4]) / (4 * rT))) / f2 +
               hb[4] * m2Ω2 / f2
    elseif i == 5
        return -((db[5] * fpT) / fT) + (hb[5] * ((2M) / r3 + ll1 / r2)) / fT +
               (4 * (-(dt * hb[4] * fpT) / 4 + (db[5] * fT * fpT) / 4 -
                (f2 * hb[7]) / (2 * r2) -
                (fT * (hb[1] - fT * hb[3] - fT * hb[6]) * ll1) / (2 * r2) +
                (fT * hb[5] * (1 - (7M) / (2 * rT))) / r2)) / f2 +
               hb[5] * m2Ω2 / f2
    elseif i == 6
        return -((db[6] * fpT) / fT) + (hb[6] * ((2M) / r3 + ll1 / r2)) / fT +
               hb[6] * m2Ω2 / f2 -
               (2 * (hb[1] - hb[5] - (hb[3] + hb[6]) * (1 - (4M) / rT))) / (fT * r2)
    elseif i == 7
        return -((db[7] * fpT) / fT) + (hb[7] * ((2M) / r3 + ll1 / r2)) / fT +
               hb[7] * m2Ω2 / f2 -
               (2 * (hb[7] + hb[5] * (-1 + ld) * (2 + ld))) / (fT * r2)
    elseif i == 8
        return -((db[8] * fpT) / fT) + (hb[8] * ((2M) / r3 + ll1 / r2)) / fT +
               (4 * (-(dt * hb[9] * fpT) / 4 + (db[8] * fT * fpT) / 4 -
                (3 * fT * fpT * hb[8]) / (4 * rT))) / f2 +
               hb[8] * m2Ω2 / f2
    elseif i == 9
        return -((db[9] * fpT) / fT) + (hb[9] * ((2M) / r3 + ll1 / r2)) / fT +
               (4 * (((-dt * hb[8] + db[9] * fT) * fpT) / 4 +
                (fT * (-(fT * hb[10]) / 2 + hb[9] * (1 - (7M) / (2 * rT)))) / r2)) / f2 +
               hb[9] * m2Ω2 / f2
    elseif i == 10
        return -((db[10] * fpT) / fT) + (hb[10] * ((2M) / r3 + ll1 / r2)) / fT +
               hb[10] * m2Ω2 / f2 -
               (2 * (hb[10] + hb[9] * (-1 + ld) * (2 + ld))) / (fT * r2)
    else
        throw(ArgumentError("d2hdr2: basis index i must be in 1:10, got $i"))
    end
end

d2hdr2(i::Integer, l::Integer, m::Integer, r, f, fp, r0,
       hbar::AbstractVector, dhbar::AbstractVector) =
    d2hdr2(Float64, i, l, m, r, f, fp, r0, hbar, dhbar)
