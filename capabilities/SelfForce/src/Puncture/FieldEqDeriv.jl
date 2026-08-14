# Lorenz-gauge wave-operator second radial derivative ddhret[i,l,m] of the retarded field, from the
# Barack-Lousto field equations. Port of the field-equation Block in Punctures/subtraction_fieldeq.jl
# (h1S-exact-order-eps2.wl L353-368, Wardell et al. arXiv:2508.00087). Purely algebraic: no HDF5
# input or radial solver — each component's second r-derivative is built from the ten mode values
# h_k, their first r-derivatives drh_k, and the time-derivative scalar dt = -i m Ω (Ω = √(M/r0³)).
# Type-generic {T}: a Complex{BigFloat} h/dh pair runs at full precision.

"""
    EVENSTATIC_M0_COMPONENTS

Even-parity components whose m=0 (static) second r-derivative is NOT taken from these field
equations but read directly from the EvenStatic solver output (L378-402 overwrites `ddhret[i,l,0]`
for these `i` at even `l`). Used by the HDF5 stage to route m=0 even-static modes around the
field-equation operator.
"""
const EVENSTATIC_M0_COMPONENTS = (1, 2, 3, 4, 5, 6, 7)

"""
    is_evenstatic_ddh(i, l, m)

True when `ddhret[i,l,m]` is an even-static mode supplied by the solver (m=0, even l, even-parity
component) rather than computed from the field equations.
"""
is_evenstatic_ddh(i::Integer, l::Integer, m::Integer) =
    m == 0 && iseven(l) && i in EVENSTATIC_M0_COMPONENTS

"""
    ddhret(i, l, m, r, r0, M, h, dh)

Second radial derivative of the retarded field component `i` (1..10) at angular numbers `(l, m)`,
evaluated at field-point radius `r` for an orbit of radius `r0` about mass `M`. `h` and `dh` are
length-10 indexable collections of the mode values `h_1..h_10` and first r-derivatives `drh_1..drh_10`
at this field point. Returns the complex second derivative. The time derivative enters as
`dt = -i m √(M/r0³)`; couplings to the higher-l-only amplitudes (`h5,h7,h10`) are gated by
`l>=1`/`l>=2` exactly as in the source.
"""
function ddhret(i::Integer, l::Integer, m::Integer, r::Real, r0::Real, M::Real, h, dh)
    T = float(promote_type(typeof(r), typeof(r0), typeof(M)))
    CT = Complex{T}
    rT, r0T, MT = T(r), T(r0), T(M)
    f = 1 - 2MT / rT
    fp = 2MT / rT^2
    dt = CT(0, -m) * sqrt(MT / r0T^3)

    h1, h2, h3, h4, h5, h6, h7, h8, h9, h10 =
        CT(h[1]), CT(h[2]), CT(h[3]), CT(h[4]), CT(h[5]),
        CT(h[6]), CT(h[7]), CT(h[8]), CT(h[9]), CT(h[10])
    d1, d2, d3, d4, d5, d6, d7, d8, d9, d10 =
        CT(dh[1]), CT(dh[2]), CT(dh[3]), CT(dh[4]), CT(dh[5]),
        CT(dh[6]), CT(dh[7]), CT(dh[8]), CT(dh[9]), CT(dh[10])

    # l>=1 / l>=2 gated amplitudes (If[l>=k, h_j, 0] in the source).
    g1(x) = l >= 1 ? x : zero(CT)
    g2(x) = l >= 2 ? x : zero(CT)

    # Common second-order operator acting on (h_i, drh_i).
    common(hi, di) = -(fp / f) * di + dt^2 * hi / f^2 +
                     ((2MT) / rT^3 + l * (l + 1) / rT^2) * hi / f

    coup = if i == 1
        (f * fp * d1) / 2 - (fp * (dt * h2)) / 2 +
        f^2 / (2rT^2) * (h1 - f * h3 - g1(h5) - f * h6)
    elseif i == 2
        (f * fp * d2) / 2 - (fp * (dt * h1)) / 2 +
        f^2 / (2rT^2) * (h2 - g1(h4))
    elseif i == 3
        (-f / (2rT^2)) * (h1 - g1(h5) - (1 - 4MT / rT) * (h3 + h6))
    elseif i == 4
        (f * fp * d4) / 4 - (fp * (dt * h5)) / 4 - (3 // 4) * fp * f / rT * h4 -
        (l * (l + 1) * f) / (2rT^2) * h2
    elseif i == 5
        (f * fp * d5) / 4 - (fp * (dt * h4)) / 4 +
        f / rT^2 * (1 - (7 // 2) * MT / rT) * h5 -
        f / (2rT^2) * l * (l + 1) * (h1 - f * h3 - f * h6) -
        f^2 / (2rT^2) * g2(h7)
    elseif i == 6
        (-f / (2rT^2)) * (h1 - g1(h5) - (1 - 4MT / rT) * (h3 + h6))
    elseif i == 7
        (-f / (2rT^2)) * (h7 + (l - 1) * (l + 2) * h5)
    elseif i == 8
        (f * fp * d8) / 4 - (fp * (dt * h9)) / 4 - (3 // 4) * fp * f / rT * h8
    elseif i == 9
        (fp * (f * d9 - dt * h8)) / 4 +
        f / rT^2 * ((1 - (7 // 2) * MT / rT) * h9 - f / 2 * g2(h10))
    elseif i == 10
        (-f / (2rT^2)) * (h10 + (l - 1) * (l + 2) * h9)
    else
        throw(ArgumentError("ddhret: component i must be 1..10, got $i"))
    end

    hi = (h1, h2, h3, h4, h5, h6, h7, h8, h9, h10)[i]
    di = (d1, d2, d3, d4, d5, d6, d7, d8, d9, d10)[i]
    return common(hi, di) + (4 / f^2) * coup
end
