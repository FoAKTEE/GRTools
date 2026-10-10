# Independent Darwin-angle integration of the equatorial Kerr equations (M = 1).
# Uses constants_of_motion, but NO trajectory, frequency, elliptic or worldline_rates code.
module EquatorialQuadratureOracle
using Spacetime: constants_of_motion
using NumericsCore: gl_nodes_weights
export quadrature_orbit, integrate_oracle

function integrate_oracle(f, lo, hi, nodes, weights)
    mid = (lo + hi) / 2
    half = (hi - lo) / 2
    return half * sum(w * f(mid + half * z) for (z, w) in zip(nodes, weights))
end

function quadrature_orbit(a::T, p::T, e::T; n=128) where {T<:AbstractFloat}
    E, L, _ = constants_of_motion(abs(a), p, e, a < 0 ? -one(T) : one(T))
    a < 0 && (L = -L)
    # R = (1-E²)(r_apo-r)(r-r_peri)(r-r3)r; r4=0 at Q=0.
    r3 = 2 / (1-E^2) - p / (1-e) - p / (1+e)
    nodes, weights = gl_nodes_weights(T, n)
    πT = T(pi)
    radius(χ) = p / (1 + e*cos(χ))
    function dlambda(χ)
        r = radius(χ)
        return sqrt(1-e^2) / ((1+e*cos(χ))*sqrt((1-E^2)*(r-r3)*r))
    end
    function rates(χ)
        r = radius(χ)
        P = E*(r^2+a^2)-a*L
        Δ = r^2-2r+a^2
        return (t=(r^2+a^2)*P/Δ+a*(L-a*E), phi=a*P/Δ+L-a*E, tau=r^2)
    end
    quad(f, hi) = integrate_oracle(f, zero(T), hi, nodes, weights)
    # Half-period reflection is exact because these integrands depend only on cos(χ).
    halfλ = quad(dlambda, πT)
    Υr = πT/halfλ
    means = map(key -> quad(χ -> getproperty(rates(χ), key)*dlambda(χ), πT)/halfλ,
                (:t, :phi, :tau))
    function at_phase(q)
        q = mod(T(q), 2πT)
        reflected = q > πT
        qhalf = reflected ? 2πT-q : q
        target = qhalf/Υr
        χ = qhalf
        lo, hi = zero(T), πT
        for _ in 1:100
            residual = quad(dlambda, χ)-target
            abs(residual) <= 32eps(T)*halfλ && break
            residual > 0 ? (hi=χ) : (lo=χ)
            next = χ-residual/dlambda(χ)
            χ = lo < next < hi ? next : (lo+hi)/2
        end
        abs(quad(dlambda, χ)-target) <= 64eps(T)*halfλ || error("Darwin inversion did not converge")
        r = radius(χ)
        Ur = p*e*sin(χ)/(1+e*cos(χ))^2/dlambda(χ)
        deviations = map((key, mean) -> quad(χ -> (getproperty(rates(χ),key)-mean)*dlambda(χ), χ),
                         (:t, :phi, :tau), means)
        s = reflected ? -one(T) : one(T)
        return (r=r, Ur=s*Ur, Δtr=s*deviations[1], Δφr=s*deviations[2], Δτr=s*deviations[3])
    end
    return (; E, L, at_phase, Υr, Γ=means[1], Υφ=means[2], Υτ=means[3], rates, dlambda, radius, nodes, weights)
end
end
