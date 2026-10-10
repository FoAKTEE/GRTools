# ZZ22 (arXiv:2112.03311) Huygens-wavelet tracks as overplotted in FR22 (arXiv:2207.14637v2) Figs 2-4, row 3.
#
# Kinematics only (no disc response, no amplitude). Units GM = a = n = 1, T_p = 2pi; FR22 disc p = q = 0,
# h_p = 0.1  =>  Omega(R) = R^(-3/2) (alpha = 3/2) and constant c_s = 0.1 (beta = 0).
# Planet: FR22 eqs (24)-(25), mean anomaly M = t, pericentre at t = 0, phi_p(0) = 0.
# Wavelet emitted at t1, seen at t2 (ZZ22 msMNRAS.tex l.147-176):
#   (4)/(5)/(6)  R = R_pl(t1) -/+ c_s (t2 - t1)            inner (-) / outer (+); inner R <= 0 discarded
#   pattern      Omega_pat = Omega_pl(t1) = dphi_p/dt(t1)   (instantaneous; the only reading consistent with ZZ22 eq 12)
#   (8)          theta_c = phi_p(t1) - sgn(R_pl(t1) - R_c) * R_c^(-1/2)/c_s * B(R_pl(t1)/R_c),  R_c = Omega_pat^(-2/3)
#   (9)          theta   = Omega_pat (t2 - t1) + theta_c + sgn(R - R_c) * R_c^(-1/2)/c_s * B(R/R_c)
#   B(x) = 3 - x - 2 x^(-1/2): the bracket of ZZ22 eq (2) at (alpha, beta) = (3/2, 0).
# Plotted y = mod(theta - phi_p(t2), 2pi)/pi, x = log10 R.
# Emission schedule: one wavelet per T_p/800 from t = 0 to the snapshot t2 = (4 + T/T_p) T_p inclusive. This window is
# the one the published FR22 vector paths encode (oldest outer wavelet at age 4 + T/T_p orbits in all 12 panels);
# see emri-disk results/emri-eccentric/paper_EMRI-ECC/evidence/fr22-zz22-overlay/.
#
# Base Julia only. Usage: julia --startup-file=no _zz22_wavelets.jl <outdir>
# writes <outdir>/zz22_e<e>_T<T>.csv (family,k,log10r,y; k = (t2 - t1)/(T_p/800)) for e in {0.1,0.25,0.5},
# T/T_p in {0,0.25,0.5,0.75}.

module ZZ22Wavelets

using Printf

const CS = 0.1           # FR22 h_p = 0.1 at a, globally isothermal (3_numerical_procedure.tex:79)
const PER_ORBIT = 800    # FR22 4_results.tex:36, one wavelet every T_p/800

"Eccentric anomaly for mean anomaly M in [0, 2pi) (Newton)."
function kepler_E(M, e)
    E = e < 0.8 ? M : oftype(M, pi)
    for _ in 1:100
        dE = (E - e * sin(E) - M) / (1 - e * cos(E))
        E -= dE
        abs(dE) < 1e-15 && break
    end
    return E
end

"Planet (r, continuous phi, dphi/dt) at time t; pericentre at t = 0 with phi = 0."
function planet(t, e)
    nturn = floor(t / 2pi)
    E = kepler_E(t - 2pi * nturn, e)
    f = 2 * atan(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))   # in [0, 2pi] for E in [0, 2pi)
    return (r = 1 - e * cos(E), phi = f + 2pi * nturn, Omega = (1 + e * cos(f))^2 / (1 - e^2)^1.5)
end

zz22_B(x) = 3 - x - 2 / sqrt(x)

"(R, theta) at t2 of the wavelet launched at t1; dir = -1 inner, +1 outer. `nothing` if R <= 0."
function wavelet(t1, t2, e, dir; cs = CS)
    p = planet(t1, e)
    R = p.r + dir * cs * (t2 - t1)
    R > 0 || return nothing
    Rc = p.Omega^(-2 / 3)
    pref = Rc^(-1 / 2) / cs
    theta_c = p.phi - sign(p.r - Rc) * pref * zz22_B(p.r / Rc)
    return (R = R, theta = p.Omega * (t2 - t1) + theta_c + sign(R - Rc) * pref * zz22_B(R / Rc))
end

"All wavelets of one family at snapshot phase T (units of T_p), newest first: vector of (k, log10 R, y)."
function track(e, T, dir; per_orbit = PER_ORBIT, window_orbits = 4 + T)
    dt = 2pi / per_orbit
    t2 = (4 + T) * 2pi
    phi2 = planet(t2, e).phi
    out = Tuple{Int,Float64,Float64}[]
    for k in 0:round(Int, window_orbits * per_orbit)
        w = wavelet(t2 - k * dt, t2, e, dir)
        w === nothing && continue
        push!(out, (k, log10(w.R), mod(w.theta - phi2, 2pi) / pi))
    end
    return out
end

function main(outdir)
    mkpath(outdir)
    for e in (0.1, 0.25, 0.5), T in (0.0, 0.25, 0.5, 0.75)
        path = joinpath(outdir, @sprintf("zz22_e%g_T%.2f.csv", e, T))
        open(path, "w") do io
            println(io, "family,k,log10r,y")
            for (fam, dir) in (("inner", -1), ("outer", +1)), (k, x, y) in track(e, T, dir)
                @printf(io, "%s,%d,%.6f,%.6f\n", fam, k, x, y)
            end
        end
        println(path)
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ZZ22Wavelets.main(only(ARGS))
end
