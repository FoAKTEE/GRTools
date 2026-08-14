# Vendored adaptive Runge-Kutta integrators — resolves the 4× ODE vendoring (Teukolsky.radial_numint,
# H1Lorenz.integrator, DiscModes.integrator, FEW.dopr853) into one type-generic family.
#
#   - DP5  (Dormand-Prince RK4(5), PI step control, target-landing) — type-generic {T}, BigFloat-clean.
#   - GBS  (Gragg-Bulirsch-Stoer modified-midpoint + Richardson extrapolation) — the high-order
#          BigFloat-feasible path (set an explicit `kmax`; `gbs_kmax` returns 0 for Float64).
#   - DOPR853 (Dormand-Prince 8(5,3) with degree-7 dense output) — the FEW trajectory integrator,
#          Float64, ported verbatim into the nested `Dopr853` submodule.
#
# Step tolerance inherits precision via `rk_tol(T)` (precision.jl): Float64 → 1e-12, BigFloat →
# 2^(20−bits). Sources: H1Lorenz/src/integrator.jl (DP5/GBS), FEW/src/trajectory/dopr853.jl (DOPR853).

# DP5 Butcher tableau (rationals → exact at any precision)
const _DP_C = (1//5, 3//10, 4//5, 8//9, 1//1, 1//1)
const _DP_A = (
    (1//5,),
    (3//40, 9//40),
    (44//45, -56//15, 32//9),
    (19372//6561, -25360//2187, 64448//6561, -212//729),
    (9017//3168, -355//33, 46732//5247, 49//176, -5103//18656),
    (35//384, 0//1, 500//1113, 125//192, -2187//6784, 11//84),
)
const _DP_B  = (35//384, 0//1, 500//1113, 125//192, -2187//6784, 11//84, 0//1)
const _DP_BH = (5179//57600, 0//1, 7571//16695, 393//640, -92097//339200, 187//2100, 1//40)

"""
    dp5_step!(f!, xc, y, k1, h, tol) -> (err, y5, k7)

One adaptive Dormand-Prince RK4(5) step from `(xc, y)` with first-stage derivative `k1` and step `h`.
Returns the scaled RMS error estimate `err` (relative to `tol`), the 5th-order solution `y5`, and the
final-stage derivative `k7` (FSAL reuse). Type-generic; `f!(g, x, y)` writes `dy/dx` into `g`.
"""
function dp5_step!(f!, xc::T, y::Vector{T}, k1::Vector{T}, h::T, tol::T) where {T}
    ny = length(y); k = Vector{Vector{T}}(undef, 7); k[1] = k1
    g = similar(y)
    for s in 1:6
        ys = copy(y)
        @inbounds for j in 1:s
            aj = T(_DP_A[s][j]); aj == 0 && continue
            ys .+= (h * aj) .* k[j]
        end
        f!(g, xc + h * T(_DP_C[s]), ys); k[s+1] = copy(g)
    end
    y5 = copy(y); y4 = copy(y)
    @inbounds for s in 1:7
        y5 .+= (h * T(_DP_B[s])) .* k[s]; y4 .+= (h * T(_DP_BH[s])) .* k[s]
    end
    err = zero(T)
    @inbounds for j in 1:ny
        sc = tol + tol * max(abs(y[j]), abs(y5[j]))
        err += ((y5[j] - y4[j]) / sc)^2
    end
    return sqrt(err / ny), y5, k[7]
end

"""
    integrate_to_targets!(f!, y0, x0, targets, dir; tol) -> Matrix (length(targets) × ny)

Adaptive DP5 integration of `dy/dx = f!(g,x,y)` from `x0` (state `y0`) in direction `dir` (+1/−1),
capping each step to land exactly on the sorted `targets`. Row `i` of the result is the state at
`targets[i]`. Type-generic, BigFloat-clean. (H1Lorenz `_integrate_to_targets!`.)
"""
function integrate_to_targets!(f!, y0::Vector{T}, x0::T, targets::AbstractVector, dir::Int; tol::T) where {T}
    n = length(targets); ny = length(y0)
    Y = zeros(T, n, ny)
    y = copy(y0); xc = x0
    g = similar(y); f!(g, xc, y); k1 = copy(g)
    safety = T(9 // 10); minfac = T(1 // 5); maxfac = T(5); order = 5
    idx = 1
    xtarget = T(targets[1])
    h = dir > 0 ? abs(xtarget - xc) : -abs(xtarget - xc)
    h == 0 && (h = T(dir) * eps(T))
    while idx <= n
        xtarget = T(targets[idx])
        if (dir > 0 && xc >= xtarget - eps(T) * abs(xtarget) * 8) ||
           (dir < 0 && xc <= xtarget + eps(T) * abs(xtarget) * 8)
            Y[idx, :] .= y; idx += 1; continue
        end
        if (dir > 0 && xc + h > xtarget) || (dir < 0 && xc + h < xtarget)
            h = xtarget - xc
        end
        err, y5, k7 = dp5_step!(f!, xc, y, k1, h, tol)
        if err <= 1 || abs(h) <= eps(T) * abs(xc) * 16
            xc += h; y = y5; k1 = k7
            fac = err == 0 ? maxfac : min(maxfac, max(minfac, safety * err^(-one(T) / order)))
            h *= fac
        else
            h *= max(minfac, safety * err^(-one(T) / order))
        end
    end
    return Y
end

"""
    gbs_kmax(T) -> Int

Gragg-Bulirsch-Stoer max extrapolation depth: `0` for non-BigFloat `T` (use DP5), else
`clamp(cld(bits,24)+3, 6, 18)` so the tableau can reach full BigFloat precision in a few columns.
"""
function gbs_kmax(::Type{T}) where {T}
    T <: BigFloat || return 0
    bits = precision(BigFloat)
    return clamp(cld(bits, 24) + 3, 6, 18)
end

"""
    gragg_midpoint!(f!, x, y, h, n) -> Vector

Gragg modified-midpoint rule: from `(x, y)` take `n` substeps of size `h/n` and return
`y(x + h)` with the Gragg smoothing endpoint formula `(z_{n-1} + z_n + (h/n)·f(x+h, z_n))/2`.
Error expansion is even in `h` (h²-extrapolable), which is what makes it the building block
of `gbs_integrate_to_targets!`. Type-generic `{T}`; mutates only its internal work arrays.
"""
function gragg_midpoint!(f!, x::T, y::Vector{T}, h::T, n::Int) where {T}
    ny = length(y); hh = h / n
    z0 = copy(y); z1 = similar(y); g = similar(y)
    f!(g, x, z0)
    @inbounds for q in 1:ny; z1[q] = z0[q] + hh * g[q]; end
    for k in 1:(n-1)
        f!(g, x + k * hh, z1)
        @inbounds for q in 1:ny
            t = z0[q] + 2 * hh * g[q]; z0[q] = z1[q]; z1[q] = t
        end
    end
    f!(g, x + h, z1)
    yn = similar(y)
    @inbounds for q in 1:ny; yn[q] = (z0[q] + z1[q] + hh * g[q]) / 2; end
    return yn
end

"""
    gbs_integrate_to_targets!(f!, y0, x0, targets, dir; tol, kmax) -> Matrix

Gragg-Bulirsch-Stoer extrapolation integrator (explicit, RHS-only): per step a modified-midpoint
sweep with substep counts `2,4,…,2·kmax`, Richardson-extrapolated in `h²`; the diagonal differences
give the error estimate and working order. Adaptive step + exact target landing. Type-generic; the
BigFloat-feasible high-order alternative to DP5. (H1Lorenz `_gbs_integrate_to_targets!`.)
"""
function gbs_integrate_to_targets!(f!, y0::Vector{T}, x0::T, targets::AbstractVector, dir::Int;
                                  tol::T, kmax::Int) where {T}
    n = length(targets); ny = length(y0)
    Y = zeros(T, n, ny)
    y = copy(y0); xc = x0
    nseq = [2k for k in 1:kmax]
    ktarget = max(3, kmax - 1)
    span1 = abs(T(targets[1]) - xc)
    hnext = T(dir) * max(span1 / 64, eps(T) * 16 * (1 + abs(xc)))
    h = hnext
    idx = 1
    while idx <= n
        xtarget = T(targets[idx])
        if (dir > 0 && xc >= xtarget - eps(T) * abs(xtarget) * 16) ||
           (dir < 0 && xc <= xtarget + eps(T) * abs(xtarget) * 16)
            Y[idx, :] .= y; idx += 1; continue
        end
        h = hnext
        if (dir > 0 && xc + h > xtarget) || (dir < 0 && xc + h < xtarget)
            h = xtarget - xc
        end
        Ttab = Vector{Vector{T}}(undef, kmax)
        prevdiag = Vector{T}(undef, ny)
        err = T(Inf); kused = kmax; accepted = false
        for k in 1:kmax
            Ttab[k] = gragg_midpoint!(f!, xc, y, h, nseq[k])
            prev = Ttab[k]
            for j in (k-1):-1:1
                r = (T(nseq[k]) / T(nseq[j]))^2
                newv = similar(prev)
                @inbounds for q in 1:ny
                    newv[q] = prev[q] + (prev[q] - Ttab[j][q]) / (r - 1)
                end
                Ttab[j] = prev; prev = newv
            end
            best = prev
            if k >= 2
                e = zero(T)
                @inbounds for q in 1:ny
                    sc = tol + tol * max(abs(y[q]), abs(best[q]))
                    e += ((best[q] - prevdiag[q]) / sc)^2
                end
                err = sqrt(e / ny)
                @inbounds prevdiag .= best
                if err <= 1
                    y = best; kused = k; accepted = true; break
                end
            else
                @inbounds prevdiag .= best
            end
        end
        if accepted || abs(h) <= eps(T) * abs(xc) * 64
            accepted || (y = copy(prevdiag))
            xc += h
            fac = max(T(1) / 4, min(T(4), T(2)^(ktarget - kused)))
            cand = abs(h) * fac
            spanmax = abs(T(targets[end]) - x0) + span1
            hnext = T(dir) * min(cand, spanmax)
        else
            h /= 8; hnext = h
        end
    end
    return Y
end

# =====================================================================================================
# DOPR853 — Dormand-Prince 8(5,3) with degree-7 dense output (FEW trajectory integrator, Float64).
# Ported verbatim from packages/FastEMRIWaveforms/src/trajectory/dopr853.jl (the tableau / dense-output
# coefficients carry full Float64 precision from upstream src/few/trajectory/dopr853.py).
# =====================================================================================================
module Dopr853

export DOPR853, integrate, dense_eval, take_step

const c2 = 0.526001519587677318785587544488e-01
const c3 = 0.789002279381515978178381316732e-01
const c4 = 0.118350341907227396726757197510e00
const c5 = 0.281649658092772603273242802490e00
const c6 = 0.333333333333333333333333333333e00
const c7 = 0.25e00
const c8 = 0.307692307692307692307692307692e00
const c9 = 0.651282051282051282051282051282e00
const c10 = 0.6e00
const c11 = 0.857142857142857142857142857142e00
const c14 = 0.1e00
const c15 = 0.2e00
const c16 = 0.777777777777777777777777777778e00

const b1 = 5.42937341165687622380535766363e-2
const b6 = 4.45031289275240888144113950566e0
const b7 = 1.89151789931450038304281599044e0
const b8 = -5.8012039600105847814672114227e0
const b9 = 3.1116436695781989440891606237e-1
const b10 = -1.52160949662516078556178806805e-1
const b11 = 2.01365400804030348374776537501e-1
const b12 = 4.47106157277725905176885569043e-2

const bhh1 = 0.244094488188976377952755905512e00
const bhh2 = 0.733846688281611857341361741547e00
const bhh3 = 0.220588235294117647058823529412e-01

const er1 = 0.1312004499419488073250102996e-01
const er6 = -0.1225156446376204440720569753e01
const er7 = -0.4957589496572501915214079952e00
const er8 = 0.1664377182454986536961530415e01
const er9 = -0.3503288487499736816886487290e00
const er10 = 0.3341791187130174790297318841e00
const er11 = 0.8192320648511571246570742613e-01
const er12 = -0.2235530786388629525884427845e-01

const a21 = 5.26001519587677318785587544488e-2
const a31 = 1.97250569845378994544595329183e-2
const a32 = 5.91751709536136983633785987549e-2
const a41 = 2.95875854768068491816892993775e-2
const a43 = 8.87627564304205475450678981324e-2
const a51 = 2.41365134159266685502369798665e-1
const a53 = -8.84549479328286085344864962717e-1
const a54 = 9.24834003261792003115737966543e-1
const a61 = 3.7037037037037037037037037037e-2
const a64 = 1.70828608729473871279604482173e-1
const a65 = 1.25467687566822425016691814123e-1
const a71 = 3.7109375e-2
const a74 = 1.70252211019544039314978060272e-1
const a75 = 6.02165389804559606850219397283e-2
const a76 = -1.7578125e-2
const a81 = 3.70920001185047927108779319836e-2
const a84 = 1.70383925712239993810214054705e-1
const a85 = 1.07262030446373284651809199168e-1
const a86 = -1.53194377486244017527936158236e-2
const a87 = 8.27378916381402288758473766002e-3
const a91 = 6.24110958716075717114429577812e-1
const a94 = -3.36089262944694129406857109825e0
const a95 = -8.68219346841726006818189891453e-1
const a96 = 2.75920996994467083049415600797e1
const a97 = 2.01540675504778934086186788979e1
const a98 = -4.34898841810699588477366255144e1
const a101 = 4.77662536438264365890433908527e-1
const a104 = -2.48811461997166764192642586468e0
const a105 = -5.90290826836842996371446475743e-1
const a106 = 2.12300514481811942347288949897e1
const a107 = 1.52792336328824235832596922938e1
const a108 = -3.32882109689848629194453265587e1
const a109 = -2.03312017085086261358222928593e-2
const a111 = -9.3714243008598732571704021658e-1
const a114 = 5.18637242884406370830023853209e0
const a115 = 1.09143734899672957818500254654e0
const a116 = -8.14978701074692612513997267357e0
const a117 = -1.85200656599969598641566180701e1
const a118 = 2.27394870993505042818970056734e1
const a119 = 2.49360555267965238987089396762e0
const a1110 = -3.0467644718982195003823669022e0
const a121 = 2.27331014751653820792359768449e0
const a124 = -1.05344954667372501984066689879e1
const a125 = -2.00087205822486249909675718444e0
const a126 = -1.79589318631187989172765950534e1
const a127 = 2.79488845294199600508499808837e1
const a128 = -2.85899827713502369474065508674e0
const a129 = -8.87285693353062954433549289258e0
const a1210 = 1.23605671757943030647266201528e1
const a1211 = 6.43392746015763530355970484046e-1
const a141 = 5.61675022830479523392909219681e-2
const a147 = 2.53500210216624811088794765333e-1
const a148 = -2.46239037470802489917441475441e-1
const a149 = -1.24191423263816360469010140626e-1
const a1410 = 1.5329179827876569731206322685e-1
const a1411 = 8.20105229563468988491666602057e-3
const a1412 = 7.56789766054569976138603589584e-3
const a1413 = -8.298e-3
const a151 = 3.18346481635021405060768473261e-2
const a156 = 2.83009096723667755288322961402e-2
const a157 = 5.35419883074385676223797384372e-2
const a158 = -5.49237485713909884646569340306e-2
const a1511 = -1.08347328697249322858509316994e-4
const a1512 = 3.82571090835658412954920192323e-4
const a1513 = -3.40465008687404560802977114492e-4
const a1514 = 1.41312443674632500278074618366e-1
const a161 = -4.28896301583791923408573538692e-1
const a166 = -4.69762141536116384314449447206e0
const a167 = 7.68342119606259904184240953878e0
const a168 = 4.06898981839711007970213554331e0
const a169 = 3.56727187455281109270669543021e-1
const a1613 = -1.39902416515901462129418009734e-3
const a1614 = 2.9475147891527723389556272149e0
const a1615 = -9.15095847217987001081870187138e0

const d41 = -0.84289382761090128651353491142e01
const d46 = 0.56671495351937776962531783590e00
const d47 = -0.30689499459498916912797304727e01
const d48 = 0.23846676565120698287728149680e01
const d49 = 0.21170345824450282767155149946e01
const d410 = -0.87139158377797299206789907490e00
const d411 = 0.22404374302607882758541771650e01
const d412 = 0.63157877876946881815570249290e00
const d413 = -0.88990336451333310820698117400e-01
const d414 = 0.18148505520854727256656404962e02
const d415 = -0.91946323924783554000451984436e01
const d416 = -0.44360363875948939664310572000e01
const d51 = 0.10427508642579134603413151009e02
const d56 = 0.24228349177525818288430175319e03
const d57 = 0.16520045171727028198505394887e03
const d58 = -0.37454675472269020279518312152e03
const d59 = -0.22113666853125306036270938578e02
const d510 = 0.77334326684722638389603898808e01
const d511 = -0.30674084731089398182061213626e02
const d512 = -0.93321305264302278729567221706e01
const d513 = 0.15697238121770843886131091075e02
const d514 = -0.31139403219565177677282850411e02
const d515 = -0.93529243588444783865713862664e01
const d516 = 0.35816841486394083752465898540e02
const d61 = 0.19985053242002433820987653617e02
const d66 = -0.38703730874935176555105901742e03
const d67 = -0.18917813819516756882830838328e03
const d68 = 0.52780815920542364900561016686e03
const d69 = -0.11573902539959630126141871134e02
const d610 = 0.68812326946963000169666922661e01
const d611 = -0.10006050966910838403183860980e01
const d612 = 0.77771377980534432092869265740e00
const d613 = -0.27782057523535084065932004339e01
const d614 = -0.60196695231264120758267380846e02
const d615 = 0.84320405506677161018159903784e02
const d616 = 0.11992291136182789328035130030e02
const d71 = -0.25693933462703749003312586129e02
const d76 = -0.15418974869023643374053993627e03
const d77 = -0.23152937917604549567536039109e03
const d78 = 0.35763911791061412378285349910e03
const d79 = 0.93405324183624310003907691704e02
const d710 = -0.37458323136451633156875139351e02
const d711 = 0.10409964950896230045147246184e03
const d712 = 0.29840293426660503123344363579e02
const d713 = -0.43533456590011143754432175058e02
const d714 = 0.96324553959188282948394950600e02
const d715 = -0.39177261675615439165231486172e02
const d716 = -0.14972683625798562581422125276e03

const beta = 0.0
const alpha = 1.0 / 8.0 - beta * 0.2
const safe = 0.9
const minscale = 1.0 / 3.0
const maxscale = 6.0

"""
    DOPR853(ode!; abstol=1e-10)

Adaptive Dormand-Prince 8(5,3) integrator state: the user ODE `ode!(dydx, x, y, args)`, the
per-component absolute tolerance, and the PI-controller carry-over.
"""
mutable struct DOPR853{F}
    ode!::F
    abstol::Float64
    errOld::Float64
    previousReject::Bool
end
DOPR853(ode!; abstol::Real=1e-10) = DOPR853(ode!, Float64(abstol), 1e-4, false)

function _dormand_prince_steps!(ode!, x, solOld, h, args,
                                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, arg)
    ode!(k1, x, solOld, args)
    @. arg = solOld + h * (a21 * k1);                                   ode!(k2, x + c2 * h, arg, args)
    @. arg = solOld + h * (a31 * k1 + a32 * k2);                        ode!(k3, x + c3 * h, arg, args)
    @. arg = solOld + h * (a41 * k1 + a43 * k3);                        ode!(k4, x + c4 * h, arg, args)
    @. arg = solOld + h * (a51 * k1 + a53 * k3 + a54 * k4);             ode!(k5, x + c5 * h, arg, args)
    @. arg = solOld + h * (a61 * k1 + a64 * k4 + a65 * k5);             ode!(k6, x + c6 * h, arg, args)
    @. arg = solOld + h * (a71 * k1 + a74 * k4 + a75 * k5 + a76 * k6);  ode!(k7, x + c7 * h, arg, args)
    @. arg = solOld + h * (a81 * k1 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7)
    ode!(k8, x + c8 * h, arg, args)
    @. arg = solOld + h * (a91 * k1 + a94 * k4 + a95 * k5 + a96 * k6 + a97 * k7 + a98 * k8)
    ode!(k9, x + c9 * h, arg, args)
    @. arg = solOld + h * (a101 * k1 + a104 * k4 + a105 * k5 + a106 * k6 + a107 * k7 +
                           a108 * k8 + a109 * k9)
    ode!(k10, x + c10 * h, arg, args)
    @. arg = solOld + h * (a111 * k1 + a114 * k4 + a115 * k5 + a116 * k6 + a117 * k7 +
                           a118 * k8 + a119 * k9 + a1110 * k10)
    ode!(k2, x + c11 * h, arg, args)
    @. arg = solOld + h * (a121 * k1 + a124 * k4 + a125 * k5 + a126 * k6 + a127 * k7 +
                           a128 * k8 + a129 * k9 + a1210 * k10 + a1211 * k2)
    ode!(k3, x + h, arg, args)
    return nothing
end

function _error!(solNew, solOld, h, k1, k2, k3, k6, k7, k8, k9, k10, abstol)
    n = length(solOld)
    err = 0.0; err2 = 0.0; sk = 1.0 / abstol
    @inbounds for i in 1:n
        temp = b1 * k1[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] +
               b9 * k9[i] + b10 * k10[i] + b11 * k2[i] + b12 * k3[i]
        solNew[i] = solOld[i] + h * temp
        e2 = (temp - bhh1 * k1[i] - bhh2 * k9[i] - bhh3 * k3[i]) * sk
        err2 += e2 * e2
        e1 = (er1 * k1[i] + er6 * k6[i] + er7 * k7[i] + er8 * k8[i] + er9 * k9[i] +
              er10 * k10[i] + er11 * k2[i] + er12 * k3[i]) * sk
        err += e1 * e1
    end
    denom = err + 0.01 * err2
    denom = denom <= 0.0 ? 1.0 : denom
    return abs(h) * err * sqrt(1.0 / (n * denom))
end

function _controller(err, errOld, previousReject, h)
    if err <= 1.0
        scale = err == 0.0 ? maxscale :
                clamp(safe * err^(-alpha) * errOld^beta, minscale, maxscale)
        hnew = previousReject ? h : h * scale
        return (true, hnew, clamp(err, 1e-4, 1e300), false)
    else
        scale = clamp(safe * err^(-alpha), minscale, 1e300)
        return (false, h * scale, errOld, true)
    end
end

function _prep_evaluate(ode!, x0, y0, h0, x1, y1, args,
                        k1, k2, k3, k4, k5, k6, k7, k8, k9, k10)
    n = length(y0)
    newDer = zeros(n)
    ode!(newDer, x1, y1, args)
    diff = y1 .- y0
    rcont1 = copy(y0)
    rcont2 = diff
    bspl = h0 .* k1 .- diff
    rcont3 = bspl
    rcont4 = diff .- h0 .* newDer .- bspl
    rcont5 = @. d41 * k1 + d46 * k6 + d47 * k7 + d48 * k8 + d49 * k9 + d410 * k10 + d411 * k2 + d412 * k3
    rcont6 = @. d51 * k1 + d56 * k6 + d57 * k7 + d58 * k8 + d59 * k9 + d510 * k10 + d511 * k2 + d512 * k3
    rcont7 = @. d61 * k1 + d66 * k6 + d67 * k7 + d68 * k8 + d69 * k9 + d610 * k10 + d611 * k2 + d612 * k3
    rcont8 = @. d71 * k1 + d76 * k6 + d77 * k7 + d78 * k8 + d79 * k9 + d710 * k10 + d711 * k2 + d712 * k3
    arg = @. y0 + h0 * (a141 * k1 + a147 * k7 + a148 * k8 + a149 * k9 + a1410 * k10 +
                        a1411 * k2 + a1412 * k3 + a1413 * newDer)
    ode!(k10, x0 + c14 * h0, arg, args)
    arg = @. y0 + h0 * (a151 * k1 + a156 * k6 + a157 * k7 + a158 * k8 + a1511 * k2 +
                        a1512 * k3 + a1513 * newDer + a1514 * k10)
    ode!(k2, x0 + c15 * h0, arg, args)
    arg = @. y0 + h0 * (a161 * k1 + a166 * k6 + a167 * k7 + a168 * k8 + a169 * k9 +
                        a1613 * newDer + a1614 * k10 + a1615 * k2)
    ode!(k3, x0 + c16 * h0, arg, args)
    @. rcont5 = h0 * (rcont5 + d413 * newDer + d414 * k10 + d415 * k2 + d416 * k3)
    @. rcont6 = h0 * (rcont6 + d513 * newDer + d514 * k10 + d515 * k2 + d516 * k3)
    @. rcont7 = h0 * (rcont7 + d613 * newDer + d614 * k10 + d615 * k2 + d616 * k3)
    @. rcont8 = h0 * (rcont8 + d713 * newDer + d714 * k10 + d715 * k2 + d716 * k3)
    cont = Matrix{Float64}(undef, 8, n)
    @inbounds for i in 1:n
        cont[1, i] = rcont1[i]; cont[2, i] = rcont2[i]
        cont[3, i] = rcont3[i]; cont[4, i] = rcont4[i]
        cont[5, i] = rcont5[i]; cont[6, i] = rcont6[i]
        cont[7, i] = rcont7[i]; cont[8, i] = rcont8[i]
    end
    return cont
end

"""
    take_step(solver, x, y, h; args=nothing) -> (success, x1, y1, hnew, cont)

One adaptive DP8(5,3) step from `(x, y)` with trial step `h`, updating the PI-controller carry-over
inside `solver`. On acceptance returns the new `(x1, y1)`, the next trial step, and the segment's
degree-7 dense-output coefficients (`8 × n`, for [`dense_eval`](@ref)). On rejection returns the
original `(x, y)` with a reduced step.
"""
function take_step(solver::DOPR853, x::Real, y::AbstractVector, h::Real; args=nothing)
    n = length(y)
    yv = collect(float.(y))
    k = ntuple(_ -> zeros(n), 10)
    arg = zeros(n); solNew = zeros(n)
    _dormand_prince_steps!(solver.ode!, x, yv, h, args, k..., arg)
    err = _error!(solNew, yv, h, k[1], k[2], k[3], k[6], k[7], k[8], k[9], k[10], solver.abstol)
    success, hnew, errOld2, prevRej2 = _controller(err, solver.errOld, solver.previousReject, h)
    solver.errOld = errOld2; solver.previousReject = prevRej2
    if success
        x1 = x + h; y1 = copy(solNew)
        cont = _prep_evaluate(solver.ode!, x, yv, h, x1, y1, args, k...)
        return (true, x1, y1, hnew, cont)
    else
        return (false, float(x), yv, hnew, Matrix{Float64}(undef, 8, n))
    end
end

"""
    integrate(solver, y0, tmax; h0=0.01, max_step=100000, args=nothing) -> (ts, ys, coeffs)

Integrate the ODE held in `solver` from `x=0` (state `y0`) until `x` reaches `tmax`. Returns the
accepted step locations, the states there, and the per-segment dense-output coefficients.
"""
function integrate(solver::DOPR853, y0, tmax::Real; h0::Real=0.01,
                   max_step::Integer=100_000, args=nothing)
    tgoal = Float64(tmax)
    (isfinite(tgoal) && tgoal >= 0) ||
        throw(ArgumentError("tmax must be finite and nonnegative, got $tmax"))
    h = Float64(h0)
    (isfinite(h) && h > 0) ||
        throw(ArgumentError("h0 must be finite and positive, got $h0"))
    max_step > 0 || throw(ArgumentError("max_step must be positive, got $max_step"))
    y = collect(float.(y0))
    ts = Float64[0.0]
    ys = Vector{Vector{Float64}}([copy(y)])
    coeffs = Vector{Matrix{Float64}}()
    x = 0.0
    attempts = 0
    solver.errOld = 1e-4; solver.previousReject = false
    while x < tgoal
        attempts += 1
        attempts <= max_step ||
            error("integrate: max_step=$max_step exhausted before tmax=$tgoal")
        remaining = tgoal - x
        htrial = min(h, remaining)
        x + htrial > x || error("integrate: step size made no forward progress at x=$x")
        success, x1, y1, hnew, cont = take_step(solver, x, y, htrial; args=args)
        (isfinite(hnew) && hnew > 0) ||
            error("integrate: adaptive controller produced invalid step $hnew")
        if success
            all(isfinite, y1) || error("integrate: accepted a non-finite state at x=$x1")
            htrial == remaining && (x1 = tgoal)
            push!(ts, x1); push!(ys, copy(y1)); push!(coeffs, cont)
            x = x1; y .= y1
        end
        h = hnew
    end
    return ts, ys, coeffs
end

"""
    dense_eval(ts, coeffs, tnew) -> Vector

Evaluate the degree-7 dense interpolant from [`integrate`](@ref) at `tnew ∈ [ts[1], ts[end]]`.
"""
function dense_eval(ts::AbstractVector, coeffs::AbstractVector, tnew::Real)
    tmin, tmax = ts[1], ts[end]
    (tmin <= tnew <= tmax) || throw(ArgumentError("tnew=$tnew outside [$tmin, $tmax]"))
    seg = clamp(searchsortedlast(ts, tnew), 1, length(ts) - 1)
    cont = coeffs[seg]
    h = ts[seg+1] - ts[seg]
    s = (tnew - ts[seg]) / h
    s1 = 1.0 - s
    n = size(cont, 2)
    out = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        out[i] = cont[1, i] + s * (cont[2, i] + s1 * (cont[3, i] + s * (cont[4, i] +
                 s1 * (cont[5, i] + s * (cont[6, i] + s1 * (cont[7, i] + s * cont[8, i]))))))
    end
    return out
end

end # module Dopr853

using .Dopr853: DOPR853, take_step, dense_eval
import .Dopr853: integrate
