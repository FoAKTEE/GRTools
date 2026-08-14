# Characterization of the self-contained confluent hypergeometric ₁F₁(a;b;z) (Kummer's function M).
#
# GENUINELY-NEW primitive (no light legacy oracle): grounded entirely by closed-form analytic
# identities, evaluated to full precision. The incremental-Pochhammer term recurrence mirrors the
# structural pattern of Teukolsky/src/hypergeom.jl `_f21_maclaurin` (READ, heavy → not live-dep'd);
# that file only holds ₂F₁ / U, so ₁F₁ is built fresh here.
#
# Identities used (DLMF 13.6):
#   ₁F₁(a;a;z) = eᶻ                      (and a=b=1 ⇒ eᶻ)
#   ₁F₁(0;b;z) = 1
#   ₁F₁(1;2;z) = (eᶻ−1)/z
#   ₁F₁(−n;b;z) = truncating polynomial  (Σ_{k=0}^{n} (−n)_k/((b)_k k!) zᵏ)
#   Kummer transform ₁F₁(a;b;z) = eᶻ ₁F₁(b−a;b;−z)
#   series partial-sum self-consistency (hand-summed Pochhammer ratio)

using NumericsCore
using NumericsCore: kummer_1F1, kummer_U, _kummer_series, _hyper2F0_asymp, _kummer_asymp
using NumericsCore: fd1c, fd2c

# explicit truncating polynomial ₁F₁(−n;b;z) from the Pochhammer definition (reference)
function _kummer_poly(n::Int, b, z)
    s = one(z); term = one(z)
    for k in 1:n
        term *= (-n + (k - 1)) / ((b + (k - 1)) * k) * z
        s += term
    end
    s
end

@testset "kummer" begin

    @testset "₁F₁(a;a;z) = eᶻ" begin
        for a in (0.5, 1.0, 2.3, -1.7), z in (0.2, 1.0, 3.5, -0.8, -4.0)
            @test kummer_1F1(a, a, z) ≈ exp(z) rtol=1e-13
        end
        @test kummer_1F1(1.0, 1.0, 2.0) ≈ exp(2.0) rtol=1e-14
        # complex argument
        z = 1.2 + 0.7im
        @test kummer_1F1(2.0 + 0im, 2.0 + 0im, z) ≈ exp(z) rtol=1e-12
    end

    @testset "₁F₁(0;b;z) = 1 and ₁F₁(1;2;z) = (eᶻ−1)/z" begin
        for b in (1.0, 2.5, 4.0), z in (0.3, 2.0, -1.5)
            @test kummer_1F1(0.0, b, z) ≈ 1.0 rtol=1e-14
        end
        for z in (0.25, 1.0, 3.0, -0.5, -2.5)
            @test kummer_1F1(1.0, 2.0, z) ≈ (exp(z) - 1) / z rtol=1e-12
        end
    end

    @testset "₁F₁(−n;b;z) truncating polynomial" begin
        for n in 1:5, b in (1.0, 2.0, 3.5), z in (0.4, 1.3, -0.9, 2.2)
            @test kummer_1F1(-float(n), b, z) ≈ _kummer_poly(n, b, z) rtol=1e-12
        end
        # explicit low-order literals
        for z in (0.3, -1.1, 2.0)
            @test kummer_1F1(-1.0, 3.0, z) ≈ 1 - z / 3 rtol=1e-13          # 1 - z/b
            @test kummer_1F1(-2.0, 1.0, z) ≈ 1 - 2z + z^2 / 2 rtol=1e-13   # L_2(z)
        end
    end

    @testset "Kummer transform ₁F₁(a;b;z) = eᶻ ₁F₁(b−a;b;−z)" begin
        for (a, b) in ((0.7, 2.0), (1.5, 3.0), (-0.4, 1.2), (2.0, 5.0))
            for z in (0.6, 1.8, -1.0, -3.0)
                @test kummer_1F1(a, b, z) ≈ exp(z) * kummer_1F1(b - a, b, -z) rtol=1e-12
            end
        end
    end

    @testset "_kummer_series partial-sum self-consistency" begin
        # raw series == hand-summed Pochhammer ratio (no stabilizing transform)
        a, b, z = 1.3, 2.7, 0.4
        s = 1.0; term = 1.0
        for k in 1:40
            term *= (a + (k - 1)) / ((b + (k - 1)) * k) * z
            s += term
        end
        @test _kummer_series(a, b, z) ≈ s rtol=1e-13
        @test _kummer_series(a, a, 0.5) ≈ exp(0.5) rtol=1e-13
    end

    @testset "BigFloat precision inheritance (2^-200)" begin
        setprecision(BigFloat, 256) do
            a = BigFloat(3) / 2; z = BigFloat(7) / 5
            @test kummer_1F1(a, a, z) isa BigFloat
            @test abs(kummer_1F1(a, a, z) - exp(z)) < big(2.0)^-200
            # negative-z stability (transform path)
            zn = -BigFloat(11) / 5
            @test abs(kummer_1F1(one(BigFloat), BigFloat(2), zn) - (exp(zn) - 1) / zn) < big(2.0)^-200
            # complex BigFloat
            zc = BigFloat(6) / 5 + (BigFloat(2) / 5) * im
            @test abs(kummer_1F1(BigFloat(2) + 0im, BigFloat(2) + 0im, zc) - exp(zc)) < big(2.0)^-190
        end
    end

    # ================= large-|z| asymptotic layer (₂F₀ / Tricomi-U / M-asymptotic) =================

    # independent finite-sum reference for the TERMINATING ₂F₀ (α or β a non-positive integer)
    function _f20_ref(α, β, w, N)
        s = one(w); t = one(w)
        for n in 0:N-1
            t *= (α + n) * (β + n) / (n + 1) * w
            s += t
        end
        s
    end

    @testset "_hyper2F0_asymp terminating series is exact" begin
        # β = 0 ⇒ single term = 1
        @test _hyper2F0_asymp(2.3, 0.0, 0.7) ≈ 1.0 rtol=1e-14
        # α = −n ⇒ degree-n polynomial in w; compare to hand sum
        for (α, β, w) in ((-1.0, 1.4, -0.03), (-3.0, 0.9, 0.02), (-2.0, -5.0, -0.1))
            @test _hyper2F0_asymp(α, β, w) ≈ _f20_ref(α, β, w, 12) rtol=1e-13
        end
    end

    @testset "kummer_U closed forms U(a,a+1,z)=z^{-a}, U(0,b,z)=1, U(-n,b,z) polynomial" begin
        for a in (0.5, 1.0, 2.3, -1.4), z in (12.0, 40.0, 80.0)
            @test kummer_U(a, a + 1, z) ≈ z^(-a) rtol=1e-12       # DLMF 13.6.4
        end
        for b in (1.0, 2.5, 3.7), z in (15.0, 50.0)
            @test kummer_U(0.0, b, z) ≈ 1.0 rtol=1e-13            # U(0,b,z)=1
        end
        # truncating U(-n,b,z) = z^{-a} · finite ₂F₀ (exact polynomial)
        for n in 1:4, b in (1.0, 2.0), z in (10.0, 30.0)
            ref = z^n * _f20_ref(-float(n), -n - b + 1, -1 / z, n + 1)
            @test kummer_U(-float(n), b, z) ≈ ref rtol=1e-12
        end
    end

    @testset "kummer_U satisfies Kummer's ODE  z U'' + (b−z) U' − a U = 0" begin
        # the DEFINING equation (DLMF 13.2.1) — independent analytic anchor via 4th-order central FD
        h = 1e-2
        for (a, b, z) in ((0.7, 1.3, 50.0), (1.5, 2.2, 60.0), (-0.4, 0.8, 45.0))
            U(x) = kummer_U(a, b, x)
            resid = z * fd2c(U, z, h) + (b - z) * fd1c(U, z, h) - a * U(z)
            @test abs(resid) < 1e-6 * abs(a * U(z))
        end
    end

    @testset "kummer_U derivative relation U'(a,b,z) = −a U(a+1,b+1,z)" begin
        # DLMF 13.3.22, independent (RHS uses a different parameter set); LHS via FD
        h = 1e-2
        for (a, b, z) in ((0.7, 1.3, 50.0), (1.9, 2.6, 70.0))
            U(x) = kummer_U(a, b, x)
            @test fd1c(U, z, h) ≈ -a * kummer_U(a + 1, b + 1, z) rtol=1e-7
        end
    end

    @testset "kummer_U transform round-trip U(a,b,z)=z^{1-b}U(a-b+1,2-b,z)" begin
        for (a, b, z) in ((0.6, 1.4, 40.0), (1.3, 2.7, 55.0), (0.9, 0.5, 60.0))
            @test kummer_U(a, b, z) ≈ z^(1 - b) * kummer_U(a - b + 1, 2 - b, z) rtol=1e-11
        end
    end

    @testset "kummer_1F1 large-|z| asymptotic branch (Γ-ratio-1 anchors)" begin
        # b=a ⇒ eᶻ exactly (asymptotic ₂F₀ truncates to one term)
        for a in (0.5, 1.7, 3.0), z in (60.0, 120.0)
            @test kummer_1F1(a, a, z) ≈ exp(z) rtol=1e-13
            @test kummer_1F1(a, a, z) == _kummer_asymp(a, a, z)   # routed through the asymptotic
        end
        # a=1,b=2 ⇒ eᶻ/z, matches (eᶻ−1)/z to the recessive tail e^{−z}
        for z in (60.0, 100.0)
            @test kummer_1F1(1.0, 2.0, z) ≈ (exp(z) - 1) / z rtol=1e-13
            @test kummer_1F1(1.0, 2.0, z) == _kummer_asymp(1.0, 2.0, z)
        end
    end

    @testset "asymptotic ≈ series in the overlap annulus (anchor family)" begin
        # the INDEPENDENT _kummer_series (direct Maclaurin) and _kummer_asymp agree at moderate |z|
        for z in (25.0, 30.0, 35.0)
            @test _kummer_asymp(1.0, 2.0, z) ≈ _kummer_series(1.0, 2.0, z) rtol=1e-9
            @test _kummer_asymp(2.0, 2.0, z) ≈ _kummer_series(2.0, 2.0, z) rtol=1e-11
        end
    end

    @testset "asymptotic layer BigFloat stability (setprecision 256)" begin
        setprecision(BigFloat, 256) do
            # kummer_U closed form z^{-a} at BigFloat precision
            a = BigFloat(3) / 2; z = BigFloat(60)
            @test abs(kummer_U(a, a + 1, z) - z^(-a)) < big(2.0)^-160
            # M asymptotic anchor a=1,b=2 needs Re z past the BigFloat threshold (~178) to be exact
            zb = BigFloat(250)
            @test abs(kummer_1F1(one(BigFloat), BigFloat(2), zb) - (exp(zb) - 1) / zb) < big(2.0)^-160 * abs((exp(zb) - 1) / zb)
            # b=a is exact at any large z
            @test abs(kummer_1F1(a, a, zb) - exp(zb)) < big(2.0)^-160 * exp(zb)
        end
    end
end
