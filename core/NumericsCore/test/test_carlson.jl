# Characterization of the consolidated Carlson/Legendre/Jacobi primitives.
#
# Oracle 1 (complex Carlson + complete Legendre, principal branch): KerrGeodesics, BAKED into
#   fixtures/numericscore_carlson_ref.jl (Phase-7-B, captured once while KG was still resolvable):
#   ._RF/._RC/._RD/._RJ/._ellipK/._ellipE/._ellipPi/._ellipFinc/._ellipEinc/._ellipPiinc
#   /._jacobiSN/._jacobiCN/._jacobiDN  — NumericsCore re-implements the same Carlson duplication-
#   theorem algorithm independently, so we expect agreement to full Float64 precision (rtol 1e-12).
# Oracle 2 (FEW real / modulus convention, NO live dep): the heavy FastEMRIWaveforms test-dep is
#   dropped, but FEW's wrappers are LITERALLY defined as real(KerrGeodesics._R*/_ellip*)
#   (packages/FastEMRIWaveforms/src/elliptic.jl:8-73: `RF(x,y,z)=real(KerrGeodesics._RF(complex...))`,
#   `EllipK(k)=real(KerrGeodesics._ellipK(k^2))`), so the retained KG oracle reproduces FEW exactly.
# Oracle 3 (closed form / Abramowitz-Stegun identities): RF(1,1,1)=1, RC(x,x)=1/sqrt(x),
#   symmetry, ellipK(0)=ellipE(0)=pi/2, ellipE(1)=1, sn^2+cn^2=1, dn^2+m sn^2=1.

using NumericsCore
# Baked from KerrGeodesics._* (captured once, Phase-7-B) — no live legacy dep.
include("fixtures/numericscore_carlson_ref.jl")

@testset "carlson" begin
    rt = 1e-12

    @testset "Carlson vs KerrGeodesics (complex inputs) (baked from KG._R*, captured once)" begin
        cx, cy, cz, cp = 1.0+2.0im, 2.0-1.0im, 3.0+0.5im, 0.7+0.3im
        @test NumericsCore.RC(cx, cy)        ≈ REF_RC_cxy   rtol=rt
        @test NumericsCore.RF(cx, cy, cz)    ≈ REF_RF_cxyz  rtol=rt
        @test NumericsCore.RD(cx, cy, cz)    ≈ REF_RD_cxyz  rtol=rt
        @test NumericsCore.RJ(cx, cy, cz, cp)≈ REF_RJ_cxyzp rtol=rt
    end

    @testset "Carlson real inputs (FEW convention == real(KG)) (baked from KG._R*, captured once)" begin
        # FEW's real wrappers ARE real(KerrGeodesics._R*(complex(args))) (elliptic.jl:25-57), so the
        # KG oracle (complexified real args, real part) reproduces the dropped FEW oracle; baked here.
        @test NumericsCore.RC(2.0, 3.0)           ≈ REF_RC_23   rtol=rt
        @test NumericsCore.RF(1.0, 2.0, 3.0)      ≈ REF_RF_123  rtol=rt
        @test NumericsCore.RD(1.0, 2.0, 3.0)      ≈ REF_RD_123  rtol=rt
        @test NumericsCore.RJ(1.0, 2.0, 3.0, 4.0) ≈ REF_RJ_1234 rtol=rt
        # real inputs return a real value (type-following collapse of the two forks)
        @test NumericsCore.RF(1.0, 2.0, 3.0) isa Real
    end

    @testset "complete Legendre vs KerrGeodesics (baked from KG._ellip*, captured once)" begin
        for m in (0.001, 0.3, 0.7, 0.999, 2.4972)   # incl. m>1 principal branch
            @test NumericsCore.ellipK(m) ≈ REF_ellipK[m] rtol=rt
            @test NumericsCore.ellipE(m) ≈ REF_ellipE[m] rtol=rt
        end
        for (n, m) in ((0.2, 0.3), (0.6, 0.7))
            @test NumericsCore.ellipPi(n, m) ≈ REF_ellipPi[(n, m)] rtol=rt
        end
        # complex characteristic n branch
        nC = 0.5 + 0.2im
        @test NumericsCore.ellipPi(nC, 0.4) ≈ REF_ellipPi_nC rtol=rt
    end

    @testset "complete Legendre modulus convention (FEW EllipK(k)=K(k^2)) (baked from KG._ellip*, captured once)" begin
        # FEW's EllipK/EllipE take the MODULUS k (parameter m=k^2) and return real(KG._ellipK/_ellipE(k^2))
        # (elliptic.jl:66,73). Ground the parameter<->modulus convention against the baked KG oracle.
        for k in (0.3, 0.6, 0.9)
            @test real(NumericsCore.ellipK(k^2)) ≈ REF_ellipK_mod[k] rtol=rt
            @test real(NumericsCore.ellipE(k^2)) ≈ REF_ellipE_mod[k] rtol=rt
        end
    end

    @testset "incomplete Legendre vs KerrGeodesics (baked from KG._ellip*inc, captured once)" begin
        for phi in (0.7, 4.0, -1.3)          # incl. phi outside [-pi/2, pi/2] (period reduction)
            @test NumericsCore.ellipFinc(phi, 0.3) ≈ REF_ellipFinc[phi] rtol=rt
            @test NumericsCore.ellipEinc(phi, 0.3) ≈ REF_ellipEinc[phi] rtol=rt
            @test NumericsCore.ellipPiinc(0.2, phi, 0.3) ≈ REF_ellipPiinc[phi] rtol=rt
        end
        @test NumericsCore._phi_reduce(4.0) == REF_phi_reduce_4
    end

    @testset "Jacobi vs KerrGeodesics (baked from KG._jacobi*, captured once)" begin
        for (u, m) in ((0.8, 0.5), (1.7, 0.2), (2.3, 0.9))
            @test NumericsCore.jacobiAM(u, m) ≈ REF_jacobiAM[(u, m)] rtol=rt
            @test NumericsCore.jacobiSN(u, m) ≈ REF_jacobiSN[(u, m)] rtol=rt
            @test NumericsCore.jacobiCN(u, m) ≈ REF_jacobiCN[(u, m)] rtol=rt
            @test NumericsCore.jacobiDN(u, m) ≈ REF_jacobiDN[(u, m)] rtol=rt
        end
    end

    @testset "closed-form / A&S identities" begin
        @test NumericsCore.RF(1.0, 1.0, 1.0) ≈ 1.0 rtol=rt
        @test NumericsCore.RC(4.0, 4.0) ≈ 1/sqrt(4.0) rtol=rt
        # RF total symmetry, RD symmetry in first two args
        @test NumericsCore.RF(1.0, 2.0, 3.0) ≈ NumericsCore.RF(3.0, 1.0, 2.0) rtol=rt
        @test NumericsCore.RD(1.0, 2.0, 3.0) ≈ NumericsCore.RD(2.0, 1.0, 3.0) rtol=rt
        @test real(NumericsCore.ellipK(0.0)) ≈ pi/2 rtol=rt
        @test NumericsCore.ellipK(1.0) == complex(Inf, 0.0)
        @test NumericsCore.ellipK(1) == complex(Inf, 0.0)
        @test real(NumericsCore.ellipE(0.0)) ≈ pi/2 rtol=rt
        @test real(NumericsCore.ellipE(1.0)) ≈ 1.0  rtol=rt
        s = NumericsCore.jacobiSN(0.8, 0.5); c = NumericsCore.jacobiCN(0.8, 0.5)
        d = NumericsCore.jacobiDN(0.8, 0.5)
        @test s^2 + c^2 ≈ 1.0 rtol=rt
        @test d^2 + 0.5*s^2 ≈ 1.0 rtol=rt
    end

    @testset "BigFloat genericity (precision inheritance)" begin
        setprecision(BigFloat, 256) do
            m = big"0.3"; u = big"0.8"
            # ellipK(0)=pi/2 at BigFloat precision (not Float64-capped)
            @test real(NumericsCore.ellipK(big(0.0))) isa BigFloat
            s = NumericsCore.jacobiSN(u, m); c = NumericsCore.jacobiCN(u, m)
            d = NumericsCore.jacobiDN(u, m)
            @test s isa BigFloat
            @test abs(s^2 + c^2 - 1) < big(2.0)^-200
            @test abs(d^2 + m*s^2 - 1) < big(2.0)^-200
            # RF(1,1,1)=1 to BigFloat precision
            @test abs(NumericsCore.RF(big(1.0), big(1.0), big(1.0)) - 1) < big(2.0)^-200
            @test real(NumericsCore.ellipK(big(1.0))) == BigFloat(Inf)

            # A nonsymmetric Carlson value must not stop at the old fixed 1e-6 residual floor.
            got = NumericsCore.RF(big"0.2", big"1.3", big"2.7")
            ref = setprecision(BigFloat, 512) do
                NumericsCore.RF(big"0.2", big"1.3", big"2.7")
            end
            @test abs(got - ref) < big"1e-65"
        end
    end
end
