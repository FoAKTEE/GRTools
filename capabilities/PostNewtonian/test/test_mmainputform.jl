# Characterization tests for PostNewtonian.MMAInputForm (the symbolic MMA-InputForm→Num walker).
#
# GROUNDING:
#  (1) WALKER-CORRECTNESS GATE — the parse census over the 91-file SeriesData store. The legacy
#      PostNewtonianSelfForce walker deserialized only ~65 of 91 (the eccentric `e`-graded generating
#      functions + complex `ExtraFunctions` (RePsi/ImPsi/PolyGamma) / RedShift / Spin / Memory atoms
#      were its documented [HOLE]). The new unified walker CLOSES that hole: all 91 files (circular
#      AND eccentric) now deserialize — the `.m` TEXT is the byte-exact oracle for the eccentric
#      shapes (NO Mathematica). ZERO deferred.
#  (2) WOLFRAMSCRIPT BYTE-EXACT ORACLE on Kerr/Spherical/Orbit/Energy.m — the InputForm coeff lists
#      (a=3/7,x=5/11), (a=-2/5,x=9/13), a->0, captured verbatim in the legacy series_data_test.jl.
#      Two independent rational substitutions guard against accidental agreement; a->0 leading orders
#      match the published Schwarzschild circular-orbit energy series {1,0,-1/2,0,3/8}.

using Test
using PostNewtonian.MMAInputForm
using PostNewtonian.Catalog
using PostNewtonian.Series: PNSeries
import Symbolics
import SpecialFunctions

# --- (1) parse census --------------------------------------------------------------------
# The unified walker now handles the circular/spherical complex-AST atom class (PolyGamma / RePsi /
# ImPsi / Re / Im / Logp / the imaginary unit `I`), unblocking the 4 Kerr Circular Memory `h-l*m0`
# files (clean, small) over the prior 26-file [HOLE].
const _MEMORY_UNBLOCKED = Set([
    "/Kerr/Circular/Memory/h-l2m0",
    "/Kerr/Circular/Memory/h-l4m0",
    "/Kerr/Circular/Memory/h-l6m0",
    "/Kerr/Circular/Memory/h-l8m0",
])

# The 8 Kerr Spherical flux/scalar-flux giants (~35k-term per-coefficient sums). Module-3's balanced
# reduction (`_accumulate_balanced` over the iteratively-flattened `±` chain) now MATERIALIZES them
# in-census — the former StackOverflow on the ~3.4MB "Horizon" files is resolved. They parse; the
# census below includes them.
const _GIANT_MATERIALIZED = Set([
    "/Kerr/Spherical/Flux/AngularMomentum/Horizon",
    "/Kerr/Spherical/Flux/AngularMomentum/Infinity",
    "/Kerr/Spherical/Flux/Energy/Horizon",
    "/Kerr/Spherical/Flux/Energy/Infinity",
    "/Kerr/Spherical/ScalarFlux/AngularMomentum/Horizon",
    "/Kerr/Spherical/ScalarFlux/AngularMomentum/Infinity",
    "/Kerr/Spherical/ScalarFlux/Energy/Horizon",
    "/Kerr/Spherical/ScalarFlux/Energy/Infinity",
])

# Module-2a unblocked the 4 PURE eccentric files (no held generating-function): the nested
# `SeriesData[e,0,{..}]` coefficient (Spin-p), the `Rational[p,q]` head, the `ResummedSeriesData`
# resummation (EHorizon/LHorizon/Spin-y), and the InputForm implicit-multiplication `(341/16) e^2` /
# `240 (1-e^2)^(3/2)(1+..)` juxtaposition are now all deserialized.
const _MODULE2A_UNBLOCKED = Set([
    "/Schwarzschild/Eccentric/Flux/EHorizon",
    "/Schwarzschild/Eccentric/Flux/LHorizon",
    "/Schwarzschild/Eccentric/Local/Spin-p",
    "/Schwarzschild/Eccentric/Local/Spin-y",
])

# Module-2b closed the eccentric generating-FUNCTION files: the held `HoldForm[oneP5PN[e,..]]` /
# `chi3PN`/`chi4PNConv`/… gen-func head mints a FORMAL Symbolics atom (its `:=` program was stripped
# with the brace-balanced `"ExtraFunctions"` block); the Kerr `TagBox[RowBox[…],HoldForm]` box-
# language, the `mantissa`prec` / `*^` numeric forms, `O[e]^n` truncation, and the `Dot` (`.`) scalar
# product are all handled. ZERO deferred files remain — the whole 91-file store deserializes.
const _MODULE2B_UNBLOCKED = Set([
    "/Kerr/Eccentric/Local/Redshift",
    "/Schwarzschild/Eccentric/Flux/EInfinity",
    "/Schwarzschild/Eccentric/Flux/LInfinity",
    "/Schwarzschild/Eccentric/Local/Redshift-p",
    "/Schwarzschild/Eccentric/Local/Redshift-y",
])
const _DEFERRED = Set(String[])   # ZERO [HOLE]: every .m file (circular AND eccentric) parses

# Module-1 unblocked complex-atom-class files (formerly in `_DEFERRED`, now clean parse).
const _MODULE1_UNBLOCKED = Set([
    "/Kerr/Circular/Flux/Energy/Horizon",
    "/Kerr/Circular/Flux/Energy/Infinity",
    "/Kerr/Circular/Local/RedShift",
    "/Kerr/Circular/Local/Spin",
    "/Kerr/Generic/Flux/EInfinity-v",
])

@testset "MMAInputForm: 91-file parse census vs legacy walker" begin
    cat = pn_series_catalog()
    @test length(cat) == 91
    parsed = String[]
    failed = String[]
    for k in cat
        path = string(SERIES_DATA_DIR, k, ".m")
        try
            d = parse_pn_data_file(path)   # returns a Dict{String,Any} on success
            push!(parsed, k)
        catch
            push!(failed, k)
        end
    end
    # WHOLE-STORE GATE: all 91 files deserialize (65 original + 4 Memory + 5 Module-1 complex-atom +
    # 8 Module-3 giants + 4 Module-2a pure eccentric + 5 Module-2b gen-func = 91). ZERO deferred.
    @test length(parsed) == 91
    @test isempty(failed)
    @test _MEMORY_UNBLOCKED ⊆ Set(parsed)
    @test _MODULE1_UNBLOCKED ⊆ Set(parsed)
    @test _GIANT_MATERIALIZED ⊆ Set(parsed)      # Module-3 balanced reduction: no StackOverflow
    @test _MODULE2A_UNBLOCKED ⊆ Set(parsed)      # Module-2a pure eccentric (nested-e/Rational/Resummed)
    @test _MODULE2B_UNBLOCKED ⊆ Set(parsed)      # Module-2b gen-func (formal atom/TagBox/*^/Dot)
    @test Set(failed) == _DEFERRED
end

# --- (2) Energy.m oracle (verbatim InputForm coeff lists from wolframscript) ---------------
const _ORACLE_SUB1 = "1, 0, -1/2, 0, 3/8, -15/77, 161883/94864, -135/154, 3888891/758912, -11670345/3652264, 134368656705/8999178496, -80073495/7304528, 1558178436663/35996713984, -12605985592965/346468372096, 53678424615147891/426849034422272, -82013278979025/692936744192, 358301726033184477/975654935822336, -6244554566744745165/16433687825257472, 12482278089264748111635/11569316228981260288, -39736449509482043415/32867375650514944, 1032941964509374635421389/323940854411475288064, -11916820951328882379214845/3117930723710449647616, 36404154842326068245476109229/3841290651611273965862912, -1530505831754173883183325/127262478518793863168, 870050085900643643101819807455/30730325212890191726903296"

const _ORACLE_SUB2 = "1, 0, -1/2, 0, 3/8, 18/65, 116667/67600, 81/65, 2917283/540800, 5067999/1098500, 74278083393/4569760000, 35297721/2197000, 891309151207/18279040000, 4024650302019/74258600000, 4539353490598839/30891577600000, 5331013531227/29703440000, 1102241607284445603/2471326208000000, 1474592427340539003/2509940680000000, 28482207980919299424133/20882706457600000000, 9581764056581645577/5019881360000000, 112268656678699405401/26729864265728000, 2093957456007004388316747/339343979936000000000, 9186552703517718462451011789/705835478266880000000000, 13493486938426520745231963/678687959872000000000, 228982485018821858393468436551/5646683826135040000000000"

const _ORACLE_A0 = "1, 0, -1/2, 0, 3/8, 0, 27/16, 0, 675/128, 0, 3969/256, 0, 45927/1024, 0, 264627/2048, 0, 12196899/32768, 0, 70366725/65536, 0, 813439341/262144, 0, 4710988269/524288, 0, 109344517191/4194304"

function _parse_rat_list(s::AbstractString)
    map(split(s, ", ")) do tok
        parts = split(strip(tok), "/")
        length(parts) == 1 ? Rational{BigInt}(parse(BigInt, parts[1])) :
            Rational{BigInt}(parse(BigInt, parts[1])) // parse(BigInt, parts[2])
    end
end

function _to_rat(num)
    v = Symbolics.value(num)
    v isa Number && return Rational{BigInt}(v)
    v = Symbolics.value(Symbolics.simplify(num))
    v isa Number && return Rational{BigInt}(v)
    error("substituted coefficient is not numeric: $(repr(v))")
end

# --- Module-2b: eccentric generating-FUNCTION files (parse + materialize contract, ZERO [HOLE]) -----
# The 5 gen-func files now deserialize: the held gen-func head mints a FORMAL Symbolics atom, the Kerr
# `TagBox`/`RowBox` box-language is unwrapped, and the `mantissa`prec`/`*^`/`O[e]^n`/`Dot(.)` forms
# are handled. Grounding = the `.m` TEXT (NO Mathematica): closed-form leading coefficient == the
# literal; e→0 collapse to the Newtonian circular 32/5 (EInfinity/LInfinity); the formal atom lands
# in the exact SeriesData[y] slot it held; the Kerr TagBox content == the RowBox gen-func head; and
# every series materializes (Normal) to a Num. This replaces the former loud-raise [HOLE] contract.
@testset "MMAInputForm: Module-2b gen-func files (formal-atom / TagBox / materialize contract)" begin
    e = Symbolics.variable(:e)

    # (2/3) EInfinity leading SeriesData[y] == (32/5+(292 e²)/15+(37 e⁴)/15)/(1-e²)^(7/2): the .m
    # literal at e=1/3 (irrational √-carrying → Float64), and the e→0 collapse to Newtonian 32/5.
    dEi = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Schwarzschild", "Eccentric", "Flux", "EInfinity.m"))
    sEi = dEi["Series"]::PNSeries
    lit_Ei(ee) = (32 / 5 + (292 * ee^2) / 15 + (37 * ee^4) / 15) / (1 - ee^2)^(7 / 2)
    gotEi = Float64(Symbolics.value(Symbolics.simplify(Symbolics.substitute(sEi.coeffs[1], Dict(e => Symbolics.Num(1 / 3))))))
    @test isapprox(gotEi, lit_Ei(1 / 3); rtol = 1e-12)
    @test _to_rat(Symbolics.substitute(sEi.coeffs[1], Dict(e => Symbolics.Num(0)))) == 32 // 5

    # (2′) formal-atom placement: the SeriesData[y] slot (coeff 4, y^(3/2)) that held
    # HoldForm[oneP5PN[e,highestPower]] materializes to the formal atom `oneP5PN_e_highestPower`.
    @test isequal(sEi.coeffs[4], Symbolics.variable(:oneP5PN_e_highestPower))

    # (4) LInfinity likewise collapses to the Newtonian 32/5 at e→0.
    dLi = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Schwarzschild", "Eccentric", "Flux", "LInfinity.m"))
    sLi = dLi["Series"]::PNSeries
    @test _to_rat(Symbolics.substitute(sLi.coeffs[1], Dict(e => Symbolics.Num(0)))) == 32 // 5

    # (3) Redshift-y closed-form coeff == the .m literal (-2+4 e²)/(1-e²): exact rationals at e=0 and
    # e=1/2 (no √, so `_to_rat` compares exactly).
    dRy = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Schwarzschild", "Eccentric", "Local", "Redshift-y.m"))
    sRy = dRy["Series"]::PNSeries
    @test _to_rat(Symbolics.substitute(sRy.coeffs[3], Dict(e => Symbolics.Num(0)))) == -2 // 1
    @test _to_rat(Symbolics.substitute(sRy.coeffs[3], Dict(e => Symbolics.Num(big(1) // 2)))) == -4 // 3

    # (5) Kerr TagBox[RowBox[{"chi4PNConv","[",…}],HoldForm] box-language unwrap: the held head appears
    # as the formal atom `chi4PNConv_e_highestPower` / `chi5PNConv_e_highestPower` in the series.
    dK = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Kerr", "Eccentric", "Local", "Redshift.m"))
    sK = dK["Series"]::PNSeries
    @test any(c -> occursin("chi4PNConv_e_highestPower", string(c)), sK.coeffs)
    @test any(c -> occursin("chi5PNConv_e_highestPower", string(c)), sK.coeffs)

    # (contract) every gen-func file materializes (Normal) to a Num — ZERO [HOLE].
    for s in (sEi, sLi, sRy, sK)
        @test PostNewtonian.Series.Normal(s) isa Symbolics.Num
    end
end

# --- Module-2a: PURE eccentric files (nested-e / Rational / ResummedSeriesData / implicit-mult) -----
# The 4 pure eccentric files (no held generating-function) now deserialize. Grounding is the `.m`
# TEXT itself (NO Mathematica): each parses; a leading coefficient reproduces the literal `.m`
# formula; the e→0 limit collapses `(1-e^2)^k → 1` onto the known circular Newtonian flux 32/5; and a
# nested `SeriesData[e,…]` / `ResummedSeriesData` coefficient carries the literal `Rational`/`Log`/
# `EulerGamma` atoms. `_atom_rational`/`_atom_resummed` are unit-checked against their MMA-head rule.
@testset "MMAInputForm: Module-2a pure eccentric (Rational/Resummed/nested-e/implicit-mult)" begin
    e  = Symbolics.variable(:e)
    EG = Symbolics.variable(:EulerGamma); PI = Symbolics.variable(:Pi)
    L2 = Symbolics.variable(:Log2); L3 = Symbolics.variable(:Log3)
    Lp = Symbolics.variable(:Logp); Ly = Symbolics.variable(:Logy)

    # (unit) `Rational[p,q]` → exact p//q; `ResummedSeriesData[sf, SeriesData[e,..]]` → sf·Normal(sd)
    # with the prefactor collapsing to 1 at e=0 (leaving inner c0 = 5).
    @test Rational{BigInt}(Symbolics.value(read_mma_inputform_symbolic("Rational[-587831, 2880]"))) == -587831 // 2880
    r = read_mma_inputform_symbolic("ResummedSeriesData[1/(1-e^2)^4, SeriesData[e, 0, {5, 0, 7}, 0, 3, 1]]")
    @test isapprox(Float64(Symbolics.value(Symbolics.simplify(Symbolics.substitute(r, Dict(e => Symbolics.Num(0)))))), 5.0; rtol = 1e-12)

    # (2) EHorizon leading SeriesData[y] coefficient == the literal `.m` formula
    #   (32/5+(496 e^2)/5+204 e^4+74 e^6+(5 e^8)/2)/(1-e^2)^(15/2).  Two independent checks:
    #   a rational spot value e=1/3, and the e→0 circular collapse to 32/5.
    dE = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Schwarzschild", "Eccentric", "Flux", "EHorizon.m"))
    sE = dE["Series"]::PNSeries
    @test sE.var == :y
    # e=1/3 is irrational ((1-e²)^(15/2) carries a √), so compare as Float64.
    lit_lead(ee) = (32 / 5 + (496 * ee^2) / 5 + 204 * ee^4 + 74 * ee^6 + (5 * ee^8) / 2) / (1 - ee^2)^(15 / 2)
    gotE = Float64(Symbolics.value(Symbolics.simplify(Symbolics.substitute(sE.coeffs[1], Dict(e => Symbolics.Num(1 / 3))))))
    @test isapprox(gotE, lit_lead(1 / 3); rtol = 1e-12)
    @test _to_rat(Symbolics.substitute(sE.coeffs[1], Dict(e => Symbolics.Num(0)))) == 32 // 5   # e→0

    # (3) LHorizon e→0 leading coefficient likewise collapses to the Newtonian 32/5.
    dL = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Schwarzschild", "Eccentric", "Flux", "LHorizon.m"))
    sL = dL["Series"]::PNSeries
    @test _to_rat(Symbolics.substitute(sL.coeffs[1], Dict(e => Symbolics.Num(0)))) == 32 // 5

    # (2′) Spin-p nested SeriesData[e,0,{Rational[..] + Rational[..] EulerGamma + ..}] first coeff:
    # its e^0 term matches the literal Rational/EulerGamma/Pi/Log/Logp combination exactly.
    dP = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Schwarzschild", "Eccentric", "Local", "Spin-p.m"))
    sP = dP["Series"]::PNSeries
    cP0 = Symbolics.substitute(sP.coeffs[7], Dict(e => Symbolics.Num(0)))   # e^0 term = inner c0
    expP = -587831 // 2880 + (1256 // 15) * EG + (31697 // 6144) * PI^2 +
           (296 // 15) * L2 + (729 // 5) * L3 + (-628 // 15) * Lp
    @test isequal(Symbolics.simplify(cP0 - expP), Symbolics.Num(0))

    # (3′) Spin-y ResummedSeriesData coefficient: at e=0 the prefactor 1/(1-e^2)^4 → 1, so (with the
    # separate Logy log-running term zeroed) the coeff equals the inner c0 = the 5-term literal.
    dY = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Schwarzschild", "Eccentric", "Local", "Spin-y.m"))
    sY = dY["Series"]::PNSeries
    cY0 = Symbolics.substitute(sY.coeffs[7], Dict(e => Symbolics.Num(0), Ly => Symbolics.Num(0)))
    expY = -587831 // 2880 + (1256 // 15) * EG + (31697 // 6144) * PI^2 + (296 // 15) * L2 + (729 // 5) * L3
    @test isequal(Symbolics.simplify(cY0 - expY), Symbolics.Num(0))

    # (4) Symbolics round-trip: every pure-eccentric series materializes (Normal) to a Num.
    for s in (sE, sL, sP, sY)
        @test all(c -> c isa Symbolics.Num, s.coeffs)
        @test PostNewtonian.Series.Normal(s) isa Symbolics.Num
    end
end

@testset "MMAInputForm: Energy.m vs wolframscript byte-exact oracle" begin
    path = joinpath(SERIES_DATA_DIR, "Kerr", "Spherical", "Orbit", "Energy.m")
    d = parse_pn_data_file(path)
    @test d isa Dict{String,Any}
    @test d["Name"] == "Kerr Spherical Orbit Energy"

    s = d["Series"]::PNSeries
    # series structure (oracle: NMIN NMAX DEN NCOEFF = 0 25 2 25, var p, x0 Infinity)
    @test s.var == :p
    @test (s.nmin, s.nmax, s.den) == (0, 25, 2)
    @test length(s.coeffs) == 25

    a = Symbolics.variable(:a)
    x = Symbolics.variable(:x)

    got1 = [_to_rat(Symbolics.substitute(c, Dict(a => Symbolics.Num(big(3) // 7),
                                                 x => Symbolics.Num(big(5) // 11)))) for c in s.coeffs]
    @test got1 == _parse_rat_list(_ORACLE_SUB1)

    got2 = [_to_rat(Symbolics.substitute(c, Dict(a => Symbolics.Num(big(-2) // 5),
                                                 x => Symbolics.Num(big(9) // 13)))) for c in s.coeffs]
    @test got2 == _parse_rat_list(_ORACLE_SUB2)

    got0 = [_to_rat(Symbolics.substitute(c, Dict(a => Symbolics.Num(big(0) // 1)))) for c in s.coeffs]
    @test got0 == _parse_rat_list(_ORACLE_A0)
    @test got0[1:5] == Rational{BigInt}[1, 0, -1//2, 0, 3//8]   # published Schwarzschild leading orders
end

# Recursively evaluate a fully-substituted Num (concrete numeric leaves; the registered `polygamma`
# term and `sqrt` do not auto-fold under `substitute`, so we apply each operation directly).
function _fold(v)
    v = Symbolics.value(v)
    v isa Number && return v
    Symbolics.iscall(v) || error("_fold: residual symbol $v")
    op = Symbolics.operation(v)
    op(map(_fold, Symbolics.arguments(v))...)
end

@testset "MMAInputForm: complex-AST atom class (PolyGamma/RePsi/ImPsi)" begin
    EG = Base.MathConstants.eulergamma
    # PolyGamma[m,z] -> SpecialFunctions.polygamma(m,z); analytic identities ground the target.
    pg01 = read_mma_inputform_symbolic("PolyGamma[0, 1]")
    @test Symbolics.operation(Symbolics.value(pg01)) === SpecialFunctions.polygamma
    @test isapprox(SpecialFunctions.polygamma(0, 1), -EG; atol = 1e-12)   # ψ(1) = -γ
    @test isapprox(SpecialFunctions.polygamma(1, 1), pi^2 / 6; atol = 1e-12)  # ψ'(1) = π²/6

    # RePsi[m,n]/ImPsi[m,n] reproduce the stored "ExtraFunctions" definition exactly: the symbolic
    # form is ψ^(m)(1 + 𝕚 a n/√(1-a²)) ± ψ^(m)(1 - 𝕚 a n/√(1-a²)) (over 𝕚 for Im). Check that
    # substituting 𝕚=>im, a=>numeric makes RePsi real and ImPsi the imaginary combination.
    a = Symbolics.variable(:a); iu = Symbolics.variable(:𝕚)
    rp = read_mma_inputform_symbolic("RePsi[0, 1]")
    ip = read_mma_inputform_symbolic("ImPsi[0, 1]")
    sub = Dict(a => 0.3, iu => im)
    rpv = ComplexF64(_fold(Symbolics.substitute(rp, sub)))
    ipv = ComplexF64(_fold(Symbolics.substitute(ip, sub)))
    z = 1 + im * 0.3 / sqrt(1 - 0.3^2)
    @test isapprox(rpv, SpecialFunctions.polygamma(0, z) + SpecialFunctions.polygamma(0, conj(z)); atol = 1e-10)
    @test isapprox(real(rpv), 2 * real(SpecialFunctions.polygamma(0, z)); atol = 1e-10)
    @test abs(imag(rpv)) < 1e-10               # RePsi is real
    @test isapprox(ipv, 2 * imag(SpecialFunctions.polygamma(0, z)); atol = 1e-10)

    # Memory file byte-exact: leading SeriesData coeff (-2 √(10π/3))/7 at y² .. y^13 (den 2).
    d = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Kerr", "Circular", "Memory", "h-l2m0.m"))
    s = d["Series"]::PNSeries
    @test (s.var, s.nmin, s.nmax, s.den) == (:y, 2, 13, 2)
    PI = Symbolics.variable(:Pi)
    c1 = Float64(_fold(Symbolics.substitute(s.coeffs[1], Dict(PI => Float64(pi)))))
    @test isapprox(c1, (-2 * sqrt(10 * pi / 3)) / 7; rtol = 1e-12)
end

@testset "MMAInputForm: Module-1 complex-atom class (PsiA/PsiB/Psibar/lowercase-i/e,Y)" begin
    a = Symbolics.variable(:a); iu = Symbolics.variable(:𝕚); q = Symbolics.variable(:q)

    # (b) PsiA[n,m]/PsiB[n,m] reproduce the stored "ExtraFunctions" base-point-3 definition:
    #   PsiA = (ψ^(n)(3 + m a 𝕚/√(1-a²)) + ψ^(n)(3 - m a 𝕚/√(1-a²)))/2 ;  PsiB = (diff)/2/𝕚.
    # Oracle = the analytic digamma identity (same polygamma oracle as the RePsi/ImPsi test).
    aval = 0.3
    for (n, m) in ((0, 2), (1, 2), (0, 1))
        pa = read_mma_inputform_symbolic("PsiA[$n, $m]")
        pb = read_mma_inputform_symbolic("PsiB[$n, $m]")
        sub = Dict(a => aval, iu => im)
        pav = ComplexF64(_fold(Symbolics.substitute(pa, sub)))
        pbv = ComplexF64(_fold(Symbolics.substitute(pb, sub)))
        z = 3 + m * aval * im / sqrt(1 - aval^2)
        pg(k) = SpecialFunctions.polygamma(n, k)
        @test isapprox(pav, (pg(z) + pg(conj(z))) / 2; atol = 1e-10)
        @test abs(imag(pav)) < 1e-10                     # PsiA is real
        @test isapprox(pbv, (pg(z) - pg(conj(z))) / 2 / im; atol = 1e-10)
        @test abs(imag(pbv)) < 1e-10                     # PsiB is real (the diff/𝕚 is real)
    end

    # (e) lowercase `i` is the imaginary unit: the PolyGamma arg folds at the complex Kerr point.
    li = read_mma_inputform_symbolic("PolyGamma[0, (i*q)/Sqrt[1 - q^2]]")
    lival = ComplexF64(_fold(Symbolics.substitute(li, Dict(q => 0.4, iu => im))))
    @test isapprox(lival, SpecialFunctions.polygamma(0, im * 0.4 / sqrt(1 - 0.16)); atol = 1e-10)

    # `HoldForm[atom]` unwraps a held SCALAR atom to the atom (identity); a held eccentric gen-func
    # head mints the FORMAL Symbolics atom `oneP5PN_e_4` (Module-2b) rather than raising.
    @test isequal(read_mma_inputform_symbolic("HoldForm[PsiA[0, 2]]"), read_mma_inputform_symbolic("PsiA[0, 2]"))
    @test isequal(read_mma_inputform_symbolic("HoldForm[oneP5PN[e, 4]]"), Symbolics.variable(:oneP5PN_e_4))

    # `Psibar[n,m]` is a named formal PN variable `Psibar_n_m` (symbol identity only).
    pbar = read_mma_inputform_symbolic("Psibar[1, 2]")
    @test isequal(pbar, Symbolics.variable(Symbol("Psibar_1_2")))

    # (a) census spot-check + (c) coefficient literal: EInfinity-v e→0 leading coeff = 32/5
    # (Newtonian quadrupole; the (1-e²)^{3/2}·96/15 prefactor collapses to 32/5), and Spin.m leads
    # with the SeriesData literal {q, 0, 1, ...}.
    e = Symbolics.variable(:e)
    d = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Kerr", "Generic", "Flux", "EInfinity-v.m"))
    s = d["Series"]::PNSeries
    @test s.var == :v
    c0 = Float64(_fold(Symbolics.substitute(s.coeffs[1], Dict(e => 0.0))))
    @test isapprox(c0, 32 / 5; rtol = 1e-12)

    dsp = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Kerr", "Circular", "Local", "Spin.m"))
    ssp = dsp["Series"]::PNSeries
    @test isequal(ssp.coeffs[1], q)                       # SeriesData literal leads with `q`
    @test Rational{BigInt}(Symbolics.value(ssp.coeffs[3])) == 1

    # (d) Symbolics round-trip on a RedShift Psibar-bearing coefficient: parse→materialize stays a Num.
    drs = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Kerr", "Circular", "Local", "RedShift.m"))
    srs = drs["Series"]::PNSeries
    @test srs.coeffs[end] isa Symbolics.Num             # highest order carries the \[Psi]bar atoms
end

@testset "MMAInputForm: read_mma_inputform_symbolic low-level walker" begin
    # SeriesData -> PNSeries
    s = read_mma_inputform_symbolic("SeriesData[p, 0, {1, 0, -1/2}, 0, 3, 1]")
    @test s isa PNSeries
    @test (s.var, s.nmin, s.nmax, s.den) == (:p, 0, 3, 1)
    # exact rational kept (not Float64): coeff[3] == -1//2
    @test Rational{BigInt}(Symbolics.value(s.coeffs[3])) == -1 // 2
    # bare Association -> Dict
    d = read_mma_inputform_symbolic("<|\"Name\" -> \"x\", \"Series\" -> SeriesData[y, 0, {1}, 0, 1, 1]|>")
    @test d isa Dict{String,Any}
    @test d["Name"] == "x"
    @test d["Series"] isa PNSeries
end

@testset "MMAInputForm: arbitrary precision/accuracy and scientific literals" begin
    # Plain decimals retain the established machine-number behavior.
    ordinary = Symbolics.value(read_mma_inputform_symbolic("1.234567890123456789"))
    @test ordinary isa Float64
    @test ordinary == parse(Float64, "1.234567890123456789")

    precise = Symbolics.value(read_mma_inputform_symbolic("1.2345678901234567890123456789`40"))
    @test precise isa BigFloat
    @test precision(precise) >= ceil(Int, 40 * log2(10))
    setprecision(BigFloat, 256) do
        ref = parse(BigFloat, "1.2345678901234567890123456789")
        @test abs(BigFloat(precise) - ref) < big"1e-30"
    end

    accurate = Symbolics.value(read_mma_inputform_symbolic("0.0012345678901234567890123456789``30"))
    @test accurate isa BigFloat
    setprecision(BigFloat, 256) do
        ref = parse(BigFloat, "0.0012345678901234567890123456789")
        @test abs(BigFloat(accurate) - ref) < big"1e-30"
    end

    scientific = Symbolics.value(
        read_mma_inputform_symbolic("1.2345678901234567890123456789`50*^25"))
    @test scientific isa BigFloat
    setprecision(BigFloat, 256) do
        ref = parse(BigFloat, "1.2345678901234567890123456789e25")
        @test abs((BigFloat(scientific) - ref) / ref) < big"1e-29"
    end

    machine_scientific = Symbolics.value(read_mma_inputform_symbolic("1.25*^25"))
    @test machine_scientific isa Float64
    @test machine_scientific == 1.25e25

    @test_throws ArgumentError read_mma_inputform_symbolic("1.25``")
end

# --- Module-3: giant-sum balanced materialization ------------------------------------------
# The 8 Spherical flux/scalar-flux giants carry ~35k-term signed coefficient sums that Meta.parse
# renders as a deeply left-nested `±` tree (mixed `+`/`-` breaks Julia's flat-`+` folding); plain
# recursion overflowed the stack on the ~3.4MB "Horizon" files. `_accumulate_balanced` over the
# iteratively-flattened chain rebuilds each coefficient with `log2(N)` tree depth.
@testset "MMAInputForm: Module-3 balanced giant-sum materialization" begin
    MM = PostNewtonian.MMAInputForm

    # (b) balanced == left-fold invariance (the reordering is value-safe before applying to giants).
    a = Symbolics.variable(:a)
    v = Symbolics.Num[Symbolics.Num(big(i) // 7) + (i % 3) * a for i in 1:37]
    @test isequal(Symbolics.simplify(MM._accumulate_balanced(v) - sum(v)), Symbolics.Num(0))
    @test isequal(MM._accumulate_balanced(Symbolics.Num[]), Symbolics.Num(0))       # empty → 0
    @test isequal(MM._accumulate_balanced(Symbolics.Num[a]), a)                      # singleton

    # (a) all 8 giants materialize (parsed in-census above); re-parse the small Infinity giant here
    # for the numeric spot-checks.
    d = parse_pn_data_file(joinpath(SERIES_DATA_DIR, "Kerr", "Spherical", "Flux", "Energy", "Infinity.m"))
    s = d["Series"]::PNSeries
    @test (s.var, s.nmin, s.nmax, s.den) == (:p, 10, 35, 2)

    # (c) published Schwarzschild a→0 leading flux: the p^10 coefficient reduces to 32/5 (Newtonian
    # quadrupole energy flux to infinity).
    c0 = Symbolics.value(Symbolics.simplify(Symbolics.substitute(s.coeffs[1],
        Dict(a => Symbolics.Num(big(0) // 1)))))
    @test c0 == 32 // 5

    # (a′) RePsi/ImPsi digamma atoms fold correctly THROUGH the balanced sum: a substituted giant
    # coefficient (3.5PN, RePsi/ImPsi-bearing) evaluates to a finite REAL number — the imaginary
    # ψ-combinations cancel exactly, as a physical flux coefficient must.
    x = Symbolics.variable(:x); iu = Symbolics.variable(:𝕚)
    sub = Dict(a => 0.35, x => 0.6, iu => im,
               Symbolics.variable(:Pi) => Float64(pi),
               Symbolics.variable(:EulerGamma) => Float64(Base.MathConstants.eulergamma),
               Symbolics.variable(:Logp) => 0.0,
               Symbolics.variable(:Log2) => log(2.0), Symbolics.variable(:Log3) => log(3.0),
               Symbolics.variable(:Log5) => log(5.0), Symbolics.variable(:Log7) => log(7.0),
               Symbolics.variable(:Zeta3) => 1.2020569031595942)
    val = ComplexF64(_fold(Symbolics.substitute(s.coeffs[7], sub)))
    @test isfinite(real(val))
    @test abs(imag(val)) < 1e-6                       # digamma imaginary parts cancel (real flux)
end
