# C-LI-0: actual emitted eccentric kernels, with no circular metric restriction.
# Run in derivation/; optional argument is an existing/creatable evidence directory.
using Symbolics, Test, SHA
println("[START] C-LI-0 exact symbolic overlap"); flush(stdout)
include(joinpath(@__DIR__, "..", "src", "R3Derivation.jl"))
const R3 = R3Derivation
const LI = R3.LocallyIsothermal
const HERE = normpath(joinpath(@__DIR__, ".."))
const OUT = isempty(ARGS) ? mktempdir() : mkpath(ARGS[1])
@variables Theta Theta_prime
const OVERLAP = Dict(R3.p => Theta*R3.eps, R3.p_r => Theta*R3.eps_r,
                     R3.Gamma1 => 1+Theta, Theta_prime => Num(0))
const NAMES = ("C00", "C01", "C10", "C11", "S0", "S1",
               "dut", "duphi_up", "deps", "duphi_low")

# Reuse the arithmetic interpreter on the emitted files; no generated bodies
# are transcribed here. Symbolic amplitude indeterminates are never conjugated
# or projected onto real inputs, so the identity extends holomorphically to
# arbitrary complex metric jets, including htr and hrphi.
function emitted_overlap_expressions(tag)
    path = joinpath(HERE, "..", "src", "SpiralWaves", "generated", "eccentric_$tag.jl")
    artifact = joinpath(HERE, "out", "ECC_$tag.jls.txt")
    regenerated = joinpath(OUT, "eccentric_$tag.jl")
    R3.emit_eccentric(regenerated, artifact; tag)
    @test read(path) == read(regenerated)
    println("FRESH $tag sha256=", bytes2hex(sha256(read(path))))
    table = Dict{Symbol,Any}(R3.ARTIFACT_ALPHABET)
    table[:im] = im
    trees = Dict(R3.emitted_expressions(path))
    Dict(name => R3.eval_tree(trees["_eccentric_$(tag)_$name/return"], table) for name in NAMES)
end

@testset "C-LI-0 exact symbolic certificate" begin
    # Both selected equations reduce to the same Eulerian pressure law. Check
    # each separately, rather than only checking their mutual difference.
    common = R3.dp/(Theta*R3.eps) - R3.deps/R3.eps
    lagrangian_li = (R3.dp + R3.XI_R*R3.p_r -
        Theta*(R3.deps+R3.XI_R*R3.eps_r) - R3.eps*Theta_prime*R3.XI_R)/R3.p
    residuals = Pair{String,Any}[
        "adiabatic_closure" => R3.subs(R3.CLOSURE_RESIDUAL, OVERLAP)-common,
        "li_closure" => R3.subs(LI.CLOSURE_RESIDUAL, OVERLAP)-common,
        "li_lagrangian_closure" => R3.subs(lagrangian_li, OVERLAP)-common,
        "shared_energy_zero_bath" => LI.DENERGY-R3.DENERGY,
        "reduced_energy_overlap" => R3.subs(R3.EQ_E-LI.EQ_E, OVERLAP)]
    @test LI.DENERGY === R3.DENERGY
    ad, li = emitted_overlap_expressions.( ("adiabatic", "li") )
    for name in NAMES
        push!(residuals, "emitted_$name" => R3.subs(ad[name]-li[name], OVERLAP))
    end
    for (tag, closure, got) in (("adiabatic",R3.CLOSURE_RESIDUAL,ad),
                                 ("li",LI.CLOSURE_RESIDUAL,li))
        solved = R3.subs(closure,Dict(R3.dp => (R3.eps+R3.p)*R3.y0,
                         R3.deps => got["deps"],R3.dur => im*R3.y1))
        push!(residuals,"solved_$(tag)_closure" => R3.subs(solved,OVERLAP))
    end
    # Compare against the canonical derivation exports as well as freshness
    # against the serialized artifact, so stale artifacts cannot self-certify.
    for (tag, exports, got) in (("adiabatic", LI.closure_exports(R3), ad),
                                 ("li", R3.li_response_exports(), li))
        for (name, ex) in exports
            push!(residuals, "derivation_$(tag)_$name" => got[name]-ex)
        end
    end
    open(joinpath(OUT, "symbolic_residuals.tsv"), "w") do io
        println(io, "identity\treal_residual\timag_residual")
        for (name, residual) in residuals
            println("SIMPLIFY ", name); flush(stdout)
            parts = map(x -> Symbolics.simplify_fractions(Num(x)), (real(residual),imag(residual)))
            zero_parts = map(x -> R3._zero_num(Num(x)),parts)
            @test all(zero_parts)
            # Log the verifier's canonical zero, not an intermediate fraction
            # whose numerator still has unexpanded cancellations.
            parts = map((s,z) -> z ? Num(0) : s,parts,zero_parts)
            println(io, name, '\t', parts[1], '\t', parts[2]); flush(io)
        end
    end
    # Independent complex amplitudes are retained in the emitted expressions.
    used = union((LI._names(e) for e in values(ad))...)
    @test all(string(h) in used for h in R3.HSYMS[1:7])
    println("METRIC_INPUTS ",join(sort!(collect(intersect(used,Set(string.(R3.HSYMS))))),','))

    # Non-identities are rejected by a single exact counterexample, not by
    # asking a symbolic simplifier to finish on an enormous nonzero expression.
    # Omitting Theta_prime=0 while still explicitly imposing p'=Theta*eps'
    # is redundant in the emitted alphabet. The meaningful thermal mutant uses
    # p'=Theta*eps'+Theta_prime*eps consistently with p=Theta*eps.
    point = R3.complex_emission_point(R3.MersenneTwister(2053))
    point[Theta] = 2//7; point[Theta_prime] = 1//11
    mutants = [
        "drop_p_law" => filter(kv -> !isequal(first(kv), R3.p), OVERLAP),
        "drop_p_prime_law" => filter(kv -> !isequal(first(kv), R3.p_r), OVERLAP),
        "drop_Gamma1" => filter(kv -> !isequal(first(kv), R3.Gamma1), OVERLAP),
        "perturb_Gamma1" => merge(OVERLAP, Dict(R3.Gamma1 => 1+Theta+1//7)),
        "nonconstant_Theta" => merge(OVERLAP, Dict(Theta_prime => Theta_prime,
                                R3.p_r => Theta*R3.eps_r+Theta_prime*R3.eps))]
    open(joinpath(OUT, "mutant_residuals.tsv"), "w") do io
        println(io, "mutant\tclosure_residual\tnonzero_emitted_coefficients\temitted_residuals")
        for (name, substitutions) in mutants
            closure = R3.exact_eval(R3.subs(R3.CLOSURE_RESIDUAL-LI.CLOSURE_RESIDUAL, substitutions), point)
            values = [R3.exact_eval(R3.subs(ad[n]-li[n], substitutions), point) for n in NAMES]
            changed = [n for (n,v) in zip(NAMES,values) if !iszero(v)]
            @test !iszero(closure)
            @test !isempty(changed)
            println(io, name, '\t', closure, '\t', join(changed, ','), '\t', values)
            println("MUTANT ", name, " REJECTED; nonzero=", join(changed, ',')); flush(stdout)
        end
    end
end
println("[PASS] C-LI-0 37 complex symbolic residuals = 0; 5 mutants rejected")

# Retain the counterexample to the FORMER connected-Reynolds-only convention.
# R3 diagnostics.jl:193-201 defines the adiabatic Noether integrand. The LI
# accessor now includes the pressure current (owner decision 2026-10-13,
# gr-li-o-amf-definition); certify_li_amf_definition.jl proves the conserved
# baryon-current split for every closure. The old Reynolds-only term still
# differs from (72a) when the pressure laws and reconstructed solutions coincide.
@testset "Former Reynolds-only AMF gap under the exact closure overlap" begin
    ad_flux = R3.onshell_nu(R3.to_omega_d(R3.AMF_MODE))
    li_flux = (R3.eps+R3.p)*(im*R3.y1)*R3.cbar(LI.ELIMINATED_ONSHELL.duphi_low)
    pressure = R3.ULOW[R3.IPHI]*R3.DP_M*R3.cbar(R3.DUUP_M[R3.IR])
    identity = R3.subs(R3.re_part(ad_flux-li_flux-pressure),OVERLAP)
    @test R3._zero_num(Num(real(identity))) && R3._zero_num(Num(imag(identity)))
    point = R3.complex_emission_point(R3.MersenneTwister(2054))
    point[Theta] = 2//7
    for (a,b) in R3.CONJ_PAIRS
        point[b] = conj(point[a])
    end
    gap = R3.exact_eval(R3.subs(R3.re_part(ad_flux-li_flux),OVERLAP),point)
    term = R3.exact_eval(R3.subs(R3.re_part(pressure),OVERLAP),point)
    @test gap == term && !iszero(gap)
    open(joinpath(OUT,"amf_observable_gap.txt"),"w") do io
        println(io,"Noether - connected Reynolds - u_phi Re(dp conjugate(dur)) = 0 exactly.")
        println(io,"Common external factor: 4pi sqrt(-g). Physical conjugacy imposed on the witness.")
        println(io,"Nonzero exact witness (without common factor): ",gap)
        println(io,"This rejects the former Reynolds-only definition; both native accessors now include this current.")
    end
end
println("[PASS] former Reynolds-only gap identified; conserved native AMF parity has its own runtime gate")
