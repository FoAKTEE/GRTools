"General adiabatic C/S, reconstruction and diagnostics; no circular zeroing."
general_response_exports() = vcat(LocallyIsothermal.closure_exports(@__MODULE__), diagnostics_exports())

"Emit closure-tagged scalar kernels using the existing R3 artifact reader/printer."
function emit_eccentric(path, artifact; tag)
    entries = read_artifact(artifact)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"# GENERATED: closure = $tag; do not edit by hand.")
        println(io,"# generator: derivation/src/eccentric_emission.jl")
        println(io,"# input: $(basename(artifact)) sha256 $(_sha(artifact))")
        println(io,"# derivation: sha256 $(derivation_sha())")
        println(io,"# Inputs use the R3 alphabet (enu, epsi, emu are exponentials).")
        println(io,"# htr/hrphi and barred amplitudes are independent complex inputs.")
        println(io,"# Physical conjugacy is supplied by the caller; no reality projection is applied.")
        println(io,"# The background is on shell; omega_d != 0 and all algebraic pivots must be nonzero.")
        println(io,"# Scalar kernels accept a NamedTuple or a property-compatible input struct.")
        for (name, ex) in entries
            println(io,"\nfunction _eccentric_$(tag)_$(name)(inp)")
            for field in sort!(collect(symbol_names(ex)); by=string)
                println(io,"    $field = inp.$field")
            end
            println(io,"    return ",print_value(ex))
            println(io,"end")
        end
    end
    path
end

"Arithmetic tree of a general complex R3 expression (amplitudes remain formal)."
emission_tree(e) = Meta.parse(print_value((Num(real(e)),Num(imag(e)))))

"""
Write an explicit straight-line representation N/D of a cleared rational
expression, with conservative total-degree bounds. No expansion of a large
multivariate numerator is needed. Every division records its excluded locus.
Leaves v,W use the canonical rational conic; all other leaves are independent.
"""
function rational_circuit(io, tree)
    serial = Ref(0)
    function node(ex)
        # Samples sit on hydrostatic balance. Substitute that rational
        # parametrization before counting numerator/denominator degrees.
        ex === :p_r && return node(emission_tree((eps+p)*HYDRO_RHS))
        if ex isa Number || ex isa Symbol
            serial[] += 1
            id = serial[]
            n,d,dn,dd = ex === :v ? ("2*tpar","1+tpar^2",1,2) :
                ex === :W ? ("1+tpar^2","1-tpar^2",2,2) :
                (string(ex),"1",ex isa Symbol && ex !== :im ? 1 : 0,0)
            println(io,"N$id = $n; D$id = $d; degrees <= ($dn,$dd)")
            return (id,dn,dd)
        end
        ex isa Expr && ex.head === :call || error("non-rational tree: $ex")
        op = ex.args[1]
        if op === :^
            a,an,ad = node(ex.args[2])
            k = Int(eval_tree(ex.args[3],Dict{Symbol,Any}()))
            serial[] += 1; id = serial[]
            if k >= 0
                n,d,dn,dd = "N$a^$k","D$a^$k",k*an,k*ad
            else
                n,d,dn,dd = "D$a^$(-k)","N$a^$(-k)",-k*ad,-k*an
                println(io,"EXCLUDE N$a = 0")
            end
            println(io,"N$id = $n; D$id = $d; degrees <= ($dn,$dd)")
            return (id,dn,dd)
        end
        op in (:+,:-,:*,:/,://) || error("unsupported rational operation $op")
        a,an,ad = node(ex.args[2])
        if length(ex.args)==2
            op === :- || error("unsupported unary operation")
            serial[] += 1; id = serial[]
            println(io,"N$id = -N$a; D$id = D$a; degrees <= ($an,$ad)")
            return (id,an,ad)
        end
        for arg in ex.args[3:end]
            b,bn,bd = node(arg)
            serial[] += 1; id = serial[]
            if op in (:+,:-)
                n,d,dn,dd = "N$a*D$b $op N$b*D$a","D$a*D$b",max(an+bd,bn+ad),ad+bd
            elseif op === :*
                n,d,dn,dd = "N$a*N$b","D$a*D$b",an+bn,ad+bd
            else
                n,d,dn,dd = "N$a*D$b","D$a*N$b",an+bd,ad+bn
                println(io,"EXCLUDE N$b = 0")
            end
            println(io,"N$id = $n; D$id = $d; degrees <= ($dn,$dd)")
            a,an,ad = id,dn,dd
        end
        (a,an,ad)
    end
    id,dn,dd = node(tree)
    println(io,"ROOT N$id / D$id; EXCLUDE D$id = 0; numerator degree <= $dn; denominator degree <= $dd")
    (dn,dd)
end

"Independent complex amplitudes, including independent barred symbols."
function complex_emission_point(rng)
    pt = rational_point(rng; onshell=true)
    for s in vcat([y0,y1],HSYMS,[b for (_,b) in CONJ_PAIRS])
        pt[s] = complex(_rand_signed(rng),_rand_signed(rng))
    end
    pt
end

"IID grid parameters followed by rational conic/hydrostatic substitutions; no rejection."
function uniform_emission_point(rng, size::Int)
    draw() = RatB(rand(rng,1:size),size)
    pt = Dict{Num,Any}(s=>draw() for s in ALL_SYMBOLS)
    t = draw()/2
    pt[v]=2t/(1+t^2); pt[W]=(1+t^2)/(1-t^2)
    pt[m]=RatB(rand(rng,1:size))
    pt[omega]=pt[omega_d]+pt[m]*exact_eval(Ω_expr,pt)
    pt[p_r]=(pt[eps]+pt[p])*exact_eval(HYDRO_RHS,pt)
    for s in vcat([y0,y1],HSYMS,[b for (_,b) in CONJ_PAIRS])
        pt[s]=complex(draw(),draw())
    end
    pt
end

"Both a compiled-kernel input and the arithmetic interpreter's symbol table."
function emission_inputs(pt)
    names = sort!(collect(keys(pt)); by=string)
    inp = NamedTuple{Tuple(Symbol(string(s)) for s in names)}(Tuple(pt[s] for s in names))
    inp, Dict{Symbol,Any}(:im=>im, (Symbol(string(s))=>val for (s,val) in pt)...)
end

"Call a newly included generated kernel using the current Julia world."
emitted_value(mod,tag,name,inp) = Base.invokelatest(
    Base.invokelatest(getproperty,mod,Symbol("_eccentric_$(tag)_$name")),inp)

"""
Determinism, freshness, exact-rational round-trip and actual compiled complex kernels.
The nonuniform ratio samples test complex execution and precision. A separate
IID grid test carries a Schwartz-Zippel bound for the cleared numerator on the
declared rational constraint surface. A denominator zero fails the run; there
is no rejection or conditioning on accepted points.
"""
function emission_certificates(path, artifact, exports; tag, npts=12, seed=260817,
                               accounting_dir, verbose=true)
    npts > 0 || throw(ArgumentError("npts must be positive"))
    C = Certificate[]
    function check(name,ok,note; modality=:structural)
        push!(C,certify_bool("$tag $name","ECC emission",ok,note;
            verbose,modality,npoints=modality===:exact_rational ? npts : 0))
        verbose && flush(stdout)
    end
    mkpath(accounting_dir)
    mktempdir() do tmp
        art = joinpath(tmp,basename(artifact))
        write_exports(art,exports;header="# ECC closure=$tag; general complex amplitudes")
        p1 = emit_eccentric(joinpath(tmp,"one.jl"),art;tag)
        p2 = emit_eccentric(joinpath(tmp,"two.jl"),art;tag)
        check("determinism and freshness",read(p1)==read(p2)==read(path) && read(art)==read(artifact),
              "two regenerated files, shipped target, and expression artifact are byte-identical")
    end
    refs = Dict(n=>e for (n,e) in exports)
    trees = emitted_expressions(path)
    table = Dict{Symbol,Any}(ARTIFACT_ALPHABET); table[:im]=im
    exact = length(trees)==length(refs)
    degrees = Dict{String,Tuple{Int,Int}}()
    open(joinpath(accounting_dir,"$(tag)-rational-circuits.txt"),"w") do io
        println(io,"v=2t/(1+t^2), W=(1+t^2)/(1-t^2); t!=+/-1; all other alphabet leaves independent.")
        for (key,tree) in trees
            name = replace(split(key,"/")[1],"_eccentric_$(tag)_"=>"")
            got = eval(_artifact_node(tree,table))
            want = refs[name]
            same = isequal(Num(real(got)),Num(real(want))) && isequal(Num(imag(got)),Num(imag(want)))
            println(io,"\nEXPORT $name: EMITTED MINUS ORIGINAL")
            degrees[name] = rational_circuit(io,Expr(:call,:-,tree,emission_tree(want)))
            println(io,"Structural identity recognized = $same; exact residuals are in samples.txt.")
        end
    end
    check("rational alphabet and export coverage",exact,"$(length(trees)) parsed complex expressions; cleared difference numerator circuits and degree bounds recorded")
    mod = Module(gensym(:EccentricEmission))
    Base.include(mod,path)
    rng = MersenneTwister(seed)
    compiled, precision = true,true
    open(joinpath(accounting_dir,"$(tag)-samples.txt"),"w") do io
        println(io,"N=$npts seed=$seed; exact Rational{BigInt} real/imaginary parts, no rejected points.")
        println(io,"S support: signed a/b, a,b in 1:37; positive metric exponentials,p,eps,Gamma1; m in 1:9.")
        println(io,"Conic: signed t=a/b; a in 1:29, b in (a+1):30. omega=omega_d+m*Omega; hydrostatic p_r computed.")
        println(io,"Excluded: enu*epsi*emu*p*eps*(eps+p)*W*omega_d=0; all divisor numerators and root denominators in rational-circuits.txt.")
        println(io,"Schwartz-Zippel bound NOT ASSERTED: nonuniform dependent ratio sampler. Round-trip is an exact-rational identity test, not a symbolic proof.")
        for j in 1:npts
            pt = complex_emission_point(rng)
            inp,tbl = emission_inputs(pt)
            for (name,want) in exports
                got = emitted_value(mod,tag,name,inp)
                ref = eval_tree(emission_tree(want),tbl)
                residual = got-ref
                compiled &= iszero(residual)
                println(io,"point=$j export=$name residual=$residual")
                if j==1
                    for T in (Float64,BigFloat)
                        typed = map(z -> z isa Complex ? Complex{T}(z) : T(z),inp)
                        value = emitted_value(mod,tag,name,typed)
                        target = Complex{T}(ref)
                        precision &= real(value) isa T && isapprox(value,target;rtol=10000Base.eps(T),atol=10000Base.eps(T))
                    end
                end
            end
        end
    end
    check("compiled independent-complex round-trip",compiled,"$(length(exports)) kernels x $npts exact samples; seed=$seed; nonuniform S recorded; no probability claim";modality=:exact_rational)
    check("Float64 and BigFloat precision inheritance",precision,"compiled kernels, complex inputs, 10000Base.eps(T) absolute/relative comparison")
    # Over Q(i), replacing a complex amplitude by (a+i*b)/|S| is linear,
    # so it cannot increase the total degree already recorded by the circuit.
    # The rational conic and hydrostatic substitutions were counted there.
    # False acceptance implies every numerator vanishes. A zero denominator
    # instead throws and fails; hence no conditional/rejection correction is needed.
    grid_size = 2^24
    grid_seed = seed+100
    grid_rng = MersenneTwister(grid_seed)
    grid_ok = true
    bound = sum(min(RatB(1),RatB(d[1],grid_size))^npts for d in values(degrees))
    open(joinpath(accounting_dir,"$(tag)-uniform-samples.txt"),"w") do io
        println(io,"S=1:$grid_size; independent uniform integer coordinates; N=$npts; seed=$grid_seed.")
        println(io,"Scale ordinary parameters and each complex real/imaginary part by 1/|S|; m unscaled; t=s/(2|S|).")
        println(io,"v,W,p_r,omega derived as in rational-circuits.txt; no rejected points; any zero denominator fails the run.")
        println(io,"Hypotheses: characteristic-zero Q(i); fixed nonzero cleared numerator on the rational constraint surface; IID uniform grid draws (PRNG sampling model).")
        println(io,"Per-component bound <= (d/|S|)^N; union bound over $(length(degrees)) complex numerators = $bound.")
        println(io,"This is a probabilistic identity guarantee for the round-trip, not a symbolic proof or a physical error bound.")
        for j in 1:npts
            inp,tbl = emission_inputs(uniform_emission_point(grid_rng,grid_size))
            for (name,want) in exports
                residual = emitted_value(mod,tag,name,inp)-eval_tree(emission_tree(want),tbl)
                grid_ok &= iszero(residual)
                println(io,"point=$j export=$name residual=$residual")
            end
        end
    end
    check("uniform-grid rational round-trip",grid_ok,
        "N=$npts seed=$grid_seed |S|=$grid_size; d<=$(maximum(first.(values(degrees)))); union false-acceptance bound=$bound";
        modality=:exact_rational)
    code = read(path,String)
    check("runtime dependency isolation",!occursin(r"(?m)^\s*(using |import |include\()",code),"generated functions use only Base arithmetic and input properties")
    C
end
