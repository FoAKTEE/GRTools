"""Read-only FR24 v1 oracle; no mirrored program is loaded or executed."""
module FR24Cube
using SHA

export read_cube, check_cube

"""
Read the 14 non-indirect .out files. Arrays use (q,p,h[,e]) order, last axis
fastest in the text. Source: read_data.ipynb:45-63,82-103 (read as text only).
Values and a conservative serialization uncertainty travel together.
No authors' data are vendored. *_indirect tables postdate v1 and are excluded.
"""
function read_cube(root)
    hashes = Dict{String,String}()
    grids = Dict{String,Any}()
    tables = Dict{String,Any}()
    for family in ("circular","eccentric")
        rel = joinpath(family,"parameters.out")
        text = read(joinpath(root,rel),String)
        hashes[rel] = bytes2hex(sha256(text))
        grid = [parse.(Float64,split(line)) for line in split(strip(text),'\n')]
        length(grid) == (family == "circular" ? 3 : 4) || error("invalid parameter axes")
        all(g -> !isempty(g) && all(isfinite,g) && issorted(g) && allunique(g),grid) ||
            error("parameter axes must be finite, nonempty and strictly increasing")
        grids[family] = grid
        shape = Tuple(length.(grid))
        names = family == "circular" ? ["t_$(part)_m$m" for m in 1:2 for part in ("net","l","c")] :
                ["t_net","t_l","t_c","tau_L_inv","tau_a_inv","tau_e_inv"]
        for name in names
            rel = joinpath(family,name*".out")
            text = read(joinpath(root,rel),String)
            hashes[rel] = bytes2hex(sha256(text))
            tokens = split(text)
            length(tokens) == prod(shape) || error("wrong table length: $rel")
            values = parse.(Float64,tokens)
            all(isfinite,values) || error("nonfinite table entry: $rel")
            uncertainty = map(tokens) do token
                matched = match(r"^[+-]?\d+\.(\d+)[eE]([+-]?\d+)$",token)
                matched === nothing && error("expected scientific decimal: $rel")
                max(5e-8,0.5*10.0^(parse(Int,matched[2])-length(matched[1])))
            end
            # Julia's first axis is fastest; reverse before permuting back to named order.
            axes = reverse(1:length(shape))
            tables[joinpath(family,name)] = (
                value=permutedims(reshape(values,reverse(shape)),axes),
                uncertainty=permutedims(reshape(uncertainty,reverse(shape)),axes))
        end
    end
    (; grids,tables,hashes)
end

# PREREG: linear propagation of serialized precision, plus bounded Float64 slack.
function identity_check(terms...)
    residual = abs.(sum(c .* table.value for (c,table) in terms))
    bound = sum(abs.(c) .* table.uncertainty for (c,table) in terms)
    scale = sum(abs.(c .* table.value) for (c,table) in terms)
    bound .+= 32eps(Float64) .* max.(1,scale)
    worst = argmax(residual)
    Dict("count"=>length(residual), "max_absolute_residual"=>maximum(residual),
         "tolerance_at_max_residual"=>bound[worst], "max_tolerance"=>maximum(bound),
         "failure_count"=>count(residual .> bound),
         "worst_absolute_index"=>collect(Tuple(worst)),
         "worst_ratio_index"=>collect(Tuple(argmax(residual ./ bound))),
         "max_residual_over_tolerance"=>maximum(residual ./ bound),
         "passed"=>all(residual .<= bound))
end

function distribution(values)
    s = sort(values)
    quantile(p) = begin
        x = 1+(length(s)-1)*p
        lo,hi = floor(Int,x),ceil(Int,x)
        s[lo]+(x-lo)*(s[hi]-s[lo])
    end
    Dict("count"=>length(s),"min"=>first(s),"p25"=>quantile(0.25),
         "median"=>quantile(0.5),"p75"=>quantile(0.75),"p95"=>quantile(0.95),
         "max"=>last(s),"rms"=>sqrt(sum(abs2,s)/length(s)))
end

"FR24 eq39,45,48-51; 4_eccentric_torques.tex:105-160, 2_linear_framework.tex:314."
function check_cube(cube)
    tables,grids = cube.tables,cube.grids
    length(cube.hashes)==14 || error("expected 14 v1 files")
    ec = grids["eccentric"]; cc = grids["circular"]
    Tuple(length.(ec)) == (4,4,3,13) || error("expected 624 eccentric points")
    Tuple(length.(cc)) == (4,4,5) || error("expected 80 circular points")
    e = reshape(ec[4],1,1,1,:)
    zero_e = findfirst(iszero,ec[4])
    zero_e === nothing && error("missing circular slice")
    checks = Dict{String,Any}()
    for (family,suffix) in (("eccentric",""),("circular","_m1"),("circular","_m2"))
        checks[family*"_partition"*suffix] = identity_check(
            (1,tables[family*"/t_net"*suffix]),(-1,tables[family*"/t_l"*suffix]),(-1,tables[family*"/t_c"*suffix]))
    end
    tauL, taua, taue = (tables["eccentric/"*s] for s in ("tau_L_inv","tau_a_inv","tau_e_inv"))
    checks["angular_momentum"] = identity_check((1,tauL),(-1 ./ sqrt.(1 .- e.^2),tables["eccentric/t_net"]))
    checks["timescale_relation"] = identity_check((1,tauL),(-1/2,taua),(e.^2 ./ (1 .- e.^2),taue))
    slice(t) = (;value=t.value[:,:,:,zero_e],uncertainty=t.uncertainty[:,:,:,zero_e])
    checks["circular_a"] = identity_check((1,slice(taua)),(-2,slice(tauL)))
    checks["circular_e"] = Dict("count"=>length(slice(taue).value),
        "max_absolute_residual"=>maximum(abs,slice(taue).value),"tolerance"=>0.0,
        "passed"=>all(iszero,slice(taue).value))
    # Separate runs: report the distribution; never apply an algebraic equality tolerance.
    cross = Dict{String,Any}()
    mapping = [[something(findfirst(==(x),cc[k])) for x in ec[k]] for k in 1:3]
    for part in ("net","l","c")
        v = slice(tables["eccentric/t_"*part]).value
        reference = tables["circular/t_"*part*"_m1"].value[mapping...]
        absolute = abs.(v-reference)
        nonzero = .!iszero.(reference)
        relative = absolute[nonzero] ./ abs.(reference[nonzero])
        metric = fill(-Inf,size(v)); metric[nonzero] .= relative
        worst = argmax(metric)
        cross[part] = Dict("relative"=>distribution(relative),"absolute"=>distribution(vec(absolute)),
            "zero_reference_count"=>count(!,nonzero),
            "worst_relative_qph"=>[ec[k][Tuple(worst)[k]] for k in 1:3])
    end
    Dict("checks"=>checks,"cross_run"=>cross,"source_sha256"=>cube.hashes,
         "passed"=>all(c["passed"] for c in values(checks)))
end
end
