# Setup-value gate only: convergence-result verification belongs to the later figure nodes.
include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
include(joinpath(@__DIR__, "_fr22_setup.jl"))

"Compare the complete setup object with the independently extracted JSON and its tex anchors."
function check_fr22_setup(setup::FR22Setup, spec, src)
    entries = vcat(spec["entries"], spec["figure_printed"])
    ids = [Symbol(ent["id"]) for ent in entries]
    complete = length(unique(ids)) == length(ids) && Set(ids) == Set(keys(setup.values))
    println(complete ? "[PASS] complete printed-value coverage" : "[FAIL] printed-value coverage differs")
    passed = complete
    for ent in entries
        key = Symbol(ent["id"])
        value = get(setup.values, key, nothing)
        comparable = value isa Tuple ? collect(value) : value
        ok = comparable == ent["value"]
        source = ent["source"]
        if haskey(source, "file")
            lines = readlines(joinpath(src, source["file"]))
            line = source["line"]
            ok &= 1 <= line <= length(lines) && occursin(source["tex"], lines[line])
        else
            ok &= all(isfile(joinpath(src, f)) for f in source["files"])
        end
        category = ent["category"]
        scope = category == "result_statement" ? "printed target only; figure verification pending" : "setup"
        println("[", ok ? "PASS" : "FAIL", "] ", key, " paper=", ent["value"], " ours=", comparable, " (", scope, ")")
        passed &= ok
    end
    v = setup.values
    consistent = v.e_list == v.e_list_results &&
        (v.m_max, v.lm_max) == v.conv_adopted &&
        v.m_max == v.conv_lm_fixed_m[1] && v.lm_max == v.conv_m_fixed_lm[1] &&
        isapprox(v.zz22_mass[1]/v.h_p^v.M_th_def, v.zz22_mass[2]; rtol=8eps(typeof(v.h_p))) &&
        v.snapshots_caption == v.snapshots_text[1:3] && v.column_titles == v.snapshots_text
    println(consistent ? "[PASS] setup cross-identities" : "[FAIL] setup cross-identities")
    println("[OPEN] epsilon = 0.3H: choice = ", setup.softening,
            "; caption omits phase 0.75; wavelet interval T is interpreted as T_p.")
    println(length(entries), "/", length(entries), " printed values checked; ",
            passed && consistent ? "ALL PASS" : "FAIL")
    passed && consistent
end

function main()
    opts = ReproHarness.parse_args()
    ReproHarness.setup!()
    ReproHarness.require_tier(opts, :fast)
    root = get(ENV, "EMRI_DISK_ROOT", normpath(joinpath(@__DIR__, "..", "..", "..", "emri-disk")))
    specpath = joinpath(root, "results", "emri-eccentric", "paper_EMRI-ECC", "evidence",
                        "fr22-setup-parameters", "setup_parameters.json")
    src = joinpath(root, "ref-paper", "arXiv-2207.14637", "src")
    isfile(specpath) && isdir(src) || ReproHarness.hole("FR22 extracted setup/source missing; set EMRI_DISK_ROOT")
    spec = ReproHarness.json_read(read(specpath, String))
    check_fr22_setup(FR22Setup(), spec, src) ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
