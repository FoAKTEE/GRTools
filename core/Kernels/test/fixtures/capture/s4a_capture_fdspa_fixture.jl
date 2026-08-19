# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/s4a_capture_fdspa_fixture.jl (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = 4d1c62ad8d893c96abeb388c4e58bcc2b09b983af5c1ad0be7ff71c57a91b6ec
# fixture(s)  : fdspa_ref.jl, fdspa_inputs.jl
# oracle      : PRE-HOIST EMRIWaveforms.ModeSum.FEW.fd_interpolated_mode_sum (removed by the S4-a hoist; oracle no longer exists in-tree)
# run         : julia --project=capabilities/EMRIWaveforms test/fixtures/capture/s4a_capture_fdspa_fixture.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const EMRIDISK_ROOT        = get(ENV, "GRTOOLS_EMRIDISK_ROOT",        "/data/haiyangw/claude/emri-disk")
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

# s4a_capture_fdspa_fixture.jl — one-shot PRE-HOIST capture of the FD-SPA mode-sum outputs from the
# committed EMRIWaveforms.ModeSum.FEW.fd_interpolated_mode_sum (the FEW-namespace local kernel at
# EMRIWaveforms/src/ModeSum.jl, BEFORE it is hoisted into Kernels).
#
# Run:  julia --project=grtools/capabilities/EMRIWaveforms scratchpad/s4a_capture_fdspa_fixture.jl
# Writes: grtools/core/Kernels/test/fixtures/fdspa_ref.jl
#
# Inputs come from the SHARED generator file (grtools/core/Kernels/test/fixtures/fdspa_inputs.jl) that
# the Kernels test will also include — bit-identical inputs, no transcription.

using EMRIWaveforms
using NumericsCore
const FEW = EMRIWaveforms.ModeSum.FEW

const ROOT = "$(EMRIDISK_ROOT)"
include(joinpath(FIXTURE_DIR, "fdspa_inputs.jl"))

# ---- run the pre-change kernel ---------------------------------------------------------------------
c1 = fdspa_case1_inputs()
hp1, hc1 = FEW.fd_interpolated_mode_sum(c1.t, c1.teuk, c1.ylms, c1.Phi_phi, c1.Phi_r,
                                        c1.f_phi, c1.f_r, c1.m_arr, c1.n_arr; f_arr = c1.f_arr)
hp1n, hc1n = FEW.fd_interpolated_mode_sum(c1.t, c1.teuk, c1.ylms, c1.Phi_phi, c1.Phi_r,
                                          c1.f_phi, c1.f_r, c1.m_arr, c1.n_arr;
                                          f_arr = c1.f_arr, include_minus_m = false)
sep1 = FEW.fd_interpolated_mode_sum(c1.t, c1.teuk, c1.ylms, c1.Phi_phi, c1.Phi_r,
                                    c1.f_phi, c1.f_r, c1.m_arr, c1.n_arr;
                                    f_arr = c1.f_arr, separate_modes = true)
c2 = fdspa_case2_inputs()
hp2, hc2 = FEW.fd_interpolated_mode_sum(c2.t, c2.teuk, c2.ylms, c2.Phi_phi, c2.Phi_r,
                                        c2.f_phi, c2.f_r, c2.m_arr, c2.n_arr; f_arr = c2.f_arr)

# ---- coverage audit: case 2 must exercise all three ±m target-bin accumulation orders ---------------
# (replicates the legacy kernel's per-(mode,segment) visited-bin-range arithmetic on the same spline)
function coverage_counts(c)
    t = c.t; nt, K = size(c.teuk); nseg = nt - 1
    freq = collect(float.(c.f_arr)); nf = length(freq)
    zidx = findfirst(==(0.0), freq)
    f_first = freq[1]; df = freq[2] - freq[1]
    Y = Matrix{Float64}(undef, 2K + 4, nt)
    @views Y[1:K, :]    .= permutedims(real.(c.teuk))
    @views Y[K+1:2K, :] .= permutedims(imag.(c.teuk))
    @views Y[2K+1, :]   .= c.f_phi
    @views Y[2K+2, :]   .= c.f_r
    @views Y[2K+3, :]   .= c.Phi_phi
    @views Y[2K+4, :]   .= c.Phi_r
    sp = NumericsCore.CubicSplineInterpolant(collect(Float64, t), Y)
    inter = dirfirst = parfirst = 0
    for mi in 1:K, seg in 1:nseg
        m = c.m_arr[mi]; n = c.n_arr[mi]
        L = t[seg+1] - t[seg]
        f_y  = m * sp.y[2K+1, seg]  + n * sp.y[2K+2, seg]
        f_c1 = m * sp.c1[2K+1, seg] + n * sp.c1[2K+2, seg]
        f_c2 = m * sp.c2[2K+1, seg] + n * sp.c2[2K+2, seg]
        f_c3 = m * sp.c3[2K+1, seg] + n * sp.c3[2K+2, seg]
        fmin, fmax = FEW._fd_seg_freq_range(f_c3, f_c2, f_c1, f_y, L)
        bl = trunc(Int, abs((fmin - f_first) / df)) + 1
        bh = trunc(Int, abs((fmax - f_first) / df))
        m == 0 && continue
        for i in max(1, bl + 1):min(nf, bh + 1)
            j = 2 * zidx - i
            (1 <= j <= nf && bl <= j - 1 <= bh) || continue
            i == j ? (inter += 1) : (i < j ? (dirfirst += 1) : (parfirst += 1))
        end
    end
    return (inter, dirfirst, parfirst)
end
inter, dirf, parf = coverage_counts(c2)
println("case2 accumulation-order coverage: interleave(i==zidx)=$inter  direct-first=$dirf  partner-first=$parf")
@assert inter > 0 && dirf > 0 && parf > 0 "case 2 does not cover all three ±m accumulation orders"
nzlo = count(i -> !iszero(hp2[i]) && !isnan(hp2[i]), 1:161)
nzhi = count(i -> !iszero(hp2[i]) && !isnan(hp2[i]), 161:321)
println("case2 nonzero finite hp bins: below/at zidx=$nzlo  at/above zidx=$nzhi")
@assert nzlo > 0 && nzhi > 0 "case 2 hp must be nonzero on both sides of the zero bin"

# ---- serialize --------------------------------------------------------------------------------------
function entlines(v::AbstractVector)
    lits = String[]; nans = Int[]
    for (i, z) in pairs(v)
        if isnan(z)
            push!(nans, i)
        elseif !iszero(z)
            push!(lits, "    ($(i), $(repr(real(z))), $(repr(imag(z)))),")
        end
    end
    return lits, nans
end
function entlines(M::AbstractMatrix)
    lits = String[]; nans = Tuple{Int,Int}[]
    for I in CartesianIndices(M)
        z = M[I]
        if isnan(z)
            push!(nans, (I[1], I[2]))
        elseif !iszero(z)
            push!(lits, "    ($(I[1]), $(I[2]), $(repr(real(z))), $(repr(imag(z)))),")
        end
    end
    return lits, nans
end

commit = readchomp(`git -C $ROOT rev-parse --short HEAD`)
io = IOBuffer()
println(io, "# fdspa_ref.jl — BAKED pre-hoist FD-SPA mode-sum reference (characterization fixture).")
println(io, "#")
println(io, "# Captured one-shot by scratchpad/s4a_capture_fdspa_fixture.jl from the committed FEW-namespace")
println(io, "# kernel EMRIWaveforms.ModeSum.FEW.fd_interpolated_mode_sum (EMRIWaveforms/src/ModeSum.jl) at")
println(io, "# repo state $commit, BEFORE the kernel was hoisted into Kernels (S4-a). Inputs regenerated")
println(io, "# bit-identically by fdspa_inputs.jl (shared with the capture script — no transcription).")
println(io, "# Entries are (bin, re, im) [vector] / (mode, bin, re, im) [separate-modes matrix] for every")
println(io, "# non-zero finite bin (repr round-trip ⇒ exact Float64); *_NAN lists the SPA-turnover NaN bins;")
println(io, "# ALL remaining bins were exactly zero in the capture.")
for (name, val) in [("FDSPA_C1_HP", hp1), ("FDSPA_C1_HC", hc1),
                    ("FDSPA_C1_NOMM_HP", hp1n), ("FDSPA_C1_NOMM_HC", hc1n),
                    ("FDSPA_C2_HP", hp2), ("FDSPA_C2_HC", hc2)]
    lits, nans = entlines(val)
    println(io, "\nconst $name = Tuple{Int,Float64,Float64}[")
    foreach(l -> println(io, l), lits)
    println(io, "]")
    println(io, "const $(name)_NAN = Int[$(join(nans, ", "))]")
end
lits, nans = entlines(sep1)
println(io, "\nconst FDSPA_C1_SEP = Tuple{Int,Int,Float64,Float64}[")
foreach(l -> println(io, l), lits)
println(io, "]")
println(io, "const FDSPA_C1_SEP_NAN = Tuple{Int,Int}[$(join(nans, ", "))]")

out = joinpath(FIXTURE_DIR, "fdspa_ref.jl")
write(out, String(take!(io)))
println("wrote $out")
println("counts: C1_HP nz=$(length(entlines(hp1)[1])) nan=$(length(entlines(hp1)[2]))  C2_HP nz=$(length(entlines(hp2)[1])) nan=$(length(entlines(hp2)[2]))  SEP nz=$(length(entlines(sep1)[1]))")
