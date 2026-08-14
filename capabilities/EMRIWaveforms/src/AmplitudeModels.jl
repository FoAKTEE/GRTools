"""
    EMRIWaveforms.AmplitudeModels

Teukolsky-amplitude surrogates (design §6.8: FEW `amplitude/{romannet,ampinterp2d}.jl`; KerrBackwards
`romanamplitude.jl`, `interp2damplitude.jl`). Two pieces:

- `RomanAmplitude` — the ROMAN (reduced-order-model artificial-neuron) dense network. The pure GEMM
  forward pass lives in `Kernels.roman_forward` (already characterized); the HDF5 weight loading is
  `Data.read_roman_nn` (S4-b: reader bodies live in `DataHDF5Ext`); this module owns only the
  `(p,e)→[y,e]` feature map + the conjugate-partner / renormalization conventions.
  ONE struct carries both the `AmplitudeVectorNorm` renorm spline and the `ModeIndexTable`; `get_amplitudes`
  takes a `backend::Symbol` selector: `:few` (default) applies the amplitude-vector renormalization,
  `:kb` returns the raw forward pass and supports the `ModeIndexTable` `(l,m,n)` partner overload.

- `AmpInterpSchwarzEcc` — the consolidated Schwarzschild-eccentric bicubic amplitude interpolant
  (FEW `AmpInterpSchwarzEcc` + KB `Interp2DAmplitude` merged into ONE `BicubicNotAKnot` surface bank
  over the `Data.read_ampgrid` grid; AHA — the two legacy classes are the same `(w,u)=(e, y)` decision).

`schwarzecc_amp_modes` is re-exported from `Data` (the closed-form mode discovery).

- `AmpInterpKerrEqEcc` — the Kerr eccentric-equatorial amplitude interpolant over the 4.74GB `ZNAmps`
  grid (`few.amplitude.ampinterp2d.AmpInterpKerrEccEq`). The grid stores tensor-cubic B-spline
  COEFFICIENTS over clamped knots at 33 spin `z`-knots (regions A & B), read by
  [`KerrAmpSpline.read_kerrampgrid`] and evaluated by [`KerrAmpSpline.bisplev_clamped`]; `get_amplitudes`
  maps the apex via `Mappings.kerrecceq_forward_map` (kind `amplitude`), selects the region by mask, and
  linearly interpolates in `z`. Grounded by the `a → 0` agreement with the already-green
  `AmpInterpSchwarzEcc` (same physical Teukolsky `A_{lmn}`). `kerr_amp_native` stays a documented
  `[Perturbation edge]` (native eccentric-Kerr `A_{lmn}` source integration is not a grounded primitive).
"""
module AmplitudeModels

using NumericsCore: BicubicNotAKnot, blas_single_thread, ModeIndexTable, mode_index_table, p_to_y
import Kernels
import Data

include("KerrAmpSpline.jl")
using .KerrAmpSpline
import .KerrAmpSpline: KerrAmpGrid

export RomanAmplitude, AmpInterpSchwarzEcc, get_amplitudes, schwarzecc_amp_modes,
       AmpInterpKerrEqEcc, kerr_amp_native

# Schwarzschild-eccentric amplitude-grid u-axis y(p,e) = log(p − 2.1 − 2e)
# (few.utils.mappings.schwarzecc.schwarzecc_p_to_y): Phase7-D dedup — canonical `NumericsCore.p_to_y`.

"`schwarzecc_amp_modes(; lmax=10, nmax=30)` — re-export of `Data.schwarzecc_amp_modes`."
schwarzecc_amp_modes(; kwargs...) = Data.schwarzecc_amp_modes(; kwargs...)

# --- shared NN-weight loader (SchwarzschildEccentricInput.hdf5) ---------------------------
# The HDF5 reader BODY moved to Data's weakdep-ext as `Data.read_roman_nn` (S4-b, design §7.2/§10) —
# EMRIWaveforms' direct HDF5 dep is the ext trigger. It returns the NamedTuple
# (weights, biases, reduced_basis, transform_factor_inv, break_index, num_teuk_modes), positionally
# identical to the legacy in-src `_load_roman_nn` 6-tuple; weights are (n_out, n_in) and
# reduced_basis (num_teuk_modes, break_index) — exactly the layout Kernels.roman_forward consumes.

# ==========================================================================================
# Unified RomanAmplitude — ONE struct + get_amplitudes(...; backend=:few|:kb).
# Fuses the FEW renorm fork (romannet.jl) and the KB no-renorm/ModeIndexTable fork
# (romanamplitude.jl): the struct carries BOTH the AmplitudeVectorNorm renorm spline AND the
# ModeIndexTable; the backend selector picks the renorm (:few) or raw / partner-overload (:kb) path.
# ==========================================================================================

_default_nn_file() = Data.datapath("SchwarzschildEccentricInput.hdf5")
_default_norm_file() = Data.datapath("AmplitudeVectorNorm.dat")

"""
    RomanAmplitude(; lmax=10, nmax=30, nn_file=…, norm_file=…)

ROMAN Teukolsky amplitude network (`few.amplitude.romannet.RomanAmplitude` ∪ KerrBackwards
`romanamplitude.jl`, one canonical form). Loads the dense weights/biases, complex `reduced_basis`, the
`AmplitudeVectorNorm` renormalization grid (`amp_norm` spline), and the `(l,m,n)` `ModeIndexTable`. The
forward pass delegates to `Kernels.roman_forward`; [`get_amplitudes`](@ref) selects the renormalizing
(`backend=:few`) or raw (`backend=:kb`) convention. The network basis is fixed at `(lmax,nmax)=(10,30)`;
smaller supported ranges normalize that full basis and return the configured mode-table slice.
"""
struct RomanAmplitude
    weights::Vector{Matrix{Float64}}
    biases::Vector{Vector{Float64}}
    reduced_basis::Matrix{ComplexF64}
    transform_factor_inv::Float64
    break_index::Int
    num_teuk_modes::Int
    num_m0::Int
    network_num_m0::Int
    network_indices::Vector{Int}
    amp_norm::BicubicNotAKnot{Float64}
    table::ModeIndexTable
end

function RomanAmplitude(; lmax::Int = 10, nmax::Int = 30,
                        nn_file::AbstractString = _default_nn_file(),
                        norm_file::AbstractString = _default_norm_file())
    2 <= lmax <= 10 || throw(ArgumentError("lmax must lie in 2:10 for the fixed ROMAN basis, got $lmax"))
    0 <= nmax <= 30 || throw(ArgumentError("nmax must lie in 0:30 for the fixed ROMAN basis, got $nmax"))
    weights, biases, reduced_basis, tfi, bi, ntm = Data.read_roman_nn(nn_file)
    # renormalization spline: RectBivariateSpline(unique_y, unique_e, norm.reshape(num_e,num_y).T)
    raw = Data.read_dat(norm_file)
    y_in = Float64.(raw[:, 1]); e_in = Float64.(raw[:, 2]); norm = Float64.(raw[:, 3])
    uy = sort(unique(y_in)); ue = sort(unique(e_in))
    num_y = length(uy); num_e = length(ue)
    vals = Matrix{Float64}(undef, num_y, num_e)
    for iy in 1:num_y, ie in 1:num_e
        vals[iy, ie] = norm[(ie - 1) * num_y + iy]
    end
    amp_norm = BicubicNotAKnot(uy, ue, vals)
    network_table = mode_index_table(; lmax = 10, nmax = 30)
    ntm == network_table.num_modes ||
        throw(DimensionMismatch("ROMAN network has $ntm output modes; expected $(network_table.num_modes) for lmax=10,nmax=30"))
    table = mode_index_table(; lmax = lmax, nmax = nmax)
    network_indices = Int[network_table.index_map[(table.l_arr[i], table.m_arr[i], table.n_arr[i])]
                          for i in 1:table.num_modes]
    return RomanAmplitude(weights, biases, reduced_basis, tfi, bi, table.num_modes, table.num_m0,
                          network_table.num_m0, network_indices, amp_norm, table)
end

"""
    get_amplitudes(roman::RomanAmplitude, p, e; backend=:few, renormalize=true) -> Vector{ComplexF64}

All configured `num_teuk_modes` Teukolsky mode amplitudes at `(p,e)` (`m ≥ 0`, `m=0`-first). The network forward
pass on `[y,e]` (`y = log(p−2.1−2e)`) is `Kernels.roman_forward`. `backend=:few` (default) rescales the
full fixed-basis vector so its conjugate-folded L2 norm matches the `amp_norm` spline (the FEW
`RomanAmplitude` renorm, skippable with `renormalize=false`), then slices a configured smaller table;
the slice's norm is therefore not itself `amp_norm`. `backend=:kb` returns the raw configured slice
(the KB convention, whose renorm lives in the top-level generator). `renormalize` is ignored for `:kb`.
"""
function get_amplitudes(roman::RomanAmplitude, p::Real, e::Real;
                        backend::Symbol = :few, renormalize::Bool = true)
    backend in (:few, :kb) ||
        throw(ArgumentError("backend must be :few or :kb, got $backend"))
    y = p_to_y(p, e)
    full_teuk = Kernels.roman_forward(roman.weights, roman.biases, roman.reduced_basis,
                                      roman.transform_factor_inv, Float64[y, e])
    if backend === :few && renormalize
        amp_norm = roman.amp_norm(y, e)
        amp_for_norm = sqrt(sum(abs2, full_teuk) +
                            sum(abs2, @view full_teuk[(roman.network_num_m0 + 1):end]))
        full_teuk .*= amp_norm / amp_for_norm
    end
    return length(roman.network_indices) == length(full_teuk) ? full_teuk :
           full_teuk[roman.network_indices]
end

"""
    get_amplitudes(roman::RomanAmplitude, p, e, modes; backend=:kb) -> Vector{ComplexF64}

Amplitudes for the requested `(l,m,n)` `modes` at `(p,e)`; `m < 0` uses the `m > 0` partner conjugated
(`special_index_map` + `conj`, the KB `RomanAmplitude.get_amplitudes` rule). The raw forward pass is
used (`backend=:kb`, the only supported backend for the mode-subset overload).
"""
function get_amplitudes(roman::RomanAmplitude, p::Real, e::Real, modes; backend::Symbol = :kb)
    backend === :kb ||
        throw(ArgumentError("the (l,m,n) mode-subset overload requires backend=:kb, got $backend"))
    teuk = get_amplitudes(roman, p, e; backend = :kb)
    out = Vector{ComplexF64}(undef, length(modes))
    for (i, mode) in enumerate(modes)
        l, m, n = Int.(Tuple(mode))
        j = roman.table.special_index_map[(l, m, n)]
        out[i] = m < 0 ? conj(teuk[j]) : teuk[j]
    end
    return out
end

# ==========================================================================================
# Consolidated AmpInterpSchwarzEcc — BicubicNotAKnot over the Data.read_ampgrid surfaces.
# ==========================================================================================

"""
    AmpInterpSchwarzEcc(; modes=schwarzecc_amp_modes(), file=…)

Schwarzschild-eccentric Teukolsky amplitude interpolant (`few.amplitude.ampinterp2d.AmpInterpSchwarzEcc`
∪ KerrBackwards `Interp2DAmplitude`, one canonical form). Loads the grid + per-mode complex surfaces via
`Data.read_ampgrid` and builds one [`BicubicNotAKnot`] real/imag surface per `(l,m,n)` over the
interpolation coordinates `(w,u) = (e, log(p−2.1−2e))`. Evaluate with [`get_amplitudes`](@ref).
"""
struct AmpInterpSchwarzEcc
    modes::Vector{NTuple{3,Int}}
    index::Dict{NTuple{3,Int},Int}
    re::Vector{BicubicNotAKnot{Float64}}
    im::Vector{BicubicNotAKnot{Float64}}
    unique_w::Vector{Float64}
    unique_u::Vector{Float64}
end

function AmpInterpSchwarzEcc(; modes = schwarzecc_amp_modes(),
                            file::AbstractString = Data.datapath("Teuk_amps_a0.0_lmax_10_nmax_30_new.h5"))
    mode_list = NTuple{3,Int}[Tuple(Int.(m)) for m in modes]
    schema, reM, imM = Data.read_ampgrid(file; modes = mode_list)
    re = [BicubicNotAKnot(schema.w, schema.u, reM[k]) for k in eachindex(mode_list)]
    im = [BicubicNotAKnot(schema.w, schema.u, imM[k]) for k in eachindex(mode_list)]
    index = Dict(m => i for (i, m) in enumerate(mode_list))
    return AmpInterpSchwarzEcc(mode_list, index, re, im, collect(schema.w), collect(schema.u))
end

"""
    get_amplitudes(amp::AmpInterpSchwarzEcc, p, e) -> Vector{ComplexF64}

All stored mode amplitudes at `(p,e)` (Schwarzschild, `a=0`, `xI=1`), at `(w,u) = (e, log(p−2.1−2e))`.
"""
function get_amplitudes(amp::AmpInterpSchwarzEcc, p::Real, e::Real)
    w = e; u = p_to_y(p, e)
    return ComplexF64[amp.re[k](w, u) + im * amp.im[k](w, u) for k in eachindex(amp.modes)]
end

"""
    get_amplitudes(amp::AmpInterpSchwarzEcc, p, e, modes) -> Vector{ComplexF64}

The requested `(l,m,n)` mode amplitudes at `(p,e)`.
"""
function get_amplitudes(amp::AmpInterpSchwarzEcc, p::Real, e::Real, modes)
    w = e; u = p_to_y(p, e)
    out = Vector{ComplexF64}(undef, length(modes))
    for (i, mode) in enumerate(modes)
        k = amp.index[Tuple(Int.(mode))]
        out[i] = amp.re[k](w, u) + im * amp.im[k](w, u)
    end
    return out
end

# =====================================================================================
# Kerr eccentric-equatorial amplitude interpolant (FEW ampinterp2d.AmpInterpKerrEccEq).
# Per-z-knot tensor-cubic B-spline banks over regions A & B (ZNAmps grid, B-spline COEFFICIENTS read by
# KerrAmpSpline.read_kerrampgrid), region-mask select, + LINEAR z(spin) interpolation, driven by
# Mappings.kerrecceq_forward_map(kind="amplitude").
# =====================================================================================

using .KerrAmpSpline: read_kerrampgrid, bisplev_clamped, _flat_to_cmat
using ..Mappings: kerrecceq_forward_map

# dominant default mode family (the (l,l,0) circular-equatorial line); callers pass their own `modes`.
_kerr_default_modes() = NTuple{3,Int}[(l, l, 0) for l in 2:10]

"""
    AmpInterpKerrEqEcc(; modes=_kerr_default_modes(), file=…, threaded=true)

Kerr eccentric-equatorial Teukolsky-amplitude interpolant over the FEW `ZNAmps_l10_m10_n55_DS2Outer.h5`
grid (`few.amplitude.ampinterp2d.AmpInterpKerrEccEq`). The grid stores tensor-cubic B-spline
COEFFICIENTS (regions A & B) at 33 spin `z`-knots; for each requested `(l,m,n)` mode this builds a
per-`z`-knot real/imag coefficient bank for both regions. Evaluate with [`get_amplitudes`](@ref): the
apex `(a,p,e,xI)` is mapped to `(u,w,z)` + a region mask by [`kerrecceq_forward_map`](@ref) (kind
`amplitude`), the bracketing `z`-knot bicubic slices are evaluated by [`bisplev_clamped`], and the
result is linearly interpolated in `z`. The per-mode bank build is threaded under
`NumericsCore.blas_single_thread` (parallel == serial, bit-identical).
"""
struct AmpInterpKerrEqEcc
    modes::Vector{NTuple{3,Int}}
    index::Dict{NTuple{3,Int},Int}
    z_knots::Vector{Float64}
    twA::Vector{Float64}; tuA::Vector{Float64}
    twB::Vector{Float64}; tuB::Vector{Float64}
    reA::Vector{Array{Float64,3}}; imA::Vector{Array{Float64,3}}
    reB::Vector{Array{Float64,3}}; imB::Vector{Array{Float64,3}}
end

function AmpInterpKerrEqEcc(; modes = _kerr_default_modes(),
                           file::AbstractString = Data.datapath("ZNAmps_l10_m10_n55_DS2Outer.h5"),
                           threaded::Bool = true)
    mode_list = NTuple{3,Int}[Tuple(Int.(m)) for m in modes]
    g = read_kerrampgrid(file; modes = mode_list)
    nm = length(mode_list); nz = length(g.z_knots)
    nwcA = length(g.twA) - 4; nucA = length(g.tuA) - 4
    nwcB = length(g.twB) - 4; nucB = length(g.tuB) - 4
    reA = [Array{Float64,3}(undef, nwcA, nucA, nz) for _ in 1:nm]
    imA = [Array{Float64,3}(undef, nwcA, nucA, nz) for _ in 1:nm]
    reB = [Array{Float64,3}(undef, nwcB, nucB, nz) for _ in 1:nm]
    imB = [Array{Float64,3}(undef, nwcB, nucB, nz) for _ in 1:nm]
    build_mode! = function (m)
        @inbounds for z in 1:nz
            reA[m][:, :, z] = _flat_to_cmat(view(g.coeffsA, :, 1, m, z), nwcA, nucA)
            imA[m][:, :, z] = _flat_to_cmat(view(g.coeffsA, :, 2, m, z), nwcA, nucA)
            reB[m][:, :, z] = _flat_to_cmat(view(g.coeffsB, :, 1, m, z), nwcB, nucB)
            imB[m][:, :, z] = _flat_to_cmat(view(g.coeffsB, :, 2, m, z), nwcB, nucB)
        end
    end
    if threaded
        blas_single_thread() do
            Threads.@threads for m in 1:nm
                build_mode!(m)
            end
        end
    else
        for m in 1:nm
            build_mode!(m)
        end
    end
    index = Dict{NTuple{3,Int},Int}(m => i for (i, m) in enumerate(mode_list))
    return AmpInterpKerrEqEcc(mode_list, index, g.z_knots,
                              g.twA, g.tuA, g.twB, g.tuB, reA, imA, reB, imB)
end

@inline _bisplev_complex(tw, tu, reB, imB, zi, w, u) =
    complex(bisplev_clamped(tw, tu, view(reB, :, :, zi), w, u),
            bisplev_clamped(tw, tu, view(imB, :, :, zi), w, u))

# Core scalar evaluation over a set of stored-bank column indices `ks` (1-based into amp.modes).
function _kerr_eval(amp::AmpInterpKerrEqEcc, a, p, e, xI, ks)
    u, w, y, z, near = kerrecceq_forward_map(a, p, e, xI; return_mask = true, kind = "amplitude")
    (0 ≤ u ≤ 1 && 0 ≤ w ≤ 1 && 0 ≤ z ≤ 1) ||
        error("AmpInterpKerrEqEcc accessed out-of-bounds (u=$u, w=$w, z=$z).")
    tw, tu, reBank, imBank = near ? (amp.twA, amp.tuA, amp.reA, amp.imA) :
                                    (amp.twB, amp.tuB, amp.reB, amp.imB)
    zk = amp.z_knots
    out = Vector{ComplexF64}(undef, length(ks))
    exact = findfirst(==(z), zk)
    if exact !== nothing
        @inbounds for (i, k) in enumerate(ks)
            out[i] = _bisplev_complex(tw, tu, reBank[k], imBank[k], exact, w, u)
        end
    else
        ia = findfirst(>(z), zk)
        (ia !== nothing && ia > 1) ||
            error("AmpInterpKerrEqEcc: z=$z outside the z-knot range [$(zk[1]), $(zk[end])].")
        ib = ia - 1
        za = zk[ia]; zb = zk[ib]; s = (z - zb) / (za - zb)
        @inbounds for (i, k) in enumerate(ks)
            Aa = _bisplev_complex(tw, tu, reBank[k], imBank[k], ia, w, u)
            Ab = _bisplev_complex(tw, tu, reBank[k], imBank[k], ib, w, u)
            out[i] = (Aa - Ab) * s + Ab
        end
    end
    return out
end

"""
    get_amplitudes(amp::AmpInterpKerrEqEcc, a, p, e, xI=1) -> Vector{ComplexF64}

All stored mode amplitudes at apex `(a, p, e, xI)` (equatorial `|xI| = 1`). At `a → 0` these reproduce
the Schwarzschild [`AmpInterpSchwarzEcc`] amplitudes (same physical Teukolsky `A_{lmn}`).
"""
get_amplitudes(amp::AmpInterpKerrEqEcc, a::Real, p::Real, e::Real, xI::Real = 1) =
    _kerr_eval(amp, a, p, e, xI, 1:length(amp.modes))

"""
    get_amplitudes(amp::AmpInterpKerrEqEcc, a, p, e, xI, modes) -> Vector{ComplexF64}

The requested `(l,m,n)` mode amplitudes at apex `(a, p, e, xI)`.
"""
get_amplitudes(amp::AmpInterpKerrEqEcc, a::Real, p::Real, e::Real, xI::Real, modes) =
    _kerr_eval(amp, a, p, e, xI, [amp.index[Tuple(Int.(m))] for m in modes])

# --- native on-the-fly compute fallback — documented [Perturbation edge] -------------------
const _KERR_NATIVE_EDGE = "[Perturbation edge] kerr_amp_native: the native eccentric-Kerr Teukolsky \
A_lmn source-integration pipeline is not yet a grounded Perturbation primitive (Perturbation.\
AmplitudesFluxes exposes the MST radial asymptotic amplitudes Btrans/Binc, NOT the waveform A_lmn). \
Use AmpInterpKerrEqEcc (the ZNAmps grid interpolant) instead."

"""
    kerr_amp_native(a, p, e, modes; kwargs...)

Native on-the-fly Kerr-ecc-eq amplitude compute fallback — a documented `[Perturbation edge]`: the
eccentric-Kerr `A_{lmn}` source integration is not yet a grounded Perturbation primitive. The grounded
path is the [`AmpInterpKerrEqEcc`] grid interpolant; this fallback errors with the documented note.
"""
kerr_amp_native(args...; kwargs...) = error(_KERR_NATIVE_EDGE)

end # module AmplitudeModels
