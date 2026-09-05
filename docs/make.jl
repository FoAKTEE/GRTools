# GrTools Documenter site (§12 of GRTOOLS_DESIGN_CAPABILITY.md).
# One site, page-per-capability. Local build only (no deploydocs).
#
# Build:  julia --project=docs docs/make.jl   (from repo root)

using Documenter

# Core (tiers 0–1)
using NumericsCore
using BigSpecial
using Wigner
using Data
using Kernels
using TestKit
using EllipticCore

# Capabilities (tiers 2–6)
using Spacetime
using Harmonics
using Perturbation
using MetricReconstruction
using SelfForce
using PostNewtonian
using ModifiedGravity
using EMRIWaveforms
using Hybrid
using AccretionDisc
using ScalarSelfForce

# Doctests are evaluated in `Main`, so a `jldoctest` block inside a package docstring sees
# nothing of its own module unless a `DocTestSetup` is registered for it. Register `using <Pkg>`
# once for every module the site documents — inert for a package with no doctest yet, and the
# reason the first one that landed (`ScalarSelfForce.select_hold_sweep`) can run at all.
for m in (NumericsCore, BigSpecial, Wigner, Data, Kernels, TestKit, EllipticCore,
          Spacetime, Harmonics, Perturbation, MetricReconstruction,
          SelfForce, ScalarSelfForce, PostNewtonian, ModifiedGravity,
          EMRIWaveforms, Hybrid, AccretionDisc)
    DocMeta.setdocmeta!(m, :DocTestSetup, :(using $(nameof(m))); recursive = true)
end

makedocs(
    sitename = "GrTools",
    # Local/standalone site: no remote-repository detection (restore a Remotes.GitHub
    # entry here if the site is ever deployed with edit links).
    remotes = nothing,
    modules = [
        NumericsCore, BigSpecial, Wigner, Data, Kernels, TestKit, EllipticCore,
        Spacetime, Harmonics, Perturbation, MetricReconstruction,
        SelfForce, ScalarSelfForce, PostNewtonian, ModifiedGravity,
        EMRIWaveforms, Hybrid, AccretionDisc,
    ],
    format = Documenter.HTML(
        prettyurls = false,
        mathengine = Documenter.MathJax3(),
        # Large API pages (Perturbation ~230 KiB) exceed the 200 KiB default;
        # this is a local-only site, so raise the hard cap (warn stays at 100 KiB).
        size_threshold = 512 * 2^10,
    ),
    # No `jldoctest` block exists in the suite yet (0 in docs/src, 0 in any docstring), so
    # this runs zero doctests today and arms the gate for the first one that lands.
    doctest = true,
    warnonly = [:missing_docs, :cross_references, :docs_block, :autodocs_block],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Capabilities" => [
            "capabilities/spacetime.md",
            "capabilities/harmonics.md",
            "capabilities/perturbation.md",
            "capabilities/metricreconstruction.md",
            "capabilities/selfforce.md",
            "capabilities/scalarselfforce.md",
            "capabilities/postnewtonian.md",
            "capabilities/modifiedgravity.md",
            "capabilities/emriwaveforms.md",
            "capabilities/hybrid.md",
            "capabilities/accretiondisc.md",
        ],
        "Core" => [
            "core/numericscore.md",
            "core/bigspecial.md",
            "core/wigner.md",
            "core/data.md",
            "core/kernels.md",
            "core/testkit.md",
            "core/ellipticcore.md",
        ],
        "Topics" => [
            "topics/conventions.md",
            "topics/precision.md",
            "topics/gpu-parallel.md",
            "topics/data-artifacts.md",
            "topics/fixtures.md",
            "topics/reproducibility.md",
            "topics/reproducing-2606.04998.md",
        ],
        "Migration" => "migration.md",
    ],
)

# No deploydocs: this is a local site (build output in docs/build/).
