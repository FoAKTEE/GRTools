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
    doctest = false,
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
        ],
        "Migration" => "migration.md",
    ],
)

# No deploydocs: this is a local site (build output in docs/build/).
