"""Runtime-facing FR certificates; deliberately separate from the R3 environment."""
module FRDerivation

using AccretionDisc
using Symbolics
using Random

# Include the existing alphabet and verifier, without importing the R3 package or
# its relativistic derivation stages. No certificate machinery is copied here.
module R3Support
using Symbolics
include("../../derivation/src/symbols.jl")
include("../../derivation/src/certify.jl")
end
using .R3Support: Certificate, RatB, exact_eval, certify_symbolic, certify_bool
const Modes = AccretionDisc.Modes

include("runtime_identities.jl")
include("kepler.jl")
include("background.jl")
include("fr24_background.jl")
include("fr24_dimensions.jl")
include("normalisation.jl")
include("fourier.jl")
include("second_order_averages.jl")
include("linearised.jl")
include("lindblad.jl")
include("phi_factorisation.jl")

end
