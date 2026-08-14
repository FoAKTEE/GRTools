# Barack-Lousto tensor-harmonic normalization + second-order-Ricci mode-selection rules.
# `a_il` is owned by NumericsCore (re-exported through the SecondOrder front door, not re-owned);
# only `selection_modes` is defined here.

"""
    selection_modes(lmax) -> Vector{NTuple{3,Int}}

Valid target `(i,l,m)` modes of the second-order Ricci source for `l ∈ 0:lmax`, `m ∈ 0:l`,
`i ∈ 1:10`, returned in the driver push-order (`l` outermost, then `m`, then `i`). Transcription of
the mode-selection rule block in `ref-code/SecondOrderRicci/src/Ricci.cc` L86-99:

- parity: drop `i ≤ 7` when `l+m` is odd, drop `i > 7` when `l+m` is even;
- drop `i ∈ {4,5,8,9}` when `l < 1`;
- drop `i ∈ {7,10}` when `l < 2`.

The push-order is stable under raising `lmax` (the `lmax=k` result is a strict prefix of `lmax=k+1`).
"""
function selection_modes(lmax::Integer)
    modes = NTuple{3,Int}[]
    for l3 in 0:lmax, m3 in 0:l3, i3 in 1:10
        ((i3 <= 7 && isodd(l3 + m3)) || (i3 > 7 && iseven(l3 + m3))) && continue
        (i3 in (4, 5, 8, 9) && l3 < 1) && continue
        (i3 in (7, 10) && l3 < 2) && continue
        push!(modes, (i3, l3, m3))
    end
    return modes
end
