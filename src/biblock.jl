#===============================================================================

    The main object defined in this file is: BiBlock.
    It represents a single block and is the simplest composite unit
    that can be used for practical settings of inference and smoothing.
    Additionally, it contains the option to sample with a preconditioned
    Crank-Nicolson step.

===============================================================================#

"""
    mutable struct BiBlock{L,TGP,TGPl,TW,TWn,TX}
        b::Block{L,TGP,TGPl,TW,TWn,TX}
        b°::Block{L,TGP,TGPl,TW,TWn,TX}
        ρ::Float64
        accpt_history::Vector{Bool}
    end

Composite unit that allows for sampling of a single block. It provides two
`Block`s: one proposal `b°`, one accepted `b` that can be used for smoothing or
inference problems. `ρ` is a memory parameter of the preconditioned
Crank-Nicolson scheme and `accpt_history` stores the history of accept/reject
decisions (useful for MCMC).

    function BiBlock(
        sp::SamplingPair,
        range::UnitRange{Int64},
        ρ=0.0,
        last_block=false,
        ll_hist_len=0
    )

Base constructor.
"""
mutable struct BiBlock{L,TGP,TGPl,TW,TWn,TX}
    b::Block{L,TGP,TGPl,TW,TWn,TX}
    b°::Block{L,TGP,TGPl,TW,TWn,TX}
    ρ::Float64
    accpt_history::Vector{Bool}

    function BiBlock(
            sp::SamplingPair{TGP,TGPl,TW,TWn,TX},
            range::UnitRange{Int64},
            ρ=0.0,
            last_block=false,
            ll_hist_len=0
        ) where {TGP,TGPl,TW,TWn,TX}
        b = Block(sp.u, range, last_block, ll_hist_len)
        b° = Block(sp.u°, range, last_block, ll_hist_len)
        new{last_block,TGP,TGPl,TW,TWn,TX}(
            b, b°, ρ,
            Vector{Bool}(undef, ll_hist_len)
        )
    end
end

"""
    set_accepted!(bb::BiBlock, i::Int, v)

Commit the accept/reject decision `v` to acceptance history of BiBlock `b` at
the position `i`.
"""
set_accepted!(bb::BiBlock, i::Int, v) = (bb.accpt_history[i] = v)


"""
    swap_paths!(bb::BiBlock)

Swap `XX` and `WW` containers between proposal-acceptance pair.
"""
function swap_paths!(bb::BiBlock)
    for i in eachindex(bb.b.XX)
        bb.b.XX[i], bb.b°.XX[i] = bb.b°.XX[i], bb.b.XX[i]
        bb.b.WW[i], bb.b°.WW[i] = bb.b°.WW[i], bb.b.WW[i]
    end
end

"""
    swap_XX!(bb::BiBlock)

Swap `XX` containers between proposal-acceptance pair.
"""
function swap_XX!(bb::BiBlock)
    for i in eachindex(bb.b.XX)
        bb.b.XX[i], bb.b°.XX[i] = bb.b°.XX[i], bb.b.XX[i]
    end
end

"""
    swap_WW!(bb::BiBlock)

Swap `WW` containers between proposal-acceptance pair.
"""
function swap_WW!(bb::BiBlock)
    for i in eachindex(bb.b.WW)
        bb.b.WW[i], bb.b°.WW[i] = bb.b°.WW[i], bb.b.WW[i]
    end
end

"""
    swap_PP!(bb::BiBlock)

Swap `PP` containers (including PP_last) between proposal-acceptance pair.
"""
function swap_PP!(bb::BiBlock)
    for i in eachindex(bb.b.PP)
        bb.b.PP[i], bb.b°.PP[i] = bb.b°.PP[i], bb.b.PP[i]
    end
    # on a terminal block i.e. BiBlock{true} this simply makes no difference
    bb.b.P_last[1], bb.b°.P_last[1] = bb.b°.P_last[1], bb.b.P_last[1]
end
