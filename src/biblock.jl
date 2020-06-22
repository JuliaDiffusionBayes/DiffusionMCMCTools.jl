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
    draw_proposal_path!(bb::BiBlock)

Sample a proposal path, compute log-likelihood along the way. Assumes
`bb.b.XX[1].x[1]` is a starting point. Uses preconditioned Crank-Nicolson scheme
with memory parameter set as `bb.ρ`.
"""
function draw_proposal_path!(bb::BiBlock) end

function draw_proposal_path!(bb::BiBlock{false})
    success, bb.b°.ll = _draw_proposal_path!(bb)
    success || return false # anything?
    success, ll°_last = _draw_proposal_path_last_segment!(bb, y1)
    bb.b°.ll += ll°_last
    return success
end

function draw_proposal_path!(bb::BiBlock{true})
    success, bb.b°.ll = _draw_proposal_path!(bb)
    return success
end

function _draw_proposal_path!(bb::BiBlock)
    rand!(
        bb.b.PP, bb.b°.XX, bb.b°.WW, bb.b.WW, bb.ρ, _LL, bb.b.XX[1].x[1];
        Wnr=bb.b°.Wnr
    )
end

function _draw_proposal_path_last_segment!(bb::BiBlock, y1)
    success, ll° = rand!(
        bb.b.P_last[1], bb.b°.XX[end], bb.b°.WW[end], bb.b.WW[end], bb.ρ, _LL,
        y1; Wnr=bb.b°.Wnr
    )
end

"""
    accept_reject_proposal_path!(bb::BiBlock, mcmciter)

Accept/reject decision of the Metropolis-Hastings algorithm for the step of path
imputation.
"""
function accept_reject_proposal_path!(bb::BiBlock, mcmciter)
    accepted = rand(Exponential(1.0)) > -(bb.b°.ll - bb.b.ll)
    accepted && swap_paths!(bb)
    set_accepted!(bb, mcmciter, accepted)
    save_ll!(bb.b, mcmciter)
    save_ll!(bb.b°, mcmciter)
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

"""
    ll_of_accepted(bb::BiBlock, i)

Return the log-likelihood of the path that was accepted at the `i`th iteration.
"""
function ll_of_accepted(bb::BiBlock, i)
    return (bb.accpt_history[i] ? bb.b°.ll_history[i] : bb.b.ll_history[i])
end

"""
    GP.set_obs!(bb::BiBlock)

Freeze an artificial observation at the terminal point of the block. For a
terminal block nothing is done.
"""
function GP.set_obs!(bb::BiBlock) end

function GP.set_obs!(bb::BiBlock{true})
    set_obs!(bb.b.P_last, bb.b.XX[end].x[end])
    set_obs!(bb.b°.P_last, bb.b.XX[end].x[end])
end

GP.set_obs!(bb::BiBlock{false}) = nothing


"""
    GP.recompute_guiding_term!(bb::BiBlock)

Recompute the guiding terms of both the proposal and the accepted laws.
"""
function GP.recompute_guiding_term!(bb::BiBlock)
    recompute_guiding_term!(bb.b)
    recompute_guiding_term!(bb.b°)
end

"""
    find_W_for_X!(bb::BiBlock)

Find the Wiener process `bb.b.WW` that reconstructs path `bb.b.XX` under the
accepted law `bb.b.PP` (possibly including `bb.b.P_last`).
"""
find_W_for_X!(bb::BiBlock) = find_W_for_X!(bb.b)
