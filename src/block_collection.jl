"""
    struct BlockCollection{TGP,TGPl,TW,TWn,TX}
        blocks::Vector{BiBlock{_L,TGP,TGPl,TW,TWn,TX} where _L}
    end

Gathers all blocks relevant for an entire single recording.

    BlockCollection(
        sp::SamplingPair{TGP,TGPl,TW,TWn,TX},
        ranges,
        ρρ=0.0,
        ll_hist_len=0
    ) where {TGP,TGPl,TW,TWn,TX}

Base constructor.
"""
struct BlockCollection{TGP,TGPl,TW,TWn,TX}
    blocks::Vector{BiBlock{_L,TGP,TGPl,TW,TWn,TX} where _L}

    function BlockCollection(
            sp::SamplingPair{TGP,TGPl,TW,TWn,TX},
            ranges,
            ρρ=0.0,
            ll_hist_len=0
        ) where {TGP,TGPl,TW,TWn,TX}
        N = length(ranges)
        ρρ = (typeof(ρρ) <: Number ? fill(ρρ, N) : ρρ)
        blocks = map(i->BiBlock(sp, ranges[i], ρρ[i], i==N, ll_hist_len), 1:N)
        new{TGP,TGPl,TW,TWn,TX}(blocks)
    end
end


#===============================================================================

        IMPUTATION

===============================================================================#
"""
    draw_proposal_path!(bc::BlockCollection)

Sample proposal paths on each block, compute log-likelihoods along the way.
Assumes `bc.blocks[i].b.XX[1].x[1]` are starting points for each block. Uses
the preconditioned Crank-Nicolson scheme.
"""
draw_proposal_path!(bc::BlockCollection) = draw_proposal_path!.(bc.blocks)

#===============================================================================

        ACCEPT/REJECT DECISION

===============================================================================#

"""
    accept_reject_proposal_path!(bb::BlockCollection, mcmciter)

Accept/reject decision of the Metropolis-Hastings algorithm for the step of path
imputation, done separately for each block.
"""
function accept_reject_proposal_path!(bc::BlockCollection, mcmciter)
    for biblock in bc.blocks
        accept_reject_proposal_path!(biblock, mcmciter)
    end
end

#===============================================================================

        SWAPS

===============================================================================#

"""
    swap_paths!(bc::BlockCollection)

For each block in the collection swap `XX` and `WW` containers between
proposal-acceptance pair.
"""
swap_paths!(bc::BlockCollection) = swap_paths!.(bc.blocks)

"""
    swap_XX!(bc::BlockCollection)

For each block in the collection swap `XX` containers between
proposal-acceptance pair.
"""
swap_XX!(bc::BlockCollection) = swap_XX!.(bc.blocks)

"""
    swap_WW!(bc::BlockCollection)

For each block in the collection swap `WW` containers between
proposal-acceptance pair.
"""
swap_WW!(bc::BlockCollection) = swap_WW!.(bc.blocks)

"""
    swap_PP!(bc::BlockCollection)

For each block in the collection swap `PP` containers between
proposal-acceptance pair.
"""
swap_PP!(bc::BlockCollection) = swap_PP!.(bc.blocks)

"""
    swap_ll!(bc::BlockCollection)

For each block in the collection swap `ll` fields between
proposal-acceptance pair.
"""
swap_ll!(bc::BlockCollection) = swap_ll!.(bc.blocks)

#===============================================================================

        UTILITY

===============================================================================#

"""
    loglikhd!(bc::BlockCollection)

For each `BiBlock` in a collection compute the log-likelihood for the accepted
block, evaluated at a sampled path and store the result in internal fields `ll`.
"""
loglikhd!(bc::BlockCollection) = loglikhd!.(bc.blocks)

"""
    loglikhd°!(bc::BlockCollection)

For each `BiBlock` in a collection compute the log-likelihood for the proposal
block, evaluated at a sampled path and store the result in internal fields `ll`.
"""
loglikhd°!(bc::BlockCollection) = loglikhd°!.(bc.blocks)

@doc raw"""
    fetch_ll(bc::BlockCollection)

Retreive the log-likelihood for the entire accepted path.
!!! warning
    The function uses only internal fields `ll` for this computation, which
    means that for the call to this function to make sense the log-likelihood
    must have been previously computed and stored in the field `ll`. If it
    hasn't been done, then you must first call `loglikhd!(bc)`.
"""
fetch_ll(bc::BlockCollection) = mapreduce(bb->bb.b.ll, +, bc.blocks)

@doc raw"""
    fetch_ll°(bc::BlockCollection)

Retreive the log-likelihood for the entire proposed path.
!!! warning
    The function uses only internal fields `ll` for this computation, which
    means that for the call to this function to make sense the log-likelihood
    must have been previously computed and stored in the field `ll`. If it
    hasn't been done, then you must first call `loglikhd°!(bc)`.
"""
fetch_ll°(bc::BlockCollection) = mapreduce(bb->bb.b°.ll, +, bc.blocks)

"""
    save_ll!(bc::BlockCollection, i::Int)

For each block in the collection commit the current proposal and accepted
log-likelihood fields `ll` to history, at index `i`.
"""
save_ll!(bc::BlockCollection, i::Int) = save_ll!.(bc.blocks, i)

"""
    ll_of_accepted(bb::BlockCollection, i)

Return an array of log-likelihoods (one for each block) of the paths that were
accepted at the `i`th iteration.
"""
ll_of_accepted(bc::BlockCollection, i) = ll_of_accepted.(bc.blocks, i)

"""
    accpt_rate(bb::BlockCollection, range)

Compute the acceptance rate over the `range` of MCMC accept/reject history for
each block in the collection.
"""
function accpt_rate(bc::BlockCollection, range)
    map(bc.blocks) do biblock
        accpt_rate(biblock, range)
    end
end

#===============================================================================

        SETTING ARTIFICIAL OBSERVATIONS

===============================================================================#

#===============================================================================

        SETTING PARAMETERS

===============================================================================#

"""
    is_critical_update(bb::BlockCollection, pnames)

For each block in the collection verify whether the update characterized by a
list of parameter names stored in `pnames` is `critical` in a sense of prompting
for recomputation of the guiding term.
"""
GP.is_critical_update(bc::BlockCollection, pnames) = map(
    bb->GP.is_critical_update(bb.b.PP, pnames.θ_local_aux, pnames.θ_local_obs),
    bc.blocks
)

"""
    set_proposal_law!(
        bc::BlockCollection,
        θ°,
        pnames,
        critical_change=GP.is_critical_update(bc, pnames);
        skip=0
    )

For each block in the collection set the parameters in `bb.b°.PP` and
`bb.b°.P_last` to `θ°` and make sure that all other parameters are shared with
`bb.b.PP` and `bb.b.P_last`. Recompute the guiding term if needed, and then,
compute the proposal trajectory `bb.b°.XX` for the proposal point `θ°`.
"""
function set_proposal_law!(
        bc::BlockCollection,
        θ°,
        pnames,
        critical_change=GP.is_critical_update(bc, pnames);
        skip=0
    )
    for i in eachindex(bc.blocks)
        set_proposal_law!(
            bc.blocks[i], θ°, pnames.blocks[i], critical_change[i]; skip=skip
        )
    end
end
