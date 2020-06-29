"""
    struct BlockEnsemble{T}
        recordings::Vector{T}
    end

Gathers all blocks relevant for an entire single recording.

    BlockEnsemble(
        se::SamplingEnsemble,
        ranges,
        ρρ=0.0,
        ll_hist_len=0
    )

Base constructor.
"""
struct BlockEnsemble{T}
    recordings::Vector{T}

    function BlockEnsemble(
            se::SamplingEnsemble,
            ranges,
            ρρ=0.0,
            ll_hist_len=0
        )
        num_rec = num_recordings(se)
        ρρ = _vec_me(ρρ, num_rec)
        ll_hist_len = _vec_me(ll_hist_len, num_rec)
        rec = map(1:num_rec) do i
            BlockCollection(se.recordings[i], ranges[i], ρρ[i], ll_hist_len[i])
        end
        new{eltype(rec)}(rec)
    end
end

OBS.num_recordings(aob::BlockEnsemble) = length(aob.recordings)



#===============================================================================

        IMPUTATION

===============================================================================#
"""
    draw_proposal_path!(be::BlockEnsemble)

Call `draw_proposal_path!` separately for each recording
"""
draw_proposal_path!(be::BlockEnsemble) = draw_proposal_path!.(be.recordings)

#===============================================================================

        ACCEPT/REJECT DECISION

===============================================================================#

"""
    accept_reject_proposal_path!(be::BlockEnsemble, mcmciter)

Call `accept_reject_proposal_path!` separately for each recording.
"""
function accept_reject_proposal_path!(be::BlockEnsemble, mcmciter)
    for rec in be.recordings
        accept_reject_proposal_path!(rec, mcmciter)
    end
end

#===============================================================================

        SWAPS

===============================================================================#

"""
    swap_paths!(be::BlockEnsemble)

Call `swap_paths!` separately for each recording.
"""
swap_paths!(be::BlockEnsemble) = swap_paths!.(be.recordings)

"""
    swap_XX!(be::BlockEnsemble)

Call `swap_XX!` separately for each recording.
"""
swap_XX!(be::BlockEnsemble) = swap_XX!.(be.recordings)

"""
    swap_WW!(be::BlockEnsemble)

Call `swap_WW!` separately for each recording.
"""
swap_WW!(be::BlockEnsemble) = swap_WW!.(be.recordings)

"""
    swap_PP!(be::BlockEnsemble)

Call `swap_PP!` separately for each recording.
"""
swap_PP!(be::BlockEnsemble) = swap_PP!.(be.recordings)

"""
    swap_ll!(be::BlockEnsemble)

Call `swap_ll!` separately for each recording.
"""
swap_ll!(be::BlockEnsemble) = swap_ll!.(be.recordings)

#===============================================================================

        UTILITY

===============================================================================#

"""
    loglikhd!(be::BlockEnsemble)

Call `loglikhd!` separately for each recording.
"""
loglikhd!(be::BlockEnsemble) = loglikhd!.(be.recordings)

"""
    loglikhd°!(be::BlockEnsemble)

Call `loglikhd°!` separately for each recording.
"""
loglikhd°!(be::BlockEnsemble) = loglikhd°!.(be.recordings)

@doc raw"""
    fetch_ll(be::BlockEnsemble)

Retreive the log-likelihood for all accepted paths
!!! warning
    The function uses only internal fields `ll` for this computation, which
    means that for the call to this function to make sense the log-likelihood
    must have been previously computed and stored in the field `ll`. If it
    hasn't been done, then you must first call `loglikhd!(be)`.
"""
fetch_ll(be::BlockEnsemble) = mapreduce(rec->fetch_ll(rec), +, be.recordings)

@doc raw"""
    fetch_ll°(be::BlockEnsemble)

Retreive the log-likelihood for all proposed paths
!!! warning
    The function uses only internal fields `ll` for this computation, which
    means that for the call to this function to make sense the log-likelihood
    must have been previously computed and stored in the field `ll`. If it
    hasn't been done, then you must first call `loglikhd°!(be)`.
"""
fetch_ll°(be::BlockEnsemble) = mapreduce(rec->fetch_ll°(rec), +, be.recordings)

"""
    save_ll!(be::BlockEnsemble, i::Int)

Call `save_ll!` separately for each recording.
"""
save_ll!(be::BlockEnsemble, i::Int) = save_ll!.(be.recordings, i)

"""
    ll_of_accepted(be::BlockEnsemble, i)

Return an array of log-likelihoods (one for each block in each recording) of the
paths that were accepted at the `i`th iteration.
"""
ll_of_accepted(be::BlockEnsemble, i) = ll_of_accepted.(be.recordings, i)

"""
    accpt_rate(be::BlockEnsemble, range)

Compute the acceptance rate over the `range` of MCMC accept/reject history for
each block for each recording.
"""
function accpt_rate(be::BlockEnsemble, range)
    map(be.recordings) do rec
        accpt_rate(rec, range)
    end
end

#===============================================================================

        SETTING ARTIFICIAL OBSERVATIONS

===============================================================================#

"""
    GP.set_obs!(be::BlockEnsemble)

Call `save_ll!` separately for each recording.
"""
GP.set_obs!(be::BlockEnsemble) = set_obs!.(be.recordings)

"""
    GP.recompute_guiding_term!(be::BlockEnsemble, [::Val{:_only}])

Call `recompute_guiding_term!` separately for each recording.
"""
function GP.recompute_guiding_term!(be::BlockEnsemble)
    recompute_guiding_term!.(be.recordings)
end

function GP.recompute_guiding_term!(be::BlockEnsemble, v::Val)
    recompute_guiding_term!.(be.recordings, v)
end

"""
    find_W_for_X!(be::BlockEnsemble)

Call `find_W_for_X!` separately for each recording.
"""
find_W_for_X!(be::BlockEnsemble) = find_W_for_X!.(be.recordings)


#===============================================================================

        SETTING PARAMETERS

===============================================================================#

"""
    is_critical_update(be::BlockEnsemble, pnames)

Call `is_critical_update` separately for each recording.
"""
GP.is_critical_update(be::BlockEnsemble, pnames) = map(
    rec->GP.is_critical_update(rec, pnames),
    be.recordings
)

"""
    set_proposal_law!(
        be::BlockEnsemble,
        θ°,
        pnames,
        critical_change=GP.is_critical_update(bc, pnames);
        skip=0
    )

Call `set_proposal_law!` separately for each recording.
"""
function set_proposal_law!(
        be::BlockEnsemble,
        θ°,
        pnames,
        critical_change=GP.is_critical_update(bc, pnames);
        skip=0
    )
    for i in eachindex(be.recordings)
        set_proposal_law!(
            be.recordings[i], θ°, pnames.recordings[i], critical_change[i];
            skip=skip
        )
    end
end
