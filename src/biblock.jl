#===============================================================================

    The main object defined in this file is: BiBlock.
    It represents a single block and is the simplest composite unit
    that can be used for practical settings of inference and smoothing.
    Additionally, it contains the option to sample with a preconditioned
    Crank-Nicolson step.

===============================================================================#



#===============================================================================

        DEFINITION

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

#===============================================================================

        IMPUTATION

===============================================================================#

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
    y1 = bb.b°.XX[end-1].x[end]
    success, ll°_last = _draw_proposal_path_last_segment!(bb, y1)
    bb.b°.ll += ll°_last
    success
end

function draw_proposal_path!(bb::BiBlock{true})
    success, bb.b°.ll = _draw_proposal_path!(bb)
    success
end

function _draw_proposal_path!(bb::BiBlock)
    rand!(
        bb.b.PP, bb.b°.XX, bb.b°.WW, bb.b.WW, bb.ρ, _LL, bb.b.XX[1].x[1];
        Wnr=bb.b°.Wnr
    )
end

function _draw_proposal_path_last_segment!(bb::BiBlock, y1)
    rand!(
        bb.b.P_last[1], bb.b°.XX[end], bb.b°.WW[end], bb.b.WW[end], bb.ρ, _LL,
        y1; Wnr=bb.b°.Wnr
    )
end


#===============================================================================

        ACCEPT/REJECT DECISION

===============================================================================#

"""
    accept_reject_proposal_path!(bb::BiBlock, mcmciter)

Accept/reject decision of the Metropolis-Hastings algorithm for the step of path
imputation.
"""
function accept_reject_proposal_path!(bb::BiBlock, mcmciter)
    accepted = rand(Exponential(1.0)) > -(bb.b°.ll - bb.b.ll)
    accepted && swap_paths!(bb)
    set_accepted!(bb, mcmciter, accepted)
    save_ll!(bb, mcmciter)
    accepted && swap_ll!(bb)
end

"""
    set_accepted!(bb::BiBlock, i::Int, v)

Commit the accept/reject decision `v` to acceptance history of BiBlock `b` at
the position `i`.
"""
set_accepted!(bb::BiBlock, i::Int, v) = (bb.accpt_history[i] = v)

#===============================================================================

        SWAPS

===============================================================================#

"""
    swap_paths!(bb::BiBlock)

Swap `XX` and `WW` containers between proposal-acceptance pair.
"""
function swap_paths!(bb::BiBlock)
    swap_XX!(bb)
    swap_WW!(bb)
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
function swap_PP!(bb::BiBlock) end

function swap_PP!(bb::BiBlock{true})
    _swap_PP!(bb)
end

function swap_PP!(bb::BiBlock{false})
    _swap_PP!(bb)
    bb.b.P_last[1], bb.b°.P_last[1] = bb.b°.P_last[1], bb.b.P_last[1]
    bb.b.P_excl[1], bb.b°.P_excl[1] = bb.b°.P_excl[1], bb.b.P_excl[1]
    for i in eachindex(bb.b.Pb_excl)
        bb.b.Pb_excl[i], bb.b°.Pb_excl[i] = bb.b°.Pb_excl[i], bb.b.Pb_excl[i]
    end
end

function _swap_PP!(bb::BiBlock)
    for i in eachindex(bb.b.PP)
        bb.b.PP[i], bb.b°.PP[i] = bb.b°.PP[i], bb.b.PP[i]
    end
end


"""
    swap_ll!(bb::BiBlock)

Swap `ll` containers between proposal-acceptance pair.
"""
function swap_ll!(bb::BiBlock)
    bb.b.ll, bb.b°.ll = bb.b°.ll, bb.b.ll
end

#===============================================================================

        UTILITY

===============================================================================#

"""
    ll_of_accepted(bb::BiBlock, i)

Return the log-likelihood of the path that was accepted at the `i`th iteration.
"""
function ll_of_accepted(bb::BiBlock, i)
    return (bb.accpt_history[i] ? bb.b°.ll_history[i] : bb.b.ll_history[i])
end


"""
    accpt_rate(bb::BiBlock, range)

Compute the acceptance rate over the `range` of MCMC accept/reject history.
"""
accpt_rate(bb::BiBlock, range) = sum(bb.accpt_history[range])/length(range)

"""
    loglikhd!(b::BiBlock)

Compute the log-likelihood for the accepted block, evaluated at a sampled path
and store the result in an internal field `ll`.
"""
loglikhd!(bb::BiBlock) = loglikhd!(bb.b)

"""
    loglikhd°!(b::BiBlock)

Compute the log-likelihood for the proposal block, evaluated at a sampled path
and store the result in an internal field `ll`.
"""
loglikhd°!(bb::BiBlock) = loglikhd!(bb.b°)

"""
    save_ll!(bb::BiBlock, i::Int)

Commit the current proposal and accepted log-likelihood fields `ll` to history,
at index `i`.
"""
function save_ll!(bb::BiBlock, i::Int)
    save_ll!(bb.b, i)
    save_ll!(bb.b°, i)
end

#===============================================================================

        SETTING ARTIFICIAL OBSERVATIONS

===============================================================================#

"""
    GP.set_obs!(bb::BiBlock)

Freeze an artificial observation at the terminal point of the block. For a
terminal block nothing is done.
"""
function GP.set_obs!(bb::BiBlock) end

function GP.set_obs!(bb::BiBlock{false})
    set_obs!(bb.b.P_last[1], bb.b.XX[end].x[end])
    set_obs!(bb.b°.P_last[1], bb.b.XX[end].x[end])
end

GP.set_obs!(bb::BiBlock{true}) = nothing


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

#===============================================================================

        SETTING PARAMETERS

===============================================================================#

"""
    is_critical_update(bb::BiBlock, pnames)

Verify whether the update characterized by a list of parameter names stored in
`pnames` is `critical` in a sense of prompting for recomputation of the guiding
term.
"""
GP.is_critical_update(bb::BiBlock, pnames) = GP.is_critical_update(
    bb.b.PP, pnames.θ_local_aux, pnames.θ_local_obs
)


"""
    set_proposal_law!(
        bb::BiBlock,
        θ°,
        pnames,
        critical_change=GP.is_critical_update(bb, pnames);
        skip=0
    )

Set the parameters in `bb.b°.PP` and `bb.b°.P_last` to `θ°` and make sure that
all other parameters are shared with `bb.b.PP` and `bb.b.P_last`. Recompute the
guiding term if needed, and then, compute the proposal trajectory `bb.b°.XX` for
the proposal point `θ°`.
"""
function set_proposal_law!(
        bb::BiBlock,
        θ°,
        pnames,
        critical_change=GP.is_critical_update(bb, pnames);
        skip=0
    )
    critical_change = DD.set_parameters!(bb, θ°, pnames, critical_change)
    critical_change && recompute_guiding_term!(bb.b°)
    recompute_path!(bb.b°, bb.b.WW; skip=skip)
end


"""
    DD.set_parameters!(
        bb::BiBlock,
        θ°,
        pnames,
        critical_change=is_critical_update(bb, pnames)
    )

Set the parameters in `bb.b°.PP` and `bb.b°.P_last` to `θ°` and make sure that
all other parameters are shared with `bb.b.PP` and `bb.b.P_last`.
"""
function DD.set_parameters!(
        bb::BiBlock,
        θ°,
        pnames,
        critical_change=GP.is_critical_update(bb, pnames)
    )
    GP.equalize_obs_params!(bb) && (critical_change = true)
    GP.equalize_law_params!(bb, pnames) && (critical_change = true)
    DD.set_parameters!(bb.b°.PP, θ°, pnames.PP)
    DD.set_parameters!(bb.b°.P_last, θ°, pnames.P_last)
    DD.set_parameters!(bb.b°.P_excl, θ°, pnames.P_excl)
    DD.set_parameters!(bb.b°.Pb_excl, θ°, pnames.Pb_excl)
    critical_change
end

function DD.set_parameters!(PP::AbstractArray{<:GuidProp}, θ°, pnames)
    DD.set_parameters!(PP, θ°, pnames.updt, pnames.updt_aux, pnames.updt_obs)
end

"""
    GP.equalize_obs_params!(bb::BiBlock)

Make sure that in the corresponding pairs of `GuidProp` structs of both `bb.b`
and `bb.b°`, the parameters `θ` of the `obs` fields are the same. If not, then
set the ones in `bb.b°` to be the same as the ones in `bb.b`.
"""
function GP.equalize_obs_params!(bb::BiBlock)
    # equalization of observations that does not matter for this update, but
    # will once the blocks switch...
    GP.equalize_obs_params!(bb.b.P_excl, bb.b°.P_excl)

    # no equalization on P_last or Pb_excl

    # and finally equalization that matters currently, may yield critical change
    GP.equalize_obs_params!(bb.b.PP, bb.b°.PP)
end

"""
    GP.equalize_law_params!(bb::BiBlock, pnames)

Make sure that in the corresponding pairs of `GuidProp` structs of both `bb.b`
and `bb.b°`, the `variable` parameters of each law are the same. If not, then
set the ones in `bb.b°` to be the same as the ones in `bb.b`.
"""
function GP.equalize_law_params!(bb::BiBlock, pnames) end

function GP.equalize_law_params!(bb::BiBlock{true}, pnames)
    _eql_PP!(bb.b.PP, bb.b°.PP, pnames.PP.var, pnames.PP.var_aux)
end

function GP.equalize_law_params!(bb::BiBlock{false}, pnames)
    critical_change = _eql_PP!(
        bb.b.PP, bb.b°.PP, pnames.PP.var, pnames.PP.var_aux
    )

    _eql_PP!(
        bb.b.P_last, bb.b°.P_last, pnames.P_last.var, pnames.P_last.var_aux
    ) && (critical_change = true)

    # this is just not needed now, but will be for the next blocks
    _eql_PP!(
        bb.b.P_excl, bb.b°.P_excl, pnames.P_excl.var, pnames.P_excl.var_aux
    )
    _eql_PP!(
        bb.b.Pb_excl, bb.b°.Pb_excl, pnames.Pb_excl.var, pnames.Pb_excl.var_aux
    )

    critical_change
end

"""
    _eql_PP!(PP, PP°, var_p_names, var_p_aux_names)

Go through all `GuidProp` structs in  `PP` and `PP°` and make sure that the
parameters listed in `var_p_names` and `var_p_aux_names` agree. Equalize them if
they do not. `var_p_names` should list all parameter names from the target law
and `var_p_aux_names` should list all parameters from the auxiliary law.
"""
function _eql_PP!(PP, PP°, var_p_names, var_p_aux_names)
    if !DD.same_entries(PP, PP°, var_p_names)
        GP.equalize_law_params!(
            PP, PP°, var_p_names, var_p_aux_names
        ) && return true
    end
    false
end
