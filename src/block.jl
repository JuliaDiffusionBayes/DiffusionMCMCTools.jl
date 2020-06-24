#===============================================================================

    The main object defined in this file is: Block.
    It is the smallest `composite unit` needed for blocking. It provides
    a window view into `SamplingUnit` that restricts the range of view
    to some range i:j that corresponds to a chosen block.

===============================================================================#


const TVIEW{T} = SubArray{T,1,Array{T,1},Tuple{UnitRange{Int64}},true}

"""
    mutable struct Block{L,TGP,TGPl,TW,TWn,TX}
        PP::TVIEW{TGP}
        P_last::TVIEW{TGPl} # view into a single element
        WW::TVIEW{TW}
        Wnr::TWn
        XX::TVIEW{TX}
        ll::Float64
        ll_history::Vector{Float64}
    end

The smallest "containerless" unit that provides a view into `SamplingUnit`
restricted to a range `i:j` of a block. `L` is an important flag that indicates
whether it is a terminal block or not.

# Fields
---
- `PP`: a vector of views into relevant `GuidProp`
- `P_last`: a view into a single (or none) `GuidProp` corresponding to the
            terminal sub-interval.
- `WW`: a vector of views into containers for sampled Wiener process
- `Wnr`: a flag for sampling Wiener processes
- `XX`: a vector of views into containers for a sampled process
- `ll`: a placeholder for computed log-likelihood
- `ll_history`: a history of computed log-likelihoods (useful for MCMC)


    function Block(
        u::SamplingUnit,
        range::UnitRange{Int64},
        last_block=false,
        ll_hist_len=0
    )

Base constructor.
"""
mutable struct Block{L,TGP,TGPl,TW,TWn,TX}
    PP::TVIEW{TGP}
    P_last::TVIEW{TGPl} # view to a single element
    WW::TVIEW{TW}
    Wnr::TWn
    XX::TVIEW{TX}
    ll::Float64
    ll_history::Vector{Float64}

    function Block(
            u::SamplingUnit{TGP,TGPl,TW,TWn,TX},
            range::UnitRange{Int64},
            last_block=false,
            ll_hist_len=0
        ) where {TGP,TGPl,TW,TWn,TX}
        PP = view(u.PP, range[1]:(range[end]-!last_block)) # omit the last law
        P_last = view(u.PPb, range[end]:range[end])

        XX = view(u.XX, range)
        WW = view(u.WW, range)

        new{last_block,TGP,TGPl,TW,TWn,TX}(
            PP, P_last, WW, u.Wnr, XX, -Inf,
            Vector{Float64}(undef, ll_hist_len)
        )
    end
end

"""
    set_ll!(b::Block, i::Int, v)

Set the internal log-likelihood history field `b.ll_history[i]` with a value `v`
"""
set_ll!(b::Block, i::Int, v) = (b.ll_history[i] = v)

"""
    save_ll!(b::Block, i::Int)

Commit the current log-likelihood field `ll` to history `b.ll_history` at index
`i`.
"""
save_ll!(b::Block, i::Int) = set_ll!(b, i, b.ll)


"""
    GP.recompute_guiding_term!(b::Block)

Recompute the guiding term.
"""
function GP.recompute_guiding_term!(b::Block) end

function GP.recompute_guiding_term!(b::Block{false})
    recompute_guiding_term!(b.PP, b.P_last[1])
end

function GP.recompute_guiding_term!(b::Block{true})
    recompute_guiding_term!(b.PP)
end

"""
    find_W_for_X!(b::Block)

Compute the Wiener process `b.WW` that is needed for obtaining path `b.XX` under
the law stored in `b`.
"""
function find_W_for_X!(b::Block) end

function find_W_for_X!(b::Block{false})
    for i in eachindex(b.PP)
        DD.invsolve!(b.XX[i], b.WW[i], b.PP[i])
    end
    DD.invsolve!(b.XX[end], b.WW[end], b.P_last[1])
end

function find_W_for_X!(b::Block{true})
    for i in eachindex(b.PP)
        DD.invsolve!(b.XX[i], b.WW[i], b.PP[i])
    end
end

"""
    GP.loglikhd(b::Block)

Compute the log-likelihood evaluated at a sampled path.
"""
function GP.loglikhd(b::Block) end

GP.loglikhd(b::Block{false}) = (
    loglikhd(b.PP, b.XX) + loglikhd(b.P_last[1], b.XX[end])
)

GP.loglikhd(b::Block{true}) = loglikhd(b.PP, b.XX)

"""
    recompute_path!(b::Block, WW=b.WW; skip=0)

Recompute the path `b.XX` for a given wiener process `WW`.
"""
function recompute_path!(b::Block, WW=b.WW; skip=0) end

function recompute_path!(b::Block{false}, WW=b.WW; skip=0)
    success = _recompute_path!(b, WW; skip=skip)
    success || return false
    y1 = b.XX[end-1].x[end]
    success, ll = GP.solve_and_ll!(
        b.XX[end], WW[end], b.P_last[1], y1; skip=skip
    )
    b.ll += ll
    success
end

function recompute_path!(b::Block{true}, WW=b.WW; skip=0)
    _recompute_path!(b, WW; skip=skip)
end

function _recompute_path!(b::Block, WW=b.WW; skip=0)
    ll_tot = loglikhd_obs(b.PP[1], y1)
    y1 = b.XX[1].x[1]
    for i in 1:length(b.PP)
        success, ll = GP.solve_and_ll!(b.XX[i], WW[i], b.PP[i], y1; skip=skip)
        success || (b.ll = ll; return false)
        ll_tot += ll
        y1 = b.XX[i].x[end]
    end
    b.ll = ll_tot
    true
end
