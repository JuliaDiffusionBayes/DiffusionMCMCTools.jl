#===============================================================================

    The main object defined in this file is: SamplingUnit.
    It is the smallest `composite unit` (i.e. a unit consisting of at least
    two `GuidProp`/containers for paths etc.) and it can be used for sampling
    conditioned diffusions based on conditioning on multiple points.

===============================================================================#
"""
    struct SamplingUnit{TGP,TGPb,TW,TWn,TX}
        PP::Vector{TGP}
        PPb::Vector{TGPb}
        WW::Vector{TW}
        Wnr::TWn
        XX::Vector{TX}
    end

Smallest composite unit with containers needed for sampling of conditioned
diffusions via Guided Proposals.

# Fields
---
- `PP`: a vector of `GuidProp`
- `PPb`: a vector of `GuidProp` that can be used in a blocking schemes as
         `GuidProp`s on terminal subintervals
- `WW`: a vector of containers for sampled Wiener process
- `Wnr`: a flag for sampling Wiener processes
- `XX`: a vector of containers for a sampled process
"""
struct SamplingUnit{TGP,TGPb,TW,TWn,TX}
    PP::Vector{TGP}
    PPb::Vector{TGPb}
    WW::Vector{TW}
    Wnr::TWn
    XX::Vector{TX}

    function SamplingUnit(
            aux_laws, recording, x0_prior, tts, args;
            aux_laws_blocking=aux_laws, artificial_noise=1e-11,
            solver_choice_blocking=args
        ) where T
        PP = build_guid_prop(aux_laws, recording, tts, args)
        PPb = guid_prop_for_blocking(
            PP,
            aux_laws_blocking,
            artificial_noise,
            solver_choice_blocking
        )
        XX, WW = trajectory(PP)
        Wnr = Wiener(PP)

        init_paths!(P, WW, Wnr, XX, x0_prior)
        new{eltype(PP),eltype(PPb),eltype(WW),typeof(Wnr),eltype(XX)}(
            PP, PPb, WW, Wnr, XX
        )
    end
end

"""
    init_paths!(P, WW, Wnr, XX, x0_prior)

Sample paths of guided proposals without using the preconditioned Crank–Nicolson
scheme.
"""
function init_paths!(PP, WW, Wnr, XX, x0_prior)
    while true
        forward_guide!(PP, XX, WW, rand(x0_prior); Wnr=Wnr) && return
    end
end

function all_eltypes(
        ::SamplingUnit{TGP,TGPb,TW,TWn,TX}
    ) where {TGP,TGPb,TW,TWn,TX}
    TGP,TGPb,TW,TWn,TX
end

"""
    GP.recompute_guiding_term!(u::SamplingUnit)

Recompute the guiding term assuming `u.PP` is a vector with guided proposal laws
"""
function GP.recompute_guiding_term!(u::SamplingUnit)
    recompute_guiding_term!(u.PP)
end

"""
    GP.loglikhd(u::SamplingUnit)

Return log-likelihood evaluated at a sampled path
"""
GP.loglikhd(u::SamplingUnit) = loglikhd(u.PP, u.XX)

"""
    draw_proposal_path!(u::SamplingUnit)

Sample a proposal path, compute log-likelihood along the way. Assumes
`u.XX[1].x[1]` is a starting point. No preconditioned Crank-Nicolson scheme is
used.
"""
function draw_proposal_path!(u::SamplingUnit)
    rand!(u.PP, u.XX, u.WW, _LL, u.XX[1].x[1]; Wnr=u.Wnr)
end
