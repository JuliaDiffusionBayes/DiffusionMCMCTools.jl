module DiffusionMCMCTools

    using DiffusionDefinition, ObservationSchemes, GuidedProposals
    using Random, Distributions

    const DD = DiffusionDefinition
    const OBS = ObservationSchemes
    const GP = GuidedProposals
    const _LL = Val(:ll)

    #temporary
    import GuidedProposals: is_critical_update

    OBS.var_parameter_names(P::DD.DiffusionProcess) = DD.var_parameter_names(P)

    include("sampling_unit.jl")
    include("sampling_pair.jl")
    include("sampling_ensemble.jl")

    include("block.jl")
    include("biblock.jl")
    include("block_collection.jl")
    include("block_ensemble.jl")

    include("param_names_collections.jl")

    # sampling_unit.jl
    export SamplingUnit
    export draw_proposal_path!

    # sampling_pair.jl
    export SamplingPair

    # sampling_ensemble.jl
    export SamplingEnsemble

    # block.jl
    export Block
    export set_ll!, save_ll!, find_W_for_X!, recompute_path!, loglikhd!

    # biblock.jl
    export BiBlock
    export accept_reject_proposal_path!, set_accepted!
    export swap_paths!, swap_XX!, swap_WW!, swap_PP!, swap_ll!
    export ll_of_accepted, accpt_rate, loglikhd°!
    export set_proposal_law!

    # block_collection.jl
    export BlockCollection
    export fetch_ll, fetch_ll°

    # block_ensemble.jl
    export BlockEnsemble

    # param_names_collections.jl
    export ParamNamesUnit
    export ParamNamesBlock
    export ParamNamesRecording
    export ParamNamesAllObs
end
