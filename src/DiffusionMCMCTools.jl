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
    include("block.jl")
    include("biblock.jl")

    export SamplingUnit, draw_proposal_path!
    export SamplingPair
    export Block, find_W_for_X!, recompute_path!
    export BiBlock, swap_paths!, swap_XX!, swap_WW!, swap_PP!, set_proposal_law!
    export save_ll!, swap_ll!
    export accept_reject_proposal_path!, ll_of_accepted, accpt_rate
end
