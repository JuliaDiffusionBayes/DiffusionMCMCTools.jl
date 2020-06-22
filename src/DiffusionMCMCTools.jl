module DiffusionMCMCTools

    using DiffusionDefinition, ObservationSchemes, GuidedProposals

    const DD = DiffusionDefinition
    const OBS = ObservationSchemes
    const GP = GuidedProposals
    const _LL = Val(:ll)

    OBS.var_parameter_names(P::DD.DiffusionProcess) = DD.var_parameter_names(P)

    include("sampling_unit.jl")
    include("block_unit.jl")

end
