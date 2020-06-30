# [BlockCollection](@id block_collection_def)
*****
```@docs
DiffusionMCMCTools.BlockCollection
```

## Imputation
```@docs
DiffusionMCMCTools.draw_proposal_path!(bc::DiffusionMCMCTools.BlockCollection)
```

## Accept/reject decision
```@docs
DiffusionMCMCTools.accept_reject_proposal_path!(bb::DiffusionMCMCTools.BlockCollection, mcmciter)
```

## Swaps
```@docs
DiffusionMCMCTools.swap_paths!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.swap_XX!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.swap_WW!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.swap_PP!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.swap_ll!(bc::DiffusionMCMCTools.BlockCollection)
```

#### Setting up a block

```@docs
DiffusionMCMCTools.set_obs!(bc::DiffusionMCMCTools.BlockCollection)
```

## utility
```@docs
DiffusionMCMCTools.loglikhd!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.loglikhd°!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.fetch_ll(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.fetch_ll°(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.save_ll!(bc::DiffusionMCMCTools.BlockCollection, i::Int)
DiffusionMCMCTools.ll_of_accepted(bb::DiffusionMCMCTools.BlockCollection, i)
DiffusionMCMCTools.accpt_rate(bb::DiffusionMCMCTools.BlockCollection, range)
```

## Setting parameters
```@docs
DiffusionMCMCTools.set_proposal_law!(
    bc::DiffusionMCMCTools.BlockCollection,
    θ°,
    pnames,
    critical_change=GP.is_critical_update(bc, pnames);
    skip=0
)
```
