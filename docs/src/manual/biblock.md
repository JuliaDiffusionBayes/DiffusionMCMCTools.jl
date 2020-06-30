# [BiBlock](@id biblock_def)
****
A composite unit that represents a block. It can be used for sampling on a block in a setting of smoothing or inference.
```@docs
DiffusionMCMCTools.BiBlock
```

There are many functions implemented for it.

#### Imputation of paths
----
```@docs
DiffusionMCMCTools.draw_proposal_path!(bb::DiffusionMCMCTools.BiBlock)
```

#### Accept/reject decision in an MCMC setting
-----
```@docs
DiffusionMCMCTools.accept_reject_proposal_path!(bb::DiffusionMCMCTools.BiBlock, mcmciter)
```

#### Adjustments made after the accept-reject decision (regardless of what it was)
----
```@docs
DiffusionMCMCTools.set_accepted!(bb::DiffusionMCMCTools.BiBlock, i::Int, v)
```

#### Adjustments to the containers in case of **acceptance** of proposals:
----
```@docs
DiffusionMCMCTools.swap_paths!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.swap_XX!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.swap_WW!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.swap_PP!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.swap_ll!(bb::DiffusionMCMCTools.BiBlock)
```

#### Setting up a block
----
```@docs
DiffusionMCMCTools.set_obs!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.recompute_guiding_term!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.find_W_for_X!(bb::DiffusionMCMCTools.BiBlock)
```

#### Setting parameters
----
```@docs
DiffusionMCMCTools.set_proposal_law!(
    bb::DiffusionMCMCTools.BiBlock,
    θ°,
    pnames,
    critical_change=DiffusionMCMCTools.is_critical_update(bb, pnames),
    skip=0
)

```

#### Utility
----
```@docs
DiffusionMCMCTools.ll_of_accepted(bb::DiffusionMCMCTools.BiBlock, i)
DiffusionMCMCTools.accpt_rate(bb::DiffusionMCMCTools.BiBlock, range)
DiffusionMCMCTools.loglikhd!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.loglikhd°!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.save_ll!(bb::DiffusionMCMCTools.BiBlock, i::Int)
```



!!! note
    As you can see above, `BiBlock` is the fundamental building block that is used for creating inference and smoothing algorithms. However, as the complexity of these algorithms grow it is useful to use some `macro structures` that operate on or are defined for multiple `BiBlock`s. This is precisely what the remaining tools defined in this package are for. Otherwise put, they aim to facilitate writing snippets of code as above in a much more compact and convenient way.
