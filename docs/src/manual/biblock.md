# BiBlock
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
```

#### Setting up a block
----
```@docs
DiffusionMCMCTools.set_obs!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.recompute_guiding_term!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.find_W_for_X!(bb::DiffusionMCMCTools.BiBlock)
```

#### Utility
----
```@docs
DiffusionMCMCTools.ll_of_accepted(bb::DiffusionMCMCTools.BiBlock, i)
```
