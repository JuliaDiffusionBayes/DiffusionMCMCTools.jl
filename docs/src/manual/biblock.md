# BiBlock
****
A composite unit that represents a block. It can be used for sampling on a block in a setting of smoothing or inference.
```@docs
DiffusionMCMCTools.BiBlock
```

There are many functions implemented for it.


#### For adjustments made after the accept-reject decision (regardless of what it was)
----
```@docs
DiffusionMCMCTools.set_accepted!(bb::DiffusionMCMCTools.BiBlock, i::Int, v)
```

#### For adjustments to the containers in case of **acceptance** of proposals:
----
```@docs
DiffusionMCMCTools.swap_paths!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.swap_XX!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.swap_WW!(bb::DiffusionMCMCTools.BiBlock)
DiffusionMCMCTools.swap_PP!(bb::DiffusionMCMCTools.BiBlock)
```
