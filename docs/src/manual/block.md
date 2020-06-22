# Block
****
A smallest unit that is needed for blocking.
```@docs
DiffusionMCMCTools.Block
```
It is rarely used on its own. Instead, it is used mainly as a building block of other, composite units. Nevertheless, there are a couple of useful functions implemented for it:

```@docs
DiffusionMCMCTools.set_ll!(b::DiffusionMCMCTools.Block, i::Int, v)
DiffusionMCMCTools.save_ll!(b::DiffusionMCMCTools.Block, i::Int)
DiffusionMCMCTools.recompute_guiding_term!(b::DiffusionMCMCTools.Block)
DiffusionMCMCTools.find_W_for_X!(b::DiffusionMCMCTools.Block)
DiffusionMCMCTools.loglikhd(b::DiffusionMCMCTools.Block)
DiffusionMCMCTools.recompute_path!(b::DiffusionMCMCTools.Block, WW=b.WW; skip=0)
```
