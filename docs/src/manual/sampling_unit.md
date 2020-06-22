# The smallest sampling unit: `SamplingUnit`
****
The struct
```@docs
DiffusionMCMCTools.SamplingUnit
```
is the main building block of all remaining units implemented in this package. It has a couple of methods implemented for it:
```@docs
DiffusionMCMCTools.recompute_guiding_term!(u::DiffusionMCMCTools.SamplingUnit)
DiffusionMCMCTools.loglikhd(u::DiffusionMCMCTools.SamplingUnit)
DiffusionMCMCTools.draw_proposal_path!(u::DiffusionMCMCTools.SamplingUnit)
```
However, because it does not contain proposal-accepted pair it is rarely used on its own for sampling. Instead, it most often appears as a member of other struct.
