# A pairing of two `SamplingUnit`s
*****
Defines all main containers for an entire single recording for a smoothing or inference problem.
```@docs
DiffusionMCMCTools.SamplingPair
```

!!! note
    In practice all sampling is done with a `BiBlock` that looks at sections of a `SamplingPair` and never through `SamplingPair` directly. Even when no blocking is needed, it should still be done by defining a `BiBlock` that simply looks at the entire `SamplingPair`.
