# DiffusionMCMCTools
**********

This is a utility package belonging to a suite of packages in
[JuliaDiffusionBayes](https://github.com/JuliaDiffusionBayes). It provides structs and methods that make it easier to write algorithms doing inference and/or smoothing for diffusion processes.

There are three categories of objects (that have their own, associated methods) provided in this package:
- Structs that gather all containers that are needed for sampling/inference of diffusion processes
- Structs that provide windowed views into the structs above (and hence, facilitate use of the so-called *blocking* technique)
- Structs that provide a translation between parameter names:
  - parameter names known to MCMC sampler
  - parameter names known to each individual diffusion law or observation

--------------------

Depending on your intended use of this package you might choose to start at different places:

- For a quick overview of [DiffusionMCMCTools.jl](https://github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl)'s main functionality see [Get started](@ref get_started)
- For a systematic introduction to all functionality introduced in this package see the [Manual](@ref manual_start)
- For a didactic introduction to problems that can be solved using [DiffusionMCMCTools.jl](https://github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl) see the [Tutorials](@ref tutorial_start)
- If you have a problem that you think can be addressed with this package, then check out the [How-to guides](@ref how_to_guides_start) to see if the answer is already there.
