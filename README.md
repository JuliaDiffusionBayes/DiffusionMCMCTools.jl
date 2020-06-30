<h1 align="center">
  <br>
  <a href="https://juliadiffusionbayes.github.io/DiffusionMCMCTools.jl/dev/"><img src="https://raw.githubusercontent.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl/master/docs/src/assets/logo.png" alt="DiffusionMCMCTools.jl" width="200"></a>
  <br>
  DiffusionMCMCTools.jl
  <br>
</h1>

> Structures and routines that facilitate more convenient code solutions to problems of inference and smoothing for diffusion processes.

<p align="center">
  <a href="https://JuliaDiffusionBayes.github.io/DiffusionMCMCTools.jl/stable">
    <img src="https://img.shields.io/badge/docs-stable-blue.svg"
         alt="Stable">
  </a>
  <a href="https://JuliaDiffusionBayes.github.io/DiffusionMCMCTools.jl/dev"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev"></a>
  <a href="https://travis-ci.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl">
      <img src="https://travis-ci.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl.svg?branch=master" alt="Build Status">
  </a>
</p>

<p align="center">
  <a href="#key-features">Key Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#how-to-use">How To Use</a> •
  <a href="#related">Related</a> •
  <a href="#license">License</a>
</p>

## Key features
- Structs that gather all containers needed for inference and smoothing of diffusions in one place
- Structs that provide views into structs above i.e. facilitate use of the *blocking* technique
- Structs that provide a translation between
  - names of the parameters as known by the Markov chain updating parameter values
  - and names of the parameters as known by the diffusion laws or observations
- Multiple methods for all of the above

## Installation

> :warning: This package is in an early development stage. In particular, **some functionality may depend on code in repositories [DiffusionDefinition.jl](https://github.com/JuliaDiffusionBayes/DiffusionDefinition.jl), [ObservationsSchemes.jl](https://github.com/JuliaDiffusionBayes/ObservationsSchemes.jl) or [GuidedProposals.jl](https://github.com/JuliaDiffusionBayes/GuidedProposals.jl) that I haven't pushed to those repos yet**. I would suggest not to use [DiffusionMCMCTools.jl](https://github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl) until this message disappears!

The package is not yet registered. To install it, type in:
```julia
] add https://github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl
```

## How To Use

See [the documentation](https://juliadiffusionbayes.github.io/DiffusionMCMCTools.jl/dev/).

## Related

DiffusionMCMCTools.jl belongs to a suite of packages in [JuliaDiffusionBayes](https://github.com/JuliaDiffusionBayes), whose aim is to facilitate Bayesian inference for diffusion processes. Some other packages in this suite are as follows:
- [DiffusionDefinition.jl](https://github.com/JuliaDiffusionBayes/DiffusionDefinition.jl): define diffusion processes and sample from their laws
- [ObservationSchemes.jl](https://github.com/JuliaDiffusionBayes/ObservationSchemes.jl): a systematic way of encoding discrete-time observations for stochastic processes
- [GuidedProposals.jl](https://github.com/JuliaDiffusionBayes/GuidedProposals.jl): defining and sampling conditioned diffusion processes
- [ExtensibleMCMC.jl](https://github.com/JuliaDiffusionBayes/ExtensibleMCMC.jl): a modular implementation of the Markov chain Monte Carlo (MCMC) algorithms
- [DiffusionMCMC.jl](https://github.com/JuliaDiffusionBayes/DiffusionMCMC.jl): Markov chain Monte Carlo (MCMC) algorithms for doing inference for diffusion processes

## License

MIT
