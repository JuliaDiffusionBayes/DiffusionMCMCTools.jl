# [Get started](@id get_started)
****

## Installation
----

```julia
] add https://github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl
```

## Structs that neatly gather all necessary containers
-----
As the title suggests. An important point is that they don't have almost any methods implemented for them and thus cannot be used for sampling or inference. You can think of them as essentially bags that carry the objects. For sampling and inference see various blocks that are presented below.

### `SamplingPair`
```julia
AuxLaw = ... # `Type` of an auxiliary diffusion (NOT its instance)
recording = (P = ..., obs = ..., t0 = ..., x0_prior = ...)
tts = ... # time grids corresponding to `obs` in `recording`

sampling_pair = SamplingPair(AuxLaw, recording, tts)
```
`sampling_pair` will contain all containers—such as `u.PP`, `u.WW`, `u.XX` i.e. accepted guided proposals, Wiener process and path of a process, `u°.PP`, `u°.WW`, `u°.XX` i.e. their proposed versions and more—that are needed for smoothing or inference based on the `recording`.

### `SamplingEnsemble`
Same as `SamplingPair`, but for multiple recordings.
```julia
AuxLaw = ... # `Type` of an auxiliary diffusion (NOT its instance)

_all_obs = AllObservations() # a structure with multiple recordings and dependency structures
add_recording!(_all_obs, ...)
all_obs, _ = initialize(all_obs)

tts = ... # time grids corresponding to `recordings` in `all_obs.recordings`

sampling_ensemble = SamplingEnsemble(AuxLaw, all_obs.recordings, tts)
```
`sampling_ensemble` will contain an instance of `SamplingPair` for each recording.

## Structs that act as windowed views into containers above
-----
Unlike the containers above, the block-structures defined below have many methods implemented for them so that they can be used for sampling and inference.

### `BiBlock`
```julia
last_block = true
block = BiBlock(sampling_pair, start_idx:end_idx, ρ, last_block, num_steps)
```
A `block` above is a view into containers of `sampling_pair` that have indices running through `start_idx:end_idx`. `last_block` is a flag that indicates whether we want `block` to be the terminal block for a given recording. If set to `true` (as above), then `block` will essentially point to `sampling_pair.u.PP[start_idx:end_idx]`, `sampling_pair.u°.PP[start_idx:end_idx]`, `sampling_pair.u.XX[start_idx:end_idx]` etc.. On the other hand, if set to `false`, then for pointers to law, i.e. to `sampling_pair.u.PP` and `sampling_pair.u°.PP`, the indices would run through `start_idx:end_idx` instead and the terminal law would be substituted with a `P_last` that can handle an exact terminal observation.

We can now perform many different operations on `block`. For instance:
```julia
# to draw a proposal path
draw_proposal_path!(block)

# use Metropolis–Hastings accept/reject rule to decide whether proposal path should be accepted
accept_reject_proposal_path!(block, mcmc_iteration_idx)

# to compute acceptance rate based on the last 100 path imputations
accpt_rate(block, (mcmc_iteration_idx-99):mcmc_iteration_idx)

# to freeze the terminal point of a block to become an artificial observation
set_obs!(block)

# to recompute the guiding term only on the "accepted" laws `bb.b.PP`
recompute_guiding_term!(block)

# to recompute the Wiener path that reconstructs the accepted path under the accepted law
find_W_for_X!(block)

# re-evaluate the log-likelihood and store it in a local field
loglikhd!(block)

# set the parameters inside the proposal law
set_proposal_law!(block, θ°, #= to be explained =# name_struct #= to be explained =#, true)
# name_struct should be an instance of `ParamNamesBlock`
```
### `BlockCollection`
```julia
block_ranges = ... # for instance [1:15, 16:30, 31:50, 51:100]
block_collection = BlockCollection(sampling_pair, block_ranges, ρ, num_steps)
```
`block_collection` will contain multiple `BiBlock`s that together split the recording into segments called *blocks*. You can call all methods listed above on `block_collection` instead of calling them separately on each `BiBlock` and the result will be most likely as you'd expect it: the function will be recursively applied to each block that `block_collection` comprises of.

### `BlockEnsemble`
`BlockEnsemble` is simply a collection of `BlockCollection`s, where a single instance of `BlockCollection` is defined per every recording.
```julia
block_ensemble = BlockEnsemble(
    sampling_ensemble,
    block_ranges, # for instance [[1:50, 51:100], [1:20, 21:35, 36:50]]
    ρ,
    num_steps
)
```
Similarly to before, you can call methods above directly on `block_ensemble` and they will be applied recursively to constituent elements of it.

## Structs that translate between different parameter names
----
It is very convenient to have the MCMC chain refer to various parameters independently from how they are internally called inside each diffusion law or inside each observation. This makes it possible to couple parameters that don't share local names. For instance you might have an instance of one law with parameter `ϵ`, say `P1.ϵ` and an instance of another law with parameter `δ`, say `P2.δ` and you'd want to impose a restriction that these parameters are the same. Worry not, **you don't need to redefine the diffusion laws!**. Simply make sure that a relevant dependency structure is passed to an `AllObservations` object, where a parameter, let's call it `shared_ϵδ`, will be defined and will point to appropriate laws. For instance:
```julia
all_obs = ...
add_dependency!(
    all_obs,
    Dict(
        :shared_ϵδ => [(1, :ϵ), (2, :δ)]
    )
)
```
The objects `ParamNamesBlock`, `ParamNamesRecording` and `ParamNamesAllObs` will define relevant translations between parameter names for you. Simply treat `shared_ϵδ` as a parameter name known to the MCMC chain that will implicitly point to relevant parameters `ϵ` and `δ` that you've defined by adding the dependency structure.

### `ParamNamesBlock`

```julia
all_obs = ...
pname = ... # for instance [:shared_ϵδ]

name_struct = ParamNamesBlock(
    block.b,
    pname,
    first(all_obs.param_depend_rev),
    first(all_obs.obs_depend_rev)
)
```

### `ParamNamesRecording`

```julia
all_obs = ...
pname = ... # for instance [:shared_ϵδ]

name_struct = ParamNamesRecording(
    block_collection,
    pname,
    first(all_obs.param_depend_rev),
    first(all_obs.obs_depend_rev)
)
```

### `ParamNamesAllObs`

```julia
all_obs = ...
pname = ... # for instance [:shared_ϵδ]

name_struct = ParamNamesAllObs(block_ensemble, pname, all_obs)
```
