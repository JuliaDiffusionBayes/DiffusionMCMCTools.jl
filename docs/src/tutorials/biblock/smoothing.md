# Smoothing with `BiBlock`s
****

!!! note "important"
    Make sure that you read the [preamble](@id tutorials_start) before you start reading this tutorial.
  

## The algorithm
----
The smoothing algorithm is very simple. We
- initialize appropriate objects:
    - the time-grids
    - all containers needed for smoothing via `SamplingPair`
    - pointers via `BiBlock` that basically point to all containers (as we do not employ any blocking schemes)
- run MCMC by repeatedly
  - proposing a new path from a proposal diffusion measure
  - accept/reject decisions and updating the chain of trajectories

```julia
"""
This is a simple smoothing function that imputes unobserved parts of the path.
No blocking is used, only the preconditioned Crank-Nicolson scheme with memory
parameter ρ.
"""
function simple_smoothing(AuxLaw, recording, dt; ρ=0.5, num_steps=10^4)
    # set up time-grids
    tts = OBS.setup_time_grids(recording, dt, standard_guid_prop_time_transf)    
    # this object contains containers
    sp = SamplingPair(AuxLaw, recording, tts)
    # and this has pointers to containers and facilitates actual sampling
    bb = BiBlock(sp, 1:length(recording.obs), ρ, true, num_steps)
    # make sure that the log-likelihood is computed for the initialized paths
    # and stored in a local field `ll`
    loglikhd!(bb)
    # we will append imputed paths once in a while to visualize the sampler
    # and see for ourselves how we are doing
    paths = []

    # MCMC
    for i in 1:num_steps
        # impute a path
        draw_proposal_path!(bb)
        # Metropolis–Hastings accept/reject step
        accept_reject_proposal_path!(bb, i)

        # progress message
        if i % 100 == 0
            println(
                "$i. ll=$(ll_of_accepted(bb, i)), acceptance rate: ",
                "$( accpt_rate(bb, (i-99):i) )"
            )
        end

        # save intermediate path for plotting
        i % 400 == 0 && append!(paths, [deepcopy(bb.b.XX)])
    end
    paths
end
```
running the algorithm is very simple
```julia
paths = simple_smoothing(FitzHughNagumoAux, recording, 0.001; ρ=0.96, num_steps=10^4)
```
## Results
----
```julia
plot_imputed_trajectories(paths)
```
![smoothing](../../assets/tutorials/biblock/smoothing.png)
