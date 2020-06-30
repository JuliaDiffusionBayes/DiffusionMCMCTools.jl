# Smoothing with blocking with `BiBlock`s
***

!!! note "important"
    Make sure that you read the [preamble](@id tutorials_start) before you start reading this tutorial.

A first level of complication to regular smoothing algorithm is an addition of a blocking scheme. Below we will define blocking "by hand" by operating directly on `BiBlock`s.

## The algorithm
---
```julia
function simple_smoothing_with_blocking(
        AuxLaw, recording, dt, AuxLawBlocking, block_layout;
        ρ=0.5, num_steps=10^4
    )
    tts = OBS.setup_time_grids(recording, dt, standard_guid_prop_time_transf)
    # this object is not changing, it still has all relevant containers
    sp = SamplingPair(AuxLaw, recording, tts)
    # and this has pointers to containers and facilitates actual sampling,
    # it's a bit more complicated than before and contains multiple sets of blocks
    blocks = [
        [
            BiBlock(sp, br, ρ, i==length(block_ranges), num_steps)
            for (i,br) in enumerate(block_ranges)
        ] for block_ranges in block_layout
    ]

    # we will again aggregate sampled paths here
    paths = []

    # MCMC
    for i in 1:num_steps
        # iterate through all sets of blocks
        for B in blocks
            # freeze terminal points of blocks to be artificial observations
            GP.set_obs!.(B)
            # recompute the guiding term only on the "accepted" laws `bb.b.PP`
            (bb->recompute_guiding_term!(bb.b)).(B)
            # recompute the Wiener path
            find_W_for_X!.(B)
            # re-evaluate the log-likelihood
            loglikhd!.(B)
            # impute a path
            draw_proposal_path!.(B)
            # Metropolis–Hastings accept/reject step
            accept_reject_proposal_path!.(B, i)

            # progress message
            if i % 100 == 0
                println(
                    "$i. ll=$(ll_of_accepted.(B, i)), acceptance rate: ",
                    "$( map(bb->accpt_rate(bb, (i-99):i), B) )"
                )
            end
        end

        # save intermediate path for plotting
        i % 400 == 0 && append!(paths, [deepcopy(sp.u.XX)])
    end
    paths
end
```

Let's run the algorithm on two sets of blocks based on 3 artificial points.

```julia
paths = simple_smoothing_with_blocking(
    FitzHughNagumoAux, recording, 0.001, FitzHughNagumoAux,
    [[1:25,26:75,76:100],[1:50, 51:100]];
    ρ=0.96, num_steps=10^4
)
```

## Results
---
```julia
plot_imputed_trajectories(paths)
```
![smoothing_with_blocking](../../assets/tutorials/biblock/smoothing_with_blocking.png)
