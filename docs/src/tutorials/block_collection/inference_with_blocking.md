# Inference with blocking with `BlockCollection`
***

```julia
function accept_reject_proposal_param!(bc::BlockCollection, mcmciter, θ, θ°)
    accepted = rand(Exponential(1.0)) > -(fetch_ll°(bc) - fetch_ll(bc))
    accepted && swap_XX!(bc)
    accepted && swap_PP!(bc)
    save_ll!(bc, mcmciter)
    accepted && swap_ll!(bc)
    accepted, copy(accepted ? θ° : θ)
end

function simple_inference_with_blocking(
        AuxLaw, all_obs, dt, AuxLawBlocking, block_layout, _θ;
        ϵ=0.3, ρ=0.5, num_steps=10^4
    )
    # making sure that things are in order...
    _pname = collect(keys(_θ))
    # for simplicity restrict to inference for a single param
    @assert length(_pname) == 1
    pname = first(_pname)
    θ = collect(values(_θ))

    # setting the initial guess θ inside the recording
    OBS.set_parameters!(all_obs, _θ)
    @assert num_recordings(all_obs) == 1
    recording = first(all_obs.recordings)

    # setting up containers
    tts = OBS.setup_time_grids(recording, dt, standard_guid_prop_time_transf)
    sp = SamplingPair(AuxLaw, recording, tts)
    blocks = [
        BlockCollection(sp, block_ranges, ρ, num_steps)
        for block_ranges in block_layout
    ]
    name_struct = [
        ParamNamesRecording(
            bc, _pname, first(all_obs.param_depend_rev),
            first(all_obs.obs_depend_rev)
        ) for bc in blocks
    ]

    paths = []

    θθ = [θ]
    a_h = Bool[]
    crit_change = [fill(true, length(bc.blocks)) for bc in blocks]

    # MCMC
    for i in 1:num_steps
        for bc in blocks
            GP.set_obs!(bc)
            recompute_guiding_term!(bc, Val(:P_only))
            find_W_for_X!(bc)
            loglikhd!(bc)
            draw_proposal_path!(bc)
            accept_reject_proposal_path!(bc, i)

            # progress message
            if i % 100 == 0
                println(
                    "$i. ll=$(ll_of_accepted(bc, i)), acceptance rate: ",
                    "$( accpt_rate(bc, (i-99):i) )"
                )
            end
        end

        θ° = customkernel(θ, ϵ)

        bc = blocks[end]
        set_proposal_law!(bc, θ°, name_struct[end], crit_change[end])
        recompute_guiding_term!(bc, Val(:P°_only))

        accpt, θ = accept_reject_proposal_param!(bc, i, θ, θ°)
        push!(θθ, θ)
        push!(a_h, accpt)

        if i % 100 == 0
            println(
                "$i. updt a-r: ",
                "$(sum(a_h[(i-99):i])/100)."
            )
        end

        # save intermediate path for plotting
        i % 400 == 0 && append!(paths, [deepcopy(sp.u.XX)])
    end
    paths, θθ
end
```
