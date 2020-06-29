# BlockCollection
*****
```@docs
DiffusionMCMCTools.BlockCollection
```

## Imputation
```@docs
DiffusionMCMCTools.draw_proposal_path!(bc::DiffusionMCMCTools.BlockCollection)
```

## Accept/reject decision
```@docs
DiffusionMCMCTools.accept_reject_proposal_path!(bb::DiffusionMCMCTools.BlockCollection, mcmciter)
```

## Swaps
```@docs
DiffusionMCMCTools.swap_paths!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.swap_XX!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.swap_WW!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.swap_PP!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.swap_ll!(bc::DiffusionMCMCTools.BlockCollection)
```

#### Setting up a block

```@docs
DiffusionMCMCTools.set_obs!(bc::DiffusionMCMCTools.BlockCollection)
```

## utility
```@docs
DiffusionMCMCTools.loglikhd!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.loglikhd°!(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.fetch_ll(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.fetch_ll°(bc::DiffusionMCMCTools.BlockCollection)
DiffusionMCMCTools.save_ll!(bc::DiffusionMCMCTools.BlockCollection, i::Int)
DiffusionMCMCTools.ll_of_accepted(bb::DiffusionMCMCTools.BlockCollection, i)
DiffusionMCMCTools.accpt_rate(bb::DiffusionMCMCTools.BlockCollection, range)
```

## Setting parameters
```@docs
DiffusionMCMCTools.set_proposal_law!(
    bc::DiffusionMCMCTools.BlockCollection,
    θ°,
    pnames,
    critical_change=GP.is_critical_update(bc, pnames);
    skip=0
)
```

# Example

## The algorithm

```julia
all_obs = AllObservations()
add_recording!(all_obs, build_recording(P, data, 0.0, KnownStartingPt(y1)))
DD.var_parameter_names(::FitzHughNagumo) = (:γ,)
all_obs, _ = initialize(all_obs)

customkernel(θ, scale=0.1) = θ .+ 2.0*scale*(rand()-0.5)


function accept_reject_proposal_param!(bc, mcmciter, θ, θ°)
    accepted = rand(Exponential(1.0)) > -(fetch_ll°(bc)-fetch_ll(bc))
    accepted && swap_XX!(bc)
    accepted && swap_PP!(bc)
    save_ll!(bc, mcmciter)
    accepted && swap_ll!(bc)
    accepted, copy(accepted ? θ° : θ)
end

function simple_inference(AuxLaw, all_obs, dt, _θ; ϵ=0.3, ρ=0.5, num_steps=10^4)
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
    num_obs = length(recording.obs)
    tts = OBS.setup_time_grids(recording, dt, standard_guid_prop_time_transf)
    sp = SamplingPair(AuxLaw, recording, tts)
    bc = BlockCollection(sp, [1:num_obs], ρ, num_steps)
    name_struct = ParamNamesRecording(
        bc, _pname, first(all_obs.param_depend_rev),
        first(all_obs.obs_depend_rev)
    )

    loglikhd!(bc)
    paths = []

    θθ = [θ]
    a_h = Bool[]
    crit_change = [true]

    for i in 1:num_steps
        draw_proposal_path!(bc)
        accept_reject_proposal_path!(bc, i)

        θ° = customkernel(θ, ϵ)
        set_proposal_law!(bc, θ°, name_struct, crit_change)

        accpt, θ = accept_reject_proposal_param!(bc, i, θ, θ°)
        push!(θθ, θ)
        push!(a_h, accpt)

        # progress message
        if i % 100 == 0
            println(
                "$i. ll=$(ll_of_accepted(bc, i)), ",
                "imp a-r: ",
                " $(accpt_rate(bc, (i-99):i)), ",
                "updt a-r: ",
                "$(sum(a_h[(i-99):i])/100)."
            )
        end

        # save intermediate path for plotting
        i % 400 == 0 && append!(paths, [deepcopy(sp.u.XX)])
    end
    paths, θθ
end
```

## Result

```julia
using OrderedCollections

θ = OrderedDict(:REC1_γ=>1.5)

DD.var_parameter_names(::FitzHughNagumo) = (:γ,)
DD.var_parameter_names(::FitzHughNagumoAux) = (:γ,)

@load_diffusion FitzHughNagumoAux
paths, θθ = simple_inference(
    FitzHughNagumoAux, all_obs, 0.001, θ; ϵ=0.3, ρ=0.96, num_steps=10^4
)
```

# Example of inference with blocking

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
