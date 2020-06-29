# Containers for parameter names
****
```@docs
DiffusionMCMCTools.ParamNamesUnit
DiffusionMCMCTools.ParamNamesBlock
DiffusionMCMCTools.ParamNamesRecording
DiffusionMCMCTools.ParamNamesAllObs
```

# Example

## The algorithm
```julia
all_obs = AllObservations()
add_recording!(all_obs, build_recording(P, data, 0.0, KnownStartingPt(y1)))
DD.var_parameter_names(::FitzHughNagumo) = (:γ,)
all_obs, _ = initialize(all_obs)

customkernel(θ, scale=0.1) = θ .+ 2.0*scale*(rand()-0.5)

#=
    NOTE we no longer need to hand-craft the `_build_struct` and
    `simple_name_structure` functions. These will be substituted with a
    struct `ParamNamesRecording`
=#

function accept_reject_proposal_param!(bb, mcmciter, θ, θ°)
    accepted = rand(Exponential(1.0)) > -(bb.b°.ll - bb.b.ll)
    accepted && swap_XX!(bb)
    accepted && swap_PP!(bb)
    save_ll!(bb.b, mcmciter)
    save_ll!(bb.b°, mcmciter)
    accepted && swap_ll!(bb)
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
    bb = BiBlock(sp, 1:num_obs, ρ, true, num_steps)
    name_struct = ParamNamesBlock(
        bb.b, _pname, first(all_obs.param_depend_rev),
        first(all_obs.obs_depend_rev)
    )

    loglikhd!(bb)
    paths = []

    θθ = [θ]
    a_h = Bool[]

    for i in 1:num_steps
        draw_proposal_path!(bb)
        accept_reject_proposal_path!(bb, i)

        θ° = customkernel(θ, ϵ)
        set_proposal_law!(bb, θ°, name_struct, true)

        accpt, θ = accept_reject_proposal_param!(bb, i, θ, θ°)
        push!(θθ, θ)
        push!(a_h, accpt)

        # progress message
        if i % 100 == 0
            println(
                "$i. ll=$(ll_of_accepted(bb, i)), ",
                "imp a-r: ",
                " $(accpt_rate(bb, (i-99):i)), ",
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
