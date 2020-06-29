"""
    struct ParamNamesUnit{N1,N2}
        var::NTuple{N1,Symbol}
        var_aux::Vector{Tuple{Vararg{Symbol,_N} where _N}}
        updt::NTuple{N2,Pair{Int64,Symbol}}
        updt_aux::Vector{Tuple{Vararg{Pair{Int64,Symbol},_N} where _N}}
        updt_obs::Vector{Tuple{Vararg{Pair{Int64,Int64},_N} where _N}}
    end

Smallest unit storing information about the names of the parameters that are
supposed to be updated at a given step of the MCMC algorithm. It is pertinent
to a single block of a single recording only, and in particular, only to a
single collection of laws, out of:
- `PP`, i.e. those relevant to a given block
- `P_last`, i.e. the artificial law of the block associated with the last obs
- `P_excl`, i.e. the last, missed out `PP`, missed out due to use of `P_last`
- `PPb`, i.e. the collection of all remaining artificial laws that are not
         `P_last`
It stores lists of parameters that are supposed to be updated at various stages
of the call to a function `set_proposal_law!`, which sets the proposal parameter
`θ°` inside the proposal laws.

# Fields
- `var`: names of `P_target` that are checked and—if need be—equalized between
         proposal and accepted laws
- `var_aux`: names of `P_aux` that are checked and—if need be—equalized between
             proposal and accepted laws
- `updt`: list of params of `P_target` that are being updated at a given MCMC
          step; listed as pairs: `idx-of-θ-to-relevant-value` =>
          `name-of-param-inside-P_target-struct`
- `updt_aux`: as above, but for `P_aux`
- `updt_obs`: list of params of `obs` (i.e. terminal observation for each guid
              prop) that are being updated at a given MCMC step; listed as
              pairs: `idx-of-θ-to-relevant-value` => `idx-of-obs-inside-obs.θ`


    ParamNamesUnit(
        PP::AbstractArray{<:GuidProp}, θnames::Vector{Symbol}, pdep, odeps
    )

Base constructor.

# Arguments
- `PP`:
- `θnames`:
- `pdep`:
- `odeps`:
"""
struct ParamNamesUnit{N1,N2}
    var::NTuple{N1,Symbol}
    var_aux::Vector{Tuple{Vararg{Symbol,_N} where _N}}
    updt::NTuple{N2,Pair{Int64,Symbol}}
    updt_aux::Vector{Tuple{Vararg{Pair{Int64,Symbol},_N} where _N}}
    updt_obs::Vector{Tuple{Vararg{Pair{Int64,Int64},_N} where _N}}

    function ParamNamesUnit(
            PP::AbstractArray{<:GuidProp},
            θnames::Vector{Symbol},
            pdep,
            odeps,
        )
        updt = find_θ_names_for_MCMC_update(θnames, pdep)
        updt_aux = find_θ_aux_names_for_MCMC_update(updt, PP)
        updt_obs = find_θ_obs_idx_for_MCMC_update(θnames, odeps)

        var = find_var_names_not_in_MCMC_update(updt, PP)
        var_aux = find_var_aux_names_not_in_MCMC_update(updt_aux, PP)

        new{length(var), length(updt)}(var, var_aux, updt, updt_aux, updt_obs)
    end

    RecordingParamNames() = new{0,0}()
end

tuple_lengths(::ParamNamesUnit{N1,N2}) where {N1,N2} = (N1, N2)

"""
    find_θ_names_for_MCMC_update(θnames, pdep)

Parse through `θnames` i.e. a list of all parameter names that are relevant for
a given MCMC update and pick out only those that are relevant for a given
diffusion law. All information about the relevant parameters of the law should
be stored inside `pdep`. Return a list of relevant parameter names in a format:
`idx-of-θ-to-relevant-value` => `name-of-param-inside-P_target-struct`.
"""
find_θ_names_for_MCMC_update(θnames, pdep) = Tuple(
    map(
        gl_p_n->(
            findfirst(
                n->(n==gl_p_n[1]),
                θnames
            ) => gl_p_n[2]
        ), # substitute a global param name (gl_p_n[1]) with its index in θnames
        filter(
            p->(p[1] in θnames),
            pdep
        ) # remove all parameter names not featured in θnames
    )
)

"""
    find_θ_aux_names_for_MCMC_update(θ_names_for_MCMC_update, PP)

Parse through `PP` and for each `GuidProp` struct enter the auxiliary law
`P_aux` and pick out names from `θ_names_for_MCMC_update` that are relevant to
this `P_aux` law.
"""
find_θ_aux_names_for_MCMC_update(θ_names_for_MCMC_update, PP) = map(PP) do P
    filter(
        p->(p[2] in DD.var_parameter_names(P.P_aux)),
        θ_names_for_MCMC_update
    )
end

"""
    find_θ_obs_idx_for_MCMC_update(θnames, odeps)

Parse through `θnames` i.e. a list of all parameter names that are relevant for
a given MCMC update and for each observation pick out only those that are
relevant for it. All information about the relevant parameters of the
observations should be stored inside `odeps`. Return a list of relevant
parameter names in a format: `idx-of-θ-to-relevant-value` =>
`idx-of-obs-inside-obs.θ`.
"""
find_θ_obs_idx_for_MCMC_update(θnames, odeps) = map(
    odep->Tuple(
        map(
            gl_p_n->(
                findfirst(
                    n->(n==gl_p_n[1]),
                    θnames
                ) => gl_p_n[2]
            ), # substitute a global parameter name with its index in θnames
            filter(
                p->(p[1] in θnames),
                odep
            ) # remove all observation params not featured in θnames
        )
    ),
    odeps
)

"""
    find_var_names_not_in_MCMC_update(composite_names_in_updt, PP)

Parse through all variable names of `P_target` and remove all those names that
already feature in `composite_names_in_updt`.
"""
function find_var_names_not_in_MCMC_update(composite_names_in_updt, PP)
    names_in_updt = getindex.(composite_names_in_updt, 2)
    length(PP) > 0 || return tuple()
    filter(pn->!(pn in names_in_updt), DD.var_parameter_names(PP))
end

"""
    find_var_aux_names_not_in_MCMC_update(composite_names_in_updt, PP)

Parse through all `PP` and for each element look through variable names of
`P_aux` and remove all those names that already feature in
`composite_names_in_updt`.
"""
function find_var_aux_names_not_in_MCMC_update(composite_names_in_updt, PP)
    names_in_updt = map(n->getindex.(n, 2), composite_names_in_updt)
    #println(names_in_updt)
    [
        filter(
            pn->!(pn in n),
            DD.var_parameter_names(P.P_aux)
        ) for (n, P) in zip(names_in_updt, PP)
    ]
end

"""
    struct ParamNamesBlock{N1,N2}
        PP::ParamNamesUnit{N1,N2}
        P_last::ParamNamesUnit{N1,N2}
        P_excl::ParamNamesUnit{N1,N2}
        Pb_excl::ParamNamesUnit{N1,N2}
    end

Stores information relevant for a `Block` about the names of the parameters that
are supposed to be updated at a given step of the MCMC algorithm. It stores
lists of parameters that are supposed to be updated at various stages of the
call to a function `set_proposal_law!`, which sets the proposal parameter `θ°`
inside the proposal laws.

# Fields
- `PP`, i.e. those relevant to a given block
- `P_last`, i.e. the artificial law of the block associated with the last obs
- `P_excl`, i.e. the last, missed out `PP`, missed out due to use of `P_last`
- `Pb_excl`, i.e. the collection of all remaining artificial laws that are not
             `P_last`

    ParamNamesBlock(b::Block, θnames, pdep, odeps)

Base constructor.

# Arguments
- `b::Block`:
- `θnames`:
- `pdep`:
- `odeps`:
"""
struct ParamNamesBlock{N1,N2}
    PP::ParamNamesUnit{N1,N2}
    P_last::ParamNamesUnit{N1,N2}
    P_excl::ParamNamesUnit{N1,N2}
    Pb_excl::ParamNamesUnit{N1,N2}

    function ParamNamesBlock(b::Block, θnames, pdep, odeps)
        i1, i2 = _idx_split(b)

        PP = ParamNamesUnit(b.PP, θnames, pdep, odeps[i1])
        P_last = ParamNamesUnit(b.P_last, θnames, pdep, [tuple() for _ in i2])
        P_excl = ParamNamesUnit(b.P_excl, θnames, pdep, odeps[i2])
        Pb_excl = ParamNamesUnit(b.Pb_excl, θnames, pdep, [tuple() for _ in i1])
        N1, N2 = tuple_lengths(PP)
        new{N1,N2}(PP, P_last, P_excl, Pb_excl)
    end
end

function _idx_split(b::Block{false})
    bidx = first(b.PP.indices) # PP is a subarray
    bidx, [bidx[end]+1]
end

function _idx_split(b::Block{true})
    bidx = first(b.PP.indices) # PP is a subarray
    bidx, Int64[]
end

tuple_lengths(::ParamNamesBlock{N1,N2}) where {N1,N2} = (N1, N2)

"""
    struct ParamNamesRecording{N1,N2}
        blocks::Vector{ParamNamesBlock{N1,N2}}
    end

Stores information relevant for an entire, single recording about the names of
the parameters that are supposed to be updated at a given step of the MCMC
algorithm. It stores lists of parameters that are supposed to be updated at
various stages of the call to a function `set_proposal_law!`, which sets the
proposal parameter `θ°` inside the proposal laws.

    ParamNamesRecording(rb::BlockCollection, θnames, pdep, odeps)

Base constructor.
"""
struct ParamNamesRecording{N1,N2}
    blocks::Vector{ParamNamesBlock{N1,N2}}

    function ParamNamesRecording(bc::BlockCollection, θnames, pdep, odeps)
        blocks = map(bb->ParamNamesBlock(bb.b, θnames, pdep, odeps), bc.blocks)
        N1, N2 = tuple_lengths(first(blocks))
        new{N1,N2}(blocks)
    end
end

"""
    struct ParamNamesAllObs{T}
        recordings::Vector{T}
    end

Stores information relevant for multiple recordings about the names of the
parameters that are supposed to be updated at a given step of the MCMC
algorithm. It stores lists of parameters that are supposed to be updated at
various stages of the call to a function `set_proposal_law!`, which sets the
proposal parameter `θ°` inside the proposal laws.

    ParamNamesAllObs(be::BlockEnsemble, θnames, all_obs)

Base constructor.
"""
struct ParamNamesAllObs{T}
    recordings::Vector{T}

    function ParamNamesAllObs(be::BlockEnsemble, θnames, all_obs)
        recordings = map(1:length(be.recordings)) do i
            ParamNamesRecording(
                be.recordings[i],
                θnames,
                all_obs.param_depend_rev[i],
                all_obs.obs_depend_rev[i]
            )
        end
        new{eltype(recordings)}(recordings)
    end
end
