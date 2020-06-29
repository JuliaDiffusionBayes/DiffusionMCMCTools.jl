"""
    struct SamplingEnsemble{T}
        recordings::Vector{T}
    end

A collection of `SamplingPair`s, can be used for smoothing or inference
problems.

    SamplingEnsemble(
        aux_laws, recordings, tts, args=tuple();
        aux_laws_blocking=aux_laws, artificial_noise=1e-11,
        solver_choice_blocking=args
    )

Base constructor.
"""
struct SamplingEnsemble{T}
    recordings::Vector{T}

    function SamplingEnsemble(
            aux_laws, recordings, tts, args=tuple();
            aux_laws_blocking=aux_laws, artificial_noise=1e-11,
            solver_choice_blocking=args
        )
        num_rec = length(recordings)
        aux_laws = _vec_me(aux_laws, num_rec)
        args = _vec_me(args, num_rec)
        aux_laws_blocking = _vec_me(aux_laws_blocking, num_rec)
        artificial_noise = _vec_me(artificial_noise, num_rec)
        solver_choice_blocking = _vec_me(solver_choice_blocking, num_rec)
        recs = map(1:num_rec) do i
            SamplingPair(
                aux_laws[i], recordings[i], tts[i], args[i];
                aux_laws_blocking = aux_laws_blocking[i],
                artificial_noise = artificial_noise[i],
                solver_choice_blocking = solver_choice_blocking[i]
            )
        end
        new{eltype(recs)}(recs)
    end
end

_vec_me(val, N) = (typeof(val) <: AbstractArray ? val : fill(val, N))


OBS.num_recordings(se::SamplingEnsemble) = length(se.recordings)
