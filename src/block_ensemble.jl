"""
    struct BlockEnsemble{T}
        recordings::Vector{T}
    end

Gathers all blocks relevant for an entire single recording.

    BlockEnsemble(
        se::SamplingEnsemble,
        ranges,
        ρρ=0.0,
        ll_hist_len=0
    )

Base constructor.
"""
struct BlockEnsemble{T}
    recordings::Vector{T}

    function BlockEnsemble(
            se::SamplingEnsemble,
            ranges,
            ρρ=0.0,
            ll_hist_len=0
        )
        num_rec = num_recordings(se)
        ρρ = _vec_me(ρρ, num_rec)
        ll_hist_len = _vec_me(ll_hist_len, num_rec)
        rec = map(1:num_rec) do i
            BlockCollection(se.recordings[i], ranges[i], ρρ[i], ll_hist_len[i])
        end
        new{eltype(rec)}(rec)
    end
end

OBS.num_recordings(aob::BlockEnsemble) = length(aob.recordings)
