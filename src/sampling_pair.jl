#===============================================================================

    The main object defined in this file is: SamplingPair.
    It is a pairing of two `SamplingUnit`s and can be used for
    smoothing or inference problems.

===============================================================================#

"""
    struct SamplingPair{TGP,TGPb,TW,TWn,TX}
        u::SamplingUnit{TGP,TGPb,TW,TWn,TX}
        u°::SamplingUnit{TGP,TGPb,TW,TWn,TX}
    end

A pairing of two `SamplingUnit`s and can be used for smoothing or inference
problems.

    SamplingPair(
        aux_laws, recording, x0_prior, tts, args;
        aux_laws_blocking=aux_laws, artificial_noise=1e-11,
        solver_choice_blocking=args
    )

Base constructor.

# Arguments
---
- `aux_laws`:
- `recording`:
- `x0_prior`:
- `tts`:
- `args`:
- `aux_laws_blocking`:
- `artificial_noise`:
- `solver_choice_blocking`:
"""
struct SamplingPair{TGP,TGPb,TW,TWn,TX}
    u::SamplingUnit{TGP,TGPb,TW,TWn,TX}
    u°::SamplingUnit{TGP,TGPb,TW,TWn,TX}

    function SamplingPair(
            aux_laws, recording, x0_prior, tts, args;
            aux_laws_blocking=aux_laws, artificial_noise=1e-11,
            solver_choice_blocking=args
        )
        u = SamplingUnit(
            aux_laws, recording, x0_prior, tts, args;
            aux_laws_blocking = aux_laws_blocking,
            artificial_noise = artificial_noise,
            solver_choice_blocking = solver_choice_blocking
        )
        u° = deepcopy(u)
        TGP,TGPb,TW,TWn,TX = all_eltypes(u)
        new{TGP,TGPb,TW,TWn,TX}(u, u°)
    end
end
