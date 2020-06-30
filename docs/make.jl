using DiffusionMCMCTools
using Documenter

makedocs(;
    modules=[DiffusionMCMCTools],
    authors="Sebastiano Grazzi, Frank van der Meulen, Marcin Mider, Moritz Schauer",
    repo="https://github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl/blob/{commit}{path}#L{line}",
    sitename="DiffusionMCMCTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaDiffusionBayes.github.io/DiffusionMCMCTools.jl",
        mathengine = Documenter.MathJax(
            Dict(
                :TeX => Dict(
                    :equationNumbers => Dict(
                        :autoNumber => "AMS"
                    ),
                    :Macros => Dict(
                        :dd => "{\\textrm d}",
                        :RR => "\\mathbb{R}",
                        :wt => ["\\widetilde{#1}", 1]
                    ),
                )
            )
        ),
        collapselevel = 1,
    ),
    pages=[
        "Home" => "index.md",
        "Get started" => joinpath("get_started", "overview.md"),
        "User manual" => Any[
            "SamplingUnit" => joinpath("manual", "sampling_unit.md"),
            "SamplingPair" => joinpath("manual", "sampling_pair.md"),
            "SamplingEnsemble" => joinpath("manual", "sampling_ensemble.md"),
            "Block" => joinpath("manual", "block.md"),
            "BiBlock" => joinpath("manual", "biblock.md"),
            "BlockCollection" => joinpath("manual", "block_collection.md"),
            "BlockEnsemble" => joinpath("manual", "block_ensemble.md"),
            "Containers for parameter names" => joinpath("manual", "param_names_collections.md"),
        ],
        "How to..." => Any[
            "..." => joinpath("how_to_guides", "under_construction.md")
        ],
        "Tutorials" => Any[
            "⚠ Preamble ⚠" => joinpath("tutorials", "preamble.md"),
            "BiBlocks" => Any[
                "Smoothing" => joinpath("tutorials", "biblock", "smoothing.md"),
                "Smoothing with blocking" => joinpath("tutorials", "biblock", "smoothing_with_blocking.md"),
                "Inference" => joinpath("tutorials", "biblock", "inference.md"),
                "Inference with blocking" => joinpath("tutorials", "biblock", "inference_with_blocking.md"),
            ],
            "Containers for parameter names" => Any[
                "Inference with BiBlocks" => joinpath("tutorials", "pnames", "inference_with_biblock.md"),
            ],
            "BlockCollection" => Any[
                "Inference" => joinpath("tutorials", "block_collection", "inference.md"),
                "Inference with blocking" => joinpath("tutorials", "block_collection", "inference_with_blocking.md"),
            ],
            "BlockEnsemble" => Any[
                "Inference" => joinpath("tutorials", "block_ensemble", "inference.md"),
            ],
        ],
        "Index" => "module_index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl",
)
