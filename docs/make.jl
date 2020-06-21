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
    ],
)

deploydocs(;
    repo="github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl",
)
