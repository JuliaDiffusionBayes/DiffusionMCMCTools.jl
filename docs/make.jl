using DiffusionMCMCTools
using Documenter

makedocs(;
    modules=[DiffusionMCMCTools],
    authors="mmider <marcin.mider@gmail.com> and contributors",
    repo="https://github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl/blob/{commit}{path}#L{line}",
    sitename="DiffusionMCMCTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaDiffusionBayes.github.io/DiffusionMCMCTools.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaDiffusionBayes/DiffusionMCMCTools.jl",
)
