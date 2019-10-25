using Documenter, SpectralElements

makedocs(;
    modules=[SpectralElements],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/pjabardo/SpectralElements.jl/blob/{commit}{path}#L{line}",
    sitename="SpectralElements.jl",
    authors="Paulo Jabardo, IPT",
    assets=String[],
)

deploydocs(;
    repo="github.com/pjabardo/SpectralElements.jl",
)
