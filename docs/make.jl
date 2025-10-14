using Documenter, DocumenterVitepress
using Microclimate

makedocs(;
    modules=[Microclimate],
    authors="Rafael Schouten, Michael Kearney",
    repo="https://github.com/BiophysicalEcology/Microclimate.jl/blob/{commit}{path}#{line}",
    sitename="Microclimate.jl",
    format=DocumenterVitepress.MarkdownVitepress(),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
    draft=false,
)

deploydocs(;
    repo="github.com/BiophysicalEcology/Microclimate.jl",
    devbranch="main",
    push_preview = true
)
