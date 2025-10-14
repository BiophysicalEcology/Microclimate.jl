using Documenter, DocumenterVitepress
using Microclimate

makedocs(;
    modules=[Microclimate],
    authors="Rafael Schouten, Michael Kearney",
    sitename="Microclimate.jl",
    clean=true,
    doctest=true,
    checkdocs=:all,
    format=DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/Microclimate.jl",
        devbranch = "main",
        devurl = "dev";
    ),
)

deploydocs(;
    repo="github.com/BiophysicalEcology/Microclimate.jl",
    branch="gh-pages",
    devbranch="main",
    push_preview = true
)
