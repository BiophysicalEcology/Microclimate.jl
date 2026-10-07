using Documenter, DocumenterVitepress
using Microclimate
using CairoMakie
using Unitful

# Don't output huge svgs for Makie plots
CairoMakie.activate!(type = "png")

makedocs(;
    modules=[Microclimate],
    authors="Michael Kearney, Rafael Schouten et al.",
    sitename="Microclimate.jl",
    clean=true,
    doctest=false,
    checkdocs=:exports,
    format=DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/Microclimate.jl",
        devbranch = "main",
        devurl = "dev";
    ),
    source = "src",
    build = "build",
    warnonly = true,
    pages=[
        "Home" => "index.md",
        "Get started" => "get_started.md",
        "Manual" => [
            joinpath("manual", "introduction.md"),
            joinpath("manual", "configuring_the_model.md"),
            joinpath("manual", "diel_time_handling.md"),
            joinpath("manual", "boundary_layer.md"),
            joinpath("manual", "radiation.md"),
            joinpath("manual", "soil_thermal_properties.md"),
            joinpath("manual", "phase_transition.md"),
            joinpath("manual", "soil_moisture.md"),
            joinpath("manual", "snow.md"),
            joinpath("manual", "canopy.md"),
            joinpath("manual", "evaporation_condensation.md"),
            joinpath("manual", "references.md"),
        ],
        "Tutorials" => [
            joinpath("tutorials", "monthly_workflow.md"),
            joinpath("tutorials", "daily_hourly_workflow.md"),
            joinpath("tutorials", "terrain_workflow.md"),
            joinpath("tutorials", "snow_demo.md"),
            joinpath("tutorials", "canopy_workflow.md"),
            joinpath("tutorials", "soil_hydraulics_algorithms.md"),
            joinpath("tutorials", "initial_conditions_and_burnin.md"),
        ],
        "API" => "api.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo="github.com/BiophysicalEcology/Microclimate.jl",
    branch="gh-pages",
    target = joinpath(@__DIR__, "build"),
    devbranch="main",
    push_preview=true,
)
