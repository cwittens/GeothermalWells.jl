using GeothermalWells
using Documenter

DocMeta.setdocmeta!(GeothermalWells, :DocTestSetup, :(using GeothermalWells); recursive=true)

makedocs(;
    modules=[GeothermalWells],
    authors="Collin Wittenstein <collin.wittenstein@gmail.com>",
    sitename="GeothermalWells.jl",
    format=Documenter.HTML(;
        canonical="https://cwittens.github.io/GeothermalWells.jl",
        edit_link="main",
        assets=String[],
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Methodology" => "methodology.md",
        "Getting Started" => "getting_started.md",
        "API Reference" => "reference.md",
        "License" => "license.md",
    ],
    warnonly=[:missing_docs],  # Don't fail on undocumented internal functions
)

deploydocs(;
    repo="github.com/cwittens/GeothermalWells.jl",
    devbranch="main",
)
