using GeothermalWells
using Documenter

DocMeta.setdocmeta!(GeothermalWells, :DocTestSetup, :(using GeothermalWells); recursive=true)

makedocs(;
    modules=[GeothermalWells],
    authors="Collin Wittenstein <collin.wittenstein@gmail.com>",
    sitename="GeothermalWells.jl",
    format=Documenter.HTML(;
        canonical="https://cwittens.github.io/GeothermalWells.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cwittens/GeothermalWells.jl",
    devbranch="master",
)
