using GeothermalWells
using Documenter
# Optional: generate a nice changelog page from NEWS.md
# Requires the `Changelog` package in the `docs` project environment
try
    using Changelog
catch
    @warn "Changelog.jl not available in docs environment; skip changelog generation"
end

DocMeta.setdocmeta!(GeothermalWells, :DocTestSetup, :(using GeothermalWells); recursive=true)

# Create changelog
Changelog.generate(Changelog.Documenter(),                           # output type
                   joinpath(@__DIR__, "..", "NEWS.md"),              # input file
                   joinpath(@__DIR__, "src", "changelog_tmp.md");    # output file
                    repo = "cwittens/GeothermalWells.jl",
                    branch = "main")
# Fix edit URL of changelog
open(joinpath(@__DIR__, "src", "changelog.md"), "w") do io
    for line in eachline(joinpath(@__DIR__, "src", "changelog_tmp.md"))
        if startswith(line, "EditURL")
            line = "EditURL = \"https://github.com/cwittens/GeothermalWells.jl/blob/main/NEWS.md\""
        end
        println(io, line)
    end
end

makedocs(;
    modules=[GeothermalWells],
    authors="Collin Wittenstein <collin.wittenstein@gmail.com>",
    sitename="GeothermalWells.jl",
    format=Documenter.HTML(;
        canonical="https://cwittens.github.io/GeothermalWells.jl",
        edit_link="main",
        assets=String["assets/favicon.ico"],
        prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Methodology" => "methodology.md",
        "Getting Started" => "getting_started.md",
        "Changelog" => "changelog.md",
        "API Reference" => "reference.md",
        "License" => "license.md",
    ],
    warnonly=[:missing_docs],  # Don't fail on undocumented internal functions
)

deploydocs(;
    repo="github.com/cwittens/GeothermalWells.jl",
    devbranch="main",
    push_preview=true,
)
