using Documenter
using Platypus

DocMeta.setdocmeta!(Platypus, :DocTestSetup, :(using Platypus); recursive=true)

makedocs(
    sitename = "Platypus.jl",
    modules = [Platypus],
    authors = "Simone Santoni",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://simoneSantoni.github.io/Platypus.jl",
        edit_link = "main",
        assets = ["assets/custom.css"],
    ),
    repo = "https://github.com/simoneSantoni/Platypus.jl/blob/{commit}{path}#{line}",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => [
            "Entities and Data" => "guide/entities.md",
            "Building Feature Pairs" => "guide/pairs.md",
            "Exploratory Analysis" => "guide/eda.md",
            "Estimation" => "guide/estimation.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Estimation" => "api/estimation.md",
            "Visualization" => "api/visualization.md",
        ],
        "References" => "references.md",
    ],
    warnonly = [:missing_docs, :docs_block],
)

deploydocs(
    repo = "github.com/simoneSantoni/Platypus.jl.git",
    devbranch = "main",
    versions = [
        "stable" => "dev",
        "dev" => "dev",
    ],
    push_preview = true,
)
