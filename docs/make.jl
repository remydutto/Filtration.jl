using Documenter, DocumenterCitations
push!(LOAD_PATH,"../src/")
using Filtration

# to add docstrings from external packages
Module = Filtration
isnothing(DocMeta.getdocmeta(Module, :DocTestSetup)) &&
        DocMeta.setdocmeta!(Module, :DocTestSetup, :(using $Module); recursive=true)

repo_url = "github.com/remydutto/Filtration.jl"

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

makedocs(;
    remotes=nothing,
    warnonly=:cross_references,
    sitename="Filtration",
    format=Documenter.HTML(;
        repolink="https://" * repo_url,
        prettyurls=false,
        size_threshold_ignore=["index.md"],
    ),
    pages=["Introduction" => "index.md",
        "API" => "api-filtration.md",
        "Developpers" => "dev-filtration.md",
        "References" => "references.md"
        ],
    plugins=[bib],
)

deploydocs(; repo=repo_url * ".git", devbranch="main")
