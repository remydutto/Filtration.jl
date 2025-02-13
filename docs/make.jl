using Documenter
push!(LOAD_PATH,"../src/")
using Filtration

repo_url = "github.com/remydutto/Filtration.jl"

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
        ],
)

deploydocs(; repo=repo_url * ".git", devbranch="main")
