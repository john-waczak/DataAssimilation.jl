using DataAssimilation
using Documenter


# NOTE: geostats.jl has good example documentation and structure in their make.jl


DocMeta.setdocmeta!(DataAssimilation, :DocTestSetup, :(using DataAssimilation); recursive=true)

makedocs(;
    modules=[DataAssimilation],
    authors="John Waczak <john.louis.waczak@gmail.com>",
    repo="https://github.com/john-waczak/DataAssimilation.jl/blob/{commit}{path}#{line}",
    sitename="DataAssimilation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://john-waczak.github.io/DataAssimilation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Extended Kalman Filter" => "ekf.md",
        "Function docs" => "fdocs.md",
    ],
)

deploydocs(;
    repo="github.com/john-waczak/DataAssimilation.jl",
    devbranch="main",
)
