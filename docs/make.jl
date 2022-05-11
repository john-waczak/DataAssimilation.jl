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
        mathengine = KaTeX(Dict(
            :macros => Dict(
                "\\cov" => "\\text{cov}",
                "\\R"   => "\\mathbb{R}",
                "\\E"   => "\\mathbb{E}",
            )
        ))
    ),
    pages=[
        "Home" => "index.md",
        "Data Assimilation" => [
            "Data Assimilation Oveview" => "dataAssimilation/dataassim.md",
            "Kalman Filter" => "dataAssimilation/kalman.md",
            "Extended Kalman Filter" => "dataAssimilation/ekf.md",
            "3d var" => "dataAssimilation/3dvar.md",
            "4d var" => "dataAssimilation/4dvar.md",
        ],
        "Function docs" => "fdocs.md",
    ],
)

deploydocs(;
    repo="github.com/john-waczak/DataAssimilation.jl",
    devbranch="main",
)
