using MasterChemicalMechanism
using Documenter

DocMeta.setdocmeta!(MasterChemicalMechanism, :DocTestSetup, :(using MasterChemicalMechanism); recursive=true)

makedocs(;
    modules=[MasterChemicalMechanism],
    authors="John Waczak <john.louis.waczak@gmail.com>",
    repo="https://github.com/john-waczak/MasterChemicalMechanism.jl/blob/{commit}{path}#{line}",
    sitename="MasterChemicalMechanism.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://john-waczak.github.io/MasterChemicalMechanism.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/john-waczak/MasterChemicalMechanism.jl",
    devbranch="main",
)
