using BasicIdealGas
using Documenter

DocMeta.setdocmeta!(BasicIdealGas, :DocTestSetup, :(using BasicIdealGas); recursive=true)

makedocs(;
    modules=[BasicIdealGas],
    authors="C. Naaktgeboren",
    sitename="BasicIdealGas.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
