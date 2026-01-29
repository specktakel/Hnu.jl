push!(LOAD_PATH, "../src/")

using Hnu
using Documenter

makedocs(
    sitename = "Hnu.jl",
    modules = [Hnu],
    pages = [
        "Home" => "index.md"
    ]
)

#deploydocs(;
#    repo="github.com/specktakel/Hnu.jl",
#)