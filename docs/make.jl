include("../src/iLP-for-jamming.jl")
using .CALiPPSO
using Documenter

makedocs(sitename="CALiPPSO",
pages = [
    "Home" => "index.md",
    "Theory" => "theory.md"]
    )