const precompile_main_function = false
include("../src/iLP-for-jamming.jl")
using .CALiPPSO
using Documenter

makedocs(sitename="CALiPPSO",
    pages = [
        "Home" => "index.md",
        "Introduction: theory and few details" => "theory.md",
        "Using CALiPPSO" => "usage.md",
        "Examples" => "tests.md",
        "`Structs` of CALiPPSO" => "types.md",
        "Possible problems and other issues" => "issues.md",
        "API Reference" => "api.md",
        "To Do's"=> "todos.md"],
        format = Documenter.HTML(prettyurls = false)
    )