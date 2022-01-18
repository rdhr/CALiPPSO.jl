const precompile_main_function = false
include("../src/iLP-for-jamming.jl")
using .CALiPPSO
using Documenter

makedocs(sitename="CALiPPSO",
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Usage" => "usage.md",
        "Tests" => "tests.md",
        "`Structs` of CALiPPSO" => "types.md",
        "Possible problems and other issues" => "issues.md",
        "To Do's"=> "todos.md"],
        format = Documenter.HTML(prettyurls = false)
    )