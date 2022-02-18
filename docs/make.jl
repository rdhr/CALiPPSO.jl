const precompile_main_function = false
include("../src/CALiPPSO.jl")
using .CALiPPSO
include("../src/random_initial_conditions.jl")
# using .RandomIC
using Documenter

makedocs(sitename="CALiPPSO",
    # modules = [CALiPPSO, RandomIC],
    authors = "Rafael Diaz Hernandez Rojas",
    pages = [
        "Home" => "index.md",
        "Introduction: theory and few details" => "theory.md",
        "User Guide" => Any["Installation" => "installation.md", "Basic usage" =>"basic_usage.md", "How the main function works" => "mainfunction.md", "Changing the default behaviour" =>"changing_default.md"],
        "Examples" => "tests.md",
        "`Structs` of CALiPPSO" => "types.md",
        "Problem solving and possible issues" => "issues.md",
        "To Do's"=> "todos.md",
        "API Reference" => "api.md"],
        format = Documenter.HTML(prettyurls = false)
    )

deploydocs(
    repo = "github.com/rdhr/CALiPPSO.git",
)