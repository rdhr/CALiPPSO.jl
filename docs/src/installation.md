# [Installation](@id Installation)

Naturally, the first thing is to have [Julia installed in your system](https://julialang.org/downloads/). Afterwards, you can add [CALiPPSO.jl package](https://github.com/rdhr/CALiPPSO.jl) using the Julia's package manager. However, CALiPPSO is not yet registered in Julia Registries (but we are working to make it so). Therefore, you'll need to add it from the GitHub repo:
```julia
using Pkg
Pkg.add(url="https://github.com/rdhr/CALiPPSO.jl.git")
```

This should also install the required dependencies, of which a very important one is [`GLPK`](https://github.com/jump-dev/GLPK.jl), because it is the default solver of CALiPPSO. In any case, once added, you can simply import CALiPPSO into your current working space (*i.e.* the REPL, a Jupyter notebook, script, etc) as any other package, namely `using CALiPPSO`. 

!!! note "Default solver"
    We chose GLPK as the default solver in our package mainly because its Julia package can be easily installed, and thus is a "safe" dependency. However, we suggest using another, more performant solver for large systems. In [this section](@ref Testing-different-solvers) we explain how to choose the solver used by CALiPPSO.

If you have everything set up, you can jump to the [Basic usage](@ref) section.



## Dependencies

All the dependencies should be downloaded and installed once you add CALiPPSO. However, as mentioned above, its default solver is GNU Linear Programming Kit (a.k.a. [GLPK](https://www.gnu.org/software/glpk/)). This solver has been amply tested in the community, but is not really performant compared to more modern options, at least for the linear optimization problems involved in CALiPPSO. Thus, we suggest to use a different solver if you are going to use CALiPSSO for systems of, say, ``N >2000`` particles.

We have mostly tested our algorithm with the [Gurobi optimizer](https://www.gurobi.com/), naturally using it through its Julia package: [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl).

!!! note "Installing Gurobi"
    Gurobi requires manual installation. This means that, before adding `Gurobi.jl` via the Julia's package manager, you need to have Gurobi properly installed in your system. Once you've downloaded Gurobi's latest version, follow the [installation instructions](https://www.gurobi.com/documentation/quickstart.html) of your operative system. 

    - Note that Gurobi is a licensed solver, but a free license is available for academic use. 
    - Consider also that, of [the different solvers we tested](@ref Testing-different-solvers), we observed that Gurobi is clearly the most performant and accurate solver, so if have access to a license, we do recommend using `CALiPPSO.jl` in combination with it for more precise results.

In any case, see the [choosing another solver](@ref changing_the_solver) for instructions on how to use another optimizer. Our implementation already includes working code for these other solvers: [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl), [Clp.jl](https://github.com/jump-dev/Clp.jl), [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl), and [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl). However, you'll need to install them separately. Fortunately, any (or all) of them are readily available from Julia's package manager interface, *i.e.* `using Pkg; Pkg.add("<solver of your choice>")` for a single one, or

```julia
using Pkg
Pkg.add(["HiGHS", "GLPK", "Clp", "Hypatia", "COSMO"])
```
for adding all of them simultaneously. This latter option might be useful if you want to run the tests described using [different solvers](@ref Testing-different-solvers).

!!! warning "Issue with COSMO and Hypatia"
    In [our tests] (@ref changing_the_solver) we observed that while our code runs without any errors using both `COSMO` and `Hypatia` solvers, the packing obtained maybe either non-isostatic or the contact forces differ significantly from the ones found using the other solvers. We think that this is due to a lower precision in the solutions obtained by these solvers.

    If this is the case, it is likely that this issue can be easily solved by a correct tunning of the solvers' parameters. We were not able to find such good combination, but we did not try much. So if you have any ideas or can help us, please [contact us](@ref Getting-help)
