# Using CALiPPSO

In this section we describe in some detail our implementation of the CALiPPSO algorithm. If you're new to [Julia](https://julialang.org/) it might be useful to have [its documentation](https://docs.julialang.org/en/v1/) at hand, specially for understanding some terminology (although we tried to keep it to a minimum).
Besides, several of the functions we define rely on [JuMP](https://jump.dev/)'s functionality, so if you are unfamiliar with this package it might be useful to consult [its documentation](https://jump.dev/JuMP.jl/stable/) (specially the parts on [creating a model](https://jump.dev/JuMP.jl/stable/manual/models/#Create-a-model) and defining/accessing [constraints](https://jump.dev/JuMP.jl/stable/manual/constraints/)).


## Dependencies

First of all, you'll need to have Julia installed and running ([here](https://julialang.org/downloads/) you can find the downloads options and installation instructions). Besides, if you want to use our code without any modification, or you would like to reproduce the results we report in our paper, you'll need to have the [Gurobi optimizer](https://www.gurobi.com/) already properly installed in your system. 
- Note that Gurobi is a licensed solver, but a free license is available for academic use. In any case, see the [choosing another solver](#changing-the-solveroptimizer) for instructions on how to use another (hopefully open-sourced) optimizer.
- Consider also that, of [the different solvers we tested](#testing-different-solvers) , we observed that Gurobi is clearly the most performant and accurate solver, so we do recommend using these scripts in combination with it for more precise results.

Now, the scripts included here depend on basic modules already included in Julia's standard library (like `LinearAlgebra` and `Statistics`), but also on few others that you need to install yourself. They are the following:

1. [JuMP.jl](https://jump.dev/JuMP.jl/stable/): As mentioned in the beginning, this is a basic package. It is used for all the mathematical modelling of the LP problems required.
   - Importantly, JuMP is *not* a solver itself, but the *interface* through which a optimization problem is created.
   - The optimization is done by a *solver* (*e.g.* GLPK, Gurobi, etc), which JuMP calls when a model is optimized (using, *e.g.* the `optimize!` function).
   - Throughout this document, we try to follow the JuMP's terminology: so we call *solver* any package for optimization (*e.g.* GLPK or HiGHS) and *optimizer* the function contained in this packages that is called for actually optimizing a model. Indeed, any solver compatible with JuMP has an `Optimizer` function defined in it.
   
2. [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl): The Julia wrapper for the [Gurobi Solver](https://www.gurobi.com/). This is the solver we tested the most and the default choice in our algorithm. Besides, we use it to produce all the results reported in our paper.
   - **Note**: Gurobi requires manual installation. This means that, before adding `Gurobi.jl` via the Julia's package manager, you need to have Gurobi properly installed in your system. Once you've downloaded Gurobi's latest version, follow the [installation instructions](https://www.gurobi.com/documentation/quickstart.html) of your operative system. 

3. [StaticArrays.jl](https://juliaarrays.github.io/StaticArrays.jl/stable/): To use static (*i.e.* fixed size) vectors for particles positions, contact vectors, etc. This type of arrays considerably speed-up vector operations (such as sum, norm, etc.) and, also importantly, guarantee that all vectors have the same dimensionality. Indeed, the value of ``d`` is inherited as the *type* of a static array (*e.g.* a contact vector is of type `SVector{d,Float64}`) and this helps to write generic code for any dimensionality.


4. [SpecialFunctions.jl](https://specialfunctions.juliamath.org/stable/): To compute the packing fraction in arbitrary dimensions we use the ``\Gamma`` function implemented in this package.

5. [Distributions.jl](https://github.com/JuliaStats/Distributions.jl): Functions to generate random initial configurations rely on some of the functions defined in this package. So it is *only necessary* if the script `random_initial_conditions.jl` is included at some point.


Fortunately, all these packages (except for `Gurobi.jl`, which requires a manual installation of the solver itself) are very easy to install using Julia's package manager. For instance,

```julia
using Pkg
Pkg.add(["JuMP", "StaticArrays", "SpecialFunctions", "Distributions"])
```
will install all of them. And once you've installed Gurobi in your computer, similarly executing `Pkg.add("Gurobi")` will add the Julia wrapper.

Besides, if you want to use another optimizer other than Gurobi, you'll need to install it too. If you want to run the tests described using [different solvers](#changing-the-solveroptimizer) you also need to install: [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl), [GLPK.jl](https://github.com/jump-dev/GLPK.jl), [Clp.jl](https://github.com/jump-dev/Clp.jl), [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl), and [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl). All these packages can also be conveniently be installed with `Pkg.add("<solver of your choice>")` for a single one, or

```julia
Pkg.add(["HiGHS", "GLPK", "Clp", "Hypatia", "COSMO"])
```
for adding all of them simultaneously (also after importing the `Pkg` module).


## Getting help

We tried our best to provide complete *docstrings* of all the functions defined in this package. So most of the information should be available simply by calling the respective documentation (*i.e.* just by typing '?'). For instance, try typing `?produce_jammed_configuration` in the REPL or in a Jupyter notebook for a detailed description of this function. 

We also tried to leave clarifying comments throughout the code, so its functioning is easier to understand.

Yet not clear enough? Found a bug or an issue? Please drop us an email at: XXXXXXX





<!-- However, all the optimizations are done using [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl) (the Julia wrapper of the [Gurobi Optimizer](https://www.gurobi.com/).) Nevertheless, the script is easily customizable to use the optimizer of your choice (see below). -->


