# CALiPPSO: A Linear Programming Algorithm for Jamming Hard Spheres

This package is a pure [Julia](https://julialang.org/) implementation of the CALiPPSO algorithm for generating jammed packings of hard spheres in arbitrary dimensions. The algorithm was introduced in [this article](https://arxiv.org/abs/2203.05654) by *Artiaco, Díaz, Parisi, and Ricci-Tersenghi*. As explained there, CALiPPSO consists on a **C**hain of **A**pproximate **Li**near **P**rogramming for **P**acking **S**pherical **O**bjects.

This package is licensed under the MIT license, so please feel free to use/modify/improve this code as better suits you. We only ask you to cite our work if you find it useful.

```
@article{calippso,
  title = {{{CALiPPSO}}: {{A Linear Programming Algorithm}} for {{Jamming Hard Spheres}}},
  author = {Artiaco, Claudia and {D{\'i}az Hern{\'a}ndez Rojas}, Rafael and Parisi, Giorgio and {Ricci-Tersenghi}, Federico},
  year = {2022},
  month = mar,
  journal = {arXiv},
  eprint = {2203.05654},
  eprinttype = {arxiv},
  primaryclass = {cond-mat, physics:physics},
  url = {http://arxiv.org/abs/2203.05654},
  archiveprefix = {arXiv}
}
```

## Documentation

You can read the full documentation of our code [here](https://rdhr.github.io/CALiPPSO.jl/dev/index.html).


## Basic usage

### Installation
Our package is not yet registered in Julia Registries (but we are working to make it so). Therefore, you'll need to add it from this GitHub repo:
```julia
import Pkg
Pkg.add(url="https://github.com/rdhr/CALiPPSO.jl.git")
```

This will also automatically install the required dependencies; the main ones are [JuMP.jl](https://jump.dev/JuMP.jl/stable/), [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) and [GLPK.jl](https://github.com/jump-dev/GLPK.jl). The latter is needed because [`GLPK`](https://www.gnu.org/software/glpk/) is the default solver of CALiPPSO. In any case, once CALiPPSO is added, you can simply import it into your current working space (*i.e.* the REPL, a Jupyter notebook, script, etc) as any other package, namely `using CALiPPSO`. 

Below we show a minimal working example (MWE) and show how to [change the solver](#Changing-the-solver) used by CALiPPSO.


### Minimal example

We tried to make this package as easy to use as possible and, indeed, it consists of *a single* main function: `produce_jammed_configuration!(Xs0, r0, L=1.0)`. We also provide a function to generate a low density random initial condition so that you can use CALiPPSO right away. However, as we explain in [our paper](https://arxiv.org/abs/2203.05654) and in [the relevant part of the documentation](https://rdhr.github.io/CALiPPSO.jl/dev/theory.html#The-initial-conditions), our algorithm works best if the initial condition is already close to its jamming point. Thus, our code is not guaranteed to work with any such low density configurations. 
However, for small systems, even a low density configuration should be suitable for initializing CALiPPSO. So, for instance, to jammed a d=3 system of 512 hard-spheres, here is a MWE

```julia
using CALiPPSO  
precompile_main_function() #optional, but highly recommended. This will produce a colorful output that you can safely ignore
using Random
Random.seed!(123) # optional, but just for reproducibility sake of this MWE
# Choosing the seed of the Julia's RNG determines the random IC produces below with `generate_random_configuration`

const d, N, φ0, L = 3, 512, 0.3, 1.0
r0, Xs0 = generate_random_configuration(d, N, φ0, L) # if L is not passed, it's assumed that the systems is in a box of size 1

packing, info, Γ_vs_t, Smax_vs_t, isostatic_vs_t = produce_jammed_configuration!(Xs0, r0; 
            ℓ0=0.2*L, max_iters=500)
```
Therefore, the main arguments of `produce_jammed_configuration!` are the particles' initial position `Xs0` and their initial radius, `r0`. So far, our implementation of `produce_jammed_configuration!` assumes the system is contained in a *periodic* (hyper-) cube of size `L`. 
The value of `L` is inferred in the following way
1. If `Xs0` is a $d\times N$ matrix specifying the position of each particle (*i.e.* each of the $N$ columns is the $d$-dimensional position vector of a particle). Then `L` should be passed as a third argument to `produce_jammed_configuration!`.
     - If left unspecified, but `Xs0` is of type `Matrix{Float64}`, then it is assumed `L=1.0`.
2. If `Xs0` is of type `Vector{SVector{d, PeriodicNumber{Float64}}}`, the elements of `Xs0` are of [`PeriodicNumber` type](https://rdhr.github.io/CALiPPSO.jl/dev/types.html#The-PeriodicNumber-type), and hen `L` is automatically inferred from them.
     - For instance, this is the case when `Xs0` is generated by calling `generate_random_configuration` as in the example above.

The usage of the keyword arguments (`ℓ0`, `max_iters`, etc.) is explained in the [dedicated section of the documentation](https://rdhr.github.io/CALiPPSO.jl/dev/changing_default.html#kwargs-control), and in the docstring of the main function. Thus, simply try
```julia
?produce_jammed_configuration!
```

### Output

The output of `produce_jammed_configuration!` is the following:
1. A jammed packing (provided convergence was attained) stored as a `MonoPacking` `struct`. This object contains an array of all the particles, and for each of them, the list of contact vectors, magnitudes of contact forces, and list of neighbours.
2. Information about the termination status of CALiPPSO, the time and amount of memory allocated during the full process, list of times of each LP optimization, etc. (see the docstring of `convergence_info` for a complete list).
3. The list of values of $\sqrt{\Gamma^\star}$ obtained after each optimization.
4. The list of values of $\max_{i,\mu} \ \{s_{i,\mu}^\star \}_{i=1,\dots,N}^{\mu=1,\dots,d}$ obtained after each optimization.
5. An analogous list that specifies (with boolean variables) if isostaticity holds at the corresponding iteration.

### Changing the solver

We used the fantastic [JuMP.jl](https://jump.dev/JuMP.jl/stable/) package for model creation within our algorithm. Thus, you should be able to use any of the available solvers (that are suited for Linear Optimization). Our implementation already includes working code for the following solvers: [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl), [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl), [GLPK.jl](https://github.com/jump-dev/GLPK.jl), [Clp.jl](https://github.com/jump-dev/Clp.jl), [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl), and [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl).
Importantly, you might need to change the source code to make other solvers work correctly.

We strongly advice using [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl) (a Julia wrapper of the [Gurobi Solver](https://www.gurobi.com/)) because it's the solver we tested the most when developing our package. 

Thus, choosing a different solver (*e.g.* Gurobi), the MWE from above will look like
```julia
using CALiPPSO 
using Random
Random.seed!(123) # optional, but just for reproducibility sake of this MWE
# Choosing the seed of the Julia's RNG determines the random IC produces below with `generate_random_configuration`

using Gurobi
const grb_args = Gurobi.Env()
const grb_attributes = Dict("OutputFlag" => 0, "FeasibilityTol" => 1e-9, "OptimalityTol" => 1e-9, "Method" => 3, "Threads" => CALiPPSO.max_threads)

precompile_main_function(Gurobi, grb_attributes, grb_args) #optional, but highly recommended. This will produce a colorful output that you can safely ignore

const d, N, φ0, L = 3, 512, 0.3, 1.0
r0, Xs0 = generate_random_configuration(d, N, φ0, L) # if L is not passed, it's assumed that the systems is in a box of size 1

packing, info, Γ_vs_t, Smax_vs_t, isostatic_vs_t = produce_jammed_configuration!(Xs0, r0; 
        ℓ0=0.2*L, max_iters=500, solver=Gurobi, solver_args=grb_args, solver_attributes=grb_attributes)
```

Note that different solvers require different choices of attributes and/or arguments. Refer to the [documentation](https://rdhr.github.io/CALiPPSO.jl/dev/changing_default.html#changing_the_solver) for more options and advanced usage.

### Advanced usage

For other features, advanced usage, and more details of how `produce_jammed_configuration!` works please refer to the [documentation](https://rdhr.github.io/CALiPPSO.jl/dev/mainfunction.html).

## ToDo's

1. Write implementation of the analogous functions for polydisperse packings.
2. Register CALiPPSO in Julia's packages registry.

## Acknowledgements

This work was supported by the Simons Foundation Grant No. 454949 (G.P.) and by the he European Research Council (ERC) under the European Union’s Horizon 2020 Grant No. 101001902.










