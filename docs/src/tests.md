# Tests and examples provided

All the scripts referred here can be found in the [`Examples` folder](https://github.com/rdhr/CALiPPSO.jl/tree/master/Examples) of the Github repo, or also in local `CALiPPSO` folder (located in the Julia's package directory, *e.g.* `.julia/packages` in Linux). However, maybe it's better to just download the scripts into a location you find more convenient. The instructions below assume that *you are executing julia inside the folder where the scripts are located*.

In any case, two of the scripts are designed to show how CALiPPSO works in different dimensions and also the influence of the [initial condition](@ref The-initial-conditions). A third one is designed to test the different solvers, but to run it you'll need to have all of them [already installed](@ref Dependencies). And a final one is just a simple script to exemplify how you can analyze the performance of `CALiPPSO.jl`.


## Trying CALiPPSO in different dimensions

These tests are designed to jam **monodisperse** packings in ``d=3, 4, 5`` dimensions to show that (i) CALiPPSO is a versatile algorithm; and, (ii) that [`produce_jammed_configuration!`](@ref) seamlessly works on different dimensions without tunning a bunch of parameters. Instead, the number of dimensions is simply inferred from the *fixed* size of the particles' position array, implemented as an [`SVector`](https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#SVector).

In both cases, the default solver used is `GLPK`. But we have also included (and commented out) the necessary lines to run the tests using `Gurobi`. This latter solver is much faster and more precise, as mentioned before.

### Highly compressed initial configuration

First, the case where CALiPPSO works better is, naturally, when the initial condition is already close to its jamming point. (See the [Theory behind CALiPPSO](@ref) for more details.) Thus, for these examples, configurations of ``N=1024`` particles have been compressed up to a pressure of ``p=10^5`` using the [Lubachevsky--Stillinger algorithm](https://doi.org/10.1007/BF01025983) (see [here](@ref The-initial-conditions) for few more details).

To run this example, simply type the following line from a terminal

```
julia jamming-monodisperse-after-LS.jl
```

or, from inside the Julia's REPL, Jupyter Notebook, etc,

```julia
include("jamming-monodisperse-after-LS.jl")
```

Then, `produce_jammed_configuration!` will be used to jam the configurations described before, and some output will be printed in the screen, showing the progress of CALiPPSO, etc. ([Here](@ref Understanding-the-output-printed-in-screen/IO) you find more info about how to interpret such output.)
For each ``d``, once CALiPPSO converges, the program also prints the values of ``sqrt{\Gamma^\star}-1`` and ``max_{i,\mu} |\mathbf{s}_{i,\mu}^\star|`` obtained at each iteration, and some statistics of the times needed to solve the LOP instances. Furthermore [`network_of_contacts`](@ref) is also called on the *jammed packing* in every dimension. Thus, you can test that CALiPPSO works properly, not only reaching the jamming point, but also extracting the full network of contacts of the jammed configuration.

!!! note
    Note that to run this script you also need the files from where the initial condition is read. The relevant `dat` files are also contained in the `Examples` folder of the repo.


### Low density initial configuration

The second example in which we test our algorithm is very similar to the previous, but now using an initial condition with low density. More precisely, this case jams configurations of ``N=512`` particles, initialized with the particles placed at random, and with initial densities ``\varphi_0 = \{0.3, 0.15, 0.1\}``, respectively for each dimension. Such initial conditions are produced using `generate_random_configuration`, which assigns uniformly random positions to the particles' centers(avoiding overlaps, of course).


To run this example, simply type the following line from a terminal

```
julia jamming-monodisperse-random-init_conf.jl 
```

or, from inside the Julia's REPL, Jupyter Notebook, etc,

```julia
include("jamming-monodisperse-random-init_conf.jl")
```

The script works essentially as the previous one, also calling  [`network_of_contacts`](@ref) once a packing has been jammed. Note however that it takes much longer to finish because much more LP optimizations are needed, as expected.
Besides, the seed for the Julia's random number generator is initialized at the beginning of the script, for reproducibility. But it is not strictly necessary for the program to run.

Consider that running this script and similar programs might produce errors. Therefore, we do **not** guarantee that CALiPPSO can be used without errors with these types of initial conditions. It will always be better to initialize our algorithm with a density close to $\varphi_J$. In fact, note the following

!!! warning
    Given that the initial condition is very far from the jamming point (*i.e.* ``(\varphi_J - \varphi_0)/\varphi_J \sim \mathcal{O}(1)``), the first instances of the jamming LOP might lead to very large displacements, or to particles being block by ``s_{bound}`` (*i.e.* the bound on ``|\mathbf{s}_{i,\mu}|``). And even though overlaps are unlikely to occur (but are certainly possible in this scenario), the main issue we observed is that at some point the solvers fail to find any solution to the LOP. Moreover, some other types of issues might occur, specially in high dimensions.


## Testing different solvers

To show that it is easy to use CALiPPSO with different solvers we provide the script `testing-different-solvers.jl`. Here we test the `Gurobi`, `HiGHS`, `Clp`, `GLPK`, and `COSMO` solvers on a same configuration, namely the ``N=1024``, highly compressed configuration in ``d=3`` (i.e. one of the ones used in [a previous example](@ref Highly-compressed-initial-configuration)). Reading through this script should provide a good idea on how to choose and tune different solvers using keywords, as [described before](@ref changing_the_solver). For instance each solver's attributes and arguments are explicitly defined; so this script can be used as a rough guideline on how to set the appropriate parameters for a given solver.

!!! Note
    Only the GLPK solver is installed alongside `CALiPPSO.jl`. Thus, **before** running this script, be sure to have installed all these solvers. As described [in this section](@ref Dependencies), all of them (except Gurobi) are straightforward to install. Only `Gurobi` requires having a license and a manual installation of the solver itself. But it's also rather easy. 


Once you have all the solvers installed, you can run the script from a terminal in the following way

```
julia testing-different-solvers.jl
```

or, from inside the Julia's REPL, Jupyter Notebook, etc,

```julia
include("testing-different-solvers.jl")
```

The main function, [`produce_jammed_configuration!`](@ref), is first precompiled with each solver. Then, it is called with the initial ``d=3`` condition mentioned above, once for each solver, in the order: `Gurobi`, `HiGHS`, `Clp`, `GLPK`, and `COSMO`. For each of them, the program outputs some info of the CALiPPSO progress, since `verbose=true`.

When a packing has been created with all the solvers, some information about how different they are is provided. More precisely, they are compared by computing the difference in radii and the mean squared distance between the particles' positions, using the packing jammed with `Gurobi` as reference.

Besides, once CALiPPSO converges with a given solver, the times required to solve each LOP instance with said solver are also printed out. Thus, it's easy to compare their performance.


!!! warning
    The  pure-Julia optimizer, `Hypatia.jl` can also be added to the list of solvers to test (naturally, after installing it). However, we did not include it by default because in our tests we observed that it consumed about 8GB of memory, compared to less than 1GB using other solvers. Besides, it was a rather slow solver, so we did not included it to avoid waiting too long. However, as we already mentioned, it is likely that we did not choose correctly its parameters and better results can be obtained.

!!! note
    Using the `COSMO.jl` solver we observed that the LOP instances are solved with a much smaller accuracy, usually leading to non-isostatic packings. In fact, in our experience `produce_jammed_configuration!` rarely converges when using `COSMO` as solver (at least not for this system size). But once again, this might be caused by our poor choice of parameters for this solver.

## Simple example of Benchmark

In `benchmark-different-ds.jl` we provide a script to benchmark the performance of CALiPPSO using the same highly compressed initial conditions [mentioned above](@ref Highly-compressed-initial-configuration). In fact, both scripts are essentially the same, except that in this one a benchmark is performed, using the `@benchmarkable` (from the [`Benchmarktools.jl`](https://juliaci.github.io/BenchmarkTools.jl/dev/manual/) package) macro and by calling [`produce_jammed_configuration!`](@ref) on these configurations 20 different times (for each value of ``d``), and setting `verbose=false` to avoid unnecessary, repeated output in screen. 

!!! note "Dependencies for this test"
    To run this test you also need to install [`Benchmarktools.jl`](https://github.com/JuliaCI/BenchmarkTools.jl) --in order to be able to call `@benchmarkable`--  and [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl). We chose to use the latter solver so the 20 repetitions can be done much more rapidly than only using GLPK.

    Naturally, if you don't want to install HiGHS and prefer the default behaviour of `produce_jammed_configuration!`, you can do so by calling this main function without specifying the `solver`, `solver_attributes`, and `solver_args` keywords.


Once you have all installed these dependencies, you can run the script from a terminal in the following way

```
julia benchmark-different-ds.jl
```

or, from inside the Julia's REPL, Jupyter Notebook, etc,

```julia
include("benchmark-different-ds.jl")
```
