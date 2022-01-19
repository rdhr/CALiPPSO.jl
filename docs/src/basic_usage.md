# Basic usage

We tried to make this package as easy to use as possible and, indeed, it consists of *a single* main function: [`produce_jammed_configuration`](@ref). (But before using it, be sure to have installed all the [dependencies](@ref Dependencies).) This function is defined in the `iLP-for-jamming.jl` file, which should be loaded (through the `include` function) whenever you want to make use of it. In other words, to include this function in your scope or Julia script, you only need to add the following lines:

```julia
include("src/iLP-for-jamming.jl")
using .CALiPPSO  
```
(Beware of the dot (.) when [importing](https://docs.julialang.org/en/v1/manual/code-loading/) the `CALiPPSO` module, which is loaded via `include`.) 

In this way, [`produce_jammed_configuration`](@ref) as well as other functions and `struct`'s will be loaded into your environment. The additional functions that are exported (such as `network_of_contacts`, `check_for_overlaps`, `get_non_rattlers`, etc.) are *not* needed for a basic usage, but might be useful for analysing the packings produced by the main function. 

Once `CALiPPSO` has been loaded, you just need to call
```julia
packing, info, Γ_vs_t, isostatic_vs_t = produce_jammed_configuration(Xs, R, L)
```
Naturally, this function should generate a jammed packing from an initial configuration of particles with initial positions `Xs`, same radius `R`, and contained in a periodic (hyper-) cube of size `L` (if this argument is left unspecified, it's assumed its value is 1). `Xs` should be a ``d\times N`` matrix specifying the position of each particle (*i.e.* each of the ``N`` columns is the ``d``-dimensional position vector of a particle). Clearly, this matrix can be constructed from importing data from a `csv` or `dat` file (or any other suitable file format for that matter).

Alternatively, `Xs` might be an array of ``N`` `StaticVector`'s, each of size ``d`` and elements of [`PeriodicNumber`](@ref) type [(see this section on the types defined here)](#types-defined-in-this-package). In such case, no value of `L` should be given when calling `produce_jammed_configuration` since it is automatically inferred from the field of `PeriodicNumber`. For convenience, if `Xs` is a matrix of `Float64` elements, it can be easily converted to an array of `StaticVector`'s with `PeriodicNumber` elements using the `PeriodicVectors` function. Thus,
```julia
Xs = PeriodicVectors(Xs, L)
produce_jammed_configuration(Xs, R)
```





See below for [more specific examples](#some-tests-included).

The output of `produce_jammed_configuration` is the following:
1. A jammed packing (provided convergence was attained) stored as a [`MonoPacking`](@ref) type (see details [below](#the-packing-type)). This object contains an array of all the particles, and for each of them, the list of contact vectors, magnitudes of contact forces, and list of neighbours.
2. Information about the termination status of CALiPPSO, the time and amount of memory allocated during the full process, list of times of each LP optimization, etc. (see the docstring of [`convergence_info`](@ref) for a complete list).
3. The list of values of ``\sqrt{\Gamma^\star}`` obtained after each iteration.
4. An analogous list that specifies (with boolean variables) if isostaticity holds at the corresponding iteration.


With its default parameters and assuming you're using Gurobi, `produce_jammed_configuration` should work in most of the cases of ``d>2`` (at least we didn't experienced any in the [tests](#some-tests-included) we ran), specially when ``\varphi_0\simeq \varphi_J``. In general situations however, be warned that even if the function terminates without throwing any error, it may happen that the output is *not* a jammed packing: it can be the case that ILP terminated because the maximal number of iterations was reached, or too many consecutive non-isostatic configurations were obtained (you can tune these parameters and several others using keyword arguments, as we explain [below](#kwargs-for-controlling-convergence-and-termination-criteria-of-ilp). This latter situation is much more likely to occur if ``\varphi_0 \ll \varphi_J`` and (maybe) with [other solvers](#changing-the-solveroptimizer) instead of Gurobi. In any case [see this section](#try-it-yourself) for some examples of its usage and [the section on keyword arguments](#controlling-produce_jammed_configuration-with-keyword-arguments) to fine tune the behaviour of `produce_jammed_configuration` to better suit your needs.

!!! note "Precompilation"
    Note that executing `include("src/iLP-for-jamming.jl")` also *pre-compiles* `produce_jammed_configuration` by calling it using a small, predefined system. This causes that whenever this file is included, several things will be printed in screen. You can safely ignore all of them. Precompilation can be avoided by including a line 
    ```julia
    const precompile_main_function=false
    ``` 
    *before* loading the `iLP-for-jamming.jl` file. 
    
    However, this is discouraged if you actually want to use such function because model creation and optimization in JuMP suffers from a ["time-to-first-solve" issue](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/performance_tips/#The-%22time-to-first-solve%22-issue), which is an analogous version of the "time-to-first-plot" one. Essentially, when a function is first called, it needs to be compiled. (This is actually the general behaviour of Julia, not only of JuMP.) Thus, by calling `produce_jammed_configuration` in a small system, the function gets pre-compiled and ready to be used in much larger systems. If `precompile_main_function` is set to `false` and then you call `produce_jammed_configuration` directly into the (presumably large) configuration you want to jam it could take much longer.

    On the other hand, avoiding precompilation is useful if you just want to import the types, functions, etc. of CALiPPSO, but will not be making use of `produce_jammed_configuration`; for instance if you are only doing data analysis.

    In any case, if you do not define `precompile_main_function`, all the functions are precompiled by default.


## Understanding the screen output


## Some other features
- The dimensionality of the system is inferred from `Xs` and the periodic boundary conditions are automatically implemented through the usage of [`PeriodicNumber` types](#the-periodicnumber-type). Of course: **no overlap should be present in the initial configuration** for ILP to run properly. 

- You can (or at least should be able to) use as input any valid HS configuration generated from your favourite method (for instance, the Lubachevsky─Stillinger (LS) compression protocol as described [above](#the-initial-conditions)).
- Alternatively, you can also use the function `generate_random_configuration(d, N, ϕ)` provided here to generate a random *low-density* initial configuration of `N` particles in `d` dimensions with density `ϕ`. See however [the caveat](#the-initial-conditions) of initializing ILP with a configuration far from jamming.
- As ILP progresses, checks for overlaps are implemented automatically.
- Stability and force balance checks are implemented, to track possible numerical issues after the LP optimizations are carried out; see the section on [understanding the output](#a-simple-example-and-understanding-the-output).
- The radius of influence, ``\ell``, is automatically adjusted in such a way that only nearby pairs of particles are considered when building the set of non-overlapping constraints.
- The behaviour and other parameters of the main function can be easily controlled through [keyword arguments](#controlling-produce_jammed_configuration-with-keyword-arguments).


