# CALiPPSO: A Julia implementation of the ILP algorithm for jamming hard-spheres packing

This package is a pure [Julia](https://julialang.org/) implemantation of the CALiPPSO algorithm for generating jammed packings of hard spheres. This algorithm was developed in **XXXXXX** by *Artiaco, Díaz, Parisi, and Ricci-Tersenghi*. As explained there, our method consists on a **C**hain of **A**pproximate **Li**near **P**rogramming for **P**acking **S**pherical **O**bjects.

Please feel free to use/modify/improve this code as better suits you. We only ask you to give credit to our work.

```
@article{...}
```

## Documentation

You can read the full documentation of our code [here](XXXXXXXXXXXXXXX)


## Basic usage

We are working on making this library available through the standard Julia's package manager. But until then, you'll need to import the CALiPSSO's module as explained next.

### Loading the main module

We tried to make this package as easy to use as possible and, indeed, it consists of *a single* main function: `produce_jammed_configuration(Xs, R, L=1.0)`. (But before using it, be sure to have all the required dependencies already installed.) This function is defined in the `iLP-for-jamming.jl` file, which should be loaded (through the `include` function) whenever you want to make use of it. In other words, to include this function in your scope or Julia script, you only need to add the following lines:

```julia
include("src/iLP-for-jamming.jl")
using .CALiPPSO  
```
(Beware of the dot (`.`) when calling the `CALiPPSO` module, beacuse it was loaded through the file in `include`.) 

Once `CALiPPSO` has been loaded, you just need to call
```julia
packing, info, Γ_vs_t, isostatic_vs_t = produce_jammed_configuration(Xs, R, L)
```
Naturally, this function should generate a jammed packing from an initial configuration of particles with initial positions `Xs`, same radius `R`, and contained in a *periodic* (hyper-) cube of size `L` (if this argument is left unspecified, it's assumed its value is 1). `Xs` should be a $d\times N$ matrix specifying the position of each particle (*i.e.* each of the $N$ columns is the $d$-dimensional position vector of a particle).


### Output

The output of `produce_jammed_configuration` is the following:
1. A jammed packing (provided convergence was attained) stored as a `MonoPacking` `struct`. This object contains an array of all the particles, and for each of them, the list of contact vectors, magnitudes of contact forces, and list of neighbours.
2. Information about the termination status of CALiPPSO, the time and amount of memory allocated during the full process, list of times of each LP optimization, etc. (see the docstring of `convergence_info` for a complete list).
3. The list of values of $\sqrt{\Gamma^\star}$ obtained after each optimization.
4. An analogous list that specifies (with boolean variables) if isostaticity holds at the corresponding iteration.


### Advanced usage

For more advanced usage and more detailes of how `produce_jammed_configuration` please refer to the [documentation](XXXXXX).

## Disclaimer













