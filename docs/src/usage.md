# Using CALiPPSO

Here we describe in some detail our implementation of the ILP algorithm. If you're new to [Julia](https://julialang.org/) it might be useful to have [its documentation](https://docs.julialang.org/en/v1/) at hand, specially for understanding some terminology, (although we tried to keep it to a minimum).
Besides, several of the functions we define rely on JuMP's functionality, so if you are unfamiliar with this package it might be useful to consult [its documentation](https://jump.dev/JuMP.jl/stable/); specially the parts on [creating a model](https://jump.dev/JuMP.jl/stable/manual/models/#Create-a-model) and defining/accessing [constraints](https://jump.dev/JuMP.jl/stable/manual/constraints/).)


## Basic usage

We tried to make this package as easy to use as possible and, indeed, it consists of *a single* main function: `produce_jammed_configuration(Xs, R, L=1.0)`. (But before using it, be sure to have all the [dependencies](#dependencies) installed.) This function is defined in the `iLP-for-jamming.jl` file, which should be loaded (through the `include` function) whenever you want to make use of it. In other words, to include this function in your scope or Julia script, you only need to add the following lines:

```julia
include("src/iLP-for-jamming.jl")
using .CALiPPSO  
```
(Beware of the dot (.) when calling the `CALiPPSO` module, since it was loaded through the file in `include`.) 

In this way, `produce_jammed_configuration` as well as other functions and `struct`'s will be loaded into your environment. The additional functions (such as `network_of_contacts`, `check_for_overlaps`, `get_non_rattlers`, etc.) are *not* needed for a basic usage, but might be useful for later analysing the packing produced by the main function. 

!!! note "Precompilation"
    Note that executing `include("src/iLP-for-jamming.jl")` also *pre-compiles* `produce_jammed_configuration` by applying it to a small system. This causes that whenever this file is included, many things will be printed in screen. You can safely ignore all of them. Printing this output can be disabled by setting `precompile_main_function=false` in the first lines of `iLP-for-jamming.jl`. However, this is discourage because model creation and optimization in JuMP suffers from a ["time-to-first-solve" issue](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/performance_tips/#The-%22time-to-first-solve%22-issue), which is an analogous version of the "time-to-first-plot" one. Essentially, when a function is called first, it needs to be compiled. (This is actually the general behaviour of Julia, not only of JuMP.) Thus, by calling `produce_jammed_configuration` in a small system, the function gets pre-compiled and ready to be used in much larger systems. If `precompile_main_function` is set to `false` and then you call `produce_jammed_configuration` directly into the configuration you want to jam it could take much longer.
 

Once `CALiPPSO` has been loaded, you just need to call
```julia
packing, info, Γ_vs_t, isostatic_vs_t = produce_jammed_configuration(Xs, R, L)
```
Naturally, this function should generate a jammed packing from an initial configuration of particles with initial positions `Xs`, same radius `R`, and contained in a (hyper-) cube of size `L` (if this argument is left unspecified, it's assumed its value is 1). `Xs` should be a $d\times N$ matrix specifying the position of each particle (*i.e.* each of the $N$ columns is the $d$-dimensional position vector of a particle). Alternatively, `Xs` might be an array of $N$ `StaticVector`'s, each of size $d$ and elements of `PeriodicNumber` type [(see this section on the types defined here)](#types-defined-in-this-package). See below for [more specific examples](#some-tests-included).

The output of `produce_jammed_configuration` is the following:
1. A jammed packing (provided convergence was attained) stored as a `MonoPacking` type (see details [below](#the-packing-type)). This object contains an array of all the particles, and for each of them, the list of contact vectors, magnitudes of contact forces, and list of neighbours.
2. Information about the termination status of CALiPPSO, the time and amount of memory allocated during the full process, list of times of each LP optimization, etc. (see the docstring of `convergence_info` for a complete list).
3. The list of values of $\sqrt{\Gamma^\star}$ obtained after each iteration.
4. An analogous list that specifies (with boolean variables) if isostaticity holds at the corresponding iteration.


With its default parameters and assuming you're using Gurobi, `produce_jammed_configuration` should work in most of the cases of $d>2$ (at least we didn't experienced any in the [tests](#some-tests-included) we ran), specially when $\varphi_0\simeq \varphi_J$. In general situations however, be warned that even if the function terminates without throwing any error, it may happen that the output is *not* a jammed packing: it can be the case that ILP terminated because the maximal number of iterations was reached, or too many consecutive non-isostatic configurations were obtained (you can tune these parameters and several others using keyword arguments, as we explain [below](#kwargs-for-controlling-convergence-and-termination-criteria-of-ilp). This latter situation is much more likely to occur if $\varphi_0 \ll \varphi_J$ and (maybe) with [other solvers](#changing-the-solveroptimizer) instead of Gurobi. In any case [see this section](#try-it-yourself) for some examples of its usage and [the section on keyword arguments](#controlling-produce_jammed_configuration-with-keyword-arguments) to fine tune the behaviour of `produce_jammed_configuration` to better suit your needs.

### Understanding the screen output


### Some other features
- The dimensionality of the system is inferred from `Xs` and the periodic boundary conditions are automatically implemented through the usage of [`PeriodicNumber` types](#the-periodicnumber-type). Of course: **no overlap should be present in the initial configuration** for ILP to run properly. 

- You can (or at least should be able to) use as input any valid HS configuration generated from your favourite method (for instance, the Lubachevsky─Stillinger (LS) compression protocol as described [above](#the-initial-conditions)).
- Alternatively, you can also use the function `generate_random_configuration(d, N, ϕ)` provided here to generate a random *low-density* initial configuration of `N` particles in `d` dimensions with density `ϕ`. See however [the caveat](#the-initial-conditions) of initializing ILP with a configuration far from jamming.
- As ILP progresses, checks for overlaps are implemented automatically.
- Stability and force balance checks are implemented, to track possible numerical issues after the LP optimizations are carried out; see the section on [understanding the output](#a-simple-example-and-understanding-the-output).
- The radius of influence, $\ell$, is automatically adjusted in such a way that only nearby pairs of particles are considered when building the set of non-overlapping constraints.
- The behaviour and other parameters of the main function can be easily controlled through [keyword arguments](#controlling-produce_jammed_configuration-with-keyword-arguments).


### Dependencies

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

3. [StaticArrays.jl](https://juliaarrays.github.io/StaticArrays.jl/stable/): To use static (*i.e.* fixed size) vectors for particles positions, contact vectors, etc. This type of arrays considerably speed-up vector operations (such as sum, norm, etc.) and, also importantly, guarantee that all vectors have the same dimensionality. Indeed, the value of $d$ is inherited as the *type* of a static array (*e.g.* a contact vector is of type `SVector{d,Float64}`) and this helps to write generic code for any dimensionality.


4. [SpecialFunctions.jl](https://specialfunctions.juliamath.org/stable/): To compute the packing fraction in arbitrary dimensions we use the $\Gamma$ function implemented in this package.

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


### Getting help

We tried our best to provide complete *docstrings* of all the functions defined in this package. So most of the information should be available simply by calling the respective documentation (*i.e.* just by typing '?'). For instance, try typing `?produce_jammed_configuration` in the REPL or in a Jupyter notebook for a detailed description of this function. 

We also tried to leave clarifying comments throughout the code, so its functioning is easier to understand.

Yet not clear enough? Found a bug or an issue? Please drop us an email at: XXXXXXX

___

## How `produce_jammed_configuration` works


Our main function consists of two essentially independent parts: (1) the main ILP loop; and (2) the packing creation from the quantities obtained after using ILP. We now describe each of them.

___

### The main Linear Optimization loop


From the initial particles' position and size (*i.e.* the input of `produce_jammed_configuration`), a `while` loop is initialized until the *convergence criteria* defined [above](#introduction-how-ilp-works-and-some-terminology) are reached. More precisely, the loop continues as long as: (1) $\sqrt{\Gamma^\star}-1 \geq$` tol_Γ_convergence`; or (2) $|s^\star_{i,\mu}| \geq$ `tol_S_convergence` for $i=1, \dots, N$ and $\mu=1,\dots, d$ (although see step 4 below); or the number of iterations (*i.e* the number of LP optimizations) does not exceed `max_iters`. The default values of these 3 quantities are [given below](#list-of-default-values) and can be easily changes as [keyword arguments](#controlling-produce_jammed_configuration-with-keyword-arguments) of `produce_jammed_configuration`.


This main loop consists of the following steps:

1. The LP model creation and optimization. (Expectedly, this is done using `JuMP`)
   1. Thus, given the the particles' position and radii, the linear optimization problem of Eqs. (1-2) above is defined using the JuMP's interface and assigned to an object call `LP_model`.
      - `LP_model` includes the relevant design variables (*i.e.* the inflation factor, $\Gamma$, and particles' displacements, $\vec{\mathbf{s}}$), as well as the set of non-overlapping (*linear*) constraints (Eq. (2) above). The constraints are added to `LP_model` using the `add_non_overlapping_constraints!` function.
      - Importantly, not all pair of particles are considered for the constraints, but only those whose distance is smaller than a cut-off, $\ell$, whose value is obtained by calling `bounds_and_cutoff`.
      - Besides, the periodic boundary conditions are automatically considered, using the so called *Minimum Image Convention*. That is, the vector differences, like $\mathbf{r}_{ij}=\mathbf{r}_i - \mathbf{r}_j$ are always computed using the virtual image of the system that corresponds to the smallest value of $|\mathbf{r}_{ij}|$. See the docstring of `MIC_vector` and `MIC_distance` for more information.  

   2. The optimization is carried out simply by calling `optimize!(LP_model)`. 
      - Provided the optimizer was able to solve the LP instance, at this point we have obtained the optimal displacements ($\vec{\mathbf{s}}^\star$) and inflation factor ($\Gamma^\star$). 
      - **Note that** both of these steps are implemented in a single function: `solve_LP_instance`.
    
2. The force balance of the current packing is assessed. To do so, a *preliminary* network of contacts is constructed from the active constraints obtained in the previous step. 
    - To do so `network_of_contacts` is applied on the particles positions and list of constraints introduced in the step 1.1. This can be done because the list of constraints of each particle is stored as a `Vector{ConstraintRef}`. (See [here](https://jump.dev/JuMP.jl/stable/manual/constraints/#Constraint-containers) for more info about `ConstraintRef` in `JuMP`.) 
    - Note that, as we mentioned [above](#introduction-how-ilp-works-and-some-terminology) and showed in our paper, even if the jamming point has not been reached, the dual variables should fulfill a force-balance equation. Thus, verifying that this is the case is a convenient way of assessing whether the optimal solution found is good or not.
    - Note that this check should be performed **before** the configuration is updated, otherwise the *wrong* contact vectors would be used.

3. The configuration is updated: $\mathbf{r}_i \to \mathbf{r}_i + \mathbf{s}_i^\star$ and $\sigma_i \to \sqrt{\Gamma^\star}\sigma_i$ for $i=1,\dots, N$. 
   - These updated values will be used to formulate the next LP instance in the next iteration of ILP. 

4. A set of *preliminary* stable particles is obtained using `obtain_non_rattlers`. Rattlers are also obtained as the complement of such set.
   - This step is important in order to check if the configuration is isostatic or not. In the latter case, the isostaticity gap (*i.e.* the difference of the number of contacts, $N_c$, and the number of degrees of freedom, $N_{dof}$) may provide insight about numerical issues when determining the contact forces. Thus, even though this step is (apparently) not strictly required in order for ILP to work, it always provides very useful information. See [examples below](#a-simple-example-and-understanding-the-output).
   - Besides, rattlers should be (almost always) excluded when testing convergence related to the magnitude of $|s^\star_{i,\mu}|$. That is, because rattlers are not blocked by their neighbours, their associated optimal displacements are notably larger than those of the stable particles, and therefore we don't consider them for checking when the main loop should terminate. So, this step is needed *in practice* for the correct functioning of ILP.

5. If the `verbose` option is set to `true`, the following information is printed (see [the dedicated section](#a-simple-example-and-understanding-the-output) on how to interpret such info):
   1. The number of LP iteration and the state of the optimization thrown after `optimize!(LP_model)` finishes, *i.e.* the value of `termination_status(LP_model)`. This latter could be, *e.g.* `OPTIMAL`, `INFEASIBLE`, `TIME_LIMIT` , etc.
   2. The value of $\sqrt{\Gamma^\star} -1$ and $\max |s^\star_{i,\mu}|$, where $i \in \text{non-rattlers}$.
   3. The time required to solve the LP problem.
   4. Few statistics about the number of constraints included in `LP_model`
   5. Information about whether the preliminary configuration is isostatic or not, as well as the isostaticity gap, $N_c - N_{dof}$, in the latter case.
   6. Few statistics of the coordination number.
   7. Number of non-rattlers (or stable particles).
   8. Maximum mismatch in the force balance (per particle).
   9. Sample of the 10 smallest forces in the system. This information is also useful for evaluating if numerical issues are present.

6. Call the `check_for_overlaps` function to check if there are any overlaps once the configuration has been updated.
   - Of course, **there shouldn't be!**
   -  ... but given that we live in a world of *finite precision* and that we actually aim for a condition in which some of the **constraints are saturated**, it can happen that the LP instance was not solved within the required accuracy. See [this section](#setting-the-required-precision) to learn how to control the overall precision of ILP, and [how to tune the options](#changing-the-default-options) for setting the tolerance with which an overlap is identified.
   -  When an overlap *does* occur, an error is thrown an ILP terminates, also terminating the main process since `error` is called. Nevertheless, some other information is shown, that can be used, hopefully, to trace back what happened.
   -  If you think that the problem is the related to numerical issues, be sure to understand [how the precision of ILP is determined](#setting-the-required-precision).
   -  Note also that a *real* overlap can also occur (*i.e*. once in which a pair of particles is overlapping by an amount much larger than the accuracy with which a solver fulfills the constraints). 
      -  When this occurs, it is likely that the value of `sbound` is too large with respect of the value of $\ell$ used. So try using a larger value of `ℓ0` (if this occurred when $\Gamma^\star$ was still relatively large).
      -  Try using a smaller value of `sbound`.
      -  Or try redefining the bounds with which $\ell$ is adjusted as ILP progresses.
      -  These three options can be set as keywords arguments of `produce_jammed_configuration` as explained [below](#changing-the-default-options). See the docstring of `bounds_and_cutoff` for more info.
7. Check if convergence criteria are fulfilled. Otherwise, go back to step 1.





Note that steps 5 and 6, by default, are only performed during the first few iterations (10) and at given intervals (also 10). To change how often the ILP info is printed out (respectively how often overlaps checks are performed) set the keyword argument `monitor_step` (respectively `interval_overlaps_check`) to the desired value. Instead, to select in how many initial iterations to include these steps, use `initial_monitor` (for printing info) and `initial_overlaps_check` for overlaps checks. More details can be found [here](#kwargs-for-controlling-output-printing-on-terminal) and [here](#kwargs-for-controlling-precision-of-overlaps-and-force-balance-tests).

__________________________

### Creating the final packing

Clearly, a lot of data is contained in a single packing, like the set of all particles position, the network of contacts, etc. Moreover, the information related to the ILP algorithm itself (*e.g.* termination status, number of iterations, etc.). To efficiently store, access, and manipulate all of them we defined few *composite types* or `struct`'s (aka *objects* in other languages). They are:

- `PeriodicNumber`
- two particle types (either `MonoParticle` or `Particle`)
- two packing types: `MonoPacking` (composed of `MonoParticle`'s) and `PolyPacking` (composed `Particle`'s), and used to model mono- or polydisperse systems, respectively.

The specific properties of each of these types are explained [below](#types-defined-in-this-package), but the idea is that each of them contains the necessary information of the object it represents. For instance:
- A `Particle` is assigned a position (as an `SVector` of `PeriodicNumber` elements), a radius, a set of neighbours, and the corresponding set of contact vectors and forces.
  - `MonoParticle` is essentially the same type of object, except that no radius is assigned to it (since it is stored as a field of the packing itself.)
- Analogously, a `PolyPacking` contains an array of `Particle` objects, and also include information about its isostaticity, mechanical stability, and whether it is jammed or not.
  - Similarly, `MonoPacking` has the same fields and an extra one to store the size of all `MonoParticles` it is composed of.


Now, once the main ILP loop has terminated (possibly producing a jammed packing), the following steps are carried out:

1. It is assessed whether, 1. ILP reached convergence (in which case we create a flag `jammed=true`); or 2. if the loop ended because `max_iters` was exceeded or too many consecutive non-isostatic solutions were obtained (in which case `jammed=false`).

2. Using the values of `Xs` and the particles' size after the last LP optimization, as well as the relevant constraints, `final_packing` is created.
   - Clearly, this is an essential step. It is done by calling the *function* `MonoPacking` or `PolyPacking`, for monodisperse or polydisperse configurations, respectively.
   - Thus, for instance, the *function* `MonoPacking` uses the set of constraints and the particles' position (as well as some other secondary arguments) and constructs the set of all `MonoParticle`'s objects (each with a position, list of contacts, etc.). This set of `MonoParticle`'s is then assigned to a `MonoPacking`, along with the particles' radius and the value of `jammed`.
   - See the docstring's of `MonoPacking` and `PolyPacking` for more info.

3. The isostaticity of `final_packing` is assessed calling `is_isostatic(final_packing)`.

4. The sum of forces on each particle in `final_packing` is computed. 
   - If any of them is greater than a  tolerance value (fixed by the kwarg `tol_mechanical_equilibrium`), then `fine_tune_forces!` is called. This function is very similar to `solve_LP_instance` defined above, but it also updates the state of `final_packing`. More precisely
     - An additional LP problem is created and optimized. An important difference with `solve_LP_instance` is that in this LP problem $|s_{i,\mu}|$ is *un*bounded.
     - The forces magnitudes and contact vectors of each particle in `final_packing` are updated by calling `update_packing_forces!`. This function essentially uses `network_of_contacts` (described above) to construct all the real contacts from the constraints of this additional LP optimization.
     - *Nb*: in this additional optimization **none of the positions are updated**.
     - The isostaticity of the updated `final_packing` is reassessed.
   - If each particle is in mechanical equilibrium, within the same tolerance value, the algorithm jumps to the next step.

5. A final overlaps check is performed on `final_packing`. Note that this is done also by the `check_for_overlaps` function, but *using the method for* `MonoPacking` type.





<!-- However, all the optimizations are done using [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl) (the Julia wrapper of the [Gurobi Optimizer](https://www.gurobi.com/).) Nevertheless, the script is easily customizable to use the optimizer of your choice (see below). -->


____
## Changing the default options

For ease of use, several default values have been defined in our code. Some of them are specific or related to the solver we used, while others specified convergence criteria, etc. In any case, we stress that all of them have been extensively tested *only with the Gurobi solver*, and when the ILP's [initial configuration was obtained after LS compression](#the-initial-conditions). Therefore, if you want to [use a different solver](#changing-the-solveroptimizer) or initialize ILP from a different type of configuration it is likely that you'll need to make some small changes to these default parameters.

### List of default values

The main default values, all of them defined as *global* variables (using `const`), are the following:


| Variable      | default value |  Role  | defined in file |
| :------------ | :-----------: |:--- | :--------------------- |
| `default_tol_overlap` | 1e-8          | Tolerance for identifying overlaps (*e.g.* in functions like `check_for_overlaps`) |`iLP-for-jamming.jl` |
| `default_tol_optimality` | 1e-9       | General precision of optimal solutions  |`iLP-for-jamming.jl` |
| `default_tol_Γ_convergence` | `eps()` (*i.e.* `2.22e-16`)|  Tolerance for testing convergence of $\sqrt{\Gamma^*}-1$ |`iLP-for-jamming.jl` |
| `default_tol_displacements` | 1e-10 | Tolerance for testing convergence of $s_{i,\mu}^\star$   | `iLP-for-jamming.jl` |
| `default_tol_zero_forces` | 1e-10   | Tolerance for considering a force different from 0 | `iLP-for-jamming.jl` |
| `default_tol_force_equilibrium` | 1e-12 | Tolerance for testing force balance (per particle) | `Particles-types-and-functions.jl` |
| `max_threads` | `Int(round(Sys.CPU_THREADS/2))` (*i.e.* half of the max threads available) | Number of threads to be used by solvers that allow parallelization | `iLP-for-jamming.jl` |

On the other hand, to use Gurobi as the default solver, the following lines are included in the `iLP-for-jamming.jl` file. They define the attributes and other options passed to the solver when a model is created, modified, or optimized. Analogous lines are also included (but have been commented out) for [other solvers](#changing-the-solveroptimizer).

```julia
const default_solver = :Gurobi
const default_args = Gurobi.Env()
const default_solver_attributes = Dict("OutputFlag" => 0, "FeasibilityTol" => tol_optimality, "OptimalityTol" => tol_optimality, "Method" => 3, "Threads" => max_threads)

```


**Remark**: The value of `tol_optimality` corresponds to the most precise one allowed by Gurobi (both for 'OptimalityTol' and 'FeasibilityTol'); go [here](https://www.gurobi.com/documentation/9.1/refman/optimalitytol.html) for more information. The value of `tol_overlap` was then chosen accordingly; that is, several times (10) larger, because the optimal value of each degree of freedom is determined with an accuracy of `tol_optimality`.

____
### Controlling `produce_jammed_configuration` with keyword arguments

The full list of keyword arguments (kwargs) of `produce_jammed_configuration` can be readily accessed from the documentation. Here we provide the same list (with their default values, defined above), with a more detailed description when needed. Thus, the value of any of them can be conveniently tunned to your needs by calling `produce_jammed_configuration(Xs0, R, L; kwarg=<your chosen value>)`.

#### Kwargs for controlling how constraints are assigned to `LP_model` and setting bounds on displacements

1. `ℓ0=3.4*R`: Initial value and upper-bound of the radius of influence ($\ell$) for assigning constraints. Used when calling `bounds_and_cutoff` from within `solve_LP_instance` and `fine_tune_forces!`. See item 4 below.
2. `sqrΓ0=1.01`: Initialization value of $\sqrt{\Gamma}$; it is used to provide a guess of the value of bound of $|s_{i,\mu}|$ when $\Gamma$ is relatively large. Used when calling `bounds_and_cutoff` from within `solve_LP_instance` and `fine_tune_forces!`. See item 4 below.
3. `sbound=0.01`: Fraction of R used as bound for $|s_{i,\mu}|$ during the last optimizations, *i.e.* when $\Gamma^\star$ has very small values. Used when calling `bounds_and_cutoff` from within `solve_LP_instance` and `fine_tune_forces!`. See next item.
4. `thresholds_bounds=(5e-4, 1e-5)`: These two value control when the different behaviours of `bounds_and_cutoff` are triggered. Thus, 
   1. when $\sqrt{\Gamma^\star}-1$ (or `sqrΓ-1` in the main function definition) is larger than the first value, $\ell=$`ℓ0` and $|s_{i,\mu}|=$ `(ℓ0 - sqrt(sqrΓ)*R)/sqrt(d)`
   2. when $\sqrt{\Gamma^\star}-1$ is between the two values, $\ell=$`min(4R, ℓ0)` and $|s_{i,\mu}|=0.1R$
   3. when $\sqrt{\Gamma^\star}-1$ is smaller than the second value, $\ell=$`min(2.7*R, ℓ0)` and $|s_{i,\mu}|=$ `sbound*R`



#### Kwargs for controlling convergence and termination criteria of ILP

1. `tol_Γ_convergence=default_tol_Γ_convergence`: Determines the value below which $\sqrt{\Gamma^\star}-1$ is considered zero (so convergence in the inflation factor has been reached).
2. `tol_S_convergence=default_tol_displacements`: Determines the value below which $|s_{i,\mu}^\star|$ is considered zero (so convergence in particles displacements ─restricted to rattlers─ has been reached).
3. `max_iters=1000`: Maximum number of iterations (*i.e.* LP optimizations) allowed before stopping the main ILP loop.
4. `non_iso_break=10`: Maximum number of preliminary configurations that can be obtained *consecutively* before the main ILP loop is terminated.

#### Kwargs for controlling precision of overlaps and force balance tests

1. `tol_mechanical_equilibrium=default_tol_force_equilibrium`: When the norm of the total force acting on any particle is smaller than this quantity, said particle is considered to be in equilibrium.
2. `zero_force=default_tol_zero_forces`: This is the threshold for determining when a force, or dual variable, is *active*. In other words, whenever `shadow_price(constraint)`, where `constraint` is a `ConstraintRef`, outputs a value larger than `zero_force`, we consider that such constraint is active, and its value is precisely `shadow_price(constraint)`. 
3. `tol_overlap=default_tol_overlap`: Maximum value of overlap that can occur between particles. That is, if $|\mathbf{r}_{ij}| - \sigma_{ij}\leq -$ `tol_overlap`, an overlap is said to have occurred.
4. `initial_overlaps_check=initial_monitor`: After each LP optimization `check_for_overlaps` is called, after the configuration has been updated, for these many *initial* iterations. (`initial_monitor` is described below.)
5. `interval_overlaps_check=10`: after the configuration has been updated, `check_for_overlaps` is also called every `interval_overlaps_check` iterations.

#### Kwargs for controlling output printing on terminal

1. `verbose=true`: turns on/off the printing of information during the main ILP loop.
2. `monitor_step=10`: The info about the ILP process is printed out after these many steps (besides other criteria).
3. `initial_monitor=monitor_step`: `verbose` is set to true for these many *initial* iterations.


In our experience, most of the problems (*e.g.* a real overlap or an optimization that didn't return `OPTIMAL` as termination status) occurred during the initial steps of ILP. Thus, obtaining as much information as possible, as well as being extra-cautious by verifying that no overlaps are present after *each* iteration, becomes important during such initial stage.

_____

### Changing the solver/optimizer

#### Solvers we tried
Given that we have used the JuMP interface for solving each LP instance, *in principle*, you can use *any* of the [JuMP compatible solvers](https://jump.dev/JuMP.jl/stable/installation/). Of course, maybe not all of them are suitable for the type of LP optimization our algorithm requires, but there should be ample choice. Indeed, besides Gurobi, we tested the following solvers:

1. [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl): the Julia wrapper of [HiGHS](https://www.maths.ed.ac.uk/hall/HiGHS/). This is also a very performant and precise solver (yet noticeably slower than Gurobi).
2. [GLPK.jl](https://github.com/jump-dev/GLPK.jl): the Julia wrapper of the famous [GNU Linear Programming Kit](https://www.gnu.org/software/glpk/). GLPK is a very well known and amply used library, so it can be used reliably.
3. [Clp.jl](https://github.com/jump-dev/Clp.jl): the Julia wrapper of the [COIN-OR Linear Programming Interface](https://projects.coin-or.org/Clp). This solver also works very well, although it is the slowest of these three options.

All these solvers produced very good results in the sense that all contacts, rattlers, etc. were correctly identified and consistent between them; no overlaps occurred, and isostaticity was always achieved. Besides, all of them are free and open-source, so they are a great alternative to Gurobi if you are unable to obtain a license. All solvers are installed automatically when wrapper is installed within julia. That is, you can use `Pkg.add("HiGHS")`, `Pkg.add("Clp")`, or `Pkg.add("GLPK")` from Julia, *without* the need of installing neither solver beforehand.

We also tested some *pure-Julia* solvers, like [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl) and [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl). We were able to execute our code with both of them without any error, although in several cases the final packing was non-isostatic due to a lack of precision in identifying the contact forces. This issue is possibly related to our poor choice of the [solver parameters](#selecting-solver-and-specifying-its-attributes), but we didn't perform any additional tests. (So if you know how to improve this *please let us know*.)

#### Selecting solver and specifying its attributes

Just as with the other options of `produce_jammed_configuration`, choosing a solver is conveniently done through [keyword arguments](#controlling-produce_jammed_configuration-with-keyword-arguments):

1. `solver::Symbol=default_solver`: This kwarg must correspond to the name of the solver that you want to use, as *a Symbol*. For example, as [mentioned above](#list-of-default-values), we defined `default_solver = :Gurobi`; or if you want to use, *e.g.* the HiGHS solver, you should pass `solver=:HiGHS` as argument. (Note the **colon** (`:`) before the name of the solver; this is what makes it of `Symbol` type, and it is very important that is included). Thus, other possible values of the `solver` kwarg are: `:GLP`, `:Clp`, `:Hypatia`, etc.
   - If you want to use a different solver, be sure to load the corresponding package *inside* the `CALiPPSO` module (see line `23` of `iLP-for-jamming.jl`). 
   - Note that `solver` is not used by the main function itself, but instead is passed, also as a kwarg, to the function where an actual optimization is performed, *i.e.* `solve_LP_instance` and `fine_tune_forces!`. 
   - In such functions, the optimizer of the given `solver` is obtained by evaluating `eval(solver).Optimizer()`. Thus, `solver` should be a symbol that yields a valid optimizer for a JuMP `Model`. For instance, if `solver=:GLPK`, `Model(eval(solver).Optimizer)` is equivalent to `Model(GLPK.Optimizer)` (note there's *no* colon before "GLPK" in the last expression). [Here](https://docs.julialang.org/en/v1/manual/metaprogramming/) you can find more information on how symbols are parsed into expression in Julia, and [here](https://jump.dev/JuMP.jl/stable/manual/models/) and [here](https://jump.dev/JuMP.jl/stable/reference/models/) the docs for JuMP model creation with different optimizers.
2. `solver_attributes::Dict=default_solver_attributes`: These are the set of options or parameters that control the behaviour of your solver. It should be passed as a `Dict` type (see [the example above](#list-of-default-values) for the default one). The idea is that this `Dict` should contain any of the parameters or options that can be set using the function [`set_optimizer_attributes`](https://jump.dev/JuMP.jl/stable/reference/models/#JuMP.set_optimizer_attributes) or other similar functions, as [explained here](https://jump.dev/JuMP.jl/stable/reference/models/#Working-with-attributes).
   - Thus, for instance, with these parameters you can specify the accuracy of the solver, the method to use, iterations limit, etc.
   - In the first lines of the file `iLP-for-jamming.jl` we have included some possible choices for each of [the solvers we tried](#solvers-we-tried).
3. `solver_args=default_args`: These are supposed to be *arguments* (and *not* attributes) that are passed as arguments to `Optimizer()` of the chosen solver. Note that the behaviour of *each solver* is different. Thus, **if you are testing a different solver to the ones we mention here, we suggest you first set** `solver_args=nothing` to be sure the function will execute properly. In fact, if `solver_args===nothing` results in `false` and you're using a non-preconfigured solver, an error will be thrown. Of course, you can modify this, as explained in the last "important" point.
   - Creating a model with arguments in JuMP has a rather strange syntax: `Model(()-> eval(solver).Optimizer(solver_args))`

   - However, given that *each* solver has different ways of passing arguments to its optimizer, we had to resort to a somewhat silly implementation for each of the solvers:
     - For the *HiGHS*, *Clp*, and *COSMO* solvers: `solver_args=nothing` because they are completely controlled by attributes (so you'll only need to specify `solver_attributes`). Hence, a model is created as `Model(()-> HiGHS.Optimizer())` for the HiGHS solver, and analogously for the other ones.
     - For *Gurobi*: `solver_args=Gurobi.Env()` and the constructor is `Model(()-> eval(solver).Optimizer(solver_args))`. This is needed to avoid Gurobi to retrieve your license every time a model is created (*i.e.* for each ILP iteration), and printing out an annoying line like `Academic license - for non-commercial use only - expires XXXXXXXXX`. If you don't care about this, you can also set `solver_args=nothing`.
     - *GLPK* and *Hypatia* passes the arguments as *keywords arguments*, an therefore `solver_args` should be a `NamedTuple`. For instance, for *GLPK*, `solver_args=(want_infeasibility_certificates=false, method=GLPK.MethodEnum(0) )` and a model is created like `Model(()->GLPK.Optimizer(;solver_args...))`. The same principle is applied for *Hypatia*. 
       - Note however that this is only because we didn't manage to fully control these solvers through attributes. If you know how to do so, please let us know.
    - **Important**: Because of the variety of ways to pass an optimizer's arguments, in our code we had to resort to a (dirty) implementation such that, whenever `solver_args` is different to `nothing`, we use a series of `if-elseif-else` statements to properly set the optimizer of the *specific* solver. These means that if you want to use a solver not considered here and pass attributes to it when the model is created, **you'll need to specify the proper argument passing syntax**. To do so, you'll have to modify the `if-elseif-else` statements that defines `optimizer` in the definition of `solve_LP_instance` and `fine_tune_forces!`.



#### Testing different solvers


#### Setting the required precision


________
