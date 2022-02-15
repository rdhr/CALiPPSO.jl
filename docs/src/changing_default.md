# Changing the default options

For ease of use, several default values have been defined in our code. Some of them are specific or related to the solver we used, while others specified convergence criteria, etc. In any case, we stress that all of them have been extensively tested *only with the Gurobi solver*, and when the CALiPPSO's [initial configuration was obtained after LS compression](@ref The-initial-conditions). Therefore, if you want to [use a different solver](@ref changing_the_solver) or initialize CALiPPSO from a different type of configuration it is likely that you'll need to make some small changes to these default parameters.

## [List of default values](@id list-defaults)

The main default values, all of them defined as *global* variables (using `const`), are the following:


| Variable      | default value |  Role  | defined in file |
| :------------ | :-----------: |:--- | :--------------------- |
| `default_tol_overlap` | 1e-8          | Tolerance for identifying overlaps (*e.g.* in functions like `check_for_overlaps`) |`CALiPPSO.jl` |
| `default_tol_optimality` | 1e-9       | General precision of optimal solutions  |`CALiPPSO.jl` |
| `default_tol_Γ_convergence` | `eps()` (*i.e.* `2.22e-16`)|  Tolerance for testing convergence of ``\sqrt{\Gamma^*}-1`` |`CALiPPSO.jl` |
| `default_tol_displacements` | 1e-10 | Tolerance for testing convergence of ``s_{i,\mu}^\star``   | `CALiPPSO.jl` |
| `default_tol_zero_forces` | 1e-10   | Tolerance for considering a force different from 0 | `CALiPPSO.jl` |
| `default_tol_force_equilibrium` | 1e-12 | Tolerance for testing force balance (per particle) | `Particles-types-and-functions.jl` |
| `max_threads` | `Int(round(Sys.CPU_THREADS/2))` (*i.e.* half of the max threads available) | Number of threads to be used by solvers that allow parallelization | `CALiPPSO.jl` |

On the other hand, to use Gurobi as the default solver, the following lines are included in the `CALiPPSO.jl` file. They define the attributes and other options passed to the solver when a model is created, modified, or optimized. Analogous lines are also included (but have been commented out) for [other solvers](@ref changing_the_solver).

```julia
const default_solver = :Gurobi
const default_args = Gurobi.Env()
const default_solver_attributes = Dict("OutputFlag" => 0, "FeasibilityTol" => default_tol_optimality, "OptimalityTol" => default_tol_optimality, "Method" => 3, "Threads" => max_threads)
```


!!! note "Default values for Gurobi"
    The value of `tol_optimality` corresponds to the most precise one allowed by Gurobi (both for 'OptimalityTol' and 'FeasibilityTol'); go [here](https://www.gurobi.com/documentation/9.1/refman/optimalitytol.html) for more information. The value of `tol_overlap` was then chosen accordingly; that is, several times (10) larger, because the optimal value of each degree of freedom is determined with an accuracy of `tol_optimality`.

    On the other hand, the values that determine the convergence criteria are clearly smaller than such tolerance value allowed by Gurobi. Nevertheless, in our tests we observed that after enough iterations these more stringent conditions can actually be met. In fact, achieving ``\sqrt{\Gamma^\star} -1 \leq `` `eps()` is relatively simple. The more complicated part is to reach a configuration with ``\max |\mathbf{s}_{i,\mu}^\star|_{i=1,\dots,N}^{\mu=1,\dots, d} \leq `` `<=tol_S_conv`. Specially for relatively large systems, *e.g.* ``N\geq 5000``.

    Thus, when dealing with large systems you might want to try something like
    ```julia
    const default_tol_Γ_convergence = default_tol_optimality
    const default_tol_displacements = default_tol_optimality
    ```
    in order to speed up convergence.

    Nevertheless, in all the configurations we tested, we did *non* find anyone for which the force balance condition couldn't be met within the tolerance defined by `default_tol_force_equilibrium`, despite this value being much smaller than Gurobi's highest accuracy.

---
## [Controlling `produce_jammed_configuration` with keyword arguments](@id kwargs-control)

The full list of keyword arguments (kwargs) of [`produce_jammed_configuration`](@ref) can be readily accessed from its docstring. Here we provide the same list (with their default values, defined above), and a more detailed description when needed. Thus, the value of any of them can be conveniently tunned to your needs by calling `produce_jammed_configuration(Xs0, R, L; kwarg=<your chosen value>)`.

### Kwargs for controlling how constraints are assigned to `LP_model` and setting bounds on displacements

As explained in our paper, the constraints a particle is subject to are assigned according to a neighbours-list approach. Thus, the associated neighbours o a given sphere are all the particles within a certain distance ``\ell``. In other words, ``\ell`` determines the radius of influence on a particle. This and other quantities are computed using [`bounds_and_cutoff`](@ref CALiPPSO.bounds_and_cutoff). This function is called within the [`solve_LP_instance`](@ref CALiPPSO.solve_LP_instance) and [`fine_tune_forces!`](@ref CALiPPSO.fine_tune_forces!) functions

1. `ℓ0=3.4*R`: Initial value of the radius of influence (``\ell``) for assigning constraints. `ℓ0` is also used as upperbound for ``\ell`` in subsequent steps. See item 4 below for more information.
2. `sqrΓ0=1.01`: Initialization value of ``\sqrt{\Gamma}``; it is used to provide a guess of the value of ``s_{bond}`` --*i.e.* and upper bound of ``|s_{i,\mu}|``-- when ``\Gamma`` is relatively large. See item 4 below for more information.
3. `sbound=0.01`: Fraction of ``R`` used as bound for ``|s_{i,\mu}|`` during the last optimizations, *i.e.* when ``\Gamma^\star`` has very small values. 
4. `thresholds_bounds=(5e-4, 1e-5)`: These two value control when the different behaviours of `bounds_and_cutoff` are triggered. Thus, 
   1. when ``\sqrt{\Gamma^\star}-1`` (or `sqrΓ-1` in the main function definition) is larger than the first value, ``\ell=```ℓ0` and ``|s_{i,\mu}|=\frac{1}{2\sqrt{d}} (\ell - \sqrt{\Gamma_0^\star}\sigma``, where ``\Gamma_0^\star`` is the optimal value of ``\Gamma`` from the previous iteration, or `sqrΓ0` squared in the first one.
   1. when ``\sqrt{\Gamma^\star}-1`` is between the two values, ``\ell=```min(4R, ℓ0)` and ``|s_{i,\mu}|=0.1R``
   2. when ``\sqrt{\Gamma^\star}-1`` is smaller than the second value, ``\ell=```min(2.7*R, ℓ0)` and ``|s_{i,\mu}|=`` `sbound*` ``R``



### Kwargs for controlling convergence and termination criteria of [the main loop](@ref mainloop)

1. `tol_Γ_convergence=default_tol_Γ_convergence`: Determines the value below which ``\sqrt{\Gamma^\star}-1`` is considered zero (so convergence in the inflation factor has been reached).
2. `tol_S_convergence=default_tol_displacements`: Determines the value below which ``|s_{i,\mu}^\star|`` is considered zero (so convergence in particles displacements ─restricted to *non*-rattlers─ has been reached).
3. `max_iters=1000`: Maximum number of iterations (*i.e.* LP optimizations) allowed before stopping the main CALiPPSO's loop.
4. `non_iso_break=10`: Maximum number of non-isostatic configurations that can be obtained *consecutively* before the main CALiPPSO's loop is terminated.

### Kwargs for controlling precision of overlaps and force balance tests

1. `tol_mechanical_equilibrium=default_tol_force_equilibrium`: When the norm of the total force acting on any particle is smaller than this quantity, said particle is considered to be in equilibrium.
2. `zero_force=default_tol_zero_forces`: This is the threshold for determining when a force, or dual variable, is *active*. In other words, whenever [`shadow_price`](https://jump.dev/JuMP.jl/stable/reference/solutions/#JuMP.shadow_price)`(constraint)` --with `constraint` being a [`ConstraintRef`](https://jump.dev/JuMP.jl/stable/reference/constraints/#JuMP.ConstraintRef)-- outputs a value larger than `zero_force`, we consider that such constraint is active, and its value is precisely `shadow_price(constraint)`. 
3. `tol_overlap=default_tol_overlap`: Maximum value of overlap that can occur between particles. That is, if ``\sigma_{ij} - |\mathbf{r}_{ij}| \geq `` `tol_overlap`, an overlap is said to have occurred.
4. `initial_overlaps_check=initial_monitor`: During each of these many *initial* iterations [`check_for_overlaps`](@ref) is called, after the configuration has been updated following an LP optimization. (`initial_monitor` is described below.)
5. `interval_overlaps_check=10`: after the configuration has been updated, [`check_for_overlaps`](@ref) is also called every `interval_overlaps_check` iterations.

### Kwargs for controlling output printing on terminal

1. `verbose=true`: turns on/off the printing of information during the main CALiPPSO's loop.
2. `monitor_step=10`: The info about the progress of CALiPPSO is printed out after these many steps (besides other criteria).
3. `initial_monitor=monitor_step`: `verbose` is set to true for these many *initial* iterations.


In our experience, most of the problems (*e.g.* a real overlap or an optimization that didn't return `OPTIMAL` as termination status) occurred during the initial steps of CALiPPSO. Thus, obtaining as much information as possible, as well as being extra-cautious by verifying that no overlaps are present after *each* iteration, becomes important during such initial stage.

---

## [Changing the solver/optimizer](@id changing_the_solver)

### Solvers we tried
Given that we have used JuMP's interface for solving each LP instance, *in principle*, you can use *any* of the [JuMP compatible solvers](https://jump.dev/JuMP.jl/stable/installation/). Of course, maybe not all of them are suitable for the type of LP optimization our algorithm requires, but there should be ample choice. Indeed, besides Gurobi, we tested the following solvers:

1. [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl): the Julia wrapper of [HiGHS](https://www.maths.ed.ac.uk/hall/HiGHS/). This is also a very performant and precise solver (yet noticeably slower than Gurobi).
2. [GLPK.jl](https://github.com/jump-dev/GLPK.jl): the Julia wrapper of the famous [GNU Linear Programming Kit](https://www.gnu.org/software/glpk/). GLPK is a very well known and amply used library, so it can be used reliably.
3. [Clp.jl](https://github.com/jump-dev/Clp.jl): the Julia wrapper of the [COIN-OR Linear Programming Interface](https://projects.coin-or.org/Clp). This solver also works very well, although it is the slowest of these three options.

All these solvers produced very good results in the sense that all contacts, rattlers, etc. were correctly identified and consistent between them; no overlaps occurred, and isostaticity was always achieved. Besides, all of them are free and open-source, so they are a great alternative to Gurobi if you are unable to obtain a license. Besides, all of these solvers are installed automatically when their wrapper is installed within julia. That is, you can use `Pkg.add("HiGHS")`, `Pkg.add("Clp")`, or `Pkg.add("GLPK")` from Julia, *without* the need of installing neither solver beforehand.

We also tested some *pure-Julia* solvers, like [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl) and [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl). We were able to execute our code with both of them without any error, although in several cases the final packing was non-isostatic due to a lack of precision in identifying the contact forces. This issue is possibly related to our poor choice of the [solver parameters](@ref Selecting-a-solver-and-specifying-its-attributes), but we didn't perform any additional tests. (So if you know how to improve this [*please contact us*](mailto:rafael.diazhernandezrojas@uniroma1.it).)

### Selecting a solver and specifying its attributes

Just as with the other options of [`produce_jammed_configuration`](@ref), choosing a solver is conveniently done through [keyword arguments](@ref kwargs-control):

1. `solver::Symbol=default_solver`: This kwarg must correspond to the name of the solver that you want to use, as *a Symbol*. For example, as [mentioned above](@ref list-defaults), we defined `default_solver = :Gurobi`; or if you want to use, *e.g.* the HiGHS solver, you should pass `solver=:HiGHS` as argument. (Note the **colon** (`:`) before the name of the solver; this is what makes it of `Symbol` type, and it is very important that is included). Thus, other possible values of the `solver` kwarg are: `:GLP`, `:Clp`, `:Hypatia`, etc.
   - If you want to use a different solver, be sure to load the corresponding package *inside* the `CALiPPSO` module (see line `23` of `CALiPPSO.jl`). 
   - Note that `solver` is not used by the main function itself, but instead is passed, also as a kwarg, to the function where an actual optimization is performed, *i.e.* [`solve_LP_instance`](@ref CALiPPSO.solve_LP_instance) and [`fine_tune_forces!`](@ref CALiPPSO.fine_tune_forces!). 
   - In such functions, the optimizer of the given `solver` is obtained by evaluating `eval(solver).Optimizer()`. Thus, `solver` should be a symbol that yields a valid optimizer for a JuMP `Model`. For instance, if `solver=:GLPK`, `Model(eval(solver).Optimizer)` is equivalent to `Model(GLPK.Optimizer)` (note there's *no* colon before "GLPK" in the last expression). [Here](https://docs.julialang.org/en/v1/manual/metaprogramming/) you can find more information on how symbols are parsed into expression in Julia, and [here](https://jump.dev/JuMP.jl/stable/manual/models/) and [here](https://jump.dev/JuMP.jl/stable/reference/models/) you can consult the relevant docs for model creation with different optimizers in JuMP.
2. `solver_attributes::Dict=default_solver_attributes`: These are the set of options or parameters that control the behaviour of your solver. It should be passed as a `Dict` type (see [the example above](@ref list-defaults) for the default one). The idea is that this `Dict` should contain any of the parameters or options that can be set using the function [`set_optimizer_attributes`](https://jump.dev/JuMP.jl/stable/reference/models/#JuMP.set_optimizer_attributes) or other similar functions, as [explained here](https://jump.dev/JuMP.jl/stable/reference/models/#Working-with-attributes).
   - Thus, for instance, with these parameters you can specify the accuracy of the solver, the method to use, iterations limit, etc.
   - In the first lines of the file `CALiPPSO.jl` we have included some possible choices for each of [the solvers we tried](@ref Solvers-we-tried).
3. `solver_args=default_args`: These are supposed to be *arguments* (and *not* attributes) that are passed as arguments to [`Optimizer()`](https://jump.dev/JuMP.jl/stable/moi/tutorials/implementing/#The-Optimizer-object) of the chosen solver. Note that the behaviour of *each solver* is different. Thus, **if you are testing a different solver to the ones we mention here, we suggest you first set** `solver_args=nothing` to be sure the function will execute properly. In fact, if `solver_args===nothing` results in `false` and you're using a solver that has not been preconfigured, an error will be thrown. Of course, you can modify this, as explained at the end.
   - Creating a model with arguments in JuMP has a rather strange syntax: `Model(()-> eval(solver).Optimizer(solver_args))`

   - However, given that *each* solver has different ways of passing arguments to its optimizer, we had to resort to a somewhat silly implementation for each of the solvers:
     - For the *HiGHS*, *Clp*, and *COSMO* solvers: `solver_args=nothing` because they are completely controlled by attributes (so you'll only need to specify `solver_attributes`). Hence, a model is created as `Model(()-> HiGHS.Optimizer())` for the HiGHS solver, and analogously for the other ones.
     - For *Gurobi*: `solver_args=Gurobi.Env()` and the constructor is `Model(()-> eval(solver).Optimizer(solver_args))`. This is needed to [avoid Gurobi to retrieve your license](https://github.com/jump-dev/Gurobi.jl#reusing-the-same-gurobi-environment-for-multiple-solves) every time a model is created (*i.e.* for each iteration of CALiPPSO), and printing out an annoying line like `Academic license - for non-commercial use only - expires XXXXXXXXX`. If you don't care about this, you can also set `solver_args=nothing`.
     - *GLPK* and *Hypatia* passes the arguments as *keywords arguments*, an therefore `solver_args` should be a `NamedTuple`. For instance, for *GLPK*, `solver_args=(want_infeasibility_certificates=false, method=GLPK.MethodEnum(0) )` and a model is created like `Model(()->GLPK.Optimizer(;solver_args...))`. The same principle is applied for *Hypatia*. 
       - Note however that this is only because we didn't manage to fully control these solvers through attributes. If you know how to do so, [please let us know](mailto:rafael.diazhernandezrojas@uniroma1.it).


!!! warning "Using a solver not listed here"
    We provide working implementations for the set of solvers [mentioned above](@ref Solvers-we-tried). And while our code should work with several other available libraries and wrappers, keep in mind that there is *no* universal way to control the behaviour of all of them. That is, while JuMP allows to seemly change between them, there is still a considerable variety of ways to pass an optimizer's arguments.
    
    Therefore, in our code we had to resort to a (dirty) implementation such that, whenever `solver_args` is different to `nothing`, we use a series of `if-elseif-else` statements to properly set the optimizer of the *specific* solver. All of this means that if you want to use a solver not considered here and pass attributes to it when the model is created, **you'll need to specify the proper argument passing syntax**. To do so, you'll have to modify the `if-elseif-else` statements that defines `optimizer` in the definition of `solve_LP_instance` (lines 526-536) and `fine_tune_forces!` (lines 753-763), and adapt it to the solver of your choice.




