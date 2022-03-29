# Changing the default options

For ease of use, several default values have been defined in our code. Some of them are specific or related to the solvers we testes, while others specified convergence criteria, etc. In any case, we stress that all of them have been extensively tested *only with the Gurobi solver*, and when the CALiPPSO's [initial configuration was obtained after LS compression](@ref The-initial-conditions). Therefore, if you want to [use a different solver](@ref changing_the_solver) or initialize CALiPPSO from a different type of configuration it is likely that you'll need to make some small changes to these default parameters.

## [List of default values](@id list-defaults)

The main default values, all of them defined as *global* variables (using `const`), can be accessed by calling `default_parameters` once CALiPPSO has been loaded. `default_parameters` is a Julia's Dictionary [(*i.e.* a `Dict` type)](https://docs.julialang.org/en/v1/base/collections/#Base.Dict) whose values are `Dict`s themselves:

```@example
using CALiPPSO # hide
default_parameters
```

Thus, you can access, say, all the default values associated with the precision of CALiPPSO by calling
```@example
using CALiPPSO # hide
default_parameters["Precision parameters"]
```

while the default value that determines, *e.g.* the convergence criterion of ``\Gamma^\star`` can be obtained as:
```@example
using CALiPPSO # hide
default_parameters["Convergence parameters"]["default_tol_Γ_convergence"]
```

!!! warning
    Note that changing the values of `default_parameters`, or any of its entries, will **not** change the default behaviour of `produce_jammed_configuration!`. To do so, the associated kwarg should be specified when calling  this function, as explained [below](@ref kwargs-control)


For completeness, we also show here a list with all the default values of such global variables. Note however that these default values have been mostly tested with Gurobi, so you might need to pass different values to `produce_jammed_configuration!`, depending on which solver you choose to use; see the end of this section for more info.

### Parameters related to the precision; accessed through `default_parameters["Precision parameters"]`

| Variable      | default value |  Role  |
| :------------ | :-----------: |:--- | 
| `default_tol_overlap` | ``10^{-8}``         | Tolerance for identifying overlaps (*e.g.* in functions like `check_for_overlaps`) |
| `default_tol_optimality` | `0.1*default_tol_overlap` (``10^{-9}``)       | General precision of optimal solutions. This value determines different features for each solver. (See below)  |
| `default_tol_zero_forces` | `0.1*default_tol_optimality` (``10^{-10}``)   | If a force's magnitude is smaller than this value, it is considered 0. | 
| `default_tol_force_equilibrium` | ``10^{-12}`` | Force balance condition (per particle) should be satisfied within this precision. |

**Note:** The value of `default_tol_optimality` was chosen having in mind that you can use Gurobi or HiGHS solvers. But it also seems to work fine with GLPK. However, we have not performed extensive tests.

### Parameters related to convergence criteria; accessed through `default_parameters["Convergence parameters"]`

| Variable      | default value |  Role  |
| :------------ | :-----------: |:--- | 
| `default_tol_Γ_convergence` | ``10^{-12}``|  Tolerance for testing convergence of ``\sqrt{\Gamma^*}-1`` |
| `default_tol_displacements_convergence` | ``10^{-9}`` | Tolerance for testing convergence of ``s_{i,\mu}^\star``, restricted to *non*-rattlers. | 
| `default_max_iterations`  | ``1000`` | Maximum number of LP optimizations performed in the [main loop](@ref mainloop) before it's terminated.|



### [Parameters related to the default solver and its behaviour; accessed through `default_parameters["Solver parameters"]`](@id def_solv_params)

| Variable      | default value |  Role  |
| :------------ | :-----------: |:--- | 
| `default_solver` | `GLPK` | Solver used by `produce_jammed_configuration!` |
| `default_solver_args` | `(want_infeasibility_certificates = false, method = GLPK.SIMPLEX)` | Specify method used by GLPK; avoid infeasibility certificates to improve speed|
| `default_solver_attributes` | `Dict("tol_dj"  => default_tol_optimality, "msg_lev" => 0, "tol_bnd" => default_tol_optimality)`| Precision of GLPK solver; avoid printing output. See GLPK's documentation for more info.|



### Other parameters defined; accessed through `default_parameters["Other parameters"]`

| Variable      | default value |  Role  |
| :------------ | :-----------: |:--- | 
| `default_threads` | `Int(round(Sys.CPU_THREADS/2))` (*i.e.* half of the max threads available) | Number of threads to be used by solvers that allow parallelization | 




!!! tip "Default values for Gurobi"
    If you have access to a Gurobi license and want to use it as solver, we suggest using the following values when calling `produce_jammed_configuration!`
    ```julia
    const default_tol_optimality = CALiPPSO.default_tol_optimality
    const max_threads = CALiPPSO.max_threads
    const grb_args = Gurobi.Env()
    const grb_attributes = Dict("OutputFlag" => 0, "FeasibilityTol" => default_tol_optimality, "OptimalityTol" => default_tol_optimality, "Method" => 3, "Threads" => max_threads)
    
    precompile_main_function(Gurobi, grb_attributes, grb_args)
    produce_jammed_configuration!(Xs0, r0; solver=Gurobi, solver_args=grb_args, solver_attributes=grb_attributes, <other kwargs>)
    ```
    
    In fact, the value of `default_tol_optimality` reported above was chosen based on highest precision allowed by Gurobi (both for 'OptimalityTol' and 'FeasibilityTol'); go [here](https://www.gurobi.com/documentation/9.1/refman/optimalitytol.html) for more information. The value of `default_tol_overlap` was then chosen accordingly; that is, several times (10) larger, because the optimal value of each degree of freedom is determined with an accuracy of `default_tol_optimality`.

    On the other hand, the value that determines the convergence criterion of ``\Gamma`` is clearly smaller than such precision allowed by Gurobi. Nevertheless, in our tests we observed that after enough iterations this more stringent condition can actually be met. In fact, achieving ``\sqrt{\Gamma^\star} -1 \leq 10^{-12}`` is relatively simple. The more complicated part is to reach a configuration with ``\max |\mathbf{s}_{i,\mu}^\star|_{i=1,\dots,N}^{\mu=1,\dots, d} \leq 10^{-9}`` . Specially for relatively large systems, *e.g.*, ``dN\geq 5000``.

    Thus, when dealing with very large systems (about ``dN>30,000``) you might want to try something like
    ```julia
    const tol_Γ_convergence = default_tol_optimality
    const tol_displacements = 10*default_tol_optimality

    produce_jammed_configuration!(Xs0, r0; solver=Gurobi, solver_args=grb_args, solver_attributes=grb_attributes, tol_Γ_convergence=tol_Γ_convergence, tol_S_convergence=tol_displacements, <other kwargs>)
    ```
    in order to speed up convergence, at the cost of a rather smaller precision.

    In any case, in all the configurations we tested, we did *not* find anyone for which the force balance condition couldn't be met within the tolerance defined by `default_tol_force_equilibrium`, despite this value being much smaller than Gurobi's highest accuracy.


!!! tip "Default values for Other Solvers"
    We do not list here the analogous values and variables for the [other solvers](@ref changing_the_solver) we [tested](@ref Testing-different-solvers). However, you can find some of the ones we found useful in the first few lines of the `CALiPPSO.jl` file and in the exampled described in the [Testing different solvers](@ref) section. Note that if you want to modify the source code of `CALiPPSO` so that other solver is used by default, uncommenting the relevant lines would be useful.

    Note also that the default values defined in the file `CALiPPSO.jl` for the `Hypathia` and `COSMO` solvers do not really produce very good results. If you find a better choice of arguments, [*please let us know*](mailto:rafael.diazhernandezrojas@uniroma1.it).


---
## [Controlling `produce_jammed_configuration!` with keyword arguments (a.k.a. kwargs)](@id kwargs-control)

The full list of keyword arguments (kwargs) of [`produce_jammed_configuration!`](@ref) can be readily accessed from its docstring. Here we provide the same list (with their default values, defined above), and a more detailed description when needed. Thus, the value of any of them can be conveniently tunned to your needs by calling `produce_jammed_configuration!(Xs0, R, L; kwarg=<your chosen value>)`.

### Kwargs for controlling how constraints are assigned to `LP_model` and setting bounds on displacements

As explained in our paper, the constraints a particle is subject to are assigned according to a neighbours-list approach. Thus, the associated neighbours o a given sphere are all the particles within a certain distance ``\ell``. In other words, ``\ell`` determines the radius of influence on a particle. This and other quantities are computed using [`bounds_and_cutoff`](@ref CALiPPSO.bounds_and_cutoff). This function is called within the [`solve_LP_instance`](@ref CALiPPSO.solve_LP_instance) and [`fine_tune_forces!`](@ref CALiPPSO.fine_tune_forces!) functions

1. `ℓ0=4*R`: Initial value of the radius of influence (``\ell``) for assigning constraints. `ℓ0` is also used as upper-bound for ``\ell`` in subsequent steps. See item 4 below for more information.
2. `sqrΓ0=1.01`: Initialization value of ``\sqrt{\Gamma}``; it is used to provide a guess of the value of ``s_{bound}`` --*i.e.* and upper bound of ``|s_{i,\mu}|``-- when ``\Gamma`` is relatively large. See item 4 below for more information.
3. `sbound=0.01`: Fraction of ``R`` used as bound for ``|s_{i,\mu}|`` during the last optimizations, *i.e.* when ``\Gamma^\star`` has very small values. 
4. `thresholds_bounds=(5e-4, 1e-5)`: These two value control when the different behaviours of `bounds_and_cutoff` are triggered. Thus, 
   1. when ``\sqrt{\Gamma^\star}-1`` (or `sqrΓ-1` in the main function definition) is larger than the first value, ``\ell=```ℓ0` and ``s_{bound}=\frac{ (\ell/2 - \sqrt{\Gamma_0^\star}\sigma)}{\sqrt{d}}``, where ``\Gamma_0^\star`` is the optimal value of ``\Gamma`` from the previous iteration, or `sqrΓ0` squared in the first one.
   1. when ``\sqrt{\Gamma^\star}-1`` is between the two values, ``\ell=```min(4R, ℓ0)` and ``s_{bound}=0.1R``
   2. when ``\sqrt{\Gamma^\star}-1`` is smaller than the second value, ``\ell=```min(2.7*R, ℓ0)` and ``s_{bound}=`` `sbound*` ``R``



### Kwargs for controlling convergence and termination criteria of [the main loop](@ref mainloop)

1. `tol_Γ_convergence=default_tol_Γ_convergence` (``10^{-12}``): Determines the value below which ``\sqrt{\Gamma^\star}-1`` is considered zero (so convergence in the inflation factor has been reached).
2. `tol_S_convergence=default_tol_displacements_convergence` (``10^{-9}``): Determines the value below which ``|s_{i,\mu}^\star|`` is considered zero (so convergence in particles displacements ─restricted to *non*-rattlers─ has been reached).
3. `max_iters=default_max_iterations` (``1000``): Maximum number of iterations (*i.e.* LP optimizations) allowed before stopping the main CALiPPSO's loop.
4. `non_iso_break=max_iters`: Maximum number of non-isostatic configurations that can be obtained *consecutively* before the main CALiPPSO's loop is terminated. Only to be able to stop the main function when dealing with strange, atypical cases (*e.g.* a configuration for which CALiPPSO repeatedly fails to converge because a stable particle is surrounded by several rattlers); in general it shouldn't be necessary to tune this parameter.

### Kwargs for controlling precision of overlaps, testing force balance, and updating lists of distances

1. `tol_mechanical_equilibrium=default_tol_force_equilibrium` (``10^{-12}``): When the norm of the total force acting on a particle is smaller than this quantity, said particle is considered to be in equilibrium.
2. `zero_force=default_tol_zero_forces` (``10^{-10}``): This is the threshold for determining when a force, or dual variable, is different from zero (*active*). In other words, whenever [`shadow_price`](https://jump.dev/JuMP.jl/stable/reference/solutions/#JuMP.shadow_price)`(constraint)` --with `constraint` being a [`ConstraintRef`](https://jump.dev/JuMP.jl/stable/reference/constraints/#JuMP.ConstraintRef)-- outputs a value larger than `zero_force`, we consider that such constraint is active, and its value is precisely `shadow_price(constraint)`. 
3. `tol_overlap=default_tol_overlap` (``10^{-8}``): Maximum value of overlap that can occur between particles. That is, if ``\sigma_{ij} - |\mathbf{r}_{ij}| \geq `` `tol_overlap`, an overlap is said to have occurred.
4. `initial_overlaps_check=initial_monitor` (``10``): During each of these many *initial* iterations [`check_for_overlaps`](@ref) is called, after the configuration has been updated following an LP optimization. (`initial_monitor` is described below.)
5. `interval_overlaps_check=10`: [`check_for_overlaps`](@ref) is also called every `interval_overlaps_check` iterations, after the configuration has been updated.
6. `ratio_sℓ_update::T=0.1`: This is the fraction of the cutoff distance ``\ell`` which is allowed for each particle to move, before updating the list of its distances with the rest of spheres. The threshold for the displacement, `s_update`, is determined as `s_update=ratio_sℓ_update*ℓ/sqrt(d)`. Setting `s_update=0.0` is equivalent to updating the lists of distances after each LP optmization. Of course, this hinders performance for large samples, but in such way it can be guaranteed that all the relevant constraints are considered. See [`update_distances!`](@ref CALiPPSO.update_distances!) for more information about how the update of the lists of distances is implemented.

### Kwargs for controlling output printing on terminal

1. `verbose=true`: turns on/off the printing of information during the main CALiPPSO's loop.
2. `monitor_step=10`: The info about the progress of CALiPPSO is printed out after these many steps (besides other criteria).
3. `initial_monitor=monitor_step`: `verbose` is set to true for these many *initial* iterations.


In our experience, most of the problems (*e.g.* a real overlap or an optimization that didn't return `OPTIMAL` as termination status) occurred during the initial steps of CALiPPSO. Thus, obtaining as much information as possible, as well as being extra-cautious by verifying that no overlaps are present after *each* iteration, becomes important during such initial stage.

---

## [Changing the solver/optimizer](@id changing_the_solver)

### Solvers we tried
Given that we have used JuMP's interface for solving each LP instance, *in principle*, you can use *any* of the [JuMP compatible solvers](https://jump.dev/JuMP.jl/stable/installation/). Of course, maybe not all of them are suitable for the type of LP optimization our algorithm requires, but there should be ample choice. Indeed, besides GLPK, we tested the following solvers:


1. [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl): the Julia wrapper of the [Gurobi solver](https://www.gurobi.com/). We used this solver in most of [our tests](@ref Tests-and-examples-provided) and to produce all the figures of [our paper](https://arxiv.org/abs/2203.05654). We also observed that it is the fastest and most precise of the others we tried, so we suggest running `CALiPPSO` with it (as shown in  the green box above). However, you'll need to install `Gurobi` manually on your computer, as described [before](#Dependencies); keep in mind there are free academic licenses.
2. [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl): the Julia wrapper of [HiGHS](https://www.maths.ed.ac.uk/hall/HiGHS/). This is also a very performant and precise solver (yet noticeably slower than Gurobi).
3. [GLPK.jl](https://github.com/jump-dev/GLPK.jl): the Julia wrapper of the famous [GNU Linear Programming Kit](https://www.gnu.org/software/glpk/). GLPK is a very well known and amply used library, so it can be used reliably. This is the default solver of `CALiPPSO`, so it's installed when this package is added.
4. [Clp.jl](https://github.com/jump-dev/Clp.jl): the Julia wrapper of the [COIN-OR Linear Programming Interface](https://projects.coin-or.org/Clp). This solver also works very well, although it is the slowest of these options.

All these solvers produced very good results in the sense that all contacts, rattlers, etc. were correctly identified and consistent between them; no overlaps occurred, and isostaticity was always achieved. Besides, except for Gurobi, all of them are free and open-source, and are installed automatically when their wrapper is installed within julia. That is, you can use `Pkg.add("HiGHS")` or `Pkg.add("Clp")` from Julia, *without* the need of installing any other package or library beforehand.

We also tested some *pure-Julia* solvers, like [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl) and [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl). We were able to execute our code with both of them without any error, although in several cases the final packing was non-isostatic due to a lack of precision in identifying the contact forces. This issue is possibly related to our poor choice of the [solver parameters](@ref Selecting-a-solver-and-specifying-its-attributes), but we didn't perform any additional tests. (So if you know how to improve this [*please contact us*](mailto:rafael.diazhernandezrojas@uniroma1.it).)

### Selecting a solver and specifying its attributes

Just as with the other options of [`produce_jammed_configuration!`](@ref), choosing a solver is conveniently done through [keyword arguments](@ref kwargs-control):

1. `solver::Module=default_solver` (`GLPK`): This kwarg must correspond to the module of a solver, *i.e.* just its name in most cases. For example, as [mentioned above](@ref list-defaults), we defined `default_solver = GLPK`; or if you want to use, *e.g.*, the HiGHS solver, you should pass `solver=HiGHS` as argument. Thus, other possible values of the `solver` kwarg are: `Gurobi`, `Clp`, `Hypatia`, etc.
   - If you want to use a different solver, be sure to load the corresponding package before calling the main function. We also suggest that you run `precompile_main_function` with the solver of your choice before calling [`produce_jammed_configuration!`](@ref). See the examples in the [basic usage section](@ref Basic-usage)
   - Note that `solver` is not used by the main function itself, but instead is passed, also as a kwarg, to the functions where the optimization is actually performed, *i.e.* [`solve_LP_instance`](@ref CALiPPSO.solve_LP_instance) and [`fine_tune_forces!`](@ref CALiPPSO.fine_tune_forces!). 
2. `solver_attributes::Dict=default_solver_attributes`: These are the set of options or parameters that control the behaviour of your solver. It should be passed as a `Dict` type (see [here](@ref def_solv_params) for the definition of `default_solver_attributes`). The idea is that this `Dict` should contain any of the parameters or options that can be set using the function [`set_optimizer_attributes`](https://jump.dev/JuMP.jl/stable/reference/models/#JuMP.set_optimizer_attributes) or other similar functions, as [explained here](https://jump.dev/JuMP.jl/stable/reference/models/#Working-with-attributes).
   - Thus, for instance, with these parameters you can specify the accuracy of the solver, the method to use, iterations limit, etc.
   - In the first lines of the file `CALiPPSO.jl` we have included some possible choices for each of [the solvers we tried](@ref Solvers-we-tried). All these lines, except the ones for `GLPK` have been commented out.
3. `solver_args=default_args`: These are supposed to be *arguments* (and *not* attributes) that are passed as arguments to [`Optimizer()`](https://jump.dev/JuMP.jl/stable/moi/tutorials/implementing/#The-Optimizer-object) of the chosen solver. Note that the behaviour of *each solver* is different. Thus, **if you are testing a different solver to the ones we mention here, we suggest you first set** `solver_args=nothing` to be sure the function will execute properly. In fact, if `solver_args===nothing` results in `false` and you're using a solver that has not been preconfigured, an error will be thrown. Of course, you can modify this, as explained at the end.
   - Creating a model with arguments in JuMP has a rather strange syntax: `Model(()-> solver.Optimizer(solver_args))`

   - However, given that *each* solver has different ways of passing arguments to its optimizer, we had to resort to a somewhat silly implementation for each of the solvers:
     - For the *HiGHS*, *Clp*, and *COSMO* solvers: `solver_args=nothing` because they are completely controlled by attributes (so you'll only need to specify `solver_attributes`). Hence, a model is created as `Model(()-> HiGHS.Optimizer())` for the HiGHS solver, and analogously for the other ones.
     - For *Gurobi*: `solver_args=Gurobi.Env()` and the constructor is `Model(()-> Gurobi.Optimizer(solver_args))`. This is needed to [avoid Gurobi to retrieve your license](https://github.com/jump-dev/Gurobi.jl#reusing-the-same-gurobi-environment-for-multiple-solves) every time a model is created (*i.e.* for each iteration of CALiPPSO), and printing out an annoying line like `Academic license - for non-commercial use only - expires XX-XX-XX`. If you don't care about this, you can also set `solver_args=nothing`.
     - *GLPK* and *Hypatia* pass the arguments as *keywords arguments*, an therefore `solver_args` should be a `NamedTuple`. For instance, for *GLPK*, `solver_args=(want_infeasibility_certificates=false, method=GLPK.MethodEnum(0) )` and a model is created like `Model(()->GLPK.Optimizer(;solver_args...))`. The same principle applies for *Hypatia*. 
       - Note however that this is only because we didn't manage to fully control these solvers through attributes. If you know how to do so, [please let us know](mailto:rafael.diazhernandezrojas@uniroma1.it).


!!! warning "Using a solver not listed here"
    We provide working implementations for the set of solvers [mentioned above](@ref Solvers-we-tried). And while our code should work with several other available libraries and wrappers, keep in mind that there is *no* universal way to control the behaviour of all of them. That is, while JuMP allows to seemly change between them, there is still a considerable variety of ways to pass an optimizer's arguments to tune its properties.
    
    Therefore, in our code we had to resort to a (dirty) implementation such that, whenever `solver_args` is different to `nothing`, we use a series of `if-elseif-else` statements to properly set the optimizer of the *specific* solver. All of this means that if you want to use a solver not considered here and pass attributes to it when the model is created, **you'll need to specify the proper argument passing syntax**. To do so, you'll have to modify the `if-elseif-else` statements that defines `optimizer` in the definition of `solve_LP_instance` (lines 526-536) and `fine_tune_forces!` (lines 753-763), and adapt it to the solver of your choice.




