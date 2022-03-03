# [How `produce_jammed_configuration!` works](@id main_function)


Our main function consists of two essentially independent parts: (1) the [main CALiPPSO loop](@ref Theory-behind-CALiPPSO); and (2) the packing creation from the quantities obtained after the main loop converged. We now describe each of them.

---

## [The main CALiPPSO's loop (a.k.a. ILP)](@id mainloop)


From the initial particles' position and size (*i.e.* the input of [`produce_jammed_configuration!`](@ref)), a `while` loop is initialized until the *convergence criteria* defined [before](@ref Theory-behind-CALiPPSO) are reached. More precisely, the loop continues until: (1) ``\sqrt{\Gamma^\star}-1 <``` tol_Γ_convergence` **and** (2) ``|s^\star_{i,\mu}| <`` `tol_S_convergence` for ``i=1, \dots, N`` and ``\mu=1,\dots, d`` (although see step 4 below); **or** the number of iterations (*i.e* the number of LP optimizations) exceeds `max_iters`. The default values of these 3 quantities are [specified later](@ref list-defaults) and can be easily changed through [keyword arguments](@ref kwargs-control) of [`produce_jammed_configuration!`](@ref).


This main loop consists of the following steps:

1. The LP model creation and optimization. (Expectedly, this is done using [`JuMP`'s funcionality](https://jump.dev/JuMP.jl/stable/manual/models/#Create-a-model))
   1. Thus, given the the particles' position and radii, the linear optimization problem of Eqs. (2) [in the theory section](@ref Theory-behind-CALiPPSO) is defined using the JuMP's API and assigned to an object call `LP_model`.
      - `LP_model` includes the relevant design variables (*i.e.* the inflation factor, ``\Gamma``, and particles' displacements, ``\vec{\mathbf{s}}``), as well as the set of non-overlapping (*linear*) constraints. The constraints are added to `LP_model` using the [`add_non_overlapping_constraints!`](@ref CALiPPSO.add_non_overlapping_constraints!) function.
      - Importantly, not all pair of particles are considered for the constraints, but only those whose distance is smaller than a cut-off, ``\ell``, whose value is obtained by calling [`bounds_and_cutoff`](@ref CALiPPSO.bounds_and_cutoff). This function also outputs the value of `s_bound` to be imposed as constraint of the displacements magnitudes when solving the LP instance.
      - Besides, the periodic boundary conditions are automatically considered, using the so called *Minimum Image Convention*. That is, the vector differences, like ``\mathbf{r}_{ij}=\mathbf{r}_i - \mathbf{r}_j`` are always computed using the virtual image of the system that corresponds to the smallest value of ``|\mathbf{r}_{ij}|``. See the docstrings of [`MIC_vector`](@ref CALiPPSO.MIC_vector) and [`MIC_distance`](@ref CALiPPSO.MIC_distance) for more information.  

   2. The optimization is carried out simply by calling [`optimize!`](https://jump.dev/JuMP.jl/stable/reference/solutions/#JuMP.optimize!)`(LP_model)`. 
      - Provided the optimizer was able to solve the LP instance, at this point we have obtained the optimal displacements (``\vec{\mathbf{s}}^\star``) and inflation factor (``\Gamma^\star``). 
      - **Note that** both of these steps are implemented in a single function: [`solve_LP_instance`](@ref CALiPPSO.solve_LP_instance).
    
2. The force balance of the current packing is assessed. To do so, a *preliminary* network of contacts is constructed from the active constraints obtained in the previous step. 
    - To do so [`network_of_contacts`](@ref) is applied on the particles positions and list of constraints introduced in the step 1.1. This can be done because the list of constraints of each particle is stored as a `Vector{ConstraintRef}`. (See [here](https://jump.dev/JuMP.jl/stable/manual/constraints/#Constraint-containers) for more info about `ConstraintRef` in `JuMP`.) 
    - As we mentioned [above](@ref Contact-forces) and showed in our paper, even if the jamming point has not been reached, the dual variables should fulfill a force-balance type of equation. Thus, verifying that this is the case is a convenient way of assessing whether the optimal solution of the LP instance found is good or not.
    - **Note that** this check should be performed **before** the configuration is updated, otherwise the *wrong* contact vectors would be used.

3. The configuration is updated: ``\mathbf{r}_i \to \mathbf{r}_i + \mathbf{s}_i^\star`` and ``\sigma_i \to \sqrt{\Gamma^\star}\sigma_i`` for ``i=1,\dots, N``. 
   - These updated values will be used to formulate the next LP instance in the next iteration of the main loop. 

4. A set of *preliminary* stable particles is obtained using [`obtain_non_rattlers`](@ref CALiPPSO.obtain_non_rattlers). Rattlers are also obtained as the complement of such set.
   - This step is important in order to check if the configuration is isostatic or not. In the latter case, the isostaticity gap (*i.e.* the difference of the number of contacts, ``N_c``, and the number of degrees of freedom, ``N_{dof}``) may provide insight about numerical issues when determining the contact forces. Thus, even though this step is (apparently) not strictly required in order for CALiPPSO to work, it usually provides very useful information.
   - Besides, rattlers should be (almost always) excluded when testing convergence related to the magnitude of ``|s^\star_{i,\mu}|``. That is, because rattlers are not blocked by their neighbours, their associated optimal displacements are notably larger than those of the stable particles, and therefore we don't consider them for checking when the main loop should terminate. For instance, compare the value of `max |sᵢ|` of *all* particles with `bound on |sᵢ|` in the [example output](@ref output-process) of before. Note that when `max |sᵢ|` of *stable* particles is considered instead, a much smaller value is obtained
   - So, this step is needed *in practice* for the correct functioning of CALiPPSO. Otherwise the convergence criterion of ``|s^\star_{i,\mu}| <`` `tol_S_convergence` would never be met due to the presence of rattlers.

5. If the kwarg `verbose=true`, some information about the progress of CALiPPSO is printed out. This is explained in detail in [the dedicated section](@ref output-process).

6. Call the [`check_for_overlaps`](@ref) function to check if there are any overlaps once the configuration has been updated.
   - Of course, **there shouldn't be!**
   -  ... but given that we live in a world of *finite precision* and that we actually aim for a condition in which some of the **constraints are saturated**, it can happen that the LP instance was not solved within the required accuracy. See [this section](@ref Setting-the-required-precision) to learn how to control the overall precision of `CALiPPSO`, and [how to tune the options](@ref kwargs-control) for setting the tolerance with which an overlap is identified.
   -  When an overlap *does* occur, an error is thrown an `produce_jammed_configuration!` terminates, also terminating the main process since `error` is called. Nevertheless, some other information is shown, that can be used, hopefully, to trace back what happened.
   -  If you think that the problem is the related to [numerical issues](@ref Possible-issues), be sure to understand [how the precision of `produce_jammed_configuration!` is determined](@ref Setting-the-required-precision).
   -  Note also that a *real* overlap can also occur (*i.e*. once in which a pair of particles is overlapping by an amount much larger than the accuracy with which a solver fulfills the constraints). If this happens, try some of the solutions [mentioned here](@ref REAL-overlaps)
      
7. Check if convergence criteria are fulfilled. If this is the case, the main loop terminates. Otherwise, go back to step 1.



Note that steps 5 and 6, by default, are only performed during the first few iterations (10) and at given intervals (also 10). To change how often information about the main loop progress is printed out (respectively how often overlaps checks are performed) set the keyword argument `monitor_step` (respectively `interval_overlaps_check`) to the desired value. Instead, to select in how many initial iterations to include these steps, use `initial_monitor` (for printing info) and `initial_overlaps_check` for overlaps checks. More details can be found [here](@ref kwargs-control) and in the docstring of [`produce_jammed_configuration!`](@ref).

---

## Creating the final packing

Clearly, a lot of data is contained in a single packing, like the set of all particles position, the network of contacts, etc. Moreover, the information related to the algorithm itself (*e.g.* termination status, number of iterations, etc.). To efficiently store, access, and manipulate all of them, `CALiPPSO` relies on few [*composite types*](https://docs.julialang.org/en/v1/manual/types/#Composite-Types) or `struct`'s (aka *objects* in other languages). In the [types section](@ref types) we describe all of them in detail, but for the purposes of this section, the most important ones are [`MonoParticle`](@ref) and [`MonoPacking`](@ref). Very briefly:

- A `MonoParticle{d,T}` is assigned a position (as an [`SVector`](https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#SVector) of size `d` and with [`PeriodicNumber`](@ref) elements of type `T`: almost surely `Float64`), a set of neighbours, and the corresponding set of contact vectors and forces. (Note that the particle's size is *not* specified.)
- A `MonoPacking{d,T}` is composed of an array of ``N`` `MonoParticle{d,T}`, their common radius, ``R``; and also includes information about whether mechanical equilibrium holds and whether the packings is jammed or not.
  
Instead, the information about convergence time, number of LP iterations, etc. are stored in a [`convergence_info`](@ref) object.

Now, once the main CALiPPSO's loop has finished (possibly producing a jammed packing), the following steps are carried out:

1. It is assessed whether, (i) the process of the main loop converged (in which case we create a flag `jammed=true`); or (ii) if the loop ended because `max_iters` was exceeded or too many non-isostatic solutions were obtained consecutively (in which case `jammed=false`).

2. Using the values of `Xs` and the particles' size after the last LP optimization, as well as the relevant constraints, `final_packing` is created.
   - Clearly, this is an essential step. It is done by calling the *constructor* [`MonoPacking`](@ref MonoPacking(::Vector, ::Vector, ::Vector, ::Real)).
   - This [method](https://docs.julialang.org/en/v1/manual/methods/) of the `MonoPacking` function uses the set of constraints and the particles' position (as well as some other secondary arguments) and constructs the set of all `MonoParticle`s objects (each with a position, list of contacts, etc.). This set of `MonoParticle`'s is then assigned to a `MonoPacking`, along with the particles' radius and the value of `jammed`. When the packing is created, it is assessed whether it is in mechanical equilibrium or not.
   - See the docstring's of [`MonoPacking`](@ref MonoPacking(::Vector, ::Vector, ::Vector, ::Real)) for more info.

3. The isostaticity of `final_packing` is assessed calling [`is_isostatic`](@ref)`(final_packing)`.

4. The sum of forces on each particle in `final_packing` is computed. 
   - If any of them is greater than a  tolerance value (fixed by the kwarg [`tol_mechanical_equilibrium`](@ref list-defaults)), then [`fine_tune_forces!`](@ref CALiPPSO.fine_tune_forces!) is called. This function is very similar to [`solve_LP_instance`](@ref CALiPPSO.solve_LP_instance) described [above](@ref mainloop), but it also updates the state of `final_packing`. More precisely
     - An additional LP problem is created and solved. An important difference with `solve_LP_instance` is that in this LP problem ``|s_{i,\mu}|`` is *un*bounded.
     - The forces magnitudes and contact vectors of each particle in `final_packing` are updated by calling [`update_packing_forces!`](@ref CALiPPSO.update_packing_forces!). This function essentially uses [`network_of_contacts`](@ref) to construct all the *real* contacts from the constraints of this additional LP optimization.
     - *Note*: in this additional optimization **none of the positions _nor_ the radius are updated**.
     - The isostaticity of the updated `final_packing` is reassessed.
   - If each particle is in mechanical equilibrium, within the tolerance value, the algorithm jumps to the next step.

5. A final overlaps check is performed on `final_packing`. Note that this is done also by the [`check_for_overlaps`](@ref check_for_overlaps(::MonoPacking)), but *using the [method](https://docs.julialang.org/en/v1/manual/methods/) for* `MonoPacking` type.


