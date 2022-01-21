# [How `produce_jammed_configuration` works](@id main_function)


Our main function consists of two essentially independent parts: (1) the main ILP loop; and (2) the packing creation from the quantities obtained after using ILP. We now describe each of them.

---

### [The main Linear Optimization loop](@id mainloop)


From the initial particles' position and size (*i.e.* the input of `produce_jammed_configuration`), a `while` loop is initialized until the *convergence criteria* defined [above](#introduction-how-ilp-works-and-some-terminology) are reached. More precisely, the loop continues as long as: (1) ``\sqrt{\Gamma^\star}-1 \geq``` tol_Γ_convergence`; or (2) ``|s^\star_{i,\mu}| \geq`` `tol_S_convergence` for ``i=1, \dots, N`` and ``\mu=1,\dots, d`` (although see step 4 below); or the number of iterations (*i.e* the number of LP optimizations) does not exceed `max_iters`. The default values of these 3 quantities are [given below](#list-of-default-values) and can be easily changes as [keyword arguments](#controlling-produce_jammed_configuration-with-keyword-arguments) of `produce_jammed_configuration`.


This main loop consists of the following steps:

1. The LP model creation and optimization. (Expectedly, this is done using `JuMP`)
   1. Thus, given the the particles' position and radii, the linear optimization problem of Eqs. (1-2) above is defined using the JuMP's interface and assigned to an object call `LP_model`.
      - `LP_model` includes the relevant design variables (*i.e.* the inflation factor, ``\Gamma``, and particles' displacements, ``\vec{\mathbf{s}}``), as well as the set of non-overlapping (*linear*) constraints (Eq. (2) above). The constraints are added to `LP_model` using the `add_non_overlapping_constraints!` function.
      - Importantly, not all pair of particles are considered for the constraints, but only those whose distance is smaller than a cut-off, ``\ell``, whose value is obtained by calling `bounds_and_cutoff`.
      - Besides, the periodic boundary conditions are automatically considered, using the so called *Minimum Image Convention*. That is, the vector differences, like ``\mathbf{r}_{ij}=\mathbf{r}_i - \mathbf{r}_j`` are always computed using the virtual image of the system that corresponds to the smallest value of ``|\mathbf{r}_{ij}|``. See the docstring of `MIC_vector` and `MIC_distance` for more information.  

   2. The optimization is carried out simply by calling `optimize!(LP_model)`. 
      - Provided the optimizer was able to solve the LP instance, at this point we have obtained the optimal displacements (``\vec{\mathbf{s}}^\star``) and inflation factor (``\Gamma^\star``). 
      - **Note that** both of these steps are implemented in a single function: `solve_LP_instance`.
    
2. The force balance of the current packing is assessed. To do so, a *preliminary* network of contacts is constructed from the active constraints obtained in the previous step. 
    - To do so `network_of_contacts` is applied on the particles positions and list of constraints introduced in the step 1.1. This can be done because the list of constraints of each particle is stored as a `Vector{ConstraintRef}`. (See [here](https://jump.dev/JuMP.jl/stable/manual/constraints/#Constraint-containers) for more info about `ConstraintRef` in `JuMP`.) 
    - Note that, as we mentioned [above](#introduction-how-ilp-works-and-some-terminology) and showed in our paper, even if the jamming point has not been reached, the dual variables should fulfill a force-balance equation. Thus, verifying that this is the case is a convenient way of assessing whether the optimal solution found is good or not.
    - Note that this check should be performed **before** the configuration is updated, otherwise the *wrong* contact vectors would be used.

3. The configuration is updated: ``\mathbf{r}_i \to \mathbf{r}_i + \mathbf{s}_i^\star`` and ``\sigma_i \to \sqrt{\Gamma^\star}\sigma_i`` for ``i=1,\dots, N``. 
   - These updated values will be used to formulate the next LP instance in the next iteration of ILP. 

4. A set of *preliminary* stable particles is obtained using `obtain_non_rattlers`. Rattlers are also obtained as the complement of such set.
   - This step is important in order to check if the configuration is isostatic or not. In the latter case, the isostaticity gap (*i.e.* the difference of the number of contacts, ``N_c``, and the number of degrees of freedom, ``N_{dof}``) may provide insight about numerical issues when determining the contact forces. Thus, even though this step is (apparently) not strictly required in order for ILP to work, it always provides very useful information. See [examples below](#a-simple-example-and-understanding-the-output).
   - Besides, rattlers should be (almost always) excluded when testing convergence related to the magnitude of ``|s^\star_{i,\mu}|``. That is, because rattlers are not blocked by their neighbours, their associated optimal displacements are notably larger than those of the stable particles, and therefore we don't consider them for checking when the main loop should terminate. So, this step is needed *in practice* for the correct functioning of ILP.

5. If the `verbose` option is set to `true`, the following information is printed (see [the dedicated section](#a-simple-example-and-understanding-the-output) on how to interpret such info):
   1. The number of LP iteration and the state of the optimization thrown after `optimize!(LP_model)` finishes, *i.e.* the value of `termination_status(LP_model)`. This latter could be, *e.g.* `OPTIMAL`, `INFEASIBLE`, `TIME_LIMIT` , etc.
   2. The value of ``\sqrt{\Gamma^\star} -1`` and ``\max |s^\star_{i,\mu}|``, where ``i \in \text{non-rattlers}``.
   3. The time required to solve the LP problem.
   4. Few statistics about the number of constraints included in `LP_model`
   5. Information about whether the preliminary configuration is isostatic or not, as well as the isostaticity gap, ``N_c - N_{dof}``, in the latter case.
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
      -  When this occurs, it is likely that the value of `sbound` is too large with respect of the value of ``\ell`` used. So try using a larger value of `ℓ0` (if this occurred when ``\Gamma^\star`` was still relatively large).
      -  Try using a smaller value of `sbound`.
      -  Or try redefining the bounds with which ``\ell`` is adjusted as ILP progresses.
      -  These three options can be set as keywords arguments of `produce_jammed_configuration` as explained [below](#changing-the-default-options). See the docstring of `bounds_and_cutoff` for more info.
7. Check if convergence criteria are fulfilled. Otherwise, go back to step 1.





Note that steps 5 and 6, by default, are only performed during the first few iterations (10) and at given intervals (also 10). To change how often the ILP info is printed out (respectively how often overlaps checks are performed) set the keyword argument `monitor_step` (respectively `interval_overlaps_check`) to the desired value. Instead, to select in how many initial iterations to include these steps, use `initial_monitor` (for printing info) and `initial_overlaps_check` for overlaps checks. More details can be found [here](#kwargs-for-controlling-output-printing-on-terminal) and [here](#kwargs-for-controlling-precision-of-overlaps-and-force-balance-tests).

---

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
     - An additional LP problem is created and optimized. An important difference with `solve_LP_instance` is that in this LP problem ``|s_{i,\mu}|`` is *un*bounded.
     - The forces magnitudes and contact vectors of each particle in `final_packing` are updated by calling `update_packing_forces!`. This function essentially uses `network_of_contacts` (described above) to construct all the real contacts from the constraints of this additional LP optimization.
     - *Nb*: in this additional optimization **none of the positions are updated**.
     - The isostaticity of the updated `final_packing` is reassessed.
   - If each particle is in mechanical equilibrium, within the same tolerance value, the algorithm jumps to the next step.

5. A final overlaps check is performed on `final_packing`. Note that this is done also by the `check_for_overlaps` function, but *using the method for* `MonoPacking` type.


