# Basic usage

## The main function: `produce_jammed_configuration`

We tried to make this package as easy to use as possible and, indeed, it consists of *a single* main function: [`produce_jammed_configuration`](@ref). (But before using it, be sure to have installed all the [dependencies](@ref Dependencies).) This function is defined in the `CALiPPSO` module, contained in the `iLP-for-jamming.jl` file. This means that this file should be loaded (through the `include` function) in any script making use of such function. In other words, to include `produce_jammed_configuration` in your scope or Julia script, you need to add the following lines to your script, REPL, Jupyter Notebook, etc.:

```julia
include("src/iLP-for-jamming.jl")
using .CALiPPSO  
```
(Beware of the dot (.) when [importing](https://docs.julialang.org/en/v1/manual/code-loading/) the `CALiPPSO` module.) 

In this way, [`produce_jammed_configuration`](@ref) as well as [other functions and `struct`'s will be loaded](@ref Main-Exported-Functions) into your environment. The additional functions that are exported when loading `CALiPPSO` (such as [`network_of_contacts`](@ref), [`check_for_overlaps`](@ref), [`get_non_rattlers`](@ref), etc.) are *not* needed for a basic usage, but might be useful for analysing the packings produced by the main function. 

Once `CALiPPSO` has been loaded, you just need to call
```julia
packing, info, Γ_vs_t, isostatic_vs_t = produce_jammed_configuration(Xs, R, L)
```
Naturally, this function should generate a jammed packing from an initial (monodisperse) configuration of particles with positions `Xs`, radius `R`, and contained in a periodic (hyper-) cube of size `L` (if this argument is left unspecified, it's assumed its value is 1). `Xs` should be a `Matrix{Float64}` of size ``d\times N``, thus specifying the position of each particle (*i.e.* each of the ``N`` columns is the ``d``-dimensional position vector of a particle). Clearly, this matrix can be constructed from importing data from a `csv` or `dat` file (or any other suitable file format for that matter).

Alternatively, `Xs` might be an array of ``N`` `StaticVector`s, each of size ``d`` and elements of [`PeriodicNumber`](@ref) type [(see here for the dedicated section on types)](@ref types). In such case, `L` should be given as an argument when calling [`produce_jammed_configuration`](@ref), because its value  automatically inferred from the corresponding field of `PeriodicNumber`. For convenience, if `Xs` is a matrix of `Float64` elements (*i.e.* a `Matrix{Float64}` type), it can be easily converted to an array of `StaticVector`s with `PeriodicNumber` elements using the `PeriodicVectors` function. Thus, the analogous version of the code above reads
```julia
Xs = PeriodicVectors(Xs, L)
packing, info, Γ_vs_t, isostatic_vs_t = produce_jammed_configuration(Xs, R)
```
Note that here we are making use of the (fantastic) [method dispatch feature](https://docs.julialang.org/en/v1/manual/methods/#Methods) of Julia.

!!! note "Precompilation"
    Note that executing `include("src/iLP-for-jamming.jl")` also *pre-compiles* `produce_jammed_configuration` by calling it using a small, predefined system. This causes that whenever this file is included, several things will be printed in screen. You can safely ignore all of them. It also makes loading `CALiPPSO` somewhat time consuming.
    
    In any case, precompilation can be avoided (hence making loading `CALiPPSO` much faster) by including a line 
    ```julia
    const precompile_main_function=false
    ``` 
    *__before__* loading the `iLP-for-jamming.jl` file. Nevertheless, **keep in mind** that if you do *not* define `precompile_main_function`, all the functions are precompiled by default.
    
    Avoiding precompilation is discouraged if you actually want to use such function because model creation and optimization in JuMP suffers from a ["time-to-first-solve" issue](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/performance_tips/#The-%22time-to-first-solve%22-issue), which is an analogous version of the "time-to-first-plot" one. Essentially, when a function is first called, it needs to be compiled. (This is actually the general behaviour of Julia, not only of JuMP.) Thus, by calling `produce_jammed_configuration` in a small system, the function gets pre-compiled and ready to be used in much larger systems. If `precompile_main_function` is set to `false` and then you call `produce_jammed_configuration` directly into the (presumably large) configuration you want to jam it could take much longer.

    On the other hand, avoiding precompilation is useful if you just want to import the types, functions, etc. of `CALiPPSO`, but you will not be making use of `produce_jammed_configuration`; for instance if you are only doing data analysis.


Using the variables introduced in the snippets above, the output of `produce_jammed_configuration` is the following:
1. `packing`: A jammed packing (provided convergence was attained) stored as a [`MonoPacking`](@ref) object (actually a [`struct`](https://docs.julialang.org/en/v1/manual/types/#Composite-Types)). 
   + This object contains an array of the ``N`` particles in the packing.
     + Each particle is stored as a [`MonoParticle`](@ref) object that contains: (i) the position of the centre; (ii) the list of all contact vectors; (iii) the list of contact forces magnitudes; and (iv) the list of neighbours in contact.
   + The radius of all the (hyper-) spheres in the packing.
   + Information of whether (i) the packing satisfies mechanical equilibrium (stored as a boolean in the `mechanical_equilibrium` field), and (ii) whether the packings reached jamming or not (also as a boolean specified using the `jammed` field).
2. `info`: Information about the process and termination status of CALiPPSO, *e.g.*, number of iterations, the time and amount of memory allocated during the full process, list of times of each LP optimization, etc. All of this is stored using a [`convergence_info`](@ref) object.
3. `Γ_vs_t`: The list of values of ``\sqrt{\Gamma^\star}`` obtained after each iteration; see [the theory section](@ref Theory-behind-CALiPPSO) for more information.
4. `isostatic_vs_t`: An analogous list that specifies (with boolean variables) if isostaticity holds at after each iteration.


With its [default parameters](@ref List-of-default-values) and assuming you're using Gurobi, [`produce_jammed_configuration`](@ref) should work in most of the cases with ``d>2`` (at least we didn't experienced any error in the [tests](@ref Some-examples-included) we ran), specially when ``p>1000`` (*i.e.* for high pressures, or ``\varphi_0\lesssim \varphi_J``). 

!!! tip "Initial condition at low pressure"
    If the input configuration is not sufficiently compressed, a simple fix is to increase the cutoff distance used to build the neighbours-list, ``\ell``. This parameter is conveniently fixed using [a keyword argument](@ref kwargs-control), namely, `ℓ0`. For instance, try
    ```julia
    produce_jammed_configuration(Xs, R, L; ℓ0=0.2*L)
    ```
    Of course, you should choose the initial value of ``\ell``, here specified by `ℓ0`, taking into account the initial size of your particles and initial density of the system. 
    See section about [The initial conditions](@ref) for more information about how to choose the initial configuration of CALiPPSO.


!!! warning "Check the packing produced by `produce_jammed_configuration`"
    In general situations, however, be warned that even if [`produce_jammed_configuration`](@ref) terminates without throwing any error, it may happen that the output is *not* a jammed packing. For instance, it might be the case that the function terminates because the maximal number of iterations was reached, or too many consecutive non-isostatic configurations were obtained (you can tune these parameters and several others using keyword arguments, as we explain in detail in the [dedicated section](@ref kwargs-control)). 
    
    This problematic situation is much likelier when ``\varphi_0 \ll \varphi_J``. We also observe that it may also occur using [other solvers](@ref changing_the_solver) instead of Gurobi.

    Therefore, be sure to verify that the `mechanical_equilibrium` and `jammed` fields of the [`MonoPacking`](@ref) obtained are both `true`.


Finally, [see this section](@ref Tests-and-examples-provided) for some examples/tests of its usage and [the section on keyword arguments](@ref kwargs-control) to fine tune the behaviour of `produce_jammed_configuration` to better suit your needs.



## Understanding the screen/IO printed output

What follows assumes that the optional argument `verbose` is set to `true` when calling `produce_jammed_configuration`; this is the default behaviour. Besides, if Julia is executed with the `color` flag and the IO supports it, the screen output is colorized and formatted to improve readability.

### [Output while `produce_jammed_configuration` is being executed](@id output-process)

When [`produce_jammed_configuration`](@ref) is running every few iterations an output as the following one will be printed in the IO. (See [here](@ref Kwargs-for-controlling-output-printing-on-terminal) for information about how often such output appears and how to change it.)

```
L1: Iteration: 3
L2: LP instance generated with (ℓ, ℓ/R) = [0.056951385144839745, 2.7]	 and bound on (|s_i|, |s_i|/R) = [0.00021093105609199902, 0.01]
L3: Status: OPTIMAL; 	 time solve LP instance (s) = 246.3279
L4: 	 Optimal values: 	 √Γ-1 = 2.4037438706159264e-11,	 max |sᵢ| (all particles) = 0.0002109310561
L5: 	 Max |sᵢ| (stable particles) = 1.0872115e-6
L6: New  values of (φ, R) = [0.644065639, 0.021093106]
L7: (Max, mean±std) constraints per particle per d: [15.0, 6.26, 3.93]	 Isostatic: true	 (Mean, Max, total in rattlers) z = [5.99975, 10.0, 0.0]	Non_rattlers= 16014
L8: ----------------------------------
L9: Force mismatch = 1.2233477227728374e-13
L10: Sample of smallest forces = 	[7.760517298819059e-7, 2.234495573171305e-6, 2.6928245176732605e-6, 2.7135404241734017e-6, 2.9159812742335627e-6, 3.1959098781771576e-6, 3.49862463912623e-6, 4.6744303036469655e-6, 7.2161329676828645e-6, 1.0484921226842317e-5]
```
**Note:** The line numbers (on the left) are only included here for reference, but are not printed.

The contents are:
+ Line 1: The number of CALiPPSO iteration (i.e. how many LP optimizations have been carried out).
+ Lines 2-4: Some information about how the LP instance (Eq. (2) in the [Theory-behind-CALiPPSO](@ref) section) was generated. This output is actually produced when calling [`solve_LP_instance`](@ref Main.CALiPPSO.solve_LP_instance) in the [main loop](@ref mainloop) of `produce_jammed_configuration`. More specifically:
  + Line 2: Values of ``\ell`` and bound on the displacements, ``s_{bound}``, and their comparison with the particles' radius, `R`. Note that the value of `R` is taken from the previous iteration, *i.e.* the one used to formulate the LP instance.
  + Line 3: Status of the solution of the LP problem obtained with the solver chosen, and time (in seconds) required to solve such problem. The status informs on whether the LP instance was solved to optimality or not. In other words, it the optimizer converged to a solution. 
    
    One reason optimality might not be attained, for instance, is that the iteration or time limits of a given solver are reached. Besides, depending on which solver you are using, you can also get information about the infeasibility of the problem, presence of duality gap, etc.

  + Line 4: The optimal value of the inflation factor, ``\Gamma^\star`` (actually ``\sqrt{\Gamma^\star}-1``), and the largest magnitude of the optimal displacements, ``\max \{ | \mathbf{s}_i^\star|\}_{i=1}^N``. 
  
    **Note:** when computing this second quantity, *all* the particles are considered. Thus, it is very likely that the value reported correspond to ``|\mathbf{s}^\star|`` of **a rattler**. Because this particles are mostly unconstrained by their neighbours, specially during the first iterations, their value of ``|\mathbf{s}^\star|`` might *saturate* `s_{bound}`. Indeed, compare these two values in the example output above.

+ Line 5: Largest magnitude of optimal displacements, bot now *only stable particles* are considered.
+ Line 6: Values of the packing fraction, ``\varphi``, and particles' radius ``R``, updated after ``\Gamma^\star`` has been found.
+ Line 7: 
  + Statistics (maximum value, mean and standard deviation) of the number of constraints induced per particle.
  + Whether the solution of the LP instance yielded an isostatic configuration or not.
  + Mean and maximum number of (possibly *linear*) contacts in stable particles; as well as the total number of contacts associated to rattlers.
  + Amount of *non*-rattlers.
+ Line 8: separator
+ Line 9: The mismatch of the force balance condition, measured as the largest total force vector. That is, the sum of forces acting on each particle is computed for all of them, and the reported value is the vector of largest magnitude.
+ Line 10: Sample of the 10 smallest forces magnitudes present in the configuration. (Note that before convergence is reached, these are not true *contact* forces.)

Lines 5-10 are the output of [`print_monitor_progress`](@ref Main.CALiPPSO.print_monitor_progress), called during the [main loop](@ref mainloop) of `produce_jammed_configuration`, at fixed intervals. Importantly, if `verbose=false` is passed to `produce_jammed_configuration`, *none* of this output will be printed.


As explained [later](@ref Possible-issues), the information of this output is very useful for problem solving, or at least identifying when an issue arises in a given configuration.


### Output when `produce_jammed_configuration` converges

When `produce_jammed_configuration` terminates, because the convergence criteria are met, the following output is given.

```
L1:	    CALiPPSO converged!
L2: Iterations to convergence = 61,	√Γ-1 = -2.4424906541753444e-15,	 Max displacement = 4.5795e-11,	 (ϕ, R) = [0.64414472, 0.02109397]
L3: 	Non_rattlers= 16000	% of rattlers = 2.344	(Max, mean±std) constraints per particle per d: [15.0, 6.26, 3.94]
L4: ┌ Warning: Force balance condition is NOT satisfied! Max force mismatch = 3.714346582469644e-12 ; 	Creating the packing anyway.
L5: └ @ Main /media/storage/Dropbox/Work/Roma-PostDoc/iLP-for-HS-jamming/Scripts/iLP/packing-type.jl:185
L6:
L7: Force balance condition not met yet; max force mismatch = 3.714346582469644e-12
L8: Performing a last LP optimization to improve force balance.
L9: LP instance generated with (ℓ, ℓ/R) = [0.056953715966077116, 2.7]	 and UNbounded displacements
L10: Max optimal displacement = 3.515092848501867e-10	 √Γ-1 = 7.327471962526033e-15
L11: Force mismatch after final optimization = 2.3502518894751413e-13
L12:
L13:
L14: Isostaticity achieved: true	 Non-rattlers = 16000	% of rattlers = 2.344	N_contacts = 47998.0	 z in rattlers = 0
L15: Maximum force equilibrium mismatch = 2.3502518894751413e-13
L16: Time to finish = 299.49 minutes;	 Memory allocated (GB): 293.17
L17: Checking for overlaps after convergence...
L18: No overlaps found! :D
```
**Note:** The line numbers (on the left) are only included here for reference, but are not printed.

Lines 1-3, and 17-18 are always printed (because [`print_converged`](@ref Main.CALiPPSO.print_converged) and [`check_for_overlaps`](@ref) are always called), while the rest of the text is only provided if `verbose=true` is used as argument of the main function. They contain the following information:

+ Line 1: Text indicating that CALiPPSO did converge.
+ Line 2: Information about convergence; specifically:
  + Number of LP optimizations required for convergence
  + Values of ``\sqrt{\Gamma^\star}-1`` and ``\max \{ | \mathbf{s}_i^\star|\}_{i=1}^N`` (only of *stable* particles) obtained in the last iteration. These values are important because they define the convergence criteria of `produce_jammed_configuration`.
  + Final values of the packing fraction and radius, ``\varphi`` and ``R``.
+ Line 3: Amount of non-rattlers (*i.e.* stable particles) and their percentage out of ``N``. Besides, provides the same statistics of the constraints per particle [as Line 7 above](@ref output-process).
+ Lines 4-13: These lines are not present in general. Instead, they will only be printed if mechanical equilibrium of the final configuration is *not* met (or is satisfied with a very low accuracy).
  + By far the most common reason this happens is because the bounds imposed on ``{|\mathbf{s}_i|}_{i=1}^N``, when calling [`solve_LP_instance`](@ref Main.CALiPPSO.solve_LP_instance). As we explain in our paper, a simple solution is to use such final configuration to generate a final LP instance and solve it, *without* any bounds on the displacements. Thus,
  + Lines 4-5: A [warning](https://docs.julialang.org/en/v1/stdlib/Logging/) thrown by the [`MonoPacking`](@ref) [*constructor*](@ref Constructors) when creating the packing. They indicate, precisely, that the force balance condition is not met and give the maximal violation to it.
  + Lines 7-11: The result of solving final LP instance to attain force balance with a higher accuracy. They print some info, such as the new values of ``\Gamma^\star``,  ``\max \{ | \mathbf{s}_i^\star|\}_{i=1}^N``, new value of force balance mismatch, etc. They are printed when [`fine_tune_forces!`](@ref CALiPPSO.fine_tune_forces!) is called. More details about this can be found in the documentation about [Creating the final packing](@ref).
+ Lines 14-16: They provide some extra information about the packing obtained after convergence, such as whether isostaticity was reached or not; number of contacts, etc. They are printed through [`print_info_convergence`](@ref CALiPPSO.print_info_convergence)
+ Lines 17-18: Check for overlaps in the final packing; printed when calling [`check_for_overlaps`](@ref).


!!! warning "Mechanical equilibrium not satisfied"
    As mentioned above when explaining Lines 4-13, the likeliest reason a packing that fails to satisfy mechanical equilibrium is small mismatch (about ``10^{-10}`` or smaller) in foce balance. However, it may also happen that a configuration indeed fails to fulfil force balance due to a poor solution of an LP instance. When this is the case, the error is much bigger (about ``10^{-2}``or even larger) and even calling `fine_tune_forces!` is not enough to correct for it. In fact, in this situation, force balance is violated at many iterations and not only the last one. (See Lines 9-10 of [the previous section](@ref output-process).)

    As described in the [Problem solving](@ref) section, there are a variety of reason why this may happen. But the most probable one, at least according to our experience, is that the optimizer used lacks the precision to solve the linear optimization problems.


### Output when `produce_jammed_configuration` *fails* to converge

TO BE ADDED

## Some other features
- The dimensionality of the system is inferred from `Xs` and the periodic boundary conditions are automatically implemented through the usage of [`PeriodicNumber`](@ref). Of course: **no overlap should be present in the initial configuration** for CALiPPSO to run properly. 

- You can (or at least should be able to) use as input any valid hard-sphere configuration generated from your favourite method (for instance, the Lubachevsky─Stillinger (LS) compression protocol as described [before](@ref The-initial-conditions)).
- Alternatively, you can also use the function [`generate_random_configuration`](@ref)`(d, N, ϕ)` provided here to generate a random *low-density* initial configuration of `N` particles in `d` dimensions with density `ϕ`. See however [the possible caveat](@ref The-initial-conditions) of initializing CALiPPSO with a configuration of low density (*i.e. far from jamming).
- As CALiPPSO progresses, checks of the absence of overlaps are implemented automatically.
- Stability and force balance checks are implemented. They are useful to track possible numerical issues after each of the LP optimizations are carried out; the details are given [above](@ref output-process) and how to solve numerical issues is discussed [here](@ref Problem-solving).
- The cutoff or radius of influence, ``\ell``, is automatically adjusted in such a way that only nearby pairs of particles are considered when building the set of non-overlapping constraints.
- The behaviour and other parameters of the main function can be easily controlled through [keyword arguments](@ref kwargs-control).




## Getting help

We tried our best to provide complete *docstrings* of all the functions defined in this package. So most of the information should be available simply by calling the respective documentation (*i.e.* just by typing '?'). For instance, try typing `?produce_jammed_configuration` in the REPL or in a Jupyter notebook for a detailed description of this function. 

We also tried to leave clarifying comments throughout the code, so its functioning is easier to understand.

If you find a problem (*e.g.* lack of convergence of the main function, low precision, non-isostaticity of the packings), please consult the [Problem solving](@ref) section.

Still there is something confusing ot not clear enough? Found a bug or an issue? Please [drop us an email](mailto:rafael.diazhernandezrojas@uniroma1.it).