# ToDo's:

- [ ] Add implementation for polydisperse packigns
- [ ] Try to implement model creation with `direct_model` instead of `Model` (to avoid copying the model).
  - Currently, when `direct` model is used it throws the error `The solver does not support an objective function of type MathOptInterface.SingleVariable.` when using HiGHS (and possibly with other solvers similar errors occur). 
- [X] Add section (or at least explanation) about rattlers
- [ ] Add file for installing the required dependencies
- [ ] Make a consistent use of "solver" and "optimizer"; follow JuMP's terminology.
- [X] Remove `monitor_force_balance` from kwargs of main function
- [X] Clarify that `precompile_main_function` can be set to `false` to avoid pre-compiling
- [X] Add precompiling behaviour conditioned on whether `precompile_main_function` has been defined or not; use `@isdefined`.
- [X] Once the point above has been resolved, removed all the output whenever `precompile_main_function=false`. This means avoid creating the first packings that throw warnings due to the absence of force balance.
- [ ] Finish the documentation
- [ ] Check all links