# Some ToDo's

- [X] Add implementation for polydisperse packings. :heavy_check_mark:
- [X] Add support for Mosek solver. :heavy_check_mark:
- [ ] Add implementation with *closed* boundary conditions
- [ ] Try to implement model creation with `direct_model` instead of `Model` (to avoid copying the model).
  - Currently, when `direct` model is used it throws the error `The solver does not support an objective function of type MathOptInterface.SingleVariable.` when using HiGHS (and possibly with other solvers similar errors occur). 
- [ ] Finish the *Issues* section of the documentation
- [ ] Try adding examples using [Literate.jl](https://fredrikekre.github.io/Literate.jl/v2/)
- [ ] Add functions to analyse packings (e.g. compute gaps, extract forces --distinguishing bucklers--, etc.)
