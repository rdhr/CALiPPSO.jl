# Tests and examples provided

## Try it yourself! 

We have included a couple of scripts so you can test that CALiPPSO is working properly. All of them are contained in the `Examples` folder, so be sure to `cd` into it before executing the commands below. Besides, be sure to have installed all the [required dependencies](#dependencies) and as well as the [solvers](#changing-the-solveroptimizer) you want to test.

Besides, consider that the [main function](#basic-usage) is defined in the file `CALiPPSO.jl`, so this file should be included in any script where `produce_jammed_configuration!` is called (with the usual syntax `include("CALiPPSO.jl")`). Now, when this script is first called through `include`, you'll see a rather long output printed out. As its first line explains, this is because a first compilation is being carried out, but you can safely ignore all of such output.

### Some examples included

We provide the following examples to test our algorithm under different circumstances

1. Testing CALiPPSO from a highly compressed configuration of $N=1204$ particles (obtained from the [LS protocol](#the-initial-conditions)) in $d=3, 4, 5$.
   - To execute this test, simply execute from the terminal
      ```
      julia test-monodisperse-after-LS.jl
      ```
   - For running this script, you'll need to have the (included) files `Centers-after-LS--N-1024--d-3.dat`, etc. in the same folder.

2. Testing CALiPPSO from a random, [low density](#the-initial-conditions) initial configuration of $N=512$ particles, in $d=2,3,4,5$.
   - Similarly, you'll need to type
       ```
      julia test-monodisperse-random_init_conf.jl
      ```
   - This script does not depend on any input file as initial condition. But be sure to have the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package already installed.

### Testing different solvers

3. Testing CALiPPSO with [different solvers](#changing-the-solveroptimizer).
   - For running this test, type:  
      ```
      julia tests-different-solvers.jl
      ``` 
   - For simplicity the same initial condition is given in all cases, namely `Centers-after-LS--N-1024--d-3.dat`, but if you want to use the analogous files corresponding to different dimensions, you just need to change the value of `d` in the `tests-different-solvers.jl` file.
   - As it is, this test is meant to work with the [solvers mentioned above](#changing-the-solveroptimizer), *i.e.* `Gurobi`, `HiGHS`, `GLPK`, `Clp`, `Hypatia`, and `COSMO`. So be sure to have installed all of them or, in any case, delete the corresponding entries of the missing ones from  the arrays `solvers`, `solvers_attributes`, and `solvers_args` defined in `tests-different-solvers.jl`.
   - If you want to add another solver, or change the default options, be sure to read the [changing the solver/optimizer section](#changing-the-solveroptimizer) above.

In all cases, `verbose` has been set to true, so you can see al the provided info about the progress and convergence status of CALiPPSO directly from the terminal.


## A simple example and understanding the output
