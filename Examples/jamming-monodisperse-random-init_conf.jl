include("../src/CALiPPSO.jl")
using .CALiPPSO  
using DelimitedFiles, Statistics, Random

include("../src/random_initial_conditions.jl")
Random.seed!(456)

#= UNcomment the following lines if you want to use Gurobi Solver, for running the tests.
    Or adjust them to the solver of your choice
=#
using Gurobi
const solver = Gurobi
const tol_overlap=1e-8 # tolerance for identifying an Overlap. Given that Gurobi's precision is 10^-9, a larger value is needed.
const tol_optimality = 1e-9; # optimality tolerances. This is the most precise value allowed by Gurobi
const solver_args = Gurobi.Env()
const solver_attributes = Dict("OutputFlag" => 0, "FeasibilityTol" => tol_optimality, "OptimalityTol" => tol_optimality, "Method" => 3, "Threads" =>  CALiPPSO.max_threads)


#= UNcomment the following lines if you want to use CALiPPSO default's solver and behaviour, *i.e.* GLPK =#
# const solver = CALiPPSO.default_solver
# const tol_overlap = CALiPPSO.default_tol_overlap
# const tol_optimality = CALiPPSO.default_tol_optimality
# const solver_args = CALiPPSO.default_args
# const solver_attributes = CALiPPSO.default_solver_attributes    


precompile_main_function(solver, solver_attributes, solver_args)

const Nrand = 512 # size of configurations to test
const L = 1.0 # size of each side of the system's volume
# const ds = [2, 3, 4, 5] # dimensions to use
# const ϕs_ds = [0.4, 0.3, 0.15, 0.1]
const ds = [3, 4, 5] # dimensions to use
const ϕs_ds = [0.3, 0.15, 0.1]



# Testing CALiPPSO in systems of different dimensions.
# The initial conditions are random configurations of low density, i.e. very far from φ_J. 
# The values of φ I'm using have been chosen in such a way that the desired configuration can be created with ease
for (nd, d) in enumerate(ds)
    ϕ = ϕs_ds[nd]
    
    r0, Xs0 = generate_random_configuration(d, Nrand, ϕ); Xsc = copy(Xs0); rc = copy(r0)
    σ0 = 2*r0

    ℓ0 = 6*r0

    printstyled("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)
    printstyled("Using CALiPPSO to jam a system of N = ", Nrand, " in d = ", d, "\t Initial (φ, R) = ", [ϕ, r0], "\n", bold=:true, color=:blue)
    printstyled("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)

    if d==2 || d==5
        max_iters =2000
    else
        max_iters = 1000
    end

    if d==2
        zero_force = 1e-6
    else
        zero_force = CALiPPSO.default_tol_zero_forces
    end

    @time jammed_packing, info_convergence, Γs_vs_t, smax_vs_t, iso_vs_t = produce_jammed_configuration!(Xs0, r0; ℓ0=ℓ0, sqrΓ0=1.5, non_iso_break=50, max_iters=max_iters, initial_monitor=20, zero_force=zero_force, solver=solver, solver_attributes=solver_attributes, solver_args=solver_args) # Allow a longer streak of non-isostatic solutions because this could be a common situation during the initial LP steps when using low density configurations
    println("_______________________________________________________________________________________\n\n")
    times = info_convergence.times_LP_optim
    println("Info about LP times")
    # println(times)
    println("(min, avg ± std, max) LP times:\t", minimum(times), ", \t", mean(times), " ± ", std(times), ", \t", maximum(times))
    println("_______________________________________________________________________________________\n\n\n")


    network_of_contacts(jammed_packing)
end