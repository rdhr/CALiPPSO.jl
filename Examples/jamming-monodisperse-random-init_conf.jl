using CALiPPSO  
using Statistics, Random

Random.seed!(456) # Just for reproducibility

#= UNcomment the following lines if you want to use Gurobi Solver, for running the tests.
    Or adjust them to the solver of your choice
=#
# using Gurobi
# const tol_overlap=1e-8 # tolerance for identifying an Overlap. Given that Gurobi's precision is 10^-9, a larger value is needed.
# const tol_optimality = 1e-9; # optimality tolerances. This is the most precise value allowed by Gurobi
# const grb_env = Gurobi.Env()
# const optimizer = Gurobi.Optimizer(grb_env)
# const solver_attributes = Dict("OutputFlag" => 0, "FeasibilityTol" => tol_optimality, "OptimalityTol" => tol_optimality, "Method" => 3, "Threads" =>  CALiPPSO.max_threads)


#= Comment the following lines if you want to use CALiPPSO with a different solver than GLPK; e.g., Gurobi (see above)=#
const optimizer = CALiPPSO.default_solver.Optimizer()
const tol_overlap = CALiPPSO.default_tol_overlap
const tol_optimality = CALiPPSO.default_tol_optimality
const solver_attributes = CALiPPSO.default_solver_attributes    


precompile_main_function(optimizer, solver_attributes)

const Nrand = 512 # size of configurations to test
const L = 1.0 # size of each side of the system's volume
const ds = [3, 4, 5] # dimensions to use
const ϕs_ds = [0.3, 0.15, 0.1]



# Testing CALiPPSO in monodisperse systems of different dimensions.
# The initial conditions are random configurations of low density, i.e. very far from φ_J. 
# The initial values of φ we're using have been chosen in such a way that the desired configuration can be created with ease
for (nd, d) in enumerate(ds)
    ϕ = ϕs_ds[nd]

    r0, Xs0 = generate_random_configuration(d, Nrand, ϕ, L); 
    ℓ0 = 6*r0
    
    printstyled("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)
    printstyled("Using CALiPPSO to jam a system of N = ", Nrand, " in d = ", d, "\t Initial (φ, R) = ", [ϕ, r0], "\n", bold=:true, color=:blue)
    printstyled("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)

    if d==5
        max_iters =5000
    else
        max_iters = 1000
    end

    @time jammed_packing, info_convergence, Γs_vs_t, smax_vs_t, iso_vs_t = produce_jammed_configuration!(Xs0, r0; ℓ0=ℓ0, sqrΓ0=1.5,  max_iters=max_iters, initial_monitor=20, optimizer=optimizer, solver_attributes=solver_attributes) 
    println("_______________________________________________________________________________________\n\n")
    times = info_convergence.times_LP_optim
    println("Info about LP times")
    # println(times)
    println("(min, avg ± std, max) LP times:\t", minimum(times), ", \t", mean(times), " ± ", std(times), ", \t", maximum(times))
    println("_______________________________________________________________________________________\n\n\n")


    network_of_contacts(jammed_packing)
end