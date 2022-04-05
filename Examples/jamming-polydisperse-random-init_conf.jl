using CALiPPSO  
using Statistics, Random
using CALiPPSO.Distributions

Random.seed!(456) # Just for reproducibility

#= UNcomment the following lines if you want to use Gurobi Solver, for running the tests.
    Or adjust them to the solver of your choice
=#
# using Gurobi
# const solver = Gurobi
# const tol_overlap=1e-8 # tolerance for identifying an Overlap. Given that Gurobi's precision is 10^-9, a larger value is needed.
# const tol_optimality = 1e-9; # optimality tolerances. This is the most precise value allowed by Gurobi
# const solver_args = Gurobi.Env()
# const solver_attributes = Dict("OutputFlag" => 0, "FeasibilityTol" => tol_optimality, "OptimalityTol" => tol_optimality, "Method" => 3, "Threads" =>  CALiPPSO.max_threads)


#= Comment the following lines if you want to use CALiPPSO with a different solver than GLPK; e.g., Gurobi (see above)=#
const solver = CALiPPSO.default_solver
const tol_overlap = CALiPPSO.default_tol_overlap
const tol_optimality = CALiPPSO.default_tol_optimality
const solver_args = CALiPPSO.default_args
const solver_attributes = CALiPPSO.default_solver_attributes



precompile_main_function(solver, solver_attributes, solver_args)

const Nrand = 512 # size of configurations to test
const L = 1.0 # size of each side of the system's volume


##################################################
##################################################
# First a test on a 2d BIdisperse configuration
##################################################
##################################################
const d_bi = 2 # dimension for testing BIdisperse packings
const N1 = ceil(Int, Nrand/2)
const ϕ_bi = 0.4
r1, r2, Xs0_bi = generate_random_configuration(d_bi, N1, N1, ϕ_bi, 1.4, L); # this creates a bidisperse configuration  with radii ratio of r2/r1= 1.4
ℓ0_bi = min(L, 6*r2)

printstyled("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)
printstyled("Using CALiPPSO to jam a BIdisperse system of N = ", 2N1, " in d = ", d_bi, "\t Initial (φ, R1, R2) = ", [ϕ_bi, r1, r2], "\n", bold=:true, color=:blue)
printstyled("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)

Rs_bi = [r1*ones(N1); r2*ones(N1)]

@time jammed_packing_bi, info_convergence_bi, Γs_vs_t_bi, smax_vs_t_bi, iso_vs_t_bi = produce_jammed_configuration!(Xs0_bi, Rs_bi; ℓ0=ℓ0_bi, sqrΓ0=1.5,  max_iters=1000, initial_monitor=20, solver=solver, solver_attributes=solver_attributes, solver_args=solver_args) 
println("_______________________________________________________________________________________\n\n")
times = info_convergence_bi.times_LP_optim
println("Info about LP times")
println("(min, avg ± std, max) LP times:\t", minimum(times), ", \t", mean(times), " ± ", std(times), ", \t", maximum(times))
println("_______________________________________________________________________________________\n\n\n")

network_of_contacts(jammed_packing_bi)


##################################################
##################################################
# Now fully polydisperse packings with the radii drawn from a Log-normal distribution
##################################################
##################################################
const ds = [2, 3, 4] # dimensions to use with POLYdisperse packings
const σRs = sqrt(log(1+0.2^2)) # this parameter determines the width of the radii distribution
const distr_Rs = LogNormal(0, σRs) # radii will be drawn from a LogNormal distribution
const scl_facs = [0.017, 0.05, 0.09]; # scale factors to obtain a feasible and good initial value of φ

for (nd, d) in enumerate(ds)
    # Random.seed!(456) # Just for reproducibility
    
    Rs = sort!(scl_facs[nd].*rand(distr_Rs, Nrand), rev=true) # by sorting by the largest values, it's easier to obtain a higher density, initial configuration. However, this could introduce unwanted correlations. So use it with care
    ϕ, Xs0 =generate_random_configuration(d, Rs, L)

    μR = mean(Rs); Rmin = minimum(Rs); Rmax = maximum(Rs)
    ℓ0 = min(L, 6*Rmax)

    printstyled("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)
    printstyled("Using CALiPPSO to jam a POLYdisperse system of N = ", Nrand, " in d = ", d, "\t Initial (φ, Rₘᵢₙ, ⟨R⟩, Rₘₐₓ) = ", [ϕ, Rmin, μR, Rmax], "\n", bold=:true, color=:blue)    
    printstyled("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", color=:blue)


    @time jammed_packing, info_convergence, Γs_vs_t, smax_vs_t, iso_vs_t = produce_jammed_configuration!(Xs0, Rs; ℓ0=ℓ0, sqrΓ0=1.5,  max_iters=3000, initial_monitor=20, solver=solver, solver_attributes=solver_attributes, solver_args=solver_args) 
    println("_______________________________________________________________________________________\n\n")
    local times = info_convergence.times_LP_optim
    println("Info about LP times")
    # println(times)
    println("(min, avg ± std, max) LP times:\t", minimum(times), ", \t", mean(times), " ± ", std(times), ", \t", maximum(times))
    println("_______________________________________________________________________________________\n\n\n")

end